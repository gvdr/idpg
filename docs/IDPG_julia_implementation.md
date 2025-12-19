# IDPG Julia Implementation Specification

## Overview

This document specifies the Julia implementation needed to support the IDPG (Intensity Dot Product Graph) paper. The implementation will:

1. Sample from Poisson point processes on the simplex
2. Generate graphs under both node-centric and edge-centric interpretations
3. Evolve intensities via PDEs (diffusion, advection, reaction-diffusion)
4. Visualize intensity landscapes, sampled graphs, and statistics over time
5. Validate theoretical formulas empirically

---

## Recommended Julia Packages (2024-2025)

### Core Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| **Graphs.jl** | ≥1.9 | Graph data structures, basic graph algorithms |
| **JumpProcesses.jl** | ≥9.17 | Point process simulation (PPP, Hawkes) via SciML |
| **PointProcesses.jl** | ≥0.5 | Abstract interface for temporal point processes |
| **DifferentialEquations.jl** | ≥7.x | ODE/PDE solving ecosystem |
| **MethodOfLines.jl** | latest | Automated finite-difference PDE discretization |
| **ModelingToolkit.jl** | latest | Symbolic PDE specification |
| **Distributions.jl** | ≥0.25 | Probability distributions for sampling |

### Visualization

| Package | Purpose |
|---------|---------|
| **Makie.jl** / **CairoMakie.jl** | High-quality static plots for paper figures |
| **GLMakie.jl** | Interactive 3D visualization of simplex and intensity |
| **GraphMakie.jl** | Graph visualization integrated with Makie |

### Utilities

| Package | Purpose |
|---------|---------|
| **StaticArrays.jl** | Fast small vectors for simplex coordinates |
| **LinearAlgebra** | Standard library for SVD, dot products |
| **StatsBase.jl** | Statistical summaries |
| **Random** | RNG control for reproducibility |

---

## Module Structure

```
IDPG/
├── src/
│   ├── IDPG.jl              # Main module
│   ├── simplex.jl           # Simplex geometry utilities
│   ├── intensity.jl         # Intensity function representations
│   ├── sampling.jl          # PPP sampling on simplex
│   ├── graph_generation.jl  # Node-centric and edge-centric graph generation
│   ├── pde_evolution.jl     # PDE dynamics on intensity
│   ├── inference.jl         # RDPG inference + density estimation
│   └── visualization.jl     # Plotting utilities
├── test/
│   ├── test_simplex.jl
│   ├── test_sampling.jl
│   ├── test_formulas.jl     # Validate E[N], E[|E|] formulas
│   └── test_pde.jl
├── examples/
│   ├── basic_idpg.jl
│   ├── product_case.jl
│   ├── diffusion_example.jl
│   └── ecological_example.jl
└── Project.toml
```

---

## Detailed Component Specifications

### 1. Simplex Geometry (`simplex.jl`)

```julia
"""
The (d-1)-dimensional probability simplex in R^d.
"""
module SimplexGeometry

using StaticArrays, LinearAlgebra

# Type alias for points on simplex
const SimplexPoint{d} = SVector{d, Float64}

"""
Check if point lies on the simplex (within tolerance).
"""
function on_simplex(x::AbstractVector; tol=1e-10)
    all(x .>= -tol) && abs(sum(x) - 1.0) < tol
end

"""
Project a point in R^d onto the simplex.
Uses algorithm from "Efficient Projections onto the ℓ1-Ball" (Duchi et al.)
"""
function project_to_simplex(x::AbstractVector)
    # Implementation of simplex projection
end

"""
Sample uniformly from the (d-1)-simplex using Dirichlet(1,1,...,1).
"""
function uniform_simplex_sample(d::Int)
    # Exponential trick: sample d exponentials, normalize
    e = randexp(d)
    return SVector{d}(e ./ sum(e))
end

"""
Barycentric coordinates utilities.
"""
function barycentric_to_cartesian(bary::AbstractVector, vertices::AbstractMatrix)
    # Convert barycentric to Cartesian for plotting
end

"""
Compute dot product (connection probability) for two simplex points.
"""
connection_probability(g::SimplexPoint, r::SimplexPoint) = dot(g, r)

end # module
```

### 2. Intensity Functions (`intensity.jl`)

```julia
"""
Representations of intensity functions ρ on G × R.
"""
module IntensityFunctions

using Distributions, LinearAlgebra

"""
Abstract type for intensity functions on Ω = G × R ⊆ Δ^{d-1} × Δ^{d-1}.
"""
abstract type AbstractIntensity{d} end

"""
Product intensity: ρ(g,r) = ρ_G(g) · ρ_R(r)
"""
struct ProductIntensity{d, FG, FR} <: AbstractIntensity{d}
    ρ_G::FG  # Function or distribution on G
    ρ_R::FR  # Function or distribution on R
end

"""
Evaluate intensity at a point.
"""
function (ρ::ProductIntensity)(g, r)
    return ρ.ρ_G(g) * ρ.ρ_R(r)
end

"""
Gaussian mixture on simplex (common parametric family).
Each component is a Dirichlet distribution.
"""
struct DirichletMixture{d}
    weights::Vector{Float64}
    components::Vector{Dirichlet{Float64}}
    scale::Float64  # Total intensity (not probability!)
end

function (dm::DirichletMixture)(x)
    p = sum(w * pdf(comp, x) for (w, comp) in zip(dm.weights, dm.components))
    return dm.scale * p
end

"""
Compute marginal quantities for product intensity.
"""
function marginal_stats(ρ::ProductIntensity{d}; n_samples=10000) where d
    # Monte Carlo integration for:
    # c_G = ∫ ρ_G(g) dg
    # c_R = ∫ ρ_R(r) dr  
    # μ_G = ∫ g · ρ_G(g) dg
    # μ_R = ∫ r · ρ_R(r) dr
end

"""
Time-varying intensity for PDE evolution.
"""
struct TimeVaryingIntensity{d, F} <: AbstractIntensity{d}
    ρ::F  # Function (g, r, t) -> intensity
end

end # module
```

### 3. Poisson Point Process Sampling (`sampling.jl`)

```julia
"""
Sample interaction sites from PPP(λ) on Ω.
"""
module PPPSampling

using ..SimplexGeometry, ..IntensityFunctions
using Distributions, Random

"""
Sample from inhomogeneous PPP on simplex using thinning.

Algorithm:
1. Find λ_max = sup ρ(g,r) over Ω
2. Sample N ~ Poisson(λ_max · |Ω|) candidate points uniformly on Ω
3. Accept each point (g,r) with probability ρ(g,r) / λ_max
"""
function sample_ppp(ρ::AbstractIntensity{d}; 
                    λ_max=nothing, 
                    rng=Random.default_rng()) where d
    
    # Estimate λ_max if not provided
    if isnothing(λ_max)
        λ_max = estimate_max_intensity(ρ)
    end
    
    # Volume of Ω = Δ^{d-1} × Δ^{d-1}
    # Volume of (d-1)-simplex is 1/((d-1)!)
    vol_simplex = 1.0 / factorial(d-1)
    vol_Ω = vol_simplex^2
    
    # Sample number of candidates
    n_candidates = rand(rng, Poisson(λ_max * vol_Ω))
    
    # Sample candidates uniformly and thin
    accepted = Vector{Tuple{SVector{d,Float64}, SVector{d,Float64}}}()
    
    for _ in 1:n_candidates
        g = uniform_simplex_sample(d)
        r = uniform_simplex_sample(d)
        
        # Accept with probability ρ(g,r) / λ_max
        if rand(rng) < ρ(g, r) / λ_max
            push!(accepted, (g, r))
        end
    end
    
    return accepted
end

"""
Sample from product intensity more efficiently.
"""
function sample_ppp_product(ρ::ProductIntensity{d}; rng=Random.default_rng()) where d
    # Sample from ρ_G and ρ_R independently
    # More efficient than joint thinning
end

"""
For temporal extension: sample from time-varying intensity.
"""
function sample_ppp_temporal(ρ::TimeVaryingIntensity, t_start, t_end; dt=0.01)
    # Return list of (g, r, t) tuples
end

end # module
```

### 4. Graph Generation (`graph_generation.jl`)

```julia
"""
Generate graphs from sampled interaction sites.
"""
module GraphGeneration

using Graphs, ..SimplexGeometry
using Random

"""
Node-centric graph generation.

Given sampled sites [(g₁,r₁), (g₂,r₂), ...], create graph where
edge i→j exists with probability gᵢ · rⱼ.
"""
function generate_node_centric(sites::Vector{Tuple{SVector{d,Float64}, SVector{d,Float64}}};
                                rng=Random.default_rng()) where d
    n = length(sites)
    g = SimpleDiGraph(n)
    
    for i in 1:n
        for j in 1:n
            g_i, _ = sites[i]
            _, r_j = sites[j]
            p = dot(g_i, r_j)  # Connection probability
            
            if rand(rng) < p
                add_edge!(g, i, j)
            end
        end
    end
    
    return g, sites  # Return graph and node positions
end

"""
Edge-centric interpretation.

Each sampled site (g,r) IS an edge from source-at-g to target-at-r.
The "graph" is a collection of directed edges in continuous space.

For visualization/analysis, we may discretize or cluster.
"""
struct EdgeCentricSample{d}
    sources::Vector{SVector{d, Float64}}  # g coordinates
    targets::Vector{SVector{d, Float64}}  # r coordinates
end

function generate_edge_centric(sites::Vector{Tuple{SVector{d,Float64}, SVector{d,Float64}}};
                               rng=Random.default_rng()) where d
    # Each site (g,r) becomes an edge with probability g·r
    sources = SVector{d, Float64}[]
    targets = SVector{d, Float64}[]
    
    for (g, r) in sites
        p = dot(g, r)
        if rand(rng) < p
            push!(sources, g)
            push!(targets, r)
        end
    end
    
    return EdgeCentricSample{d}(sources, targets)
end

"""
Discretize edge-centric sample into a graph by clustering.
"""
function discretize_edge_centric(sample::EdgeCentricSample{d}, n_clusters::Int) where d
    # K-means or other clustering on source and target spaces
    # Return weighted graph between clusters
end

end # module
```

### 5. PDE Evolution (`pde_evolution.jl`)

```julia
"""
PDE dynamics on intensity functions.

Key insight: Work on R^d with boundary condition ρ = 0 outside simplex.
"""
module PDEEvolution

using DifferentialEquations, MethodOfLines, ModelingToolkit
using DomainSets

"""
Set up diffusion equation on the simplex.

∂ρ/∂t = D ∇²ρ

with ρ = 0 on boundary of simplex.
"""
function setup_diffusion(d::Int, D::Float64; resolution=50)
    # For d=2 (1-simplex = line segment), this is straightforward
    # For d=3 (2-simplex = triangle), need triangular domain
    
    @parameters t x[1:d-1]  # Use d-1 independent coordinates on simplex
    @variables ρ(..)
    
    Dt = Differential(t)
    Dx = [Differential(x[i]) for i in 1:d-1]
    
    # Laplacian in simplex coordinates
    # (needs care with metric tensor for non-Cartesian coords)
    
    # For now, work in embedding space and project
end

"""
Alternative: Finite-difference on regular grid, mask outside simplex.
"""
struct SimplexGrid{d}
    resolution::Int
    points::Vector{SVector{d, Float64}}
    inside_simplex::BitVector
    neighbors::Vector{Vector{Int}}
end

function create_simplex_grid(d::Int, resolution::Int)
    # Create regular grid in bounding box, mark points inside simplex
end

"""
Evolve intensity via explicit finite differences.
"""
function evolve_diffusion!(ρ_values::Vector{Float64}, 
                           grid::SimplexGrid, 
                           D::Float64, 
                           dt::Float64, 
                           n_steps::Int)
    # Simple explicit Euler with Laplacian stencil
    # Boundary condition: ρ = 0 outside simplex (already zero, stays zero)
end

"""
Advection equation: ∂ρ/∂t = -v · ∇ρ
"""
function evolve_advection!(ρ_values::Vector{Float64},
                           grid::SimplexGrid,
                           v::SVector,  # Velocity field
                           dt::Float64,
                           n_steps::Int)
    # Upwind scheme for stability
end

"""
Track graph statistics as intensity evolves.
"""
function evolve_and_track(ρ_initial, grid, pde_params, t_final; 
                          sample_interval=0.1, n_graph_samples=100)
    # Returns time series of:
    # - E[N](t) (analytical from intensity)
    # - E[|E|](t) (analytical from intensity)  
    # - Empirical N, |E| from sampled graphs
end

end # module
```

### 6. Visualization (`visualization.jl`)

```julia
"""
Visualization utilities for IDPG.
"""
module IDPGVisualization

using Makie, CairoMakie, GLMakie
using GraphMakie, Graphs
using ..SimplexGeometry, ..IntensityFunctions

"""
Plot intensity function on 2-simplex (triangle).
"""
function plot_intensity_2simplex(ρ; resolution=100, colormap=:viridis)
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1,1], aspect=DataAspect())
    
    # Triangle vertices in 2D
    v1 = [0.0, 0.0]
    v2 = [1.0, 0.0]
    v3 = [0.5, sqrt(3)/2]
    
    # Sample intensity on triangle
    # Plot as heatmap or contour
    
    # Draw simplex boundary
    lines!(ax, [v1, v2, v3, v1], color=:black, linewidth=2)
    
    return fig
end

"""
Plot sampled nodes on simplex with graph edges.
"""
function plot_node_centric_graph(sites, graph; ax=nothing)
    # Project simplex points to 2D for visualization
    # Draw nodes colored by position
    # Draw edges
end

"""
Plot edge-centric sample as arrows in simplex × simplex space.
"""
function plot_edge_centric(sample::EdgeCentricSample)
    # 2D: plot source simplex and target simplex side by side
    # Draw arrows from source points to target points
end

"""
Animate intensity evolution.
"""
function animate_evolution(ρ_history, times; filename="evolution.mp4")
    # Use Makie Observables for smooth animation
end

"""
Plot theoretical vs empirical statistics.
"""
function plot_formula_validation(theoretical, empirical; 
                                  quantity=:edges, n_trials=100)
    fig = Figure()
    ax = Axis(fig[1,1], 
              xlabel="Theoretical E[$quantity]",
              ylabel="Empirical mean")
    
    scatter!(ax, theoretical, empirical)
    lines!(ax, [0, maximum(theoretical)], [0, maximum(theoretical)], 
           color=:red, linestyle=:dash, label="y=x")
    
    return fig
end

end # module
```

---

## Validation Tests

### Formula Validation (`test/test_formulas.jl`)

```julia
using Test, Statistics
using IDPG

@testset "Expected node count" begin
    # Product intensity with known c_G, c_R
    d = 3
    ρ = ProductIntensity(DirichletMixture(...), DirichletMixture(...))
    
    # Theoretical E[N]
    c_G, c_R, _, _ = marginal_stats(ρ)
    E_N_theory = c_G * c_R
    
    # Empirical E[N]
    n_trials = 1000
    N_samples = [length(sample_ppp(ρ)) for _ in 1:n_trials]
    E_N_empirical = mean(N_samples)
    
    @test isapprox(E_N_theory, E_N_empirical, rtol=0.05)
end

@testset "Expected edge count - node centric" begin
    # Similar structure
    # E[|E|] = E[N]² · (μ̃_G · μ̃_R)
end

@testset "Expected edge count - edge centric" begin  
    # E[L] = E[N] · (μ̃_G · μ̃_R)
end

@testset "Scaling relationship" begin
    # E[|E|] / E[L] ≈ E[N]
end
```

---

## Example Scripts

### `examples/basic_idpg.jl`

```julia
using IDPG
using CairoMakie

# Set up 3-dimensional latent space (2-simplex)
d = 3

# Define product intensity using Dirichlet mixtures
ρ_G = DirichletMixture(
    weights = [0.6, 0.4],
    components = [Dirichlet([5.0, 1.0, 1.0]), Dirichlet([1.0, 5.0, 1.0])],
    scale = 50.0
)

ρ_R = DirichletMixture(
    weights = [0.5, 0.5],
    components = [Dirichlet([1.0, 1.0, 5.0]), Dirichlet([2.0, 2.0, 2.0])],
    scale = 50.0
)

ρ = ProductIntensity(ρ_G, ρ_R)

# Sample interaction sites
sites = sample_ppp(ρ)
println("Sampled $(length(sites)) sites")

# Generate node-centric graph
graph, positions = generate_node_centric(sites)
println("Generated graph with $(nv(graph)) nodes and $(ne(graph)) edges")

# Visualize
fig = Figure(size=(1200, 500))
ax1 = Axis(fig[1,1], title="Intensity ρ_G")
ax2 = Axis(fig[1,2], title="Intensity ρ_R")  
ax3 = Axis(fig[1,3], title="Sampled Graph")

plot_intensity_2simplex!(ax1, ρ_G)
plot_intensity_2simplex!(ax2, ρ_R)
plot_node_centric_graph!(ax3, positions, graph)

save("basic_idpg.png", fig)
```

### `examples/diffusion_example.jl`

```julia
using IDPG
using GLMakie  # For animation

d = 3

# Initial intensity: concentrated Gaussian blob
ρ_initial = DirichletMixture(
    weights = [1.0],
    components = [Dirichlet([10.0, 2.0, 2.0])],  # Concentrated near vertex 1
    scale = 100.0
)

# Set up diffusion
D = 0.01
grid = create_simplex_grid(d, resolution=50)
ρ_values = [ρ_initial(p) for p in grid.points]

# Evolve and track
results = evolve_and_track(
    ρ_values, grid,
    pde_type = :diffusion,
    D = D,
    t_final = 10.0,
    sample_interval = 0.5,
    n_graph_samples = 50
)

# Animate intensity evolution
animate_evolution(results.ρ_history, results.times, 
                  filename="diffusion.mp4")

# Plot statistics over time
fig = Figure()
ax1 = Axis(fig[1,1], xlabel="Time", ylabel="E[N]")
ax2 = Axis(fig[2,1], xlabel="Time", ylabel="E[|E|]")

lines!(ax1, results.times, results.E_N_theory, label="Theory")
scatter!(ax1, results.times, results.E_N_empirical, label="Empirical")

lines!(ax2, results.times, results.E_E_theory, label="Theory")
scatter!(ax2, results.times, results.E_E_empirical, label="Empirical")

save("diffusion_stats.png", fig)
```

---

## Implementation Priorities

### Phase 1: Core Functionality (Week 1-2)
1. ✅ Simplex geometry utilities
2. ✅ Product intensity with Dirichlet mixtures
3. ✅ PPP sampling via thinning
4. ✅ Node-centric graph generation
5. ✅ Basic visualization on 2-simplex

### Phase 2: Validation (Week 2-3)
1. ✅ Formula validation tests (E[N], E[|E|], E[L])
2. ✅ Scaling relationship tests
3. ✅ Monte Carlo integration for marginal stats

### Phase 3: Edge-Centric & Dynamics (Week 3-4)
1. ✅ Edge-centric generation
2. ✅ Discretization/clustering for edge-centric graphs
3. ✅ Simple diffusion on simplex grid
4. ✅ Time-evolving statistics

### Phase 4: Advanced Features (Week 4+)
1. ⬜ Hawkes process extension (using JumpProcesses.jl)
2. ⬜ Advection and reaction-diffusion
3. ⬜ Inference from observed graphs
4. ⬜ Publication-quality figures

---

## Notes on Boundary Conditions

The key insight for PDE evolution: embed the simplex Δ^{d-1} in ℝ^d and impose **ρ = 0 outside the simplex**. This means:

1. Standard PDE operators (∇², ∇) work on the embedding space
2. No nodes with invalid probabilities are ever sampled (zero intensity there)
3. The simplex boundary acts as an **absorbing boundary** for diffusion
4. For **reflecting boundaries**, need to project flux back into simplex

For the product case, we can evolve ρ_G and ρ_R separately on their respective simplices.

---

## References

- JumpProcesses.jl: https://docs.sciml.ai/JumpProcesses/stable/
- MethodOfLines.jl: https://github.com/SciML/MethodOfLines.jl
- Graphs.jl: https://juliagraphs.org/Graphs.jl/
- Makie.jl: https://docs.makie.org/
- Zagatti et al. (2024) "Extending JumpProcesses.jl for fast point process simulation"
