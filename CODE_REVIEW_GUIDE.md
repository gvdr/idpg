# IDPG Code Review Guide

This document provides a structured guide for reviewing the IDPG (Intensity Dot Product Graph) codebase. It explains the architecture, what each component does, and how the pieces fit together.

---

## Table of Contents

1. [Repository Overview](#repository-overview)
2. [Library Architecture](#library-architecture)
3. [Type Hierarchy](#type-hierarchy)
4. [Data Flow](#data-flow)
5. [Module Reference](#module-reference)
6. [Scripts](#scripts)
7. [Tests](#tests)
8. [Review Checklists](#review-checklists)

---

## Repository Overview

```
idpg/
├── src/                      # Core library
│   ├── IDPG.jl              # Main module (imports, includes, exports)
│   └── modules/             # 7 implementation files
├── scripts/
│   ├── scaling_laws/        # Validation simulations (fig1, fig2, fig6)
│   └── heat_maps/           # Paper figure generation (A-E)
├── test/                    # Unit tests
├── output/
│   ├── scaling_laws/        # Simulation outputs
│   └── heat_maps/           # Figure outputs
└── archive/                 # Old/superseded material
```

### Design Principles

1. **Single module**: All code lives in the `IDPG` namespace (files are `include`d, not submodules)
2. **Dimension-generic**: Types parameterized by `{d}` for 2D, 4D, or higher
3. **Callable structs**: Intensities are evaluated as `rho(g, r)`
4. **Dual representations**: `EdgeCentricSample` (minimal) vs `FullEdgeCentricSample` (complete)
5. **Library-first**: Generic functionality in `src/modules/`, scripts only orchestrate

---

## Library Architecture

### Module Dependency Graph

Files must be included in this order (defined in `IDPG.jl`):

```
IDPG.jl
│
├── LatentSpace.jl       [no dependencies - foundation]
│
├── Intensity.jl         [depends on: LatentSpace]
│
├── PDEEvolution.jl      [depends on: LatentSpace]
│
├── Sampling.jl          [depends on: LatentSpace, Intensity, PDEEvolution]
│
├── GraphGeneration.jl   [depends on: LatentSpace, Sampling]
│
├── Visualization.jl     [depends on: LatentSpace, GraphGeneration]
│
└── EcologicalUtils.jl   [depends on: LatentSpace, GraphGeneration]
```

**Critical**: `PDEEvolution` must precede `Sampling` because `Sampling` uses `BdPlusGrid`.

### External Dependencies

| Module | Key External Packages |
|--------|----------------------|
| LatentSpace | StaticArrays, LinearAlgebra |
| Intensity | Distributions, StatsBase |
| PDEEvolution | LinearAlgebra |
| Sampling | Distributions, Random |
| GraphGeneration | Graphs, Clustering |
| Visualization | CairoMakie, GraphMakie |
| EcologicalUtils | LinearAlgebra |

---

## Type Hierarchy

### Latent Space

```
LatentPoint{d} = SVector{d, Float64}    # Point in B^d_+
```

### Intensity Types

```
AbstractIntensity{d}                     # Base type (abstract)
├── BdPlusMixture{d}                    # Gaussian mixture on B^d_+
│   fields: weights, means, concentrations, scale
│
├── ProductIntensity{d}                 # Factorized: rho(g,r) = rho_G(g) * rho_R(r)
│   fields: rho_G, rho_R (both BdPlusMixture)
│
├── TimeVaryingIntensity{d}             # Time-dependent intensity
│
└── MixtureOfProductIntensities{d}      # Multi-species: sum_m rho_{G,m} * rho_{R,m}
    fields: species (Vector of ProductIntensity)


AbstractEdgeIntensity{d}                 # Base type for edge-centric (abstract)
├── ScaledProductEdgeIntensity{d}       # Scaled for E[N]/2 semantics
│   fields: rho_source, rho_target, C_edge
│
└── SymmetricEdgeIntensity{d}           # Symmetric case
```

### Sampling Types

```
InteractionSite{d}                       # A sampled (g, r) pair
  fields: g::LatentPoint{d}, r::LatentPoint{d}
```

### Graph Types

```
EdgeCentricSample{d}                     # Minimal: only connection-relevant coords
  fields: sources (g coords), targets (r coords)

FullEdgeCentricSample{d}                 # Complete: full site for each endpoint
  fields: source_sites, target_sites (both Vector{InteractionSite})
```

### PDE Types

```
BdPlusGrid{d}                            # Finite difference grid on B^d_+
  fields: points, inside_mask, grid_to_Bd_plus, Bd_plus_to_grid, h
```

---

## Data Flow

### Node-Centric Pipeline

```
Intensity Definition
        │
        ▼
┌─────────────────────────────────────┐
│  ProductIntensity{d}                │
│  or MixtureOfProductIntensities{d}  │
└─────────────────────────────────────┘
        │
        │ sample_ppp_product(rho)
        │ or sample_ppp_mixture(rho)
        ▼
┌─────────────────────────────────────┐
│  Vector{InteractionSite{d}}         │
│  (N sites, where N ~ Poisson)       │
└─────────────────────────────────────┘
        │
        │ generate_node_centric(sites)
        │ [tests all N^2 pairs]
        ▼
┌─────────────────────────────────────┐
│  SimpleDiGraph                      │
│  E[edges] = E[N]^2 * E[g . r]       │
└─────────────────────────────────────┘
```

### Edge-Centric Pipeline

```
Intensity Definition
        │
        ▼
┌─────────────────────────────────────┐
│  ScaledProductEdgeIntensity{d}      │
│  (C_edge = E[N]/2)                  │
└─────────────────────────────────────┘
        │
        │ sample_edge_centric(ei)
        │ [samples E[N]/2 opportunities]
        ▼
┌─────────────────────────────────────┐
│  FullEdgeCentricSample{d}           │
│  E[edges] = (E[N]/2) * E[g . r]     │
└─────────────────────────────────────┘
        │
        │ discretize_edge_centric(sample, n_clusters)
        ▼
┌─────────────────────────────────────┐
│  SimpleDiGraph                      │
│  (clustered representation)         │
└─────────────────────────────────────┘
```

### PDE Evolution Pipeline

```
Initial Intensity (callable)
        │
        │ create_Bd_plus_grid(d, resolution)
        ▼
┌─────────────────────────────────────┐
│  BdPlusGrid{d}                      │
│  + rho_values (on grid points)      │
└─────────────────────────────────────┘
        │
        │ evolve_diffusion!(rho, grid, D, dt, n_steps)
        │ or evolve_advection!(rho, grid, v, dt, n_steps)
        ▼
┌─────────────────────────────────────┐
│  Evolved rho_values on grid         │
└─────────────────────────────────────┘
        │
        │ sample_from_grid(rho_G, rho_R, grid)
        ▼
┌─────────────────────────────────────┐
│  EdgeCentricSample{d}               │
│  or FullEdgeCentricSample{d}        │
└─────────────────────────────────────┘
```

### Key Insight: Edge Accounting

Each edge opportunity "consumes" 2 node-equivalents (1 source + 1 target):

| Model | Pairs/Opportunities | Expected Edges |
|-------|---------------------|----------------|
| Node-centric | N^2 | E[N]^2 * E[g . r] |
| Edge-centric | E[N]/2 | (E[N]/2) * E[g . r] |
| **Ratio** | | **2 * Lambda** |

---

## Module Reference

### 1. LatentSpace.jl

**Purpose**: Geometry of B^d_+ (non-negative unit ball).

**No internal dependencies** - pure geometry foundation.

**Key Exports**:
| Function | Purpose |
|----------|---------|
| `in_Bd_plus(x)` | Check if point is in B^d_+ |
| `on_Bd_plus_boundary(x)` | Check if on boundary |
| `project_to_Bd_plus(x)` | Project onto B^d_+ |
| `Bd_plus_volume(d)` | Volume: pi^(d/2) / (Gamma(d/2+1) * 2^d) |
| `uniform_Bd_plus_sample(d)` | Sample uniformly (type-unstable) |
| `uniform_Bd_plus_sample(Val(d))` | Sample uniformly (type-stable) |
| `connection_probability(g, r)` | Compute g . r |
| `hyperspherical_to_cartesian(r, angles)` | Coordinate transform (type-unstable) |
| `hyperspherical_to_cartesian(Val(d), r, angles)` | Coordinate transform (type-stable) |
| `cartesian_to_hyperspherical` | Inverse coordinate transform |
| `hyperspherical_jacobian(r, angles)` | Jacobian for volume integration |

---

### 2. Intensity.jl

**Purpose**: Intensity functions rho(g, r) - the heart of IDPG theory.

**Depends on**: LatentSpace

**Key Exports**:

*Types*:
- `AbstractIntensity{d}`, `BdPlusMixture{d}`, `ProductIntensity{d}`
- `MixtureOfProductIntensities{d}`, `TimeVaryingIntensity{d}`
- `AbstractEdgeIntensity{d}`, `ScaledProductEdgeIntensity{d}`

*Evaluation (callable structs)*:
```julia
(rho::BdPlusMixture)(x)           # Evaluate at point
(rho::ProductIntensity)(g, r)     # Evaluate product
(mop::MixtureOfProductIntensities)(g, r)  # Evaluate mixture
```

*Statistics*:
| Function | Returns |
|----------|---------|
| `total_intensity(rho)` | E[N] |
| `marginal_total_intensity(rho)` | (c_G, c_R) |
| `intensity_weighted_mean(rho)` | mu = integral of x * rho(x) |
| `normalized_mean(rho)` | mu_tilde = mu / c |
| `marginal_stats(rho)` | Full tuple: (c_G, c_R, mu_G, mu_R, E_N, E_edges_node, E_edges_edge, avg_conn_prob) |

*Desire Operator (spectral structure)*:
| Function | Returns |
|----------|---------|
| `second_moment_matrix(rho)` | Σ where Σ_{jk} = E[x_j x_k] under normalized intensity |
| `desire_operator_singular_values(rho)` | σ_k(D̃) = √(λ_k(Σ_G Σ_R)), sorted descending |
| `desire_stats(rho)` | Named tuple: (Σ_G, Σ_R, σ, μ̃_G, μ̃_R) |

*Construction*:
| Function | Purpose |
|----------|---------|
| `create_calibrated_product_intensity(lambda)` | Build ProductIntensity with target E[N] |
| `sample_from_mixture(bm)` | Sample point from BdPlusMixture |

*1D Utilities*:
| Function | Purpose |
|----------|---------|
| `truncated_gaussian_sigma(x, mu, sigma)` | 1D truncated Gaussian |
| `truncated_gaussian_kappa(x, mu, kappa)` | 1D truncated Gaussian (concentration) |
| `compute_1d_marginals(grid, mu_G, mu_R, sigma)` | Compute marginals on grid |
| `compute_bound_heat_matrix(grid, rho_G, rho_R)` | Compute h_bar(g, r) |

---

### 3. PDEEvolution.jl

**Purpose**: Finite difference PDE solver on B^d_+.

**Depends on**: LatentSpace

**Key Exports**:

*Grid*:
| Function | Purpose |
|----------|---------|
| `create_Bd_plus_grid(d, resolution)` | Create masked FD grid |
| `get_neighbors(grid, idx)` | 2d-neighbor stencil |

*Evolution (explicit Euler, in-place)*:
| Function | PDE |
|----------|-----|
| `evolve_diffusion!(rho, grid, D, dt, n_steps)` | d_rho/dt = D * nabla^2 rho |
| `evolve_advection!(rho, grid, v, dt, n_steps)` | d_rho/dt = -v . nabla rho |
| `evolve_advection_field!(rho, grid, v_field, dt, n_steps)` | Space-dependent v(x) |
| `evolve_reaction_diffusion!(rho, grid, D, f, dt, n_steps)` | d_rho/dt = D * nabla^2 rho + f(rho) |

*Tracking*:
| Function | Purpose |
|----------|---------|
| `evolve_and_track(rho_initial, grid; ...)` | Evolve and record statistics |
| `compute_mean_position(rho_values, grid)` | Intensity-weighted centroid |

---

### 4. Sampling.jl

**Purpose**: Sample sites from Poisson point process.

**Depends on**: LatentSpace, Intensity, PDEEvolution

**Key Exports**:

*Generic PPP*:
| Function | Purpose |
|----------|---------|
| `sample_ppp(rho)` | Thinning algorithm (slow but general) |
| `estimate_max_intensity(rho)` | Find supremum for thinning |

*Specialized PPP (fast)*:
| Function | Purpose |
|----------|---------|
| `sample_ppp_product(rho)` | For ProductIntensity |
| `sample_ppp_mixture(mop)` | For MixtureOfProductIntensities |
| `sample_ppp_mixture_sites_only(mop)` | Discard species labels |

*Grid-based (for evolved intensities)*:
| Function | Returns |
|----------|---------|
| `sample_from_grid(rho_G, rho_R, grid)` | EdgeCentricSample |
| `sample_from_grid_full(rho_G, rho_R, grid)` | FullEdgeCentricSample |
| `initialize_grid_from_mixture(grid, ...)` | Initialize grid values |

*True edge-centric (E[N]/2 semantics)*:
| Function | Returns |
|----------|---------|
| `sample_edge_centric(ei)` | FullEdgeCentricSample |
| `sample_edge_centric_ppp(rho_source, rho_target)` | Convenience wrapper |
| `sample_edge_centric_symmetric(rho)` | Symmetric case |

---

### 5. GraphGeneration.jl

**Purpose**: Generate graphs from sampled sites.

**Depends on**: LatentSpace, Sampling

**Key Exports**:

*Node-centric (N^2 pairs)*:
| Function | Returns |
|----------|---------|
| `generate_node_centric(sites)` | (SimpleDiGraph, sites) |
| `generate_node_centric_full(sites)` | FullEdgeCentricSample |

*Edge-centric (linear)*:
| Function | Returns |
|----------|---------|
| `generate_edge_centric(sites)` | EdgeCentricSample |

*Discretization (clustering)*:
| Function | Purpose |
|----------|---------|
| `discretize_edge_centric(sample, n_clusters)` | K-means clustering |
| `discretize_edge_centric_joint(sample; eps, min_samples)` | DBSCAN |
| `discretize_with_weights(sample, n_clusters)` | Weighted edges |

*Accessors for FullEdgeCentricSample*:
| Function | Returns |
|----------|---------|
| `source_g(sample)` | g coordinates of sources |
| `source_r(sample)` | r coordinates of sources |
| `target_g(sample)` | g coordinates of targets |
| `target_r(sample)` | r coordinates of targets |
| `to_edge_centric(sample)` | Convert to EdgeCentricSample |

---

### 6. Visualization.jl

**Purpose**: Plotting and animation.

**Depends on**: LatentSpace, GraphGeneration, CairoMakie, GraphMakie

**Key Exports**:

| Function | Purpose |
|----------|---------|
| `Bd_plus_to_2d(point)` | Project d-dim to 2D for plotting |
| `draw_Bd_plus_boundary!(ax)` | Draw B^2_+ boundary |
| `plot_intensity_Bd_plus(rho)` | Visualize intensity |
| `plot_intensity_Bd_plus!(ax, rho)` | Mutating version |
| `plot_node_centric_graph(sites, graph)` | Two panels (G, R) |
| `plot_edge_centric(sample)` | Sources and targets |
| `animate_evolution(...)` | PDE animation |
| `latex_figure_theme(; fontsize)` | LaTeX-style theme |
| `with_latex_theme(f; fontsize)` | Theme wrapper |

---

### 7. EcologicalUtils.jl

**Purpose**: Guild structure and food web utilities.

**Depends on**: LatentSpace, GraphGeneration

**Key Exports**:

*Guild assignment*:
| Function | Purpose |
|----------|---------|
| `assign_site_to_guild(site, guild_means)` | Nearest guild in 2d-space |
| `assign_point_to_guild(point, guild_means)` | Nearest guild in d-space |
| `build_full_guild_means(means_G, means_R)` | Concatenate [g; r] |

*Food web*:
| Function | Purpose |
|----------|---------|
| `compute_foodweb_matrix(sample, guild_means)` | Count guild-guild edges |
| `normalize_foodweb_matrix(M; mode)` | Row/col/total normalization |
| `default_trophic_affinity()` | 5x5 trophic matrix |

*Centroid construction*:
| Function | Purpose |
|----------|---------|
| `construct_food_web_centroids()` | SVD + alternating optimization |
| `svd_initialize_centroids(K_star, d)` | SVD initialization |
| `project_rows_to_Bd_plus!(M)` | Project rows to B^d_+ |
| `verify_centroids(means_G, means_R, K_star)` | Check reconstruction |

*Layout*:
| Function | Purpose |
|----------|---------|
| `trophic_layout(n_guilds)` | Positions for food web plots |

---

## Scripts

### Scaling Laws (`scripts/scaling_laws/`)

Validate theoretical predictions about quadratic vs linear scaling.

| Script | Purpose | Key Result |
|--------|---------|------------|
| `fig1_intensity_scaling.jl` | Node vs edge comparison | Ratio = 2 * Lambda |
| `fig2_entity_lifetime.jl` | Lifetime interpolation | p_overlap formula |
| `fig6_temporal_evolution.jl` | PDE-driven evolution | Ratio holds at all t |

### Heat Maps (`scripts/heat_maps/`)

Generate paper figures.

| Script | Content |
|--------|---------|
| `construct_centroids.jl` | Food web centroid optimization |
| `figure_A_anatomy.jl` | Product intensity basics |
| `figure_A_prime_asymmetric.jl` | Asymmetric case |
| `figure_A_double_prime_nonproduct.jl` | Non-product 4D |
| `figure_B_pde_absorbing.jl` | Absorbing BC |
| `figure_B_pde_reflecting.jl` | Reflecting BC |
| `figure_C_spectral_1d.jl` | 1D spectral decomposition |
| `figure_C_spectral_4d.jl` | 4D spectral decomposition |
| `figure_D_foodweb_static.jl` | Static food web |
| `figure_E_foodweb_dynamic.jl` | Dynamic food web |

### Desire Validation (`scripts/desire_validation/`)

Validate the Desire operator spectral consistency theorem: σ_k(A)/N → σ_k(D̃) as N → ∞.

| Script | Purpose | Key Result |
|--------|---------|------------|
| `figure_desire_spectral.jl` | Spectral consistency validation | 3-panel figure: convergence, bias O(1/√Λ), variance O(1/√m) |
| `validate_desire_spectral.jl` | Text-based validation tests | Console output verification |

**Output files** (`output/desire_validation/`):
- `figure_desire_spectral.png/pdf` - Main validation figure
- `figure_desire_spectral_notes.md` - Detailed notes for manuscript integration
- `desire_spectral_data.jld2` - Cached computation results

**Script pattern**:
```julia
function compute_figure_X_data()  # Compute and save
function plot_figure_X(data)       # Load and plot
function (@main)(args)             # Entry point, supports --plot
```

---

## Tests

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

| File | Coverage |
|------|----------|
| `test_latent_space.jl` | B^d_+ geometry, projections, volume |
| `test_intensity.jl` | Evaluation, marginal_stats, products |
| `test_sampling.jl` | PPP sampling, thinning, max intensity |
| `test_formulas.jl` | E[N], E[|E|], E[L] validation |
| `test_graph_generation.jl` | Node/edge-centric, discretization |
| `test_pde.jl` | PDE solver (OrdinaryDiffEq) |

---

## Review Checklists

### Mathematical Correctness

- [ ] Connection probability: P(edge) = g . r
- [ ] Node-centric edges: E[|E|] = E[N]^2 * E[g . r]
- [ ] Edge-centric edges: E[L] = (E[N]/2) * E[g . r]
- [ ] Ratio: E[|E|]_node / E[L]_edge = 2 * Lambda
- [ ] Product intensity: rho(g, r) = rho_G(g) * rho_R(r)
- [ ] Charges: c_G = integral rho_G, c_R = integral rho_R, Lambda = c_G * c_R
- [ ] Centroids: mu_tilde = (integral x * rho(x)) / c
- [ ] Second moment: Sigma_{jk} = E[x_j x_k] under normalized intensity
- [ ] Desire singular values: sigma_k(D̃) = sqrt(lambda_k(Sigma_G * Sigma_R))
- [ ] Spectral consistency: sigma_k(A)/N → sigma_k(D̃) as N → infinity

### Code Quality

- [ ] Types parameterized by dimension `{d}`
- [ ] Callable structs for intensity evaluation
- [ ] String concatenation (not interpolation) per project convention
- [ ] In-place functions end with `!`
- [ ] RNG parameter always last: `rng=Random.default_rng()`

### Numerical

- [ ] CFL condition for PDE: dt <= dx^2 / (2 * D * d)
- [ ] Mass conservation under reflecting BC
- [ ] Rejection sampling: proposal >= true intensity everywhere
- [ ] Monte Carlo: sufficient replications for convergence

---

## Changelog

### 2026-01-23

- **Desire Operator Spectral Validation**: New functions and validation scripts
  - Added `second_moment_matrix(rho)` - computes E[x_j x_k] under normalized intensity
  - Added `desire_operator_singular_values(rho)` - computes σ_k(D̃) = √(λ_k(Σ_G Σ_R))
  - Added `desire_stats(rho)` - convenience wrapper returning (Σ_G, Σ_R, σ, μ̃_G, μ̃_R)
  - New exports in IDPG.jl: `second_moment_matrix`, `desire_operator_singular_values`, `desire_stats`
- New validation scripts in `scripts/desire_validation/`:
  - `figure_desire_spectral.jl` - generates 3-panel validation figure
  - `validate_desire_spectral.jl` - text-based validation tests
- Validates spectral consistency theorem: σ_k(A)/N → σ_k(D̃) as Λ → ∞
  - Panel (a): Single graph convergence
  - Panel (b): Bias scaling O(1/√Λ)
  - Panel (c): Multi-graph variance reduction O(1/√m)

### 2026-01-12

- Archived PDESciML.jl (MethodOfLines-based) to `archive/deprecated_pde/`
- Type stability improvements: added `Val{d}` overloads for hot-path functions
  - `uniform_Bd_plus_sample(Val(d))` - type-stable sampling
  - `hyperspherical_to_cartesian(Val(d), r, angles)` - type-stable conversion
- Julia Performance Tips compliance: typed empty arrays, `@view` for slices
- SciML Style Guide compliance: `eachindex` iteration, type parameters
- Verified hyperspherical Jacobian against Wikipedia N-sphere formulas
- Removed refactoring artifacts: `docs/pde_refactoring_2026.md`, `refactor_verification/`

### 2025-01-09

- Reorganized: `simulations/` -> `scripts/scaling_laws/`, `heat_map_sims/` -> `scripts/heat_maps/`
- Renamed scripts to match figures: `sim1_*.jl` -> `fig1_*.jl`
- Updated paths: all scripts use `pkgdir(IDPG)`
- Moved `docs/` out of repo (kept locally)
- Moved `CODE_REVIEW_GUIDE.md` to root

### 2025-01-08

- Extracted to library: 1D utilities, ecological utils, centroid construction
- Deleted: `simulations/sim_utils.jl`
- Archived: `examples/`
