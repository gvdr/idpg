# PDE dynamics on intensity functions over B^d_+
# Refactored January 2026 to use OrdinaryDiffEq.jl for time integration

"""
    BdPlusGrid{d}

Finite difference grid on B^d_+ (non-negative unit ball), implemented by embedding
in a regular grid and masking points outside B^d_+.

# Fields
- `d::Int`: Dimension of the space
- `resolution::Int`: Number of grid points along each axis
- `points::Vector{LatentPoint{d}}`: Grid points inside B^d_+
- `inside_mask::BitVector`: Which linear indices are inside B^d_+
- `grid_to_Bd_plus::Dict{CartesianIndex, Int}`: Map from grid coords to point index
- `Bd_plus_to_grid::Vector{CartesianIndex}`: Map from point index to grid coords
- `h::Float64`: Grid spacing
"""
struct BdPlusGrid{d}
    resolution::Int
    points::Vector{LatentPoint{d}}
    inside_mask::BitVector
    grid_to_Bd_plus::Dict{CartesianIndex{d}, Int}
    Bd_plus_to_grid::Vector{CartesianIndex{d}}
    h::Float64
end

"""
    create_Bd_plus_grid(d::Int, resolution::Int) -> BdPlusGrid{d}

Create a regular grid on B^d_+ by embedding in [0,1]^d
and retaining only points with ||x|| ≤ 1.

# Arguments
- `d`: Dimension of the space
- `resolution`: Number of grid points along each axis
"""
function create_Bd_plus_grid(d::Int, resolution::Int)
    h = 1.0 / (resolution - 1)
    points = LatentPoint{d}[]
    grid_to_Bd_plus = Dict{CartesianIndex{d}, Int}()
    Bd_plus_to_grid = CartesianIndex{d}[]

    # Total number of grid points in bounding box
    n_total = resolution^d
    inside_mask = falses(n_total)

    # Iterate over grid
    for idx in CartesianIndices(ntuple(_ -> resolution, d))
        # Convert to coordinates in [0,1]^d
        coords = [(idx[i] - 1) * h for i in 1:d]

        # Check if in B^d_+ (non-negative and ||x|| ≤ 1)
        if all(c -> c >= 0, coords) && norm(coords) <= 1.0
            point = SVector{d, Float64}(coords)

            push!(points, point)
            linear_idx = LinearIndices(ntuple(_ -> resolution, d))[idx]
            inside_mask[linear_idx] = true
            grid_to_Bd_plus[idx] = length(points)
            push!(Bd_plus_to_grid, idx)
        end
    end

    return BdPlusGrid{d}(resolution, points, inside_mask, grid_to_Bd_plus, Bd_plus_to_grid, h)
end

"""
    get_neighbors(grid::BdPlusGrid{d}, point_idx::Int) -> Vector{Int}

Get indices of neighboring points on the grid.
Uses 2d neighbors in the grid (±1 along each axis).
"""
function get_neighbors(grid::BdPlusGrid{d}, point_idx::Int) where d
    idx = grid.Bd_plus_to_grid[point_idx]
    neighbors = Int[]

    for dim in 1:d
        for δ in [-1, 1]
            # Create neighbor index
            neighbor_tuple = ntuple(i -> i == dim ? idx[i] + δ : idx[i], d)

            # Check bounds
            if all(1 <= neighbor_tuple[i] <= grid.resolution for i in 1:d)
                neighbor_idx = CartesianIndex(neighbor_tuple)
                if haskey(grid.grid_to_Bd_plus, neighbor_idx)
                    push!(neighbors, grid.grid_to_Bd_plus[neighbor_idx])
                end
            end
        end
    end

    return neighbors
end

"""
    laplacian_stencil(grid::BdPlusGrid{d}, ρ::AbstractVector, point_idx::Int) -> Float64

Compute discrete Laplacian ∇²ρ at a point using finite differences.
"""
function laplacian_stencil(grid::BdPlusGrid{d}, ρ::AbstractVector, point_idx::Int) where d
    neighbors = get_neighbors(grid, point_idx)
    n_neighbors = length(neighbors)

    if n_neighbors == 0
        return 0.0
    end

    h² = grid.h^2
    ρ_center = ρ[point_idx]

    # Standard Laplacian stencil: ∇²ρ ≈ Σ(ρ_neighbor - ρ_center) / h²
    ∇²ρ = 0.0
    for neighbor_idx in neighbors
        ∇²ρ += ρ[neighbor_idx] - ρ_center
    end

    return ∇²ρ / h²
end

"""
    gradient_component(grid::BdPlusGrid{d}, ρ::AbstractVector, point_idx::Int, dim::Int) -> Float64

Compute gradient component ∂ρ/∂xᵢ in direction `dim` using central differences.
"""
function gradient_component(grid::BdPlusGrid{d}, ρ::AbstractVector,
                           point_idx::Int, dim::Int) where d
    idx = grid.Bd_plus_to_grid[point_idx]

    # Get neighbors in this dimension
    idx₊ = ntuple(i -> i == dim ? idx[i] + 1 : idx[i], d)
    idx₋ = ntuple(i -> i == dim ? idx[i] - 1 : idx[i], d)

    # Check if neighbors exist
    has_plus = haskey(grid.grid_to_Bd_plus, CartesianIndex(idx₊))
    has_minus = haskey(grid.grid_to_Bd_plus, CartesianIndex(idx₋))

    if has_plus && has_minus
        # Central difference: (ρ₊ - ρ₋) / 2h
        ρ₊ = ρ[grid.grid_to_Bd_plus[CartesianIndex(idx₊)]]
        ρ₋ = ρ[grid.grid_to_Bd_plus[CartesianIndex(idx₋)]]
        return (ρ₊ - ρ₋) / (2 * grid.h)
    elseif has_plus
        # Forward difference: (ρ₊ - ρ_center) / h
        ρ₊ = ρ[grid.grid_to_Bd_plus[CartesianIndex(idx₊)]]
        return (ρ₊ - ρ[point_idx]) / grid.h
    elseif has_minus
        # Backward difference: (ρ_center - ρ₋) / h
        ρ₋ = ρ[grid.grid_to_Bd_plus[CartesianIndex(idx₋)]]
        return (ρ[point_idx] - ρ₋) / grid.h
    else
        return 0.0
    end
end

# =============================================================================
# ODE Right-Hand-Side Functions for OrdinaryDiffEq
# =============================================================================

"""
Create the RHS function for the diffusion equation: ∂ρ/∂t = D∇²ρ
"""
function make_diffusion_rhs(grid::BdPlusGrid{d}, D::Float64) where d
    function diffusion_rhs!(dρ, ρ, p, t)
        for i in eachindex(ρ)
            dρ[i] = D * laplacian_stencil(grid, ρ, i)
        end
    end
    return diffusion_rhs!
end

"""
Create the RHS function for advection: ∂ρ/∂t = -v⃗ · ∇ρ
"""
function make_advection_rhs(grid::BdPlusGrid{d}, v⃗::AbstractVector) where d
    function advection_rhs!(dρ, ρ, p, t)
        for i in eachindex(ρ)
            # -v⃗ · ∇ρ = -Σ vₖ ∂ρ/∂xₖ
            advection = 0.0
            for dim in 1:d
                ∂ρ_∂xₖ = gradient_component(grid, ρ, i, dim)
                advection -= v⃗[dim] * ∂ρ_∂xₖ
            end
            dρ[i] = advection
        end
    end
    return advection_rhs!
end

"""
Create the RHS function for advection with space-dependent velocity field: ∂ρ/∂t = -v⃗(x) · ∇ρ
"""
function make_advection_field_rhs(grid::BdPlusGrid{d}, v⃗_field::Function) where d
    function advection_field_rhs!(dρ, ρ, p, t)
        for i in eachindex(ρ)
            x = grid.points[i]
            v⃗ = v⃗_field(Vector(x))
            advection = 0.0
            for dim in 1:d
                ∂ρ_∂xₖ = gradient_component(grid, ρ, i, dim)
                advection -= v⃗[dim] * ∂ρ_∂xₖ
            end
            dρ[i] = advection
        end
    end
    return advection_field_rhs!
end

"""
Create the RHS function for reaction-diffusion: ∂ρ/∂t = D∇²ρ + f(ρ)
"""
function make_reaction_diffusion_rhs(grid::BdPlusGrid{d}, D::Float64, f::Function) where d
    function reaction_diffusion_rhs!(dρ, ρ, p, t)
        for i in eachindex(ρ)
            ∇²ρ = laplacian_stencil(grid, ρ, i)
            reaction = f(ρ[i])
            dρ[i] = D * ∇²ρ + reaction
        end
    end
    return reaction_diffusion_rhs!
end

# =============================================================================
# Main Evolution Functions Using OrdinaryDiffEq
# =============================================================================

"""
    evolve_diffusion(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                     D::Float64, tspan::Tuple;
                     solver=Tsit5(), saveat=nothing, kwargs...) -> ODESolution

Evolve intensity via diffusion equation: ∂ρ/∂t = D∇²ρ

Uses OrdinaryDiffEq.jl for adaptive time-stepping and error control.

# Arguments
- `ρ₀`: Initial intensity values at grid points
- `grid`: B^d_+ grid
- `D`: Diffusion coefficient
- `tspan`: Time span as (t₀, t_final)
- `solver`: ODE solver (default: Tsit5 for non-stiff; use ROCK4 or Rodas5P for stiff)
- `saveat`: Times to save solution (default: automatic)
- `kwargs...`: Additional arguments passed to solve()

# Returns
ODE solution object. Access solution at time t via sol(t).

# Example
```julia
grid = create_Bd_plus_grid(2, 20)
ρ₀ = [exp(-10 * norm(p - [0.5, 0.5])^2) for p in grid.points]
sol = evolve_diffusion(ρ₀, grid, 0.01, (0.0, 1.0))
ρ_final = sol.u[end]
```
"""
function evolve_diffusion(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                          D::Float64, tspan::Tuple;
                          solver=Tsit5(),
                          saveat=nothing,
                          kwargs...) where d
    rhs! = make_diffusion_rhs(grid, D)
    prob = ODEProblem(rhs!, copy(ρ₀), tspan)

    solve_kwargs = Dict{Symbol, Any}(kwargs)
    if !isnothing(saveat)
        solve_kwargs[:saveat] = saveat
    end

    return solve(prob, solver; solve_kwargs...)
end

"""
    evolve_diffusion!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                      D::Float64, dt::Float64, n_steps::Int;
                      boundary::Symbol=:absorbing) -> Vector{Float64}

Legacy in-place interface for backwards compatibility.
Internally uses OrdinaryDiffEq but returns final state in ρ.
"""
function evolve_diffusion!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                           D::Float64, dt::Float64, n_steps::Int;
                           boundary::Symbol=:absorbing) where d
    t_final = dt * n_steps
    sol = evolve_diffusion(ρ, grid, D, (0.0, t_final))
    copyto!(ρ, max.(sol.u[end], 0.0))
    return ρ
end

"""
    evolve_advection(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                     v⃗::AbstractVector, tspan::Tuple;
                     solver=Tsit5(), saveat=nothing, kwargs...) -> ODESolution

Evolve intensity via advection equation: ∂ρ/∂t = -v⃗ · ∇ρ

# Arguments
- `ρ₀`: Initial intensity values at grid points
- `grid`: B^d_+ grid
- `v⃗`: Velocity vector (constant in space)
- `tspan`: Time span as (t₀, t_final)
- `solver`: ODE solver
- `saveat`: Times to save solution
"""
function evolve_advection(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                          v⃗::AbstractVector, tspan::Tuple;
                          solver=Tsit5(),
                          saveat=nothing,
                          kwargs...) where d
    rhs! = make_advection_rhs(grid, v⃗)
    prob = ODEProblem(rhs!, copy(ρ₀), tspan)

    solve_kwargs = Dict{Symbol, Any}(kwargs)
    if !isnothing(saveat)
        solve_kwargs[:saveat] = saveat
    end

    return solve(prob, solver; solve_kwargs...)
end

"""
    evolve_advection!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                      v⃗::AbstractVector, dt::Float64, n_steps::Int;
                      boundary::Symbol=:absorbing) -> Vector{Float64}

Legacy in-place interface for backwards compatibility.
"""
function evolve_advection!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                           v⃗::AbstractVector, dt::Float64, n_steps::Int;
                           boundary::Symbol=:absorbing) where d
    t_final = dt * n_steps
    sol = evolve_advection(ρ, grid, v⃗, (0.0, t_final))
    copyto!(ρ, max.(sol.u[end], 0.0))
    return ρ
end

"""
    evolve_advection_field(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                           v⃗_field::Function, tspan::Tuple;
                           solver=Tsit5(), saveat=nothing, kwargs...) -> ODESolution

Evolve intensity via advection with space-dependent velocity: ∂ρ/∂t = -v⃗(x) · ∇ρ

# Arguments
- `ρ₀`: Initial intensity values at grid points
- `grid`: B^d_+ grid
- `v⃗_field`: Function v⃗_field(x::AbstractVector) → velocity vector at point x
- `tspan`: Time span as (t₀, t_final)

# Example
```julia
# Velocity field pointing toward a target
target = [0.5, 0.5]
v⃗_toward(x) = (target .- x) ./ max(norm(target .- x), 0.1)
sol = evolve_advection_field(ρ₀, grid, v⃗_toward, (0.0, 1.0))
```
"""
function evolve_advection_field(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                                v⃗_field::Function, tspan::Tuple;
                                solver=Tsit5(),
                                saveat=nothing,
                                kwargs...) where d
    rhs! = make_advection_field_rhs(grid, v⃗_field)
    prob = ODEProblem(rhs!, copy(ρ₀), tspan)

    solve_kwargs = Dict{Symbol, Any}(kwargs)
    if !isnothing(saveat)
        solve_kwargs[:saveat] = saveat
    end

    return solve(prob, solver; solve_kwargs...)
end

"""
    evolve_advection_field!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                            v⃗_field::Function, dt::Float64, n_steps::Int;
                            boundary::Symbol=:absorbing) -> Vector{Float64}

Legacy in-place interface for backwards compatibility.
"""
function evolve_advection_field!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                                  v⃗_field::Function, dt::Float64, n_steps::Int;
                                  boundary::Symbol=:absorbing) where d
    t_final = dt * n_steps
    sol = evolve_advection_field(ρ, grid, v⃗_field, (0.0, t_final))
    copyto!(ρ, max.(sol.u[end], 0.0))
    return ρ
end

"""
    evolve_reaction_diffusion(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                              D::Float64, f::Function, tspan::Tuple;
                              solver=Rodas5P(), saveat=nothing, kwargs...) -> ODESolution

Evolve intensity via reaction-diffusion equation: ∂ρ/∂t = D∇²ρ + f(ρ)

# Arguments
- `ρ₀`: Initial intensity values at grid points
- `grid`: B^d_+ grid
- `D`: Diffusion coefficient
- `f`: Reaction function f(ρ) → rate of change
- `tspan`: Time span as (t₀, t_final)
- `solver`: ODE solver (default: Rodas5P for stiff reaction-diffusion)

# Example
```julia
# Logistic growth with diffusion
f_logistic(ρ) = 0.1 * ρ * (1 - ρ / 10)
sol = evolve_reaction_diffusion(ρ₀, grid, 0.01, f_logistic, (0.0, 10.0))
```
"""
function evolve_reaction_diffusion(ρ₀::Vector{Float64}, grid::BdPlusGrid{d},
                                   D::Float64, f::Function, tspan::Tuple;
                                   solver=Rodas5P(),
                                   saveat=nothing,
                                   kwargs...) where d
    rhs! = make_reaction_diffusion_rhs(grid, D, f)
    prob = ODEProblem(rhs!, copy(ρ₀), tspan)

    solve_kwargs = Dict{Symbol, Any}(kwargs)
    if !isnothing(saveat)
        solve_kwargs[:saveat] = saveat
    end

    return solve(prob, solver; solve_kwargs...)
end

"""
    evolve_reaction_diffusion!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                               D::Float64, f::Function, dt::Float64, n_steps::Int;
                               boundary::Symbol=:absorbing) -> Vector{Float64}

Legacy in-place interface for backwards compatibility.
"""
function evolve_reaction_diffusion!(ρ::Vector{Float64}, grid::BdPlusGrid{d},
                                    D::Float64, f::Function, dt::Float64, n_steps::Int;
                                    boundary::Symbol=:absorbing) where d
    t_final = dt * n_steps
    sol = evolve_reaction_diffusion(ρ, grid, D, f, (0.0, t_final))
    copyto!(ρ, max.(sol.u[end], 0.0))
    return ρ
end

# =============================================================================
# Tracking and Analysis
# =============================================================================

"""
    evolve_and_track(ρ_initial::Function, grid::BdPlusGrid{d};
                     pde_type::Symbol=:diffusion,
                     D::Float64=0.01, v⃗=nothing, f=nothing,
                     t_final::Float64=1.0,
                     sample_times::AbstractVector=0.0:0.1:1.0,
                     solver=nothing,
                     rng=Random.default_rng())

Evolve intensity and track statistics over time using OrdinaryDiffEq.

# Arguments
- `ρ_initial`: Function to evaluate initial intensity at grid points
- `grid`: B^d_+ grid
- `pde_type`: Type of PDE (:diffusion, :advection, :reaction_diffusion)
- `D`: Diffusion coefficient
- `v⃗`: Velocity vector for advection
- `f`: Reaction function for reaction-diffusion
- `t_final`: Final time
- `sample_times`: Times at which to record statistics
- `solver`: ODE solver (auto-selected if nothing)
- `rng`: Random number generator

# Returns
Named tuple with:
- `times`: Time points
- `ρ_history`: Intensity values at each time
- `total_intensity`: Total intensity over time
- `mean_position`: Mean position over time
- `solution`: The full ODE solution object
"""
function evolve_and_track(ρ_initial::Function, grid::BdPlusGrid{d};
                          pde_type::Symbol=:diffusion,
                          D::Float64=0.01, v⃗=nothing, f=nothing,
                          t_final::Float64=1.0,
                          sample_times::AbstractVector=0.0:0.1:1.0,
                          solver=nothing,
                          rng::AbstractRNG=Random.default_rng()) where d

    # Initialize
    ρ₀ = [ρ_initial(p) for p in grid.points]
    tspan = (0.0, t_final)

    # Select solver and solve
    if pde_type == :diffusion
        default_solver = isnothing(solver) ? Tsit5() : solver
        sol = evolve_diffusion(ρ₀, grid, D, tspan; solver=default_solver, saveat=sample_times)
    elseif pde_type == :advection
        @assert !isnothing(v⃗) "Velocity v⃗ required for advection"
        default_solver = isnothing(solver) ? Tsit5() : solver
        sol = evolve_advection(ρ₀, grid, v⃗, tspan; solver=default_solver, saveat=sample_times)
    elseif pde_type == :reaction_diffusion
        @assert !isnothing(f) "Reaction function f required for reaction-diffusion"
        default_solver = isnothing(solver) ? Rodas5P() : solver
        sol = evolve_reaction_diffusion(ρ₀, grid, D, f, tspan; solver=default_solver, saveat=sample_times)
    else
        error("Unknown PDE type: " * string(pde_type))
    end

    # Extract statistics at each saved time
    times = sol.t
    ρ_history = [max.(sol.u[i], 0.0) for i in eachindex(times)]
    total_intensities = [sum(ρ) * grid.h^d for ρ in ρ_history]
    mean_positions = [[compute_mean_position(ρ, grid)] for ρ in ρ_history]

    return (
        times = times,
        ρ_history = ρ_history,
        total_intensity = total_intensities,
        mean_position = mean_positions,
        grid = grid,
        solution = sol
    )
end

"""
Compute intensity-weighted mean position μ = ∫x·ρ(x)dx / ∫ρ(x)dx from grid values.
"""
function compute_mean_position(ρ::AbstractVector, grid::BdPlusGrid{d}) where d
    total = sum(ρ)
    if total < 1e-10
        return SVector{d, Float64}(fill(0.5/sqrt(d), d))  # Center of B^d_+
    end

    μ = sum(ρᵢ * xᵢ for (ρᵢ, xᵢ) in zip(ρ, grid.points))
    return SVector{d, Float64}(μ ./ total)
end
