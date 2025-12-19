# PDE dynamics on intensity functions over B^d_+

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
and retaining only points with ||x|| <= 1.

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

        # Check if in B^d_+ (non-negative and ||x|| <= 1)
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
        for delta in [-1, 1]
            # Create neighbor index
            neighbor_tuple = ntuple(i -> i == dim ? idx[i] + delta : idx[i], d)

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
    laplacian_stencil(grid::BdPlusGrid{d}, ρ_values::Vector{Float64}, point_idx::Int) -> Float64

Compute discrete Laplacian at a point using finite differences.
"""
function laplacian_stencil(grid::BdPlusGrid{d}, ρ_values::Vector{Float64}, point_idx::Int) where d
    neighbors = get_neighbors(grid, point_idx)
    n_neighbors = length(neighbors)

    if n_neighbors == 0
        return 0.0
    end

    h² = grid.h^2
    center_val = ρ_values[point_idx]

    # Standard Laplacian stencil
    laplacian = 0.0
    for neighbor_idx in neighbors
        laplacian += ρ_values[neighbor_idx] - center_val
    end

    return laplacian / h²
end

"""
    gradient_component(grid::BdPlusGrid{d}, ρ_values::Vector{Float64}, point_idx::Int, dim::Int) -> Float64

Compute gradient component in direction `dim` using central differences.
"""
function gradient_component(grid::BdPlusGrid{d}, ρ_values::Vector{Float64},
                           point_idx::Int, dim::Int) where d
    idx = grid.Bd_plus_to_grid[point_idx]

    # Get neighbors in this dimension
    idx_plus = ntuple(i -> i == dim ? idx[i] + 1 : idx[i], d)
    idx_minus = ntuple(i -> i == dim ? idx[i] - 1 : idx[i], d)

    # Check if neighbors exist
    has_plus = haskey(grid.grid_to_Bd_plus, CartesianIndex(idx_plus))
    has_minus = haskey(grid.grid_to_Bd_plus, CartesianIndex(idx_minus))

    if has_plus && has_minus
        # Central difference
        plus_idx = grid.grid_to_Bd_plus[CartesianIndex(idx_plus)]
        minus_idx = grid.grid_to_Bd_plus[CartesianIndex(idx_minus)]
        return (ρ_values[plus_idx] - ρ_values[minus_idx]) / (2 * grid.h)
    elseif has_plus
        # Forward difference
        plus_idx = grid.grid_to_Bd_plus[CartesianIndex(idx_plus)]
        return (ρ_values[plus_idx] - ρ_values[point_idx]) / grid.h
    elseif has_minus
        # Backward difference
        minus_idx = grid.grid_to_Bd_plus[CartesianIndex(idx_minus)]
        return (ρ_values[point_idx] - ρ_values[minus_idx]) / grid.h
    else
        return 0.0
    end
end

"""
    evolve_diffusion!(ρ_values::Vector{Float64}, grid::BdPlusGrid{d},
                      D::Float64, dt::Float64, n_steps::Int;
                      boundary::Symbol=:absorbing) -> Vector{Float64}

Evolve intensity via diffusion equation: ∂ρ/∂t = D∇²ρ

Uses explicit Euler with Laplacian stencil.

# Arguments
- `ρ_values`: Current intensity values at grid points (modified in place)
- `grid`: B^d_+ grid
- `D`: Diffusion coefficient
- `dt`: Time step
- `n_steps`: Number of time steps
- `boundary`: Boundary condition (:absorbing or :reflecting)

# Returns
The evolved ρ_values vector.
"""
function evolve_diffusion!(ρ_values::Vector{Float64}, grid::BdPlusGrid{d},
                           D::Float64, dt::Float64, n_steps::Int;
                           boundary::Symbol=:absorbing) where d
    n_points = length(grid.points)
    ρ_new = similar(ρ_values)

    for _ in 1:n_steps
        for i in 1:n_points
            lap = laplacian_stencil(grid, ρ_values, i)
            ρ_new[i] = ρ_values[i] + dt * D * lap
        end

        # Enforce non-negativity
        @. ρ_new = max(ρ_new, 0.0)

        # Copy back
        copyto!(ρ_values, ρ_new)
    end

    return ρ_values
end

"""
    evolve_advection!(ρ_values::Vector{Float64}, grid::BdPlusGrid{d},
                      v::AbstractVector, dt::Float64, n_steps::Int;
                      boundary::Symbol=:absorbing) -> Vector{Float64}

Evolve intensity via advection equation: ∂ρ/∂t = -v · ∇ρ

Uses upwind scheme for stability.

# Arguments
- `ρ_values`: Current intensity values at grid points (modified in place)
- `grid`: B^d_+ grid
- `v`: Velocity vector (constant in space)
- `dt`: Time step
- `n_steps`: Number of time steps
- `boundary`: Boundary condition
"""
function evolve_advection!(ρ_values::Vector{Float64}, grid::BdPlusGrid{d},
                           v::AbstractVector, dt::Float64, n_steps::Int;
                           boundary::Symbol=:absorbing) where d
    n_points = length(grid.points)
    ρ_new = similar(ρ_values)

    for _ in 1:n_steps
        for i in 1:n_points
            # Compute -v · ∇ρ using upwind differencing
            advection = 0.0
            for dim in 1:d
                grad = gradient_component(grid, ρ_values, i, dim)
                advection -= v[dim] * grad
            end
            ρ_new[i] = ρ_values[i] + dt * advection
        end

        # Enforce non-negativity
        @. ρ_new = max(ρ_new, 0.0)

        copyto!(ρ_values, ρ_new)
    end

    return ρ_values
end

"""
    evolve_reaction_diffusion!(ρ_values::Vector{Float64}, grid::BdPlusGrid{d},
                               D::Float64, f::Function, dt::Float64, n_steps::Int;
                               boundary::Symbol=:absorbing) -> Vector{Float64}

Evolve intensity via reaction-diffusion equation: ∂ρ/∂t = D∇²ρ + f(ρ)

# Arguments
- `ρ_values`: Current intensity values at grid points (modified in place)
- `grid`: B^d_+ grid
- `D`: Diffusion coefficient
- `f`: Reaction function f(ρ) -> rate of change
- `dt`: Time step
- `n_steps`: Number of time steps
- `boundary`: Boundary condition
"""
function evolve_reaction_diffusion!(ρ_values::Vector{Float64}, grid::BdPlusGrid{d},
                                    D::Float64, f::Function, dt::Float64, n_steps::Int;
                                    boundary::Symbol=:absorbing) where d
    n_points = length(grid.points)
    ρ_new = similar(ρ_values)

    for _ in 1:n_steps
        for i in 1:n_points
            lap = laplacian_stencil(grid, ρ_values, i)
            reaction = f(ρ_values[i])
            ρ_new[i] = ρ_values[i] + dt * (D * lap + reaction)
        end

        # Enforce non-negativity
        @. ρ_new = max(ρ_new, 0.0)

        copyto!(ρ_values, ρ_new)
    end

    return ρ_values
end

"""
    evolve_and_track(ρ_initial::Function, grid::BdPlusGrid{d};
                     pde_type::Symbol=:diffusion,
                     D::Float64=0.01, v=nothing, f=nothing,
                     dt::Float64=0.001, t_final::Float64=1.0,
                     sample_interval::Float64=0.1,
                     n_graph_samples::Int=100,
                     rng=Random.default_rng())

Evolve intensity and track statistics over time.

# Arguments
- `ρ_initial`: Function to evaluate initial intensity at grid points
- `grid`: B^d_+ grid
- `pde_type`: Type of PDE (:diffusion, :advection, :reaction_diffusion)
- `D`: Diffusion coefficient
- `v`: Velocity vector for advection
- `f`: Reaction function for reaction-diffusion
- `dt`: Time step
- `t_final`: Final time
- `sample_interval`: How often to record statistics
- `n_graph_samples`: Number of graph samples for empirical statistics
- `rng`: Random number generator

# Returns
Named tuple with:
- `times`: Time points
- `ρ_history`: Intensity values at each time
- `total_intensity`: Total intensity over time
- `mean_position`: Mean position over time
"""
function evolve_and_track(ρ_initial::Function, grid::BdPlusGrid{d};
                          pde_type::Symbol=:diffusion,
                          D::Float64=0.01, v=nothing, f=nothing,
                          dt::Float64=0.001, t_final::Float64=1.0,
                          sample_interval::Float64=0.1,
                          n_graph_samples::Int=100,
                          rng::AbstractRNG=Random.default_rng()) where d

    # Initialize
    ρ_values = [ρ_initial(p) for p in grid.points]

    # Storage for tracking
    n_samples = Int(ceil(t_final / sample_interval)) + 1
    times = Float64[]
    ρ_history = Vector{Float64}[]
    total_intensities = Float64[]
    mean_positions = Vector{SVector{d, Float64}}[]

    # Record initial state
    push!(times, 0.0)
    push!(ρ_history, copy(ρ_values))
    push!(total_intensities, sum(ρ_values) * grid.h^d)  # Approximate integral

    # Compute initial mean position
    μ = compute_mean_position(ρ_values, grid)
    push!(mean_positions, [μ])

    # Evolve
    steps_per_sample = Int(round(sample_interval / dt))
    t = 0.0

    while t < t_final
        # Evolve for one sample interval
        if pde_type == :diffusion
            evolve_diffusion!(ρ_values, grid, D, dt, steps_per_sample)
        elseif pde_type == :advection
            @assert !isnothing(v) "Velocity v required for advection"
            evolve_advection!(ρ_values, grid, v, dt, steps_per_sample)
        elseif pde_type == :reaction_diffusion
            @assert !isnothing(f) "Reaction function f required for reaction-diffusion"
            evolve_reaction_diffusion!(ρ_values, grid, D, f, dt, steps_per_sample)
        else
            error("Unknown PDE type: " * string(pde_type))
        end

        t += sample_interval

        # Record state
        push!(times, t)
        push!(ρ_history, copy(ρ_values))
        push!(total_intensities, sum(ρ_values) * grid.h^d)

        μ = compute_mean_position(ρ_values, grid)
        push!(mean_positions, [μ])
    end

    return (
        times = times,
        ρ_history = ρ_history,
        total_intensity = total_intensities,
        mean_position = mean_positions,
        grid = grid
    )
end

"""
Compute intensity-weighted mean position from grid values.
"""
function compute_mean_position(ρ_values::Vector{Float64}, grid::BdPlusGrid{d}) where d
    total = sum(ρ_values)
    if total < 1e-10
        return SVector{d, Float64}(fill(0.5/sqrt(d), d))  # Center of B^d_+
    end

    μ = zeros(d)
    for (i, p) in enumerate(grid.points)
        μ .+= ρ_values[i] .* p
    end
    return SVector{d, Float64}(μ ./ total)
end
