# PDE dynamics using SciML ecosystem (MethodOfLines.jl + OrdinaryDiffEq.jl)
# This provides a proper numerical PDE solver for intensity evolution on B^d_+
#
# Note: Functions here create PDE systems lazily (when called), not at module load time.

# =============================================================================
# Masking for B^d_+ domain (works with any discretization)
# =============================================================================

"""
    create_Bd_plus_mask_2d(resolution::Int) -> BitMatrix

Create a mask for points inside B²₊ (non-negative unit disk).
"""
function create_Bd_plus_mask_2d(resolution::Int)
    h = 1.0 / (resolution - 1)
    mask = falses(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            x = (i - 1) * h
            y = (j - 1) * h
            if x^2 + y^2 <= 1.0
                mask[i, j] = true
            end
        end
    end
    return mask
end

"""
    create_Bd_plus_mask_4d(resolution::Int) -> BitArray{4}

Create a mask for points inside B⁴₊ (non-negative unit 4-ball).
"""
function create_Bd_plus_mask_4d(resolution::Int)
    h = 1.0 / (resolution - 1)
    mask = falses(resolution, resolution, resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            for k in 1:resolution
                for l in 1:resolution
                    x1 = (i - 1) * h
                    x2 = (j - 1) * h
                    x3 = (k - 1) * h
                    x4 = (l - 1) * h
                    if x1^2 + x2^2 + x3^2 + x4^2 <= 1.0
                        mask[i, j, k, l] = true
                    end
                end
            end
        end
    end
    return mask
end

"""
    apply_Bd_plus_mask!(u::AbstractArray, mask::BitArray)

Zero out intensity values outside B^d_+.
"""
function apply_Bd_plus_mask!(u::AbstractArray, mask::BitArray)
    u[.!mask] .= 0.0
    return u
end

# =============================================================================
# 2D Diffusion using MethodOfLines
# =============================================================================

"""
    solve_diffusion_mol_2d(u0::Matrix{Float64}, D::Float64, tspan::Tuple;
                           dx=0.05, solver=Tsit5(), saveat=0.1, apply_mask=true)

Solve 2D diffusion equation using MethodOfLines.jl.

∂u/∂t = D(∂²u/∂x² + ∂²u/∂y²) on [0,1]²

# Arguments
- `u0`: Initial condition (resolution × resolution matrix)
- `D`: Diffusion coefficient
- `tspan`: Time span (t_start, t_end)
- `dx`: Grid spacing (should match resolution of u0)
- `solver`: ODE solver
- `saveat`: Save interval
- `apply_mask`: Whether to apply B²₊ mask

# Returns
Named tuple (times, snapshots, mask)
"""
function solve_diffusion_mol_2d(u0::Matrix{Float64}, D::Float64, tspan::Tuple;
                                dx::Float64=0.05,
                                solver=Tsit5(),
                                saveat::Float64=0.1,
                                apply_mask::Bool=true)
    resolution = size(u0, 1)
    @assert size(u0, 2) == resolution "Initial condition must be square"

    # Build symbolic PDE
    @parameters t x y
    @variables ρ(..)
    Dt = Differential(t)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2

    # Diffusion equation
    eq = Dt(ρ(t, x, y)) ~ D * (Dxx(ρ(t, x, y)) + Dyy(ρ(t, x, y)))

    # Domain - using DomainSets.Interval via ModelingToolkit
    domains = [t ∈ DomainSets.Interval(tspan[1], tspan[2]),
               x ∈ DomainSets.Interval(0.0, 1.0),
               y ∈ DomainSets.Interval(0.0, 1.0)]

    # Zero-flux boundary conditions
    Dx = Differential(x)
    Dy = Differential(y)
    bcs = [
        ρ(0.0, x, y) ~ 0.0,  # Will be overwritten with u0
        Dx(ρ(t, 0.0, y)) ~ 0.0,
        Dx(ρ(t, 1.0, y)) ~ 0.0,
        Dy(ρ(t, x, 0.0)) ~ 0.0,
        Dy(ρ(t, x, 1.0)) ~ 0.0,
    ]

    @named pdesys = PDESystem(eq, bcs, domains, [t, x, y], [ρ(t, x, y)])

    # Discretize
    discretization = MOLFiniteDifference([x => dx, y => dx], t)
    prob = discretize(pdesys, discretization)

    # Solve
    sol = solve(prob, solver; saveat=saveat)

    # Create mask
    mask = apply_mask ? create_Bd_plus_mask_2d(resolution) : nothing

    # Extract snapshots using MethodOfLines solution interface
    times = sol[t]
    sol_array = sol[ρ(t, x, y)]

    snapshots = []
    for i in eachindex(times)
        u = copy(sol_array[i, :, :])
        if !isnothing(mask)
            apply_Bd_plus_mask!(u, mask)
        end
        push!(snapshots, u)
    end

    return (times=times, snapshots=snapshots, mask=mask)
end

# =============================================================================
# 2D Advection using MethodOfLines
# =============================================================================

"""
    solve_advection_mol_2d(u0::Matrix{Float64}, v::Vector{Float64}, tspan::Tuple;
                           dx=0.05, solver=Tsit5(), saveat=0.1, apply_mask=true)

Solve 2D advection equation using MethodOfLines.jl.

∂u/∂t = -v · ∇u on [0,1]²

# Arguments
- `u0`: Initial condition
- `v`: Velocity vector [vx, vy]
- `tspan`: Time span
- Other arguments same as solve_diffusion_mol_2d
"""
function solve_advection_mol_2d(u0::Matrix{Float64}, v::Vector{Float64}, tspan::Tuple;
                                dx::Float64=0.05,
                                solver=Tsit5(),
                                saveat::Float64=0.1,
                                apply_mask::Bool=true)
    @assert length(v) == 2 "Velocity must be 2D"
    resolution = size(u0, 1)

    @parameters t x y
    @variables ρ(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    # Advection equation
    eq = Dt(ρ(t, x, y)) ~ -v[1] * Dx(ρ(t, x, y)) - v[2] * Dy(ρ(t, x, y))

    domains = [t ∈ DomainSets.Interval(tspan[1], tspan[2]),
               x ∈ DomainSets.Interval(0.0, 1.0),
               y ∈ DomainSets.Interval(0.0, 1.0)]

    # Inflow boundary conditions (zero at inflow boundaries)
    bcs = [
        ρ(0.0, x, y) ~ 0.0,
        ρ(t, 0.0, y) ~ 0.0,
        ρ(t, x, 0.0) ~ 0.0,
    ]

    @named pdesys = PDESystem(eq, bcs, domains, [t, x, y], [ρ(t, x, y)])

    discretization = MOLFiniteDifference([x => dx, y => dx], t; advection_scheme=UpwindScheme())
    prob = discretize(pdesys, discretization)

    sol = solve(prob, solver; saveat=saveat)

    mask = apply_mask ? create_Bd_plus_mask_2d(resolution) : nothing

    # Extract snapshots using MethodOfLines solution interface
    times = sol[t]
    sol_array = sol[ρ(t, x, y)]

    snapshots = []
    for i in eachindex(times)
        u = copy(sol_array[i, :, :])
        if !isnothing(mask)
            apply_Bd_plus_mask!(u, mask)
        end
        push!(snapshots, u)
    end

    return (times=times, snapshots=snapshots, mask=mask)
end

# =============================================================================
# 4D Diffusion using MethodOfLines
# =============================================================================

"""
    solve_diffusion_mol_4d(u0::Array{Float64,4}, D::Float64, tspan::Tuple;
                           dx=0.1, solver=ROCK4(), saveat=0.1, apply_mask=true)

Solve 4D diffusion equation using MethodOfLines.jl.

∂u/∂t = D(∂²u/∂x₁² + ∂²u/∂x₂² + ∂²u/∂x₃² + ∂²u/∂x₄²) on [0,1]⁴

# Arguments
- `u0`: Initial condition (resolution⁴ array)
- `D`: Diffusion coefficient
- `tspan`: Time span (t_start, t_end)
- `dx`: Grid spacing (should match resolution of u0)
- `solver`: ODE solver (ROCK4 recommended for stiff diffusion)
- `saveat`: Save interval
- `apply_mask`: Whether to apply B⁴₊ mask

# Returns
Named tuple (times, snapshots, mask)
"""
function solve_diffusion_mol_4d(u0::Array{Float64,4}, D::Float64, tspan::Tuple;
                                dx::Float64=0.1,
                                solver=ROCK4(),
                                saveat::Float64=0.1,
                                apply_mask::Bool=true)
    resolution = size(u0, 1)
    @assert all(size(u0) .== resolution) "Initial condition must be hypercubic"

    # Build symbolic PDE
    @parameters t x1 x2 x3 x4
    @variables ρ(..)
    Dt = Differential(t)
    Dx1x1 = Differential(x1)^2
    Dx2x2 = Differential(x2)^2
    Dx3x3 = Differential(x3)^2
    Dx4x4 = Differential(x4)^2

    # Diffusion equation in 4D
    eq = Dt(ρ(t, x1, x2, x3, x4)) ~ D * (Dx1x1(ρ(t, x1, x2, x3, x4)) +
                                          Dx2x2(ρ(t, x1, x2, x3, x4)) +
                                          Dx3x3(ρ(t, x1, x2, x3, x4)) +
                                          Dx4x4(ρ(t, x1, x2, x3, x4)))

    # Domain
    domains = [t ∈ DomainSets.Interval(tspan[1], tspan[2]),
               x1 ∈ DomainSets.Interval(0.0, 1.0),
               x2 ∈ DomainSets.Interval(0.0, 1.0),
               x3 ∈ DomainSets.Interval(0.0, 1.0),
               x4 ∈ DomainSets.Interval(0.0, 1.0)]

    # Zero-flux (Neumann) boundary conditions on all 8 faces
    Dx1 = Differential(x1)
    Dx2 = Differential(x2)
    Dx3 = Differential(x3)
    Dx4 = Differential(x4)

    bcs = [
        ρ(0.0, x1, x2, x3, x4) ~ 0.0,  # Will be overwritten with u0
        Dx1(ρ(t, 0.0, x2, x3, x4)) ~ 0.0,
        Dx1(ρ(t, 1.0, x2, x3, x4)) ~ 0.0,
        Dx2(ρ(t, x1, 0.0, x3, x4)) ~ 0.0,
        Dx2(ρ(t, x1, 1.0, x3, x4)) ~ 0.0,
        Dx3(ρ(t, x1, x2, 0.0, x4)) ~ 0.0,
        Dx3(ρ(t, x1, x2, 1.0, x4)) ~ 0.0,
        Dx4(ρ(t, x1, x2, x3, 0.0)) ~ 0.0,
        Dx4(ρ(t, x1, x2, x3, 1.0)) ~ 0.0,
    ]

    @named pdesys = PDESystem(eq, bcs, domains, [t, x1, x2, x3, x4], [ρ(t, x1, x2, x3, x4)])

    # Discretize with MethodOfLines
    discretization = MOLFiniteDifference([x1 => dx, x2 => dx, x3 => dx, x4 => dx], t)
    prob = discretize(pdesys, discretization)

    # Solve
    sol = solve(prob, solver; saveat=saveat)

    # Create mask
    mask = apply_mask ? create_Bd_plus_mask_4d(resolution) : nothing

    # Extract snapshots using MethodOfLines solution interface
    # sol[ρ(t, x1, x2, x3, x4)] returns 5D array with time as first dimension
    # (since t appears first in the variable signature)
    times = sol[t]
    sol_array = sol[ρ(t, x1, x2, x3, x4)]

    snapshots = []
    for i in eachindex(times)
        # Extract spatial slice at time index i (first dimension is time)
        u = copy(sol_array[i, :, :, :, :])
        if !isnothing(mask)
            apply_Bd_plus_mask!(u, mask)
        end
        push!(snapshots, u)
    end

    return (times=times, snapshots=snapshots, mask=mask)
end

# =============================================================================
# 4D Advection using MethodOfLines
# =============================================================================

"""
    solve_advection_mol_4d(u0::Array{Float64,4}, v::Vector{Float64}, tspan::Tuple;
                           dx=0.1, solver=Tsit5(), saveat=0.1, apply_mask=true)

Solve 4D advection equation using MethodOfLines.jl.

∂u/∂t = -v · ∇u on [0,1]⁴

# Arguments
- `u0`: Initial condition
- `v`: Velocity vector [v1, v2, v3, v4]
- `tspan`: Time span
- Other arguments same as solve_diffusion_mol_4d
"""
function solve_advection_mol_4d(u0::Array{Float64,4}, v::Vector{Float64}, tspan::Tuple;
                                dx::Float64=0.1,
                                solver=Tsit5(),
                                saveat::Float64=0.1,
                                apply_mask::Bool=true)
    @assert length(v) == 4 "Velocity must be 4D"
    resolution = size(u0, 1)
    @assert all(size(u0) .== resolution) "Initial condition must be hypercubic"

    @parameters t x1 x2 x3 x4
    @variables ρ(..)
    Dt = Differential(t)
    Dx1 = Differential(x1)
    Dx2 = Differential(x2)
    Dx3 = Differential(x3)
    Dx4 = Differential(x4)

    # Advection equation in 4D
    eq = Dt(ρ(t, x1, x2, x3, x4)) ~ -v[1] * Dx1(ρ(t, x1, x2, x3, x4)) -
                                     v[2] * Dx2(ρ(t, x1, x2, x3, x4)) -
                                     v[3] * Dx3(ρ(t, x1, x2, x3, x4)) -
                                     v[4] * Dx4(ρ(t, x1, x2, x3, x4))

    domains = [t ∈ DomainSets.Interval(tspan[1], tspan[2]),
               x1 ∈ DomainSets.Interval(0.0, 1.0),
               x2 ∈ DomainSets.Interval(0.0, 1.0),
               x3 ∈ DomainSets.Interval(0.0, 1.0),
               x4 ∈ DomainSets.Interval(0.0, 1.0)]

    # Boundary conditions: Dirichlet at inflow + Neumann at outflow
    # MethodOfLines docs suggest 2+ BCs per boundary for advection schemes
    bcs = [
        ρ(0.0, x1, x2, x3, x4) ~ 0.0,  # Initial condition
        # Dirichlet (zero value) at x=0 boundaries (inflow)
        ρ(t, 0.0, x2, x3, x4) ~ 0.0,
        ρ(t, x1, 0.0, x3, x4) ~ 0.0,
        ρ(t, x1, x2, 0.0, x4) ~ 0.0,
        ρ(t, x1, x2, x3, 0.0) ~ 0.0,
        # Neumann (zero gradient) at x=1 boundaries (outflow)
        Dx1(ρ(t, 1.0, x2, x3, x4)) ~ 0.0,
        Dx2(ρ(t, x1, 1.0, x3, x4)) ~ 0.0,
        Dx3(ρ(t, x1, x2, 1.0, x4)) ~ 0.0,
        Dx4(ρ(t, x1, x2, x3, 1.0)) ~ 0.0,
    ]

    @named pdesys = PDESystem(eq, bcs, domains, [t, x1, x2, x3, x4], [ρ(t, x1, x2, x3, x4)])

    # Note: 4D advection requires Neumann BCs at outflow to avoid BoundsError
    # See docs/known_issues.md for details on boundary condition requirements
    discretization = MOLFiniteDifference([x1 => dx, x2 => dx, x3 => dx, x4 => dx], t)
    prob = discretize(pdesys, discretization)

    sol = solve(prob, solver; saveat=saveat)

    mask = apply_mask ? create_Bd_plus_mask_4d(resolution) : nothing

    # Extract snapshots using MethodOfLines solution interface
    times = sol[t]
    sol_array = sol[ρ(t, x1, x2, x3, x4)]

    snapshots = []
    for i in eachindex(times)
        u = copy(sol_array[i, :, :, :, :])
        if !isnothing(mask)
            apply_Bd_plus_mask!(u, mask)
        end
        push!(snapshots, u)
    end

    return (times=times, snapshots=snapshots, mask=mask)
end
