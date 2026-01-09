# Figure B: PDE Evolution (Reflecting Boundary)
#
# Dimension: d = 1
# Purpose: Show how different dynamics reshape the bound heat
#
# Layout: 4 rows × 2 columns
#   Row 1: Static (t=0, t=T)
#   Row 2: Diffusion (t=0, t=T)
#   Row 3: Advection (t=0, t=T)
#   Row 4: Pursuit-Evasion (t=0, t=T)
#
# Boundary: Reflecting (Neumann) — mass conserved
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_B_pde_reflecting.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_B_pde_reflecting.jl --plot   # plot from saved data

using IDPG
using CairoMakie
using LaTeXStrings
using Distributions
using Statistics
using LinearAlgebra
using JLD2

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "heat_maps")
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_B_data.jld2")
const GRID_RES = 100

# Intensity parameters (d=1)
const MU_G = 0.3
const MU_R = 0.7
const SIGMA = 0.12

# PDE parameters
const T_FINAL = 2.0         # Longer evolution time
const DT = 0.0005           # Smaller time step for stability
const D_COEFF = 0.05        # Diffusion coefficient (keep stable)
const V_G = 0.5             # Advection velocity for G (increased)
const V_R = -0.4            # Advection velocity for R (increased)
const ALPHA = 1.0           # Prey evasion speed (increased)
const BETA = 0.8            # Predator pursuit speed (increased)
const GAMMA = 0.1           # Centering strength (reduced for more movement)
const CENTER = 0.5          # Center position for elastic centering

# ============================================================================
# Intensity Functions
# ============================================================================

# Note: Uses truncated_gaussian_sigma from IDPG library
# Alias for backward compatibility
const truncated_gaussian = truncated_gaussian_sigma

function compute_initial_intensities(grid::Vector{Float64})
    ρ_G, ρ_R = compute_1d_marginals(grid, MU_G, MU_R, SIGMA)
    return ρ_G, ρ_R
end

# ============================================================================
# PDE Evolution (1D with Reflecting Boundaries)
# ============================================================================

function apply_neumann_bc!(ρ::Vector{Float64})
    ρ[1] = ρ[2]
    ρ[end] = ρ[end-1]
end

function evolve_diffusion!(ρ::Vector{Float64}, dx::Float64, dt::Float64, D::Float64)
    n = length(ρ)
    ρ_new = copy(ρ)
    for i in 2:(n-1)
        laplacian = (ρ[i+1] - 2*ρ[i] + ρ[i-1]) / dx^2
        ρ_new[i] = ρ[i] + dt * D * laplacian
    end
    ρ .= ρ_new
    apply_neumann_bc!(ρ)
    ρ .= max.(ρ, 0.0)
end

function evolve_advection!(ρ::Vector{Float64}, dx::Float64, dt::Float64, v::Float64)
    n = length(ρ)
    ρ_new = copy(ρ)
    for i in 2:(n-1)
        if v > 0
            dρ_dx = (ρ[i] - ρ[i-1]) / dx
        else
            dρ_dx = (ρ[i+1] - ρ[i]) / dx
        end
        ρ_new[i] = ρ[i] - dt * v * dρ_dx
    end
    ρ .= ρ_new
    apply_neumann_bc!(ρ)
    ρ .= max.(ρ, 0.0)
end

function compute_interaction_gradient(grid::Vector{Float64}, ρ_other::Vector{Float64}, dx::Float64)
    μ_other = sum(grid .* ρ_other) * dx / (sum(ρ_other) * dx + 1e-10)
    return μ_other
end

function evolve_pursuit_evasion!(ρ_G::Vector{Float64}, ρ_R::Vector{Float64},
                                  grid::Vector{Float64}, dx::Float64, dt::Float64,
                                  α::Float64, β::Float64, γ::Float64, center::Float64)
    n = length(grid)
    ρ_G_new = copy(ρ_G)
    ρ_R_new = copy(ρ_R)

    μ_R = compute_interaction_gradient(grid, ρ_R, dx)
    μ_G = compute_interaction_gradient(grid, ρ_G, dx)

    for i in 2:(n-1)
        g = grid[i]
        r = grid[i]

        v_G_flee = -α * (μ_R - center)
        v_G_center = -γ * (g - center)
        v_G = v_G_flee + v_G_center

        v_R_chase = β * (μ_G - center)
        v_R_center = -γ * (r - center)
        v_R = v_R_chase + v_R_center

        if v_G > 0
            dρ_G_dx = (ρ_G[i] - ρ_G[i-1]) / dx
        else
            dρ_G_dx = (ρ_G[i+1] - ρ_G[i]) / dx
        end

        if v_R > 0
            dρ_R_dx = (ρ_R[i] - ρ_R[i-1]) / dx
        else
            dρ_R_dx = (ρ_R[i+1] - ρ_R[i]) / dx
        end

        ρ_G_new[i] = ρ_G[i] - dt * v_G * dρ_G_dx
        ρ_R_new[i] = ρ_R[i] - dt * v_R * dρ_R_dx
    end

    ρ_G .= ρ_G_new
    ρ_R .= ρ_R_new
    apply_neumann_bc!(ρ_G)
    apply_neumann_bc!(ρ_R)
    ρ_G .= max.(ρ_G, 0.0)
    ρ_R .= max.(ρ_R, 0.0)
end

# ============================================================================
# Bound Heat Computation
# ============================================================================

function compute_bound_heat(grid::Vector{Float64}, ρ_G::Vector{Float64}, ρ_R::Vector{Float64})
    n = length(grid)
    dx = grid[2] - grid[1]

    c_G = sum(ρ_G) * dx
    c_R = sum(ρ_R) * dx
    ρ_G_norm = ρ_G ./ (c_G + 1e-10)
    ρ_R_norm = ρ_R ./ (c_R + 1e-10)

    bound_heat = zeros(n, n)
    for (i, g) in enumerate(grid)
        for (j, r) in enumerate(grid)
            K = g * r
            bound_heat[i, j] = K * ρ_G_norm[i] * ρ_R_norm[j]
        end
    end

    return bound_heat
end

# ============================================================================
# Full Evolution with Tracking
# ============================================================================

"""
Compute centroid of a distribution.
"""
function compute_centroid(grid::Vector{Float64}, ρ::Vector{Float64}, dx::Float64)
    mass = sum(ρ) * dx
    if mass < 1e-10
        return 0.5  # Default to center if no mass
    end
    return sum(grid .* ρ) * dx / mass
end

"""
Run evolution and track snapshots + centroid trajectory.
"""
function run_evolution_with_tracking(regime::Symbol, grid::Vector{Float64},
                                      ρ_G_init::Vector{Float64}, ρ_R_init::Vector{Float64},
                                      T::Float64, dt::Float64, snapshot_times::Vector{Float64})
    ρ_G = copy(ρ_G_init)
    ρ_R = copy(ρ_R_init)
    dx = grid[2] - grid[1]
    n_steps = Int(round(T / dt))

    # Storage for snapshots
    snapshots = Dict{Float64, NamedTuple}()

    # Storage for trajectories (sample every few steps)
    trajectory_times = Float64[]
    μ_G_trajectory = Float64[]
    μ_R_trajectory = Float64[]

    current_time = 0.0
    snapshot_idx = 1

    # Store initial state
    if snapshot_times[1] ≈ 0.0
        h = compute_bound_heat(grid, ρ_G, ρ_R)
        snapshots[0.0] = (ρ_G = copy(ρ_G), ρ_R = copy(ρ_R), heat = h)
        snapshot_idx = 2
    end

    push!(trajectory_times, 0.0)
    push!(μ_G_trajectory, compute_centroid(grid, ρ_G, dx))
    push!(μ_R_trajectory, compute_centroid(grid, ρ_R, dx))

    for step in 1:n_steps
        if regime == :static
            # No evolution
        elseif regime == :diffusion
            evolve_diffusion!(ρ_G, dx, dt, D_COEFF)
            evolve_diffusion!(ρ_R, dx, dt, D_COEFF)
        elseif regime == :advection
            evolve_advection!(ρ_G, dx, dt, V_G)
            evolve_advection!(ρ_R, dx, dt, V_R)
        elseif regime == :pursuit_evasion
            evolve_pursuit_evasion!(ρ_G, ρ_R, grid, dx, dt, ALPHA, BETA, GAMMA, CENTER)
        end

        current_time += dt

        # Record trajectory (every 10 steps to avoid too many points)
        if step % 10 == 0
            push!(trajectory_times, current_time)
            push!(μ_G_trajectory, compute_centroid(grid, ρ_G, dx))
            push!(μ_R_trajectory, compute_centroid(grid, ρ_R, dx))
        end

        # Check for snapshot
        if snapshot_idx <= length(snapshot_times) && current_time >= snapshot_times[snapshot_idx] - dt/2
            t_snap = snapshot_times[snapshot_idx]
            h = compute_bound_heat(grid, ρ_G, ρ_R)
            snapshots[t_snap] = (ρ_G = copy(ρ_G), ρ_R = copy(ρ_R), heat = h)
            snapshot_idx += 1
        end
    end

    return snapshots, trajectory_times, μ_G_trajectory, μ_R_trajectory
end

# ============================================================================
# Data Computation
# ============================================================================

const SNAPSHOT_TIMES = [0.0, 0.5, 1.0, 1.5, 2.0]

function compute_figure_B_data()
    println("=" ^ 60)
    println("Figure B: Computing PDE Evolution Data (Reflecting BC)")
    println("=" ^ 60)

    grid = collect(range(0, 1, length=GRID_RES))
    dx = grid[2] - grid[1]
    ρ_G_init, ρ_R_init = compute_initial_intensities(grid)

    println("Initial conditions:")
    println("  μ_G = " * string(MU_G) * ", μ_R = " * string(MU_R) * ", σ = " * string(SIGMA))
    println("  Total mass G: " * string(round(sum(ρ_G_init) * dx, digits=4)))
    println("  Total mass R: " * string(round(sum(ρ_R_init) * dx, digits=4)))

    regimes = [:diffusion, :advection, :pursuit_evasion]
    regime_names = ["Diffusion", "Advection", "Pursuit-Evasion"]

    # Storage for all regime data
    all_snapshots = Dict{Symbol, Dict{Float64, NamedTuple}}()
    all_trajectories = Dict{Symbol, NamedTuple}()

    for regime in regimes
        println("\nRunning " * string(regime) * " regime...")

        snapshots, traj_times, μ_G_traj, μ_R_traj = run_evolution_with_tracking(
            regime, grid, ρ_G_init, ρ_R_init, T_FINAL, DT, SNAPSHOT_TIMES)

        all_snapshots[regime] = snapshots
        all_trajectories[regime] = (times = traj_times, μ_G = μ_G_traj, μ_R = μ_R_traj)

        println("  Final μ_G: " * string(round(μ_G_traj[end], digits=4)))
        println("  Final μ_R: " * string(round(μ_R_traj[end], digits=4)))
    end

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        grid = grid,
        regimes = regimes,
        regime_names = regime_names,
        snapshot_times = SNAPSHOT_TIMES,
        all_snapshots = all_snapshots,
        all_trajectories = all_trajectories,
        T_FINAL = T_FINAL,
        params = (μ_G = MU_G, μ_R = MU_R, σ = SIGMA, D = D_COEFF, V_G = V_G, V_R = V_R),
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

function load_figure_B_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

function plot_figure_B(data)
    println("=" ^ 60)
    println("Figure B: Generating Plot")
    println("=" ^ 60)

    grid = data.grid
    regimes = data.regimes
    regime_names = data.regime_names
    snapshot_times = data.snapshot_times
    all_snapshots = data.all_snapshots
    all_trajectories = data.all_trajectories
    T_FINAL = data.T_FINAL

    # Select 3 time points to show
    display_times = [0.0, 1.0, 2.0]
    colormap = :viridis

    # Compute global color limits
    all_heats = Float64[]
    for regime in regimes
        for t in display_times
            if haskey(all_snapshots[regime], t)
                append!(all_heats, vec(all_snapshots[regime][t].heat))
            end
        end
    end
    clims = (0, maximum(all_heats))

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
        fig = Figure(size=(1400, 800))

        for (row, (regime, name)) in enumerate(zip(regimes, regime_names))
            snapshots = all_snapshots[regime]
            traj = all_trajectories[regime]

            # Columns 1-3: Heatmaps at different times
            last_hm = nothing
            for (col, t) in enumerate(display_times)
                # Build title string without interpolation inside L"..."
                t_str = string(round(t, digits=1))
                title_str = col == 1 ? name * " (t=" * t_str * ")" : "t=" * t_str
                ax = Axis(fig[row, col],
                    xlabel = row == 3 ? L"$g$" : "",
                    ylabel = col == 1 ? L"$r$" : "",
                    title = title_str,
                    aspect = 1)

                if haskey(snapshots, t)
                    last_hm = heatmap!(ax, grid, grid, snapshots[t].heat, colormap=colormap, colorrange=clims)
                end
            end

            # Column 5: Trajectory plot (position on x-axis, time on y-axis)
            ax_traj = Axis(fig[row, 5],
                xlabel = row == 3 ? "Position" : "",
                ylabel = row == 1 ? L"Time $t$" : "",
                title = "Centroid trajectory",
                aspect = 1)

            # Draw trajectory lines
            lines!(ax_traj, traj.μ_G, traj.times, color=:blue, linewidth=1.5, label=L"$\mu_G(t)$")
            lines!(ax_traj, traj.μ_R, traj.times, color=:red, linewidth=1.5, label=L"$\mu_R(t)$")

            # Show intermediate points (subsample to ~20 points for clarity)
            n_pts = length(traj.times)
            step = max(1, n_pts ÷ 20)
            idx = 1:step:n_pts
            scatter!(ax_traj, traj.μ_G[idx], traj.times[idx], color=:blue, markersize=6, alpha=0.7)
            scatter!(ax_traj, traj.μ_R[idx], traj.times[idx], color=:red, markersize=6, alpha=0.7)

            # Mark start (circle) and end (star) points
            scatter!(ax_traj, [traj.μ_G[1]], [traj.times[1]], color=:blue, markersize=12, marker=:circle)
            scatter!(ax_traj, [traj.μ_R[1]], [traj.times[1]], color=:red, markersize=12, marker=:circle)
            scatter!(ax_traj, [traj.μ_G[end]], [traj.times[end]], color=:blue, markersize=12, marker=:star5)
            scatter!(ax_traj, [traj.μ_R[end]], [traj.times[end]], color=:red, markersize=12, marker=:star5)

            xlims!(ax_traj, 0, 1)
            ylims!(ax_traj, 0, T_FINAL)

            if row == 1
                axislegend(ax_traj, position=:rt)
            end
        end

        # Add shared colorbar between heatmaps and trajectory (column 4, spanning all rows)
        Colorbar(fig[1:3, 4], colormap=colormap, colorrange=clims, label=L"$\bar{h}(g,r)$")

        # Set column widths: heatmaps, colorbar, trajectory
        colsize!(fig.layout, 1, Relative(0.2))
        colsize!(fig.layout, 2, Relative(0.2))
        colsize!(fig.layout, 3, Relative(0.2))
        colsize!(fig.layout, 4, Relative(0.05))
        colsize!(fig.layout, 5, Relative(0.2))

        Label(fig[0, 1:5], L"Bound Heat $\bar{h}$ Evolution (Reflecting BC)", fontsize=22)

        fig
    end

    output_file = joinpath(OUTPUT_DIR, "figure_B_pde_reflecting.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_B_data()
        plot_figure_B(data)
    else
        data = compute_figure_B_data()
        plot_figure_B(data)
    end
    return 0
end
