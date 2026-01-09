# Figure B': PDE Evolution (Absorbing Boundary)
#
# Dimension: d = 1
# Purpose: Show wave-like transients when mass can leave the system
#
# Layout: Time series (5 panels) + mass evolution
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_B_pde_absorbing.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_B_pde_absorbing.jl --plot   # plot from saved data

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
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_B_prime_data.jld2")
const GRID_RES = 100

const MU_G = 0.3
const MU_R = 0.7
const SIGMA = 0.12

const T_FINAL = 1.0
const DT = 0.001
const D_COEFF = 0.05

const SNAPSHOT_TIMES = [0.0, 0.25, 0.5, 0.75, 1.0]

# ============================================================================
# Intensity Functions
# ============================================================================

# Note: Uses truncated_gaussian_sigma from IDPG library
const truncated_gaussian = truncated_gaussian_sigma

function compute_initial_intensities(grid::Vector{Float64})
    ρ_G, ρ_R = compute_1d_marginals(grid, MU_G, MU_R, SIGMA)
    return ρ_G, ρ_R
end

# ============================================================================
# PDE Evolution (Absorbing BC)
# ============================================================================

function apply_dirichlet_bc!(ρ::Vector{Float64})
    ρ[1] = 0.0
    ρ[end] = 0.0
end

function evolve_diffusion_absorbing!(ρ::Vector{Float64}, dx::Float64, dt::Float64, D::Float64)
    n = length(ρ)
    ρ_new = copy(ρ)
    for i in 2:(n-1)
        laplacian = (ρ[i+1] - 2*ρ[i] + ρ[i-1]) / dx^2
        ρ_new[i] = ρ[i] + dt * D * laplacian
    end
    ρ .= ρ_new
    apply_dirichlet_bc!(ρ)
    ρ .= max.(ρ, 0.0)
end

# ============================================================================
# Bound Heat Computation
# ============================================================================

function compute_bound_heat_unnorm(grid::Vector{Float64}, ρ_G::Vector{Float64}, ρ_R::Vector{Float64})
    n = length(grid)
    bound_heat = zeros(n, n)
    for (i, g) in enumerate(grid)
        for (j, r) in enumerate(grid)
            K = g * r
            bound_heat[i, j] = K * ρ_G[i] * ρ_R[j]
        end
    end
    return bound_heat
end

# ============================================================================
# Data Computation
# ============================================================================

function compute_figure_B_prime_data()
    println("=" ^ 60)
    println("Figure B': Computing PDE Evolution Data (Absorbing BC)")
    println("=" ^ 60)

    grid = collect(range(0, 1, length=GRID_RES))
    dx = grid[2] - grid[1]
    ρ_G_init, ρ_R_init = compute_initial_intensities(grid)

    println("Initial conditions:")
    println("  μ_G = " * string(MU_G) * ", μ_R = " * string(MU_R) * ", σ = " * string(SIGMA))
    println("  Initial mass G: " * string(round(sum(ρ_G_init) * dx, digits=4)))
    println("  Initial mass R: " * string(round(sum(ρ_R_init) * dx, digits=4)))

    ρ_G = copy(ρ_G_init)
    ρ_R = copy(ρ_R_init)

    snapshots = Dict{Float64, NamedTuple}()
    mass_history = Float64[]
    time_history = Float64[]

    current_time = 0.0
    snapshot_idx = 1

    # Store initial
    if SNAPSHOT_TIMES[1] ≈ 0.0
        h = compute_bound_heat_unnorm(grid, ρ_G, ρ_R)
        snapshots[0.0] = (ρ_G = copy(ρ_G), ρ_R = copy(ρ_R), heat = h)
        snapshot_idx = 2
    end

    push!(mass_history, sum(ρ_G) * dx * sum(ρ_R) * dx)
    push!(time_history, 0.0)

    println("\nRunning diffusion with absorbing BC...")
    while current_time < T_FINAL
        evolve_diffusion_absorbing!(ρ_G, dx, DT, D_COEFF)
        evolve_diffusion_absorbing!(ρ_R, dx, DT, D_COEFF)
        current_time += DT

        if length(time_history) < 1000 || mod(length(time_history), 10) == 0
            push!(mass_history, sum(ρ_G) * dx * sum(ρ_R) * dx)
            push!(time_history, current_time)
        end

        if snapshot_idx <= length(SNAPSHOT_TIMES) && current_time >= SNAPSHOT_TIMES[snapshot_idx] - DT/2
            t_snap = SNAPSHOT_TIMES[snapshot_idx]
            h = compute_bound_heat_unnorm(grid, ρ_G, ρ_R)
            snapshots[t_snap] = (ρ_G = copy(ρ_G), ρ_R = copy(ρ_R), heat = h)
            snapshot_idx += 1
        end
    end

    println("  Final mass G: " * string(round(sum(ρ_G) * dx, digits=4)))
    println("  Final mass R: " * string(round(sum(ρ_R) * dx, digits=4)))
    println("  Mass lost: " * string(round(100 * (1 - mass_history[end] / mass_history[1]), digits=1)) * "%")

    mkpath(OUTPUT_DIR)
    data = (
        grid = grid,
        snapshots = snapshots,
        snapshot_times = SNAPSHOT_TIMES,
        time_history = time_history,
        mass_history = mass_history,
        params = (μ_G = MU_G, μ_R = MU_R, σ = SIGMA, D = D_COEFF),
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

function load_figure_B_prime_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

function plot_figure_B_prime(data)
    println("=" ^ 60)
    println("Figure B': Generating Plot")
    println("=" ^ 60)

    grid = data.grid
    dx = grid[2] - grid[1]
    snapshots = data.snapshots
    snapshot_times = data.snapshot_times
    time_history = data.time_history
    mass_history = data.mass_history
    colormap = :viridis

    all_heats = vcat([vec(snapshots[t].heat) for t in snapshot_times]...)
    clims = (0, maximum(all_heats))

    # Use a 2x3 grid for heatmaps with equal column widths
    panel_positions = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2)]

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
        fig = Figure(size=(1200, 800))

        # Store last heatmap for colorbar
        last_hm = nothing

        for (idx, t) in enumerate(snapshot_times)
            row, col = panel_positions[idx]
            h = snapshots[t].heat
            mass_G = sum(snapshots[t].ρ_G) * dx
            mass_R = sum(snapshots[t].ρ_R) * dx

            cG_str = string(round(mass_G, digits=2))
            cR_str = string(round(mass_R, digits=2))
            t_str = string(t)
            # Use latexstring to build title with proper subscripts
            title_str = latexstring("t=" * t_str * " (\$c_G\$=" * cG_str * ", \$c_R\$=" * cR_str * ")")
            ax = Axis(fig[row, col],
                xlabel = row == 2 ? L"$g$" : "",
                ylabel = col == 1 ? L"$r$" : "",
                title = title_str,
                aspect = 1)
            last_hm = heatmap!(ax, grid, grid, h, colormap=colormap, colorrange=clims)
        end

        # Put colorbar next to last heatmap, matching panel height
        Colorbar(fig[2, 3], last_hm, label=L"$\bar{h}(g,r)$", height=Relative(0.7))

        # Set equal widths for heatmap columns
        colsize!(fig.layout, 1, Relative(0.22))
        colsize!(fig.layout, 2, Relative(0.22))
        colsize!(fig.layout, 3, Relative(0.22))

        ax_mass = Axis(fig[1:2, 4],
            xlabel = L"Time $t$",
            ylabel = L"Total Heat $\iint\bar{h}$",
            title = "Mass Decay (Absorbing BC)")

        lines!(ax_mass, time_history, mass_history, color=:blue, linewidth=2)
        for t in snapshot_times
            vlines!(ax_mass, [t], color=:red, linestyle=:dash, alpha=0.5)
        end

        mass_lost_str = string(round(100 * (1 - mass_history[end] / mass_history[1]), digits=1))
        text!(ax_mass, T_FINAL * 0.5, mass_history[1] * 0.9,
              text = "Mass lost: " * mass_lost_str * "%",
              fontsize = 14)

        Label(fig[0, 1:4], "Bound Heat Evolution with Absorbing Boundary (Diffusion)", fontsize=22)

        fig
    end

    output_file = joinpath(OUTPUT_DIR, "figure_B_pde_absorbing.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_B_prime_data()
        plot_figure_B_prime(data)
    else
        data = compute_figure_B_prime_data()
        plot_figure_B_prime(data)
    end
    return 0
end
