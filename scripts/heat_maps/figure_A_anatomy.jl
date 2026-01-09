# Figure A: Basic Heat Anatomy (Product Intensity)
#
# Dimension: d = 1
# Purpose: Show that raw, bound, and edge heats have the same shape under product intensity
#
# Panels (2×2):
#   (a) Marginals: ρ_G(g) and ρ_R(r) as 1D curves
#   (b) Bound heat: h̄(g,r) heatmap
#   (c) Kernel: K(g,r) = g·r heatmap
#   (d) Summary: Equations and integrated values
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_A_anatomy.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_A_anatomy.jl --plot   # plot from saved data

using IDPG
using CairoMakie
using LaTeXStrings
using Distributions
using Statistics
using JLD2

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "heat_maps")
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_A_data.jld2")
const GRID_RES = 100

# Intensity parameters (d=1)
const G_CENTER = 0.3   # Green intensity centered at g₀ = 0.3
const R_CENTER = 0.7   # Red intensity centered at r₀ = 0.7
const KAPPA = 20.0     # Concentration parameter

# ============================================================================
# Intensity Functions (Truncated Gaussians on [0,1])
# ============================================================================

# Note: Uses truncated_gaussian_kappa and compute_1d_marginals from IDPG library
# truncated_gaussian_kappa(x, μ, κ) computes truncated Gaussian density on [0,1]

"""
Compute intensities on a grid using library functions.
"""
function compute_intensities(grid::Vector{Float64}, g_center::Float64, r_center::Float64, κ::Float64)
    σ = 1.0 / sqrt(κ)
    ρ_G, ρ_R = compute_1d_marginals(grid, g_center, r_center, σ)
    return ρ_G, ρ_R
end

# ============================================================================
# Heat Computations
# ============================================================================

"""
Compute the kernel K(g,r) = g·r on a 2D grid.
"""
function compute_kernel(grid::Vector{Float64})
    n = length(grid)
    K = zeros(n, n)
    for (i, g) in enumerate(grid)
        for (j, r) in enumerate(grid)
            K[i, j] = g * r
        end
    end
    return K
end

"""
Compute bound heat and related quantities.
"""
function compute_heats(grid::Vector{Float64}, ρ_G::Vector{Float64}, ρ_R::Vector{Float64})
    n = length(grid)
    dx = grid[2] - grid[1]

    # Total intensities (charges)
    c_G = sum(ρ_G) * dx
    c_R = sum(ρ_R) * dx
    Λ = c_G * c_R

    # Centroids (unnormalized)
    μ_G_unnorm = sum(grid .* ρ_G) * dx
    μ_R_unnorm = sum(grid .* ρ_R) * dx

    # Normalized centroids
    μ_G_tilde = μ_G_unnorm / c_G
    μ_R_tilde = μ_R_unnorm / c_R

    # Mean connection probability
    p_bar = μ_G_tilde * μ_R_tilde

    # Kernel
    K = compute_kernel(grid)

    # Bound heat: h̄(g,r) = K(g,r) · ρ_G(g) · ρ_R(r)
    bound_heat = zeros(n, n)
    for i in 1:n
        for j in 1:n
            bound_heat[i, j] = K[i, j] * ρ_G[i] * ρ_R[j]
        end
    end

    # Total bound heat
    H_bar_total = sum(bound_heat) * dx^2

    # Raw heat total: H(Ω,Ω) = Λ² · p̄
    H_total = Λ^2 * p_bar

    # Edge heat total: H⃗(Ω,Ω) = Λ · p̄ / 2
    H_edge_total = Λ * p_bar / 2

    stats = (
        c_G = c_G,
        c_R = c_R,
        Λ = Λ,
        μ_G_tilde = μ_G_tilde,
        μ_R_tilde = μ_R_tilde,
        p_bar = p_bar,
        H_bar_total = H_bar_total,
        H_total = H_total,
        H_edge_total = H_edge_total,
    )

    return K, bound_heat, stats
end

# ============================================================================
# Data Computation
# ============================================================================

function compute_figure_A_data()
    println("=" ^ 60)
    println("Figure A: Computing Heat Map Data (Product Intensity)")
    println("=" ^ 60)

    grid = collect(range(0, 1, length=GRID_RES))
    dx = grid[2] - grid[1]

    # Compute intensities
    ρ_G, ρ_R = compute_intensities(grid, G_CENTER, R_CENTER, KAPPA)

    println("Parameters:")
    println("  g₀ = " * string(G_CENTER) * ", r₀ = " * string(R_CENTER) * ", κ = " * string(KAPPA))

    # Compute heats
    K, bound_heat, stats = compute_heats(grid, ρ_G, ρ_R)

    println("\nCharges:")
    println("  c_G = " * string(round(stats.c_G, digits=4)))
    println("  c_R = " * string(round(stats.c_R, digits=4)))
    println("  Λ = c_G · c_R = " * string(round(stats.Λ, digits=4)))

    println("\nCentroids:")
    println("  μ̃_G = " * string(round(stats.μ_G_tilde, digits=4)))
    println("  μ̃_R = " * string(round(stats.μ_R_tilde, digits=4)))

    println("\nMean connection probability:")
    println("  p̄ = μ̃_G · μ̃_R = " * string(round(stats.p_bar, digits=4)))

    println("\nTotal heats:")
    println("  H̄(B,B) = " * string(round(stats.H_bar_total, digits=4)))
    println("  H(Ω,Ω) = Λ² · p̄ = " * string(round(stats.H_total, digits=4)))
    println("  H⃗(Ω,Ω) = Λ · p̄ / 2 = " * string(round(stats.H_edge_total, digits=4)))

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        grid = grid,
        ρ_G = ρ_G,
        ρ_R = ρ_R,
        K = K,
        bound_heat = bound_heat,
        stats = stats,
        params = (g_center = G_CENTER, r_center = R_CENTER, κ = KAPPA),
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

function load_figure_A_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

function plot_figure_A(data)
    println("=" ^ 60)
    println("Figure A: Generating Plot")
    println("=" ^ 60)

    grid = data.grid
    ρ_G = data.ρ_G
    ρ_R = data.ρ_R
    K = data.K
    bound_heat = data.bound_heat
    stats = data.stats
    params = data.params

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
        # Single row layout: 1×3
        fig = Figure(size=(1500, 500))

        # ---------------------------------------------------------------------
        # Panel (a): Marginal Intensities [1,1]
        # ---------------------------------------------------------------------
        ax_a = Axis(fig[1, 1],
            xlabel = "Coordinate",
            ylabel = "Intensity",
            title = "(a) Marginal Intensities")

        lines!(ax_a, grid, ρ_G, color=:green, linewidth=2.5, label=L"$\rho_G(g)$")
        lines!(ax_a, grid, ρ_R, color=:red, linewidth=2.5, label=L"$\rho_R(r)$")

        # Mark centers with vertical dashed lines
        vlines!(ax_a, [params.g_center], color=:green, linestyle=:dash, alpha=0.5)
        vlines!(ax_a, [params.r_center], color=:red, linestyle=:dash, alpha=0.5)

        axislegend(ax_a, position=:cb)

        # ---------------------------------------------------------------------
        # Panel (b): Bound Heat [1,2]
        # ---------------------------------------------------------------------
        ax_b = Axis(fig[1, 2],
            xlabel = L"$g$",
            ylabel = L"$r$",
            title = L"(b) Bound Heat $\bar{h}(g,r)$",
            aspect = 1)

        hm_b = heatmap!(ax_b, grid, grid, bound_heat, colormap=:viridis)
        Colorbar(fig[1, 2][1, 2], hm_b, label=L"$\bar{h}$")

        # ---------------------------------------------------------------------
        # Panel (c): Kernel [1,3]
        # ---------------------------------------------------------------------
        ax_c = Axis(fig[1, 3],
            xlabel = L"$g$",
            ylabel = L"$r$",
            title = L"(c) Kernel $K(g,r) = g \cdot r$",
            aspect = 1)

        hm_c = heatmap!(ax_c, grid, grid, K, colormap=:YlOrRd)
        Colorbar(fig[1, 3][1, 2], hm_c, label=L"$K$")

        # ---------------------------------------------------------------------
        # Overall title
        # ---------------------------------------------------------------------
        Label(fig[0, 1:3], "Heat Map Anatomy (Product Intensity, d=1)", fontsize=22)

        fig
    end

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_A_anatomy.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    # Save caption/summary text to file
    caption_file = joinpath(OUTPUT_DIR, "figure_A_caption.txt")
    c_G_str = string(round(stats.c_G, digits=4))
    c_R_str = string(round(stats.c_R, digits=4))
    μ_G_str = string(round(stats.μ_G_tilde, digits=4))
    μ_R_str = string(round(stats.μ_R_tilde, digits=4))
    Λ_str = string(round(stats.Λ, digits=4))
    p_str = string(round(stats.p_bar, digits=4))
    H_bar_str = string(round(stats.H_bar_total, digits=4))
    H_str = string(round(stats.H_total, digits=4))
    H_edge_str = string(round(stats.H_edge_total, digits=4))

    caption = """Figure A: Heat Map Anatomy (Product Intensity, d=1)

Parameters:
  g₀ = $(params.g_center), r₀ = $(params.r_center), κ = $(params.κ)

Charges:
  c_G = ∫ρ_G dg = $(c_G_str)
  c_R = ∫ρ_R dr = $(c_R_str)
  Λ = c_G · c_R = $(Λ_str)

Centroids:
  μ̃_G = $(μ_G_str)
  μ̃_R = $(μ_R_str)

Mean connection probability:
  p̄ = μ̃_G · μ̃_R = $(p_str)

Heat scaling (product intensity):
  H = Λ · H̄  (raw = charge × bound)
  H⃗ = H̄ / 2  (edge = bound / 2)

Total heats:
  H̄(B,B) = $(H_bar_str)
  H(Ω,Ω) = Λ² · p̄ = $(H_str)  (= E[|E|] node-centric)
  H⃗(Ω,Ω) = Λ · p̄ / 2 = $(H_edge_str)  (= E[L] edge-centric)

Key formula:
  h̄(g,r) = K(g,r) · ρ_G(g) · ρ_R(r)

Under product intensity, raw, bound, and edge heats all have the same shape
(just rescaled), so the 2D bound heat h̄(g,r) is a complete summary.
"""
    open(caption_file, "w") do f
        write(f, caption)
    end
    println("Saved caption: " * caption_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_A_data()
        plot_figure_A(data)
    else
        data = compute_figure_A_data()
        plot_figure_A(data)
    end
    return 0
end
