# Figure A'': Non-Product Intensity (4D Raw Heat)
#
# Dimension: d = 1
# Purpose: Show that bound heat only works under product intensity.
#          With non-product ρ(g,r), different slices of the 4D raw heat have different shapes.
#
# Site intensity: Two off-diagonal blobs (NOT a product)
#   ρ(g,r) = N((0.3, 0.7), σ²I) + N((0.7, 0.3), σ²I)
#
# Layout (2×3):
#   Row 1: (a) Site intensity ρ(g,r)  | (b) Slice r_s=0.3, g_t=0.3 | (c) Slice r_s=0.7, g_t=0.3
#   Row 2: (d) Marginals             | (e) Slice r_s=0.3, g_t=0.7 | (f) Slice r_s=0.7, g_t=0.7
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_A_double_prime_nonproduct.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_A_double_prime_nonproduct.jl --plot   # plot from saved data

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
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_A_double_prime_data.jld2")
const GRID_RES = 100

# Two off-diagonal blobs
const BLOB_1 = (g = 0.3, r = 0.7)  # High r, low g
const BLOB_2 = (g = 0.7, r = 0.3)  # High g, low r
const SIGMA = 0.12                 # Blob width

# Slice positions for 4D visualization
const SLICE_VALUES = [0.3, 0.7]

# ============================================================================
# Intensity Functions
# ============================================================================

"""
Non-product site intensity: two off-diagonal Gaussian blobs.
This CANNOT be written as ρ_G(g) · ρ_R(r).
"""
function site_intensity(g::Float64, r::Float64)
    if g < 0.0 || g > 1.0 || r < 0.0 || r > 1.0
        return 0.0
    end
    blob1 = exp(-((g - BLOB_1.g)^2 + (r - BLOB_1.r)^2) / (2 * SIGMA^2))
    blob2 = exp(-((g - BLOB_2.g)^2 + (r - BLOB_2.r)^2) / (2 * SIGMA^2))
    return blob1 + blob2
end

"""
Compute the kernel K(g,r) = g·r.
"""
function kernel(g_s::Float64, r_t::Float64)
    return g_s * r_t
end

# ============================================================================
# Heat Computations
# ============================================================================

"""
Compute site intensity on 2D grid.
"""
function compute_site_intensity(grid::Vector{Float64})
    n = length(grid)
    ρ = zeros(n, n)
    for (i, g) in enumerate(grid)
        for (j, r) in enumerate(grid)
            ρ[i, j] = site_intensity(g, r)
        end
    end
    return ρ
end

"""
Compute marginals of the site intensity.
"""
function compute_marginals(grid::Vector{Float64}, ρ::Matrix{Float64})
    dx = grid[2] - grid[1]
    n = length(grid)

    # ρ^{(g)}(g) = ∫ρ(g,r) dr
    ρ_g = [sum(ρ[i, :]) * dx for i in 1:n]

    # ρ^{(r)}(r) = ∫ρ(g,r) dg
    ρ_r = [sum(ρ[:, j]) * dx for j in 1:n]

    return ρ_g, ρ_r
end

"""
Compute a slice of the 4D raw heat h(g_s, r_s, g_t, r_t) at fixed (r_s, g_t).
Returns h(g_s, r_t | r_s_fixed, g_t_fixed).
"""
function compute_heat_slice(grid::Vector{Float64}, r_s_fixed::Float64, g_t_fixed::Float64)
    n = length(grid)
    h = zeros(n, n)

    for (i, g_s) in enumerate(grid)
        for (j, r_t) in enumerate(grid)
            ρ_source = site_intensity(g_s, r_s_fixed)
            ρ_target = site_intensity(g_t_fixed, r_t)
            K = kernel(g_s, r_t)
            h[i, j] = K * ρ_source * ρ_target
        end
    end

    return h
end

# ============================================================================
# Data Computation
# ============================================================================

function compute_figure_A_double_prime_data()
    println("=" ^ 60)
    println("Figure A'': Computing Non-Product Intensity Data")
    println("=" ^ 60)

    grid = collect(range(0, 1, length=GRID_RES))
    dx = grid[2] - grid[1]

    println("Blob configuration:")
    println("  Blob 1: (g, r) = (" * string(BLOB_1.g) * ", " * string(BLOB_1.r) * ")")
    println("  Blob 2: (g, r) = (" * string(BLOB_2.g) * ", " * string(BLOB_2.r) * ")")
    println("  σ = " * string(SIGMA))

    # Compute site intensity
    ρ = compute_site_intensity(grid)

    # Compute marginals
    ρ_g, ρ_r = compute_marginals(grid, ρ)

    # Compute product of marginals (for comparison)
    ρ_product = zeros(GRID_RES, GRID_RES)
    for i in 1:GRID_RES
        for j in 1:GRID_RES
            ρ_product[i, j] = ρ_g[i] * ρ_r[j]
        end
    end

    println("\nNon-product test:")
    println("  ∫∫ρ(g,r) dg dr = " * string(round(sum(ρ) * dx^2, digits=4)))
    println("  ∫ρ_g(g) dg = " * string(round(sum(ρ_g) * dx, digits=4)))
    println("  ∫ρ_r(r) dr = " * string(round(sum(ρ_r) * dx, digits=4)))

    # Compute four slices
    slices = Dict{Tuple{Float64, Float64}, Matrix{Float64}}()
    for r_s in SLICE_VALUES
        for g_t in SLICE_VALUES
            h = compute_heat_slice(grid, r_s, g_t)
            slices[(r_s, g_t)] = h
            total = sum(h) * dx^2
            println("  Slice (r_s=" * string(r_s) * ", g_t=" * string(g_t) * "): total = " * string(round(total, digits=4)))
        end
    end

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        grid = grid,
        ρ = ρ,
        ρ_g = ρ_g,
        ρ_r = ρ_r,
        ρ_product = ρ_product,
        slices = slices,
        slice_values = SLICE_VALUES,
        params = (blob_1 = BLOB_1, blob_2 = BLOB_2, σ = SIGMA),
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

function load_figure_A_double_prime_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

function plot_figure_A_double_prime(data)
    println("=" ^ 60)
    println("Figure A'': Generating Plot")
    println("=" ^ 60)

    grid = data.grid
    ρ = data.ρ
    ρ_g = data.ρ_g
    ρ_r = data.ρ_r
    ρ_product = data.ρ_product
    slices = data.slices
    slice_values = data.slice_values
    params = data.params

    # Common colorscale for all heat slices
    all_slice_values = vcat([vec(slices[(r_s, g_t)]) for r_s in slice_values for g_t in slice_values]...)
    clims = (0, maximum(all_slice_values))

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
        # 2×3 layout
        fig = Figure(size=(1500, 900))

        # ---------------------------------------------------------------------
        # Panel (a): Site Intensity ρ(g,r) [1,1]
        # ---------------------------------------------------------------------
        ax_a = Axis(fig[1, 1],
            xlabel = L"$g$",
            ylabel = L"$r$",
            title = L"(a) Site Intensity $\rho(g,r)$",
            aspect = 1)

        hm_a = heatmap!(ax_a, grid, grid, ρ, colormap=:viridis)
        Colorbar(fig[1, 1][1, 2], hm_a, label=L"$\rho$")

        # Mark blob centers
        scatter!(ax_a, [params.blob_1.g], [params.blob_1.r], color=:white, marker=:cross, markersize=12)
        scatter!(ax_a, [params.blob_2.g], [params.blob_2.r], color=:white, marker=:cross, markersize=12)

        # ---------------------------------------------------------------------
        # Panel (b): Slice r_s=0.3, g_t=0.3 [1,2]
        # ---------------------------------------------------------------------
        ax_b = Axis(fig[1, 2],
            xlabel = L"$g_s$",
            ylabel = L"$r_t$",
            title = L"(b) $r_s=0.3$, $g_t=0.3$",
            aspect = 1)

        hm_b = heatmap!(ax_b, grid, grid, slices[(0.3, 0.3)], colormap=:inferno, colorrange=clims)
        Colorbar(fig[1, 2][1, 2], hm_b)

        # ---------------------------------------------------------------------
        # Panel (c): Slice r_s=0.7, g_t=0.3 [1,3]
        # ---------------------------------------------------------------------
        ax_c = Axis(fig[1, 3],
            xlabel = L"$g_s$",
            ylabel = L"$r_t$",
            title = L"(c) $r_s=0.7$, $g_t=0.3$",
            aspect = 1)

        hm_c = heatmap!(ax_c, grid, grid, slices[(0.7, 0.3)], colormap=:inferno, colorrange=clims)
        Colorbar(fig[1, 3][1, 2], hm_c)

        # ---------------------------------------------------------------------
        # Panel (d): Marginals [2,1]
        # ---------------------------------------------------------------------
        ax_d = Axis(fig[2, 1],
            xlabel = "Coordinate",
            ylabel = "Marginal Intensity",
            title = "(d) Marginals of the site intensity")

        # Use different styles since ρ_g and ρ_r are identical (symmetric blobs)
        lines!(ax_d, grid, ρ_g, color=:green, linewidth=3, label=L"$\rho^{(g)}$")
        lines!(ax_d, grid, ρ_r, color=:red, linewidth=2, linestyle=:dash, label=L"$\rho^{(r)}$")

        axislegend(ax_d, position=:cb)

        # ---------------------------------------------------------------------
        # Panel (e): Slice r_s=0.3, g_t=0.7 [2,2]
        # ---------------------------------------------------------------------
        ax_e = Axis(fig[2, 2],
            xlabel = L"$g_s$",
            ylabel = L"$r_t$",
            title = L"(e) $r_s=0.3$, $g_t=0.7$",
            aspect = 1)

        hm_e = heatmap!(ax_e, grid, grid, slices[(0.3, 0.7)], colormap=:inferno, colorrange=clims)
        Colorbar(fig[2, 2][1, 2], hm_e)

        # ---------------------------------------------------------------------
        # Panel (f): Slice r_s=0.7, g_t=0.7 [2,3]
        # ---------------------------------------------------------------------
        ax_f = Axis(fig[2, 3],
            xlabel = L"$g_s$",
            ylabel = L"$r_t$",
            title = L"(f) $r_s=0.7$, $g_t=0.7$",
            aspect = 1)

        hm_f = heatmap!(ax_f, grid, grid, slices[(0.7, 0.7)], colormap=:inferno, colorrange=clims)
        Colorbar(fig[2, 3][1, 2], hm_f, label=L"$h$")

        # ---------------------------------------------------------------------
        # Overall title
        # ---------------------------------------------------------------------
        Label(fig[0, 1:3], "Non-Product Intensity: 4D Raw Heat Slices (d=1)", fontsize=22)

        fig
    end

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_A_double_prime_nonproduct.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    # Save caption text to file
    caption_file = joinpath(OUTPUT_DIR, "figure_A_double_prime_caption.txt")
    caption = """Figure A'': Non-Product Intensity (4D Raw Heat)

Site intensity: Two off-diagonal Gaussian blobs
  ρ(g,r) = N(($(params.blob_1.g), $(params.blob_1.r)), σ²I) + N(($(params.blob_2.g), $(params.blob_2.r)), σ²I)
  σ = $(params.σ)

This intensity is NOT a product: ρ(g,r) ≠ ρ^(g)(g) · ρ^(r)(r)

The 4D raw heat h(g_s, r_s, g_t, r_t) = K(g_s, r_t) · ρ(g_s, r_s) · ρ(g_t, r_t)
depends on all four coordinates. Panels (b,c,e,f) show slices at fixed (r_s, g_t):
  - Different slices have DIFFERENT SHAPES
  - This is because the "inactive" coordinates (r_s, g_t) don't factor out

Under product intensity ρ(g,r) = ρ_G(g)·ρ_R(r), all slices would have the same shape
(differing only by a scalar factor ρ_R(r_s)·ρ_G(g_t)), enabling the 2D bound heat summary.

With non-product intensity, there is no 2D summary — the full 4D structure matters.
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
        data = load_figure_A_double_prime_data()
        plot_figure_A_double_prime(data)
    else
        data = compute_figure_A_double_prime_data()
        plot_figure_A_double_prime(data)
    end
    return 0
end
