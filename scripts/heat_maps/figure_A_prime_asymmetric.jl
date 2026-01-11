# Figure A': Asymmetric Source-Target Intensity (d=1)
#
# Dimension: d = 1
# Purpose: Show that raw heat (node-centric) and edge-weighted heat (edge-centric)
#          have genuinely different shapes when source/target distributions differ
#
# Two species with different source/target roles:
#   - Species 1 (Prey-like): g_center=0.3, r_center=0.7, w_S=0.2, w_T=0.8
#   - Species 2 (Predator-like): g_center=0.7, r_center=0.3, w_S=0.8, w_T=0.2
#
# Panels (2×2):
#   (a) Source vs Target marginals: ρ_S and ρ_T showing asymmetry
#   (b) Raw heat: h(g,r) node-centric (combined site intensity)
#   (c) Edge-weighted heat: h_E(g,r) edge-centric (asymmetric source/target)
#   (d) Comparison: ratio or difference between (b) and (c)
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_A_prime_asymmetric.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_A_prime_asymmetric.jl --plot   # plot from saved data

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
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_A_prime_data.jld2")
const GRID_RES = 100

# Species parameters (d=1)
# Each species has separate g and r centers
const SPECIES = [
    # Species 1 (Prey-like): low proposal (g), high acceptance (r)
    (g_center = 0.3, r_center = 0.7, w_S = 0.2, w_T = 0.8),
    # Species 2 (Predator-like): high proposal (g), low acceptance (r)
    (g_center = 0.7, r_center = 0.3, w_S = 0.8, w_T = 0.2),
]

const KAPPA = 20.0  # Concentration parameter

# ============================================================================
# Intensity Functions (Truncated Gaussians on [0,1])
# ============================================================================

# Note: Uses truncated_gaussian_kappa from IDPG library
# Alias for backward compatibility
const truncated_gaussian = truncated_gaussian_kappa

# ============================================================================
# Heat Computations
# ============================================================================

"""
Compute all intensities and heats for asymmetric source-target case.
"""
function compute_asymmetric_heats(grid::Vector{Float64})
    n = length(grid)
    dx = grid[2] - grid[1]

    # Compute species-specific marginals in g and r
    # ρ_{G,m}(g) = Gaussian centered at g_center_m
    # ρ_{R,m}(r) = Gaussian centered at r_center_m
    ρ_G_species = [[truncated_gaussian(g, sp.g_center, KAPPA) for g in grid] for sp in SPECIES]
    ρ_R_species = [[truncated_gaussian(r, sp.r_center, KAPPA) for r in grid] for sp in SPECIES]

    # Combined site intensity (for node-centric): sum over species
    # ρ_G_combined(g) = Σ_m ρ_{G,m}(g)
    # ρ_R_combined(r) = Σ_m ρ_{R,m}(r)
    ρ_G_combined = sum(ρ_G_species)
    ρ_R_combined = sum(ρ_R_species)

    # Source marginal (weighted by w_S): ρ_S(g) = Σ_m w_{S,m} · ρ_{G,m}(g)
    ρ_S = sum(sp.w_S .* ρ_G_species[i] for (i, sp) in enumerate(SPECIES))

    # Target marginal (weighted by w_T): ρ_T(r) = Σ_m w_{T,m} · ρ_{R,m}(r)
    ρ_T = sum(sp.w_T .* ρ_R_species[i] for (i, sp) in enumerate(SPECIES))

    # Kernel K(g,r) = g·r (outer product)
    K = grid * grid'

    # Joint intensity for mixture of products (sum of products, NOT product of sums):
    # ρ(g,r) = Σ_m ρ_{G,m}(g) · ρ_{R,m}(r)
    ρ_joint = zeros(n, n)
    for i in 1:n
        for j in 1:n
            for m in 1:length(SPECIES)
                ρ_joint[i, j] += ρ_G_species[m][i] * ρ_R_species[m][j]
            end
        end
    end

    # Raw heat (node-centric view): h̄(g,r) = K(g,r) · ρ(g,r)
    h_node = zeros(n, n)
    for i in 1:n
        for j in 1:n
            h_node[i, j] = K[i, j] * ρ_joint[i, j]
        end
    end

    # Edge-weighted heat (edge-centric with asymmetry):
    # Source and target sampled independently from ρ_S and ρ_T
    # h_edge(g,r) = K(g,r) · ρ_S(g) · ρ_T(r)
    h_edge = zeros(n, n)
    for i in 1:n
        for j in 1:n
            h_edge[i, j] = K[i, j] * ρ_S[i] * ρ_T[j]
        end
    end

    # Compute totals
    h_node_total = sum(h_node) * dx^2
    h_edge_total = sum(h_edge) * dx^2

    # Ratio for comparison (NaN where heats are negligible)
    ratio = fill(NaN, n, n)
    for i in 1:n
        for j in 1:n
            if h_node[i, j] > 1e-10
                ratio[i, j] = h_edge[i, j] / h_node[i, j]
            end
        end
    end

    stats = (
        c_S = sum(ρ_S) * dx,
        c_T = sum(ρ_T) * dx,
        c_G = sum(ρ_G_combined) * dx,
        c_R = sum(ρ_R_combined) * dx,
        h_node_total = h_node_total,
        h_edge_total = h_edge_total,
    )

    return (
        ρ_G_species = ρ_G_species,
        ρ_R_species = ρ_R_species,
        ρ_G_combined = ρ_G_combined,
        ρ_R_combined = ρ_R_combined,
        ρ_S = ρ_S,
        ρ_T = ρ_T,
        K = K,
        h_node = h_node,
        h_edge = h_edge,
        ratio = ratio,
        stats = stats,
    )
end

# ============================================================================
# Data Computation
# ============================================================================

function compute_figure_A_prime_data()
    println("=" ^ 60)
    println("Figure A': Computing Asymmetric Source-Target Data")
    println("=" ^ 60)

    grid = collect(range(0, 1, length=GRID_RES))
    dx = grid[2] - grid[1]

    println("Species configuration:")
    for (i, sp) in enumerate(SPECIES)
        println("  Species " * string(i) * ": g_center=" * string(sp.g_center) *
                ", r_center=" * string(sp.r_center) *
                ", w_S=" * string(sp.w_S) * ", w_T=" * string(sp.w_T))
    end

    # Compute all heats
    result = compute_asymmetric_heats(grid)

    println("\nMarginal integrals:")
    println("  ∫ρ_S dg = " * string(round(result.stats.c_S, digits=3)))
    println("  ∫ρ_T dr = " * string(round(result.stats.c_T, digits=3)))
    println("  ∫ρ_G dg = " * string(round(result.stats.c_G, digits=3)))
    println("  ∫ρ_R dr = " * string(round(result.stats.c_R, digits=3)))

    println("\nHeat totals:")
    println("  h_node (node-centric): " * string(round(result.stats.h_node_total, digits=4)))
    println("  h_edge (edge-centric): " * string(round(result.stats.h_edge_total, digits=4)))

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        grid = grid,
        ρ_G_species = result.ρ_G_species,
        ρ_R_species = result.ρ_R_species,
        ρ_G_combined = result.ρ_G_combined,
        ρ_R_combined = result.ρ_R_combined,
        ρ_S = result.ρ_S,
        ρ_T = result.ρ_T,
        K = result.K,
        h_node = result.h_node,
        h_edge = result.h_edge,
        ratio = result.ratio,
        stats = result.stats,
        species = SPECIES,
        κ = KAPPA,
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

function load_figure_A_prime_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

function plot_figure_A_prime(data)
    println("=" ^ 60)
    println("Figure A': Generating Plot")
    println("=" ^ 60)

    grid = data.grid
    ρ_S = data.ρ_S
    ρ_T = data.ρ_T
    ρ_G_combined = data.ρ_G_combined
    ρ_R_combined = data.ρ_R_combined
    h_node = data.h_node
    h_edge = data.h_edge
    ratio = data.ratio
    stats = data.stats
    species = data.species

    # 2×2 grid layout
    fig = Figure(size=(1200, 1000))

    # -------------------------------------------------------------------------
    # Panel (a): Source vs Target Marginals [1,1]
    # -------------------------------------------------------------------------
    ax_a = Axis(fig[1, 1],
        xlabel = "Coordinate",
        ylabel = "Intensity",
        title = L"(a) Source $\rho_S$ vs Target $\rho_T$")

    # Plot source and target marginals (use different styles since they overlap)
    lines!(ax_a, grid, ρ_S, color=:orange, linewidth=3, label=L"$\rho_S$ (source)")
    lines!(ax_a, grid, ρ_T, color=:purple, linewidth=2, linestyle=:dash, label=L"$\rho_T$ (target)")

    axislegend(ax_a, position=:lt)

    # -------------------------------------------------------------------------
    # Panel (b): Raw Heat (Node-centric) [1,2]
    # -------------------------------------------------------------------------
    ax_b = Axis(fig[1, 2],
        xlabel = L"$g$",
        ylabel = L"$r$",
        title = L"(b) Node-centric $h(g,r)$",
        aspect = 1)

    hm_b = heatmap!(ax_b, grid, grid, h_node, colormap=:viridis)
    Colorbar(fig[1, 2][1, 2], hm_b, label=L"$h$")

    # -------------------------------------------------------------------------
    # Panel (c): Edge-weighted Heat (Edge-centric) [2,1]
    # -------------------------------------------------------------------------
    ax_c = Axis(fig[2, 1],
        xlabel = L"$g$",
        ylabel = L"$r$",
        title = L"(c) Edge-centric $h_{\mathcal{E}}(g,r)$",
        aspect = 1)

    hm_c = heatmap!(ax_c, grid, grid, h_edge, colormap=:inferno)
    Colorbar(fig[2, 1][1, 2], hm_c, label=L"$h_{\mathcal{E}}$")

    # -------------------------------------------------------------------------
    # Panel (d): Comparison (ratio) [2,2]
    # -------------------------------------------------------------------------
    ax_d = Axis(fig[2, 2],
        xlabel = L"$g$",
        ylabel = L"$r$",
        title = L"(d) Ratio $h_{\mathcal{E}}/h$",
        aspect = 1)

    hm_d = heatmap!(ax_d, grid, grid, ratio, colormap=:RdBu)
    Colorbar(fig[2, 2][1, 2], hm_d, label="Ratio")

    # -------------------------------------------------------------------------
    # Overall title
    # -------------------------------------------------------------------------
    Label(fig[0, 1:2], "Asymmetric Source-Target: Node vs Edge Centric (d=1)", fontsize=20)

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_A_prime_asymmetric.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    # Save caption text to file
    caption_file = joinpath(OUTPUT_DIR, "figure_A_prime_caption.txt")
    caption = """Figure A': Asymmetric Source-Target Intensity (d=1)

Two species with different source/target roles:
- Species 1 (Prey-like): g_center=$(species[1].g_center), r_center=$(species[1].r_center)
  w_S=$(species[1].w_S) (rarely initiates), w_T=$(species[1].w_T) (often targeted)
- Species 2 (Predator-like): g_center=$(species[2].g_center), r_center=$(species[2].r_center)
  w_S=$(species[2].w_S) (often initiates), w_T=$(species[2].w_T) (rarely targeted)

Source marginal: ρ_S(g) = Σ_m w_{S,m} · ρ_{G,m}(g) — dominated by predators
Target marginal: ρ_T(r) = Σ_m w_{T,m} · ρ_{R,m}(r) — dominated by prey

Key insight:
- Node-centric heat h(g,r) uses combined site intensity (symmetric treatment)
- Edge-centric heat h_E(g,r) uses asymmetric source/target intensities
- The ratio (panel d) shows where edge-centric over/under-weights relative to node-centric
- Edge-centric concentrates at (predator g, prey r) = (0.7, 0.7) region
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
        data = load_figure_A_prime_data()
        plot_figure_A_prime(data)
    else
        data = compute_figure_A_prime_data()
        plot_figure_A_prime(data)
    end
    return 0
end
