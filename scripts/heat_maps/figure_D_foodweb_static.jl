# Figure D: Food Web from Heat Map (Static)
#
# Dimension: d = 4
# Purpose: Show how IDPG generates realistic trophic networks
#
# Uses guild centroids from Phase 1, static intensity
#
# Panels (2×2 grid):
#   (a) Target affinity K*_ij (heatmap, 5×5)
#   (b) Achieved affinity K̂_ij = μ_i^G · μ_j^R (heatmap, 5×5)
#   (c) Expected edges H_ij = Λ² π_i π_j K̂_ij (heatmap, 5×5)
#   (d) Realized food web graph (trophic layout)
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_D_foodweb_static.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_D_foodweb_static.jl --plot   # plot from saved data

using IDPG
using CairoMakie
using LaTeXStrings
using Distributions
using Statistics
using LinearAlgebra
using JLD2
using Random
using Graphs
using GraphMakie

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "heat_maps")
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_D_data.jld2")
const CENTROIDS_FILE = joinpath(pkgdir(IDPG), "output", "heat_maps", "foodweb_centroids.jld2")

# Sampling parameters
const LAMBDA = 100.0          # Total intensity
const KAPPA = 30.0            # Concentration
const N_SAMPLES = 5           # Samples to average for realized graph

Random.seed!(42)

# ============================================================================
# Load Centroids
# ============================================================================

"""
    load_centroids()

Load guild centroids from Phase 1 output.
"""
function load_centroids()
    @load CENTROIDS_FILE result
    return result
end

# ============================================================================
# Graph Generation
# ============================================================================

# Note: Uses sample_guild_position and project_to_Bd_plus from IDPG library

"""
    generate_idpg_graph(M_G, M_R, π, Λ, κ)

Generate a graph from IDPG with given guild centroids.
Returns adjacency matrix at guild level (aggregated).
"""
function generate_idpg_graph(M_G::Matrix{Float64}, M_R::Matrix{Float64},
                              π::Vector{Float64}, Λ::Float64, κ::Float64)
    n_guilds = size(M_G, 1)

    # Sample number of entities
    N = rand(Poisson(Λ))

    if N == 0
        return zeros(Int, n_guilds, n_guilds)
    end

    # Assign each entity to a guild
    guild_counts = zeros(Int, n_guilds)
    guild_assignments = zeros(Int, N)
    for i in 1:N
        g = rand(Distributions.Categorical(π))
        guild_assignments[i] = g
        guild_counts[g] += 1
    end

    # Sample positions for each entity
    g_positions = Vector{Vector{Float64}}(undef, N)
    r_positions = Vector{Vector{Float64}}(undef, N)

    for i in 1:N
        g_idx = guild_assignments[i]
        g_positions[i] = sample_guild_position(M_G[g_idx, :], κ)
        r_positions[i] = sample_guild_position(M_R[g_idx, :], κ)
    end

    # Generate edges
    edge_counts = zeros(Int, n_guilds, n_guilds)

    for i in 1:N
        for j in 1:N
            if i != j  # No self-loops
                # Probability of edge from i to j: g_i · r_j
                p = dot(g_positions[i], r_positions[j])
                if rand() < p
                    g_i = guild_assignments[i]
                    g_j = guild_assignments[j]
                    edge_counts[g_i, g_j] += 1
                end
            end
        end
    end

    return edge_counts
end

# ============================================================================
# Theoretical Quantities
# ============================================================================

# Note: Uses compute_expected_guild_edges and trophic_layout from IDPG library

# ============================================================================
# Data Computation
# ============================================================================

"""
    compute_figure_D_data()

Compute all data for Figure D and save to JLD2.
"""
function compute_figure_D_data()
    println("=" ^ 60)
    println("Figure D: Computing Food Web Data (Static)")
    println("=" ^ 60)

    # Load centroids
    println("Loading centroids from Phase 1...")
    result = load_centroids()
    M_G = result.M_G
    M_R = result.M_R
    K_star = result.K_star
    K_approx = result.K_approx
    guild_names = result.guild_names

    n_guilds = size(M_G, 1)

    println("  Guilds: " * join(guild_names, ", "))

    # Guild weights (equal)
    π = fill(1.0 / n_guilds, n_guilds)

    # Theoretical expected edges (using library function)
    H_expected = compute_expected_guild_edges(K_approx, π, LAMBDA)

    println("\nTheoretical quantities:")
    println("  Λ = " * string(LAMBDA))
    println("  Total expected edges: " * string(round(sum(H_expected), digits=1)))

    # Sample realized graphs
    println("\nSampling " * string(N_SAMPLES) * " realized graphs...")
    total_edges = zeros(n_guilds, n_guilds)
    for s in 1:N_SAMPLES
        edges = generate_idpg_graph(M_G, M_R, π, LAMBDA, KAPPA)
        total_edges .+= edges
    end
    avg_edges = total_edges ./ N_SAMPLES

    println("  Average realized edges: " * string(round(sum(avg_edges), digits=1)))

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        K_star = K_star,
        K_approx = K_approx,
        H_expected = H_expected,
        avg_edges = avg_edges,
        guild_names = guild_names,
        n_guilds = n_guilds,
        LAMBDA = LAMBDA,
        N_SAMPLES = N_SAMPLES,
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

"""
    load_figure_D_data()

Load precomputed data for Figure D.
"""
function load_figure_D_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

"""
    plot_figure_D(data)

Generate Figure D visualization from precomputed data.
"""
function plot_figure_D(data)
    println("=" ^ 60)
    println("Figure D: Generating Plot")
    println("=" ^ 60)

    K_star = data.K_star
    K_approx = data.K_approx
    H_expected = data.H_expected
    avg_edges = data.avg_edges
    guild_names = data.guild_names
    n_guilds = data.n_guilds
    LAMBDA = data.LAMBDA
    N_SAMPLES = data.N_SAMPLES

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
    # Create figure
    fig = Figure(size=(1100, 1000))

    # Color scale for heatmaps
    colormap = :YlOrRd

    # Panel (a): Target affinity K*
    ax_a = Axis(fig[1, 1],
        xlabel = "Consumer (predator)",
        ylabel = "Resource (prey)",
        title = L"(a) Target Affinity $K^*$",
        aspect = 1,
        xticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]),
        yticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]))

    hm_a = heatmap!(ax_a, 1:n_guilds, 1:n_guilds, K_star, colormap=colormap)
    Colorbar(fig[1, 2], hm_a, label=L"$K^*_{ij}$")

    # Add text annotations
    for i in 1:n_guilds, j in 1:n_guilds
        val = K_star[i, j]
        txt_color = val > 0.5 ? :white : :black
        text!(ax_a, i, j, text=string(round(val, digits=2)), fontsize=10, color=txt_color, align=(:center, :center))
    end

    # Panel (b): Achieved affinity K̂
    ax_b = Axis(fig[1, 3],
        xlabel = "Consumer (predator)",
        ylabel = "Resource (prey)",
        title = L"(b) Achieved Affinity $\hat{K} = \mu^G \cdot \mu^R$",
        aspect = 1,
        xticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]),
        yticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]))

    hm_b = heatmap!(ax_b, 1:n_guilds, 1:n_guilds, K_approx, colormap=colormap)
    Colorbar(fig[1, 4], hm_b, label=L"$\hat{K}_{ij}$")

    for i in 1:n_guilds, j in 1:n_guilds
        val = K_approx[i, j]
        txt_color = val > 0.5 ? :white : :black
        text!(ax_b, i, j, text=string(round(val, digits=2)), fontsize=10, color=txt_color, align=(:center, :center))
    end

    # Panel (c): Expected edges H
    ax_c = Axis(fig[2, 1],
        xlabel = "Consumer (predator)",
        ylabel = "Resource (prey)",
        title = L"(c) Expected Edges $H = \Lambda^2 \pi_i \pi_j \hat{K}$",
        aspect = 1,
        xticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]),
        yticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]))

    hm_c = heatmap!(ax_c, 1:n_guilds, 1:n_guilds, H_expected, colormap=colormap)
    Colorbar(fig[2, 2], hm_c, label=L"$H_{ij}$")

    for i in 1:n_guilds, j in 1:n_guilds
        val = H_expected[i, j]
        txt = val < 10 ? string(round(val, digits=1)) : string(round(Int, val))
        txt_color = val > maximum(H_expected) / 2 ? :white : :black
        text!(ax_c, i, j, text=txt, fontsize=10, color=txt_color, align=(:center, :center))
    end

    # Panel (d): Realized food web graph using GraphMakie
    ax_d = Axis(fig[2, 3:4],
        title = "(d) Realized Food Web (averaged over " * string(N_SAMPLES) * " samples)",
        aspect = 1.2)
    hidedecorations!(ax_d)
    hidespines!(ax_d)

    # Build directed graph from edge matrix
    G = SimpleDiGraph(n_guilds)
    edge_weights = Float64[]
    max_edge = maximum(avg_edges)

    for i in 1:n_guilds
        for j in 1:n_guilds
            if avg_edges[i, j] > 0
                add_edge!(G, i, j)
                push!(edge_weights, avg_edges[i, j] / max_edge)
            end
        end
    end

    # Trophic layout positions (using library function)
    x_positions, trophic_levels = trophic_layout(n_guilds)
    layout_pos = [Point2f(x_positions[i], trophic_levels[i]) for i in 1:n_guilds]

    # Edge styling (steeper decay for less important edges)
    edge_widths = [0.5 + 4.0 * w^2 for w in edge_weights]
    edge_colors = [RGBAf(0, 0, 0, 0.15 + 0.85 * w^2) for w in edge_weights]

    # Node styling
    guild_colors = [:green, :orange, :purple, :red, :black]

    # Plot graph with curved edges
    graphplot!(ax_d, G,
        layout = _ -> layout_pos,
        node_color = guild_colors,
        node_size = 30,
        node_strokewidth = 2,
        node_strokecolor = :white,
        edge_color = edge_colors,
        edge_width = edge_widths,
        arrow_show = true,
        arrow_size = 15,
        arrow_shift = 0.8,
        curve_distance = 0.08,
        nlabels = guild_names,
        nlabels_align = (:center, :bottom),
        nlabels_distance = 15,
        nlabels_fontsize = 10)

    # Add y-axis labels for trophic levels
    text!(ax_d, -0.1, 0.0, text="TL 0", fontsize=9, align=(:right, :center))
    text!(ax_d, -0.1, 1.0, text="TL 1", fontsize=9, align=(:right, :center))
    text!(ax_d, -0.1, 2.0, text="TL 2", fontsize=9, align=(:right, :center))
    text!(ax_d, -0.1, 3.0, text="TL 3", fontsize=9, align=(:right, :center))

    xlims!(ax_d, -0.2, 1.2)
    ylims!(ax_d, -0.3, 3.5)

    # Add overall title
    lambda_str = string(Int(LAMBDA))
    Label(fig[0, 1:4], "Food Web from Heat Map (Static, Λ=" * lambda_str * ")", fontsize=22)

    fig
    end  # end with_theme

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_D_foodweb_static.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_D_data()
        plot_figure_D(data)
    else
        data = compute_figure_D_data()
        plot_figure_D(data)
    end
    return 0
end
