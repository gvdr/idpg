# Figure E: Dynamic Food Web (Pursuit-Evasion)
#
# Dimension: d = 4
# Purpose: Show how trophic structure evolves under pursuit-evasion dynamics
#
# Uses guild centroids from Phase 1, pursuit-evasion PDE with elastic centering
#
# Layout: Time series
#   Top row: H_ij(t) heatmaps at t = 0, T/2, T
#   Bottom row: Realized food webs at same times
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_E_foodweb_dynamic.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_E_foodweb_dynamic.jl --plot   # plot from saved data

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
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_E_data.jld2")
const CENTROIDS_FILE = joinpath(pkgdir(IDPG), "output", "heat_maps", "foodweb_centroids.jld2")

# PDE and sampling parameters
const LAMBDA = 100.0          # Total intensity
const KAPPA = 30.0            # Concentration
const T_FINAL = 1.0           # Final time
const DT = 0.01               # Time step
const N_SAMPLES = 3           # Samples per time point

# Pursuit-evasion parameters
const ALPHA = 0.3             # Prey evasion speed
const BETA = 0.3              # Predator pursuit speed
const GAMMA = 0.15            # Elastic centering strength

# Center point for elastic centering (interior of B^d_+)
const CENTER = 0.4

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
# Pursuit-Evasion Dynamics in 4D
# ============================================================================

"""
    project_to_Bd_plus(x)

Project point to B^d_+ (positive orthant of unit ball).
"""
function project_to_Bd_plus(x::Vector{Float64})
    x_pos = max.(x, 0.0)
    n = norm(x_pos)
    if n > 0.99  # Stay slightly interior
        x_pos .*= 0.99 / n
    end
    return x_pos
end

"""
    compute_mean_position(M, weights)

Compute weighted mean position of guilds.
"""
function compute_mean_position(M::Matrix{Float64}, weights::Vector{Float64})
    n_guilds = size(M, 1)
    d = size(M, 2)
    μ = zeros(d)
    for g in 1:n_guilds
        μ .+= weights[g] .* M[g, :]
    end
    return μ ./ sum(weights)
end

"""
    evolve_centroids_pursuit_evasion!(M_G, M_R, dt, α, β, γ, center)

One time step of pursuit-evasion dynamics for guild centroids.
Prey (G) flee from predators, predators (R) chase prey, both attracted to center.
"""
function evolve_centroids_pursuit_evasion!(M_G::Matrix{Float64}, M_R::Matrix{Float64},
                                           weights::Vector{Float64},
                                           dt::Float64, α::Float64, β::Float64, γ::Float64, center::Float64)
    n_guilds = size(M_G, 1)
    d = size(M_G, 2)

    # Mean positions
    μ_G = compute_mean_position(M_G, weights)
    μ_R = compute_mean_position(M_R, weights)

    # Center vector
    c = fill(center, d)

    # Update each guild centroid
    for g in 1:n_guilds
        # Prey velocity: flee from predator mean + centering
        flee_dir = M_G[g, :] .- μ_R  # Direction away from predators
        center_pull = c .- M_G[g, :]  # Direction toward center

        v_G = α .* flee_dir .+ γ .* center_pull
        M_G[g, :] .+= dt .* v_G

        # Predator velocity: chase prey mean + centering
        chase_dir = μ_G .- M_R[g, :]  # Direction toward prey
        center_pull_R = c .- M_R[g, :]

        v_R = β .* chase_dir .+ γ .* center_pull_R
        M_R[g, :] .+= dt .* v_R

        # Project back to B^d_+
        M_G[g, :] .= project_to_Bd_plus(M_G[g, :])
        M_R[g, :] .= project_to_Bd_plus(M_R[g, :])
    end
end

# ============================================================================
# Graph Generation
# ============================================================================

# Note: Uses sample_guild_position from IDPG library

"""
    generate_idpg_graph(M_G, M_R, π, Λ, κ)

Generate a graph from IDPG with given guild centroids.
Returns adjacency matrix at guild level.
"""
function generate_idpg_graph(M_G::Matrix{Float64}, M_R::Matrix{Float64},
                              π::Vector{Float64}, Λ::Float64, κ::Float64)
    n_guilds = size(M_G, 1)

    # Sample number of entities
    N = rand(Distributions.Poisson(Λ))

    if N == 0
        return zeros(Int, n_guilds, n_guilds)
    end

    # Assign each entity to a guild
    guild_assignments = zeros(Int, N)
    for i in 1:N
        guild_assignments[i] = rand(Distributions.Categorical(π))
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
            if i != j
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

# Note: Uses compute_guild_affinity and compute_expected_guild_edges from IDPG library

# ============================================================================
# Data Computation
# ============================================================================

"""
    compute_figure_E_data()

Compute all data for Figure E and save to JLD2.
"""
function compute_figure_E_data()
    println("=" ^ 60)
    println("Figure E: Computing Dynamic Food Web Data")
    println("=" ^ 60)

    # Load initial centroids
    println("Loading centroids from Phase 1...")
    result = load_centroids()
    M_G_init = copy(result.M_G)
    M_R_init = copy(result.M_R)
    guild_names = result.guild_names

    n_guilds = size(M_G_init, 1)
    d = size(M_G_init, 2)

    println("  Guilds: " * join(guild_names, ", "))
    println("  PDE parameters: α=" * string(ALPHA) * ", β=" * string(BETA) * ", γ=" * string(GAMMA))

    # Guild weights (equal)
    π = fill(1.0 / n_guilds, n_guilds)

    # Time snapshots
    snapshot_times = [0.0, T_FINAL / 2, T_FINAL]
    n_snapshots = length(snapshot_times)

    # Storage for results at each time
    K_hats = Dict{Float64, Matrix{Float64}}()
    H_expecteds = Dict{Float64, Matrix{Float64}}()
    avg_edges_all = Dict{Float64, Matrix{Float64}}()
    M_Gs = Dict{Float64, Matrix{Float64}}()
    M_Rs = Dict{Float64, Matrix{Float64}}()

    # Evolution
    M_G = copy(M_G_init)
    M_R = copy(M_R_init)
    current_time = 0.0
    n_steps = Int(round(T_FINAL / DT))

    # Store initial state
    K_hats[0.0] = compute_guild_affinity(M_G, M_R)
    H_expecteds[0.0] = compute_expected_guild_edges(K_hats[0.0], π, LAMBDA)
    M_Gs[0.0] = copy(M_G)
    M_Rs[0.0] = copy(M_R)

    # Sample initial graphs
    total_edges = zeros(n_guilds, n_guilds)
    for _ in 1:N_SAMPLES
        total_edges .+= generate_idpg_graph(M_G, M_R, π, LAMBDA, KAPPA)
    end
    avg_edges_all[0.0] = total_edges ./ N_SAMPLES

    println("\nRunning pursuit-evasion dynamics...")
    println("  t=0.0: Total expected edges = " * string(round(sum(H_expecteds[0.0]), digits=1)))

    for step in 1:n_steps
        evolve_centroids_pursuit_evasion!(M_G, M_R, π, DT, ALPHA, BETA, GAMMA, CENTER)
        current_time += DT

        # Check for snapshots
        for t_snap in snapshot_times[2:end]
            if abs(current_time - t_snap) < DT / 2 && !haskey(K_hats, t_snap)
                K_hats[t_snap] = compute_guild_affinity(M_G, M_R)
                H_expecteds[t_snap] = compute_expected_guild_edges(K_hats[t_snap], π, LAMBDA)
                M_Gs[t_snap] = copy(M_G)
                M_Rs[t_snap] = copy(M_R)

                # Sample graphs
                total_edges = zeros(n_guilds, n_guilds)
                for _ in 1:N_SAMPLES
                    total_edges .+= generate_idpg_graph(M_G, M_R, π, LAMBDA, KAPPA)
                end
                avg_edges_all[t_snap] = total_edges ./ N_SAMPLES

                println("  t=" * string(round(t_snap, digits=2)) * ": Total expected edges = " * string(round(sum(H_expecteds[t_snap]), digits=1)))
            end
        end
    end

    # Compute metrics
    println("\nMetrics over time:")
    for t in snapshot_times
        total_H = sum(H_expecteds[t])
        K_hat = K_hats[t]
        σ1 = svd(K_hat).S[1]
        println("  t=" * string(round(t, digits=2)) * ": Total H = " * string(round(total_H, digits=1)) * ", σ₁(K̂) = " * string(round(σ1, digits=4)))
    end

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        snapshot_times = snapshot_times,
        K_hats = K_hats,
        H_expecteds = H_expecteds,
        avg_edges_all = avg_edges_all,
        guild_names = guild_names,
        n_guilds = n_guilds,
        N_SAMPLES = N_SAMPLES,
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

"""
    load_figure_E_data()

Load precomputed data for Figure E.
"""
function load_figure_E_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

"""
    plot_figure_E(data)

Generate Figure E visualization from precomputed data.
"""
function plot_figure_E(data)
    println("=" ^ 60)
    println("Figure E: Generating Plot")
    println("=" ^ 60)

    snapshot_times = data.snapshot_times
    H_expecteds = data.H_expecteds
    avg_edges_all = data.avg_edges_all
    guild_names = data.guild_names
    n_guilds = data.n_guilds
    N_SAMPLES = data.N_SAMPLES

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
    # Create figure: 2 rows × 3 columns
    fig = Figure(size=(1200, 800))

    colormap = :YlOrRd

    # Find global color limits
    all_H = vcat([vec(H_expecteds[t]) for t in snapshot_times]...)
    clims_H = (0, maximum(all_H) * 1.1)

    guild_colors = [:green, :orange, :purple, :red, :black]

    for (col, t) in enumerate(snapshot_times)
        # Top row: Expected edges heatmap
        ax_top = Axis(fig[1, col],
            xlabel = "Consumer",
            ylabel = col == 1 ? "Resource" : "",
            title = "t = " * string(round(t, digits=2)),
            aspect = 1,
            xticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]),
            yticks = (1:n_guilds, ["P", "SH", "LH", "SP", "A"]))

        hm = heatmap!(ax_top, 1:n_guilds, 1:n_guilds, H_expecteds[t], colormap=colormap, colorrange=clims_H)

        # Add text annotations
        for i in 1:n_guilds, j in 1:n_guilds
            val = H_expecteds[t][i, j]
            txt = val < 10 ? string(round(val, digits=1)) : string(round(Int, val))
            txt_color = val > clims_H[2] / 2 ? :white : :black
            text!(ax_top, i, j, text=txt, fontsize=8, color=txt_color, align=(:center, :center))
        end

        if col == 3
            Colorbar(fig[1, 4], hm, label=L"Expected edges $H_{ij}$")
        end

        # Bottom row: Realized food web using GraphMakie
        ax_bot = Axis(fig[2, col],
            title = "Realized web (n=" * string(N_SAMPLES) * ")",
            aspect = 1.2)
        hidedecorations!(ax_bot)
        hidespines!(ax_bot)

        # Build directed graph from edge matrix
        G = SimpleDiGraph(n_guilds)
        edge_weights = Float64[]
        avg_edges = avg_edges_all[t]
        max_edge = maximum(avg_edges) + 1e-10

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

        # Plot graph with curved edges
        graphplot!(ax_bot, G,
            layout = _ -> layout_pos,
            node_color = guild_colors,
            node_size = 25,
            node_strokewidth = 2,
            node_strokecolor = :white,
            edge_color = edge_colors,
            edge_width = edge_widths,
            arrow_show = true,
            arrow_size = 12,
            arrow_shift = 0.8,
            curve_distance = 0.08)

        xlims!(ax_bot, -0.1, 1.1)
        ylims!(ax_bot, -0.3, 3.3)
    end

    # Add horizontal legend below graphs
    ax_legend = Axis(fig[3, 1:3], height=40)
    hidedecorations!(ax_legend)
    hidespines!(ax_legend)

    for (i, name) in enumerate(guild_names)
        x_pos = (i - 1) / (n_guilds - 1)  # Spread across [0, 1]
        scatter!(ax_legend, [x_pos], [0.5],
                color=guild_colors[i], markersize=15)
        text!(ax_legend, x_pos + 0.03, 0.5,
              text=name, fontsize=9, align=(:left, :center))
    end
    xlims!(ax_legend, -0.05, 1.15)
    ylims!(ax_legend, 0.0, 1.0)

    # Add overall title
    Label(fig[0, 1:3], "Dynamic Food Web (Pursuit-Evasion)", fontsize=22)

    fig
    end  # end with_theme

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_E_foodweb_dynamic.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_E_data()
        plot_figure_E(data)
    else
        data = compute_figure_E_data()
        plot_figure_E(data)
    end
    return 0
end
