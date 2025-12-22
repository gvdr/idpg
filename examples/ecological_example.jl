# Ecological Example: Food Web as Edge-Centric IDPG
# Demonstrates how trophic interactions in an ecosystem can be modeled
# as an edge-centric IDPG, following Section 5 of the paper.
#
# Key insight: A food web edge (predator eats prey) is an ephemeral interaction
# between individuals, not a permanent relationship between species.

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using GraphMakie
using NetworkLayout

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Ecological Example: Food Web Modeling")
println("=" ^ 60)

# In this model:
# - g represents the "resource role" - an individual's position when being consumed
# - r represents the "consumer role" - an individual's position when consuming
#
# Under the edge-centric interpretation, each sampled point (g, r) represents
# a trophic interaction: a consumer "at r" eating a resource "at g".
#
# Note: In edge-centric, g and r come from DIFFERENT individuals in each interaction.
# The same individual has both g and r positions, but we only see them in one role per edge.

# Define scenarios with different kernel concentrations
# Higher concentration = more specialized niches
# Lower concentration = more generalist niches
# Scales are adjusted to keep E[L] roughly constant across scenarios

# Common structure across scenarios
weights_G = [0.7, 0.2, 0.1]  # Most biomass at base
means_G = [[0.8, 0.2], [0.5, 0.5], [0.3, 0.7]]

weights_R = [0.3, 0.4, 0.3]  # More even consumer distribution
means_R = [[0.7, 0.3], [0.5, 0.5], [0.4, 0.8]]

scenarios = [
    (
        name = "Spread (generalist)",
        κ_G = [5.0, 5.0, 5.0],
        κ_R = [5.0, 5.0, 5.0],
        scale_G = 50.0,
        scale_R = 30.0,
    ),
    (
        name = "Medium",
        κ_G = [40.0, 35.0, 35.0],
        κ_R = [35.0, 30.0, 40.0],
        scale_G = 195.0,
        scale_R = 118.0,
    ),
    (
        name = "Concentrated (specialist)",
        κ_G = [120.0, 100.0, 100.0],
        κ_R = [100.0, 90.0, 120.0],
        scale_G = 560.0,
        scale_R = 320.0,
    ),
]

# Create comparison figure (larger for better visibility)
fig = Figure(size=(1800, 1800))

println("\n--- Comparing Scenarios ---")

# Store data for network plotting
scenario_data = []

for (col, scenario) in enumerate(scenarios)
    println("\n" * "=" ^ 50)
    println("Scenario: ", scenario.name)
    println("=" ^ 50)

    # Build intensities for this scenario (using per-scenario scales)
    ρ_G = BdPlusMixture(weights_G, means_G, scenario.κ_G, scenario.scale_G)
    ρ_R = BdPlusMixture(weights_R, means_R, scenario.κ_R, scenario.scale_R)
    ρ = ProductIntensity(ρ_G, ρ_R)

    # Compute statistics
    stats = marginal_stats(ρ; rng=rng)
    println("  Expected opportunities E[N]: ", round(stats.E_N, digits=1))
    println("  Avg connection prob: ", round(stats.avg_conn_prob, digits=3))
    println("  Expected interactions E[L]: ", round(stats.E_edges_edge_centric, digits=1))

    # Sample interactions
    sites = sample_ppp_product(ρ; rng=MersenneTwister(42 + col))
    trophic_interactions = generate_edge_centric(sites; rng=MersenneTwister(42 + col))
    println("  Sampled opportunities: ", length(sites))
    println("  Realized interactions: ", length(trophic_interactions))

    # Joint clustering with DBSCAN (adaptive eps based on concentration)
    # Higher concentration → tighter clusters → smaller eps needed
    avg_kappa = mean(scenario.κ_G)
    eps_val = 0.25 / sqrt(avg_kappa / 10)  # Scale eps with concentration
    eps_val = clamp(eps_val, 0.05, 0.3)

    _, _, joint_assign, n_joint = discretize_edge_centric_joint(
        trophic_interactions; eps=eps_val, min_samples=3)
    n_noise = count(==(0), joint_assign)
    println("  Joint DBSCAN (eps=", round(eps_val, digits=3), "): ", n_joint, " clusters, ", n_noise, " noise")

    # Use DBSCAN result to guide k-means clusters, with minimum of 3
    n_clusters = max(3, n_joint)
    graph, edge_weights, src_assign, tgt_assign = discretize_with_weights(
        trophic_interactions, n_clusters; rng=MersenneTwister(42 + col))
    println("  K-means clusters: ", n_clusters, ", Links: ", ne(graph))
    if any(edge_weights .> 0)
        println("  Edge weight range: ", minimum(edge_weights[edge_weights .> 0]), "-", maximum(edge_weights))
    end

    # Statistical significance test for links
    # H0: interactions are randomly distributed (no cluster preference)
    # Expected: E[i,j] = n_source_i × n_target_j / N
    # where n_source_i = count of interactions with source in cluster i
    # Test: keep edge if observed >> expected (using binomial or ratio threshold)

    N = length(trophic_interactions)
    n_source = [count(==(k), src_assign) for k in 1:n_clusters]
    n_target = [count(==(k), tgt_assign) for k in 1:n_clusters]

    expected_weights = zeros(n_clusters, n_clusters)
    significance = zeros(n_clusters, n_clusters)  # O/E ratio

    for i in 1:n_clusters
        for j in 1:n_clusters
            expected_weights[i, j] = n_source[i] * n_target[j] / N
            if expected_weights[i, j] > 0
                significance[i, j] = edge_weights[i, j] / expected_weights[i, j]
            end
        end
    end

    # Count non-self-loop edges
    n_links = 0
    for i in 1:n_clusters
        for j in 1:n_clusters
            if i != j && edge_weights[i, j] > 0
                n_links += 1
            end
        end
    end
    println("  Links (excl. self-loops): ", n_links)

    push!(scenario_data, (ρ_G=ρ_G, ρ_R=ρ_R, interactions=trophic_interactions,
                          graph=graph, edge_weights=edge_weights, stats=stats,
                          n_clusters=n_clusters, n_links=n_links,
                          src_assign=src_assign, tgt_assign=tgt_assign))

    # Row 1: Resource intensity
    ax1 = Axis(fig[1, col], aspect=DataAspect(),
        title = scenario.name * "\nρ_G (resource intensity)")
    plot_intensity_Bd_plus!(ax1, ρ_G; resolution=40)

    # Row 2: Consumer intensity
    ax2 = Axis(fig[2, col], aspect=DataAspect(),
        title = "ρ_R (consumer intensity)")
    plot_intensity_Bd_plus!(ax2, ρ_R; resolution=40)

    # Row 3: Realized interactions - resource role (g) colored by source cluster
    ax3 = Axis(fig[3, col], aspect=DataAspect(),
        title = "g positions (colored by cluster)")
    draw_Bd_plus_boundary!(ax3)
    if !isempty(trophic_interactions.sources)
        source_2d = [Bd_plus_to_2d(s) for s in trophic_interactions.sources]
        # Color by source cluster assignment
        colors = [Makie.wong_colors()[mod1(c, 7)] for c in src_assign]
        scatter!(ax3, source_2d, color=colors, markersize=7, alpha=0.7)
    end

    # Row 4: Realized interactions - consumer role (r) colored by target cluster
    ax4 = Axis(fig[4, col], aspect=DataAspect(),
        title = "r positions (colored by cluster)")
    draw_Bd_plus_boundary!(ax4)
    if !isempty(trophic_interactions.targets)
        target_2d = [Bd_plus_to_2d(t) for t in trophic_interactions.targets]
        # Color by target cluster assignment
        colors = [Makie.wong_colors()[mod1(c, 7)] for c in tgt_assign]
        scatter!(ax4, target_2d, color=colors, markersize=7, alpha=0.7)
    end

    # Row 5: Aggregated food web (excluding self-loops)
    ax5 = Axis(fig[5, col], aspect=DataAspect(),
        title = "Food web (" * string(n_links) * " links)")
    hidedecorations!(ax5)
    hidespines!(ax5)
    limits!(ax5, -1.5, 1.5, -1.5, 1.5)

    if n_clusters > 0 && n_links > 0
        # Build graph excluding self-loops
        food_web = SimpleDiGraph(n_clusters)
        edge_widths = Float64[]
        arrow_sizes = Float64[]
        max_weight = maximum(edge_weights)

        for i in 1:n_clusters
            for j in 1:n_clusters
                if i != j && edge_weights[i, j] > 0
                    add_edge!(food_web, j, i)  # Arrow: consumer j → resource i (who eats whom)
                    # Edge width proportional to interaction count
                    w_normalized = edge_weights[i, j] / max_weight
                    width = 0.1 + 8.0 * w_normalized^0.6  # thinner base, steeper scaling
                    push!(edge_widths, width)
                    push!(arrow_sizes, 3.0 + 22.0 * w_normalized^0.6)  # scale arrows with edges
                end
            end
        end

        # Node colors by cluster index
        node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_clusters]

        graphplot!(ax5, food_web,
            layout=Shell(),
            node_size=30,
            node_color=node_colors,
            nlabels=string.(1:n_clusters),
            nlabels_align=(:center, :center),
            nlabels_color=:white,
            nlabels_fontsize=12,
            edge_color=(:black, 0.8),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    end
end

save("output/applications/ecological_foodweb.png", fig)
println("\nSaved ecological_foodweb.png")

# --- Monte Carlo comparison across scenarios ---
println("\n" * "=" ^ 60)
println("Monte Carlo Comparison")
println("=" ^ 60)

n_trials = 100

fig2 = Figure(size=(1400, 450))

for (idx, scenario) in enumerate(scenarios)
    # Use per-scenario scales
    ρ_G = BdPlusMixture(weights_G, means_G, scenario.κ_G, scenario.scale_G)
    ρ_R = BdPlusMixture(weights_R, means_R, scenario.κ_R, scenario.scale_R)
    ρ = ProductIntensity(ρ_G, ρ_R)
    stats = marginal_stats(ρ; rng=rng)

    L_samples = Int[]
    for t in 1:n_trials
        sites_trial = sample_ppp_product(ρ; rng=MersenneTwister(1000*idx + t))
        interactions_trial = generate_edge_centric(sites_trial; rng=MersenneTwister(1000*idx + t))
        push!(L_samples, length(interactions_trial))
    end

    println("\n", scenario.name, ":")
    println("  E[L] theory:    ", round(stats.E_edges_edge_centric, digits=1))
    println("  E[L] empirical: ", round(mean(L_samples), digits=1), " ± ", round(std(L_samples), digits=1))

    ax = Axis(fig2[1, idx], xlabel="Number of Interactions (L)", ylabel="Count",
        title=scenario.name)
    hist!(ax, L_samples, bins=20, color=[:steelblue, :forestgreen, :purple][idx])
    vlines!(ax, [stats.E_edges_edge_centric], color=:red, linestyle=:dash, linewidth=2, label="E[L]")
    if idx == 1
        axislegend(ax, position=:rt)
    end
end

save("output/applications/ecological_mc_comparison.png", fig2)
println("\nSaved ecological_mc_comparison.png")

# --- Ecological interpretation ---
println("\n" * "=" ^ 60)
println("Ecological Interpretation:")
println("=" ^ 60)
println("""
In this food web model:

1. The intensity ρ(g, r) = ρ_G(g) · ρ_R(r) represents the 'opportunity structure'
   for trophic interactions - which resource-consumer pairs have the chance to meet.

2. The connection probability g · r represents whether an encounter results in
   predation, based on the alignment of resource and consumer positions.

3. The edge-centric interpretation is natural because each predation event involves
   ephemeral individuals - the consumed individual ceases to exist after the interaction.

4. Kernel concentration (κ) controls niche specialization:
   - Low κ (spread): generalist niches, individuals sample broadly in trait space
   - High κ (concentrated): specialist niches, tight clustering in trait space

5. Same cluster means across scenarios allows direct comparison of how
   concentration affects the realized interaction patterns.

6. Species-level food webs (as typically studied) are statistical summaries:
   they cluster individual interactions into species nodes and count link frequencies.
""")

println("\nDone!")
