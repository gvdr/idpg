# Ecological Example: Food Web in 4D Latent Space
# Demonstrates how higher-dimensional latent space allows for more
# distinct ecological niches and clearer trophic structure.
#
# With 4D instead of 2D, species have more "corners" to occupy,
# enabling finer niche differentiation and more realistic food web structure.
#
# Using hyperspherical coordinates (r, φ₁, φ₂, φ₃) to define guild positions:
# - r: how "active" in the network (0=invisible, 1=maximally present)
# - φ₁, φ₂, φ₃: angular positions defining niche location (0 to π/2)

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using GraphMakie
using NetworkLayout
using LinearAlgebra

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Ecological Example: 4D Food Web Modeling")
println("=" ^ 60)

# In 4D hyperspherical coordinates, we specify:
#   r ∈ [0, 1]: radius (network activity level)
#   φ₁, φ₂, φ₃ ∈ [0, π/2]: angles defining niche position
#
# Key insight for trophic structure:
#   - Small angle (≈0) concentrates weight in current dimension
#   - Large angle (≈π/2) passes weight to later dimensions
#
# To point to dimension k, use (k-1) large angles followed by small angles:
#   dim 1: [small, small, small]  → producers (food source)
#   dim 2: [large, small, small]  → herbivores
#   dim 3: [large, large, small]  → small predators
#   dim 4: [large, large, large]  → apex predators

guild_names = ["Producers", "Small Herb.", "Large Herb.", "Small Pred.", "Apex Pred."]

# Angle constants for clean dimensional separation
# More extreme angles for better orthogonality between trophic levels
const small_angle = π/60   # ≈ 3° - strongly concentrates weight in current dimension
const large_angle = 29π/60 # ≈ 87° - passes almost all weight to later dimensions

# Resource roles (g) - where species sit when being consumed
# Each trophic level peaks in its own dimension with minimal overlap
g_hyperspherical = [
    (0.95, [small_angle, small_angle, small_angle]),  # Producers: peak in dim 1, highly visible
    (0.90, [large_angle, small_angle, small_angle]),  # Small herbivores: peak in dim 2
    (0.85, [large_angle, small_angle, small_angle]),  # Large herbivores: peak in dim 2 (similar niche)
    (0.75, [large_angle, large_angle, small_angle]),  # Small predators: peak in dim 3
    (0.20, [large_angle, large_angle, large_angle]),  # Apex predators: peak in dim 4, LOW visibility
]

# Consumer roles (r) - each level targets the level below
# More focused targeting to reduce cross-trophic noise
r_hyperspherical = [
    (0.02, [small_angle, small_angle, small_angle]),   # Producers: nearly zero consumption
    (0.95, [small_angle, small_angle, small_angle]),   # Small herbivores: target dim 1 (producers)
    (0.90, [small_angle, small_angle, small_angle]),   # Large herbivores: target dim 1 (producers)
    (0.90, [large_angle, small_angle, small_angle]),   # Small predators: target dim 2 (herbivores)
    (0.95, [large_angle, large_angle, small_angle]),   # Apex: target dim 3 (small predators)
]

# Convert to Cartesian coordinates
means_G = [Vector(Bd_plus_from_hyperspherical(r, angles)) for (r, angles) in g_hyperspherical]
means_R = [Vector(Bd_plus_from_hyperspherical(r, angles)) for (r, angles) in r_hyperspherical]

# Weights for mixture components (abundance distribution)
weights_G = [0.35, 0.25, 0.20, 0.15, 0.05]  # Most biomass at base
weights_R = [0.05, 0.25, 0.20, 0.30, 0.20]  # Consumers more evenly distributed

println("\nGuild structure (4D positions):")
for i in 1:5
    println("  ", guild_names[i], ":")
    println("    g = ", round.(means_G[i], digits=2))
    println("    r = ", round.(means_R[i], digits=2))
    g_dot_r = dot(means_G[i], means_R[i])
    println("    self-interaction (g·r) = ", round(g_dot_r, digits=3))
end

# Print expected interactions between guilds
println("\nExpected interaction strengths (g_i · r_j):")
print("           ")
for j in 1:5
    print(rpad(guild_names[j][1:6], 10))
end
println()
for i in 1:5
    print(rpad(guild_names[i], 11))
    for j in 1:5
        interaction = dot(means_G[i], means_R[j])
        print(rpad(string(round(interaction, digits=2)), 10))
    end
    println()
end

# Compare 2D vs 4D
# Note: 4D needs higher scale because B^4_+ has smaller volume than B^2_+
# Volume ratio: B^d_+ ~ π^(d/2) / (2^d * Γ(d/2+1))
scenarios = [
    (
        name = "2D Latent Space",
        dim = 2,
        means_G = [Vector(Bd_plus_from_hyperspherical(0.85, [π/10])),
                   Vector(Bd_plus_from_hyperspherical(0.85, [π/4])),
                   Vector(Bd_plus_from_hyperspherical(0.85, [π/3]))],
        means_R = [Vector(Bd_plus_from_hyperspherical(0.85, [π/6])),
                   Vector(Bd_plus_from_hyperspherical(0.85, [π/5])),
                   Vector(Bd_plus_from_hyperspherical(0.85, [π/4]))],
        weights_G = [0.5, 0.3, 0.2],
        weights_R = [0.3, 0.4, 0.3],
        κ = 50.0,
        scale = 120.0,
    ),
    (
        name = "4D Latent Space",
        dim = 4,
        means_G = means_G,
        means_R = means_R,
        weights_G = weights_G,
        weights_R = weights_R,
        κ = 30.0,
        scale = 8000.0,  # Very high scale to compensate for low radii
    ),
]

# Helper function: assign point to nearest guild mean (defined early for use in comparison)
function assign_to_nearest_guild(point, guild_means)
    min_dist = Inf
    best_guild = 1
    for (i, mean) in enumerate(guild_means)
        dist = norm(point .- mean)
        if dist < min_dist
            min_dist = dist
            best_guild = i
        end
    end
    return best_guild
end

fig = Figure(size=(1000, 900))

println("\n" * "=" ^ 60)
println("Comparing 2D vs 4D Latent Space")
println("=" ^ 60)

# Minimum edge count for graph visualization (filters noise)
comparison_min_count = 3.0

for (col, scenario) in enumerate(scenarios)
    println("\n--- ", scenario.name, " ---")

    n_guilds = length(scenario.weights_G)
    κ_G = fill(scenario.κ, n_guilds)
    κ_R = fill(scenario.κ, n_guilds)

    ρ_G = BdPlusMixture(scenario.weights_G, scenario.means_G, κ_G, scenario.scale)
    ρ_R = BdPlusMixture(scenario.weights_R, scenario.means_R, κ_R, scenario.scale)
    ρ = ProductIntensity(ρ_G, ρ_R)

    stats = marginal_stats(ρ; rng=rng)
    println("  Dimension: ", scenario.dim)
    println("  Number of guilds: ", n_guilds)
    println("  Expected opportunities E[N]: ", round(stats.E_N, digits=1))
    println("  Avg connection prob: ", round(stats.avg_conn_prob, digits=3))
    println("  Expected interactions E[L]: ", round(stats.E_edges_edge_centric, digits=1))

    # Sample interactions
    sites = sample_ppp_product(ρ; rng=MersenneTwister(42 + col))
    interactions = generate_edge_centric(sites; rng=MersenneTwister(42 + col))
    println("  Sampled opportunities: ", length(sites))
    println("  Realized interactions: ", length(interactions))

    # Assign to guilds using nearest-mean (not k-means)
    n_clusters = n_guilds
    edge_weights = zeros(Int, n_clusters, n_clusters)
    for k in 1:length(interactions)
        src_guild = assign_to_nearest_guild(interactions.sources[k], scenario.means_G)
        tgt_guild = assign_to_nearest_guild(interactions.targets[k], scenario.means_R)
        edge_weights[src_guild, tgt_guild] += 1
    end

    # Count filtered edges (above threshold, excluding diagonal)
    n_links = 0
    for i in 1:n_clusters
        for j in 1:n_clusters
            if i != j && edge_weights[i, j] >= comparison_min_count
                n_links += 1
            end
        end
    end
    println("  Guilds: ", n_clusters, ", Links (count≥", comparison_min_count, "): ", n_links)

    # Row 1: Food web structure (filtered by count threshold)
    ax1 = Axis(fig[1, col], aspect=DataAspect(),
        title = scenario.name * "\nFood web (" * string(n_links) * " links, count≥" * string(Int(comparison_min_count)) * ")")
    hidedecorations!(ax1)
    hidespines!(ax1)
    limits!(ax1, -1.5, 1.5, -1.5, 1.5)

    if n_clusters > 0 && n_links > 0
        food_web = SimpleDiGraph(n_clusters)
        edge_widths = Float64[]
        arrow_sizes = Float64[]
        max_weight = maximum(edge_weights)

        for i in 1:n_clusters
            for j in 1:n_clusters
                if i != j && edge_weights[i, j] >= comparison_min_count
                    add_edge!(food_web, j, i)  # Arrow: consumer j → resource i
                    w_normalized = edge_weights[i, j] / max_weight
                    width = 0.5 + 6.0 * w_normalized^0.6
                    push!(edge_widths, width)
                    push!(arrow_sizes, 5.0 + 18.0 * w_normalized^0.6)
                end
            end
        end

        node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_clusters]

        graphplot!(ax1, food_web,
            layout=Shell(),
            node_size=35,
            node_color=node_colors,
            nlabels=string.(1:n_clusters),
            nlabels_align=(:center, :center),
            nlabels_color=:white,
            nlabels_fontsize=14,
            edge_color=(:black, 0.8),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    end

    # Row 2: Edge weight matrix as heatmap (with colorbar below)
    ax2 = Axis(fig[2, col],
        title = "Interaction counts",
        xlabel = "Consumer",
        ylabel = "Resource",
        xticks = 1:n_clusters,
        yticks = 1:n_clusters,
        yreversed = true)

    # Exclude diagonal for visualization
    weights_no_diag = copy(edge_weights)
    for i in 1:n_clusters
        weights_no_diag[i, i] = 0
    end

    hm = heatmap!(ax2, permutedims(weights_no_diag), colormap=:viridis)
    Colorbar(fig[3, col], hm, label="Count", vertical=false)
end

# Add title
Label(fig[0, :], "Food Web Structure: 2D vs 4D Latent Space", fontsize=18, font=:bold)

save("output/applications/ecological_4d_comparison.png", fig)
println("\nSaved ecological_4d_comparison.png")

# --- 2x2 Grid: Variability across realizations ---
println("\n" * "=" ^ 60)
println("4D Food Web Variability (2x2 grid of samples)")
println("=" ^ 60)

# Use the 4D scenario settings
scenario_4d = scenarios[2]
n_guilds_4d = length(scenario_4d.weights_G)
κ_4d = fill(scenario_4d.κ, n_guilds_4d)

ρ_G_4d = BdPlusMixture(scenario_4d.weights_G, scenario_4d.means_G, κ_4d, scenario_4d.scale)
ρ_R_4d = BdPlusMixture(scenario_4d.weights_R, scenario_4d.means_R, κ_4d, scenario_4d.scale)
ρ_4d = ProductIntensity(ρ_G_4d, ρ_R_4d)

fig_var = Figure(size=(1000, 1000))
Label(fig_var[0, :], "4D Food Web: Variability Across Samples", fontsize=18, font=:bold)

# Trophic layout for consistent visualization
trophic_levels_4d = [0, 1, 1, 2, 3]
x_offsets_4d = [0.0, -0.5, 0.5, 0.0, 0.0]
trophic_positions_4d = [Point2f(x_offsets_4d[i], trophic_levels_4d[i]) for i in 1:n_guilds_4d]
node_colors_4d = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_guilds_4d]

sample_seeds = [1001, 1002, 1003, 1004]  # Different seeds for variety

for (idx, seed) in enumerate(sample_seeds)
    row = div(idx - 1, 2) + 1
    col = mod(idx - 1, 2) + 1

    # Sample one realization
    sites = sample_ppp_product(ρ_4d; rng=MersenneTwister(seed))
    interactions = generate_edge_centric(sites; rng=MersenneTwister(seed))

    # Assign to guilds
    edge_weights = zeros(Int, n_guilds_4d, n_guilds_4d)
    for k in 1:length(interactions)
        src_guild = assign_to_nearest_guild(interactions.sources[k], scenario_4d.means_G)
        tgt_guild = assign_to_nearest_guild(interactions.targets[k], scenario_4d.means_R)
        edge_weights[src_guild, tgt_guild] += 1
    end

    # Count edges (no threshold - show raw sample)
    n_edges = 0
    for i in 1:n_guilds_4d
        for j in 1:n_guilds_4d
            if i != j && edge_weights[i, j] > 0
                n_edges += 1
            end
        end
    end

    ax = Axis(fig_var[row, col], aspect=DataAspect(),
        title = "Sample " * string(idx) * " (" * string(length(interactions)) * " interactions, " * string(n_edges) * " edges)")
    hidedecorations!(ax)
    hidespines!(ax)
    limits!(ax, -1.5, 1.5, -0.5, 3.5)

    # Build and plot graph
    sample_web = SimpleDiGraph(n_guilds_4d)
    edge_widths = Float64[]
    arrow_sizes = Float64[]
    max_w = max(1, maximum(edge_weights))

    for i in 1:n_guilds_4d
        for j in 1:n_guilds_4d
            if i != j && edge_weights[i, j] > 0
                add_edge!(sample_web, j, i)
                w_norm = edge_weights[i, j] / max_w
                push!(edge_widths, 0.5 + 4.0 * w_norm^0.6)
                push!(arrow_sizes, 6.0 + 14.0 * w_norm^0.6)
            end
        end
    end

    if ne(sample_web) > 0
        graphplot!(ax, sample_web,
            layout=(_) -> trophic_positions_4d,
            node_size=40,
            node_color=node_colors_4d,
            edge_color=(:black, 0.7),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    else
        # Just show nodes if no edges
        for (i, pos) in enumerate(trophic_positions_4d)
            scatter!(ax, [pos], color=node_colors_4d[i], markersize=40)
        end
    end

    # Add trophic level labels on leftmost column
    if col == 1
        ax.yticks = ([0, 1, 2, 3], ["Prod.", "Herb.", "Pred.", "Apex"])
        ax.yticklabelsvisible = true
        ax.yticksvisible = false
    end

    println("  Sample ", idx, ": ", length(interactions), " interactions, ", n_edges, " edges")
end

save("output/applications/ecological_4d_variability.png", fig_var)
println("Saved ecological_4d_variability.png")

# --- Detailed 4D analysis ---
println("\n" * "=" ^ 60)
println("Detailed 4D Food Web Analysis")
println("=" ^ 60)

# Use 4D scenario with more interactions
n_guilds = 5
κ_vals = [40.0, 40.0, 40.0, 40.0, 40.0]
scale_4d = 40000.0  # High scale to get sufficient interactions with orthogonal niches

ρ_G = BdPlusMixture(weights_G, means_G, κ_vals, scale_4d)
ρ_R = BdPlusMixture(weights_R, means_R, κ_vals, scale_4d)
ρ = ProductIntensity(ρ_G, ρ_R)

stats = marginal_stats(ρ; rng=rng)
println("Expected opportunities E[N]: ", round(stats.E_N, digits=1))
println("Expected interactions E[L]: ", round(stats.E_edges_edge_centric, digits=1))

# Multiple realizations - assign interactions to guilds by nearest mean
n_trials = 50

interaction_matrices = []
edge_presence = zeros(n_guilds, n_guilds)  # Count trials where each edge is present

for t in 1:n_trials
    sites = sample_ppp_product(ρ; rng=MersenneTwister(100 + t))
    interactions = generate_edge_centric(sites; rng=MersenneTwister(100 + t))

    if length(interactions) >= 5
        # Assign each source/target to nearest guild mean
        edge_weights = zeros(n_guilds, n_guilds)
        for k in 1:length(interactions)
            src_guild = assign_to_nearest_guild(interactions.sources[k], means_G)
            tgt_guild = assign_to_nearest_guild(interactions.targets[k], means_R)
            edge_weights[src_guild, tgt_guild] += 1
        end
        push!(interaction_matrices, edge_weights)

        # Track edge presence (binary: is edge present in this trial?)
        edge_presence .+= (edge_weights .> 0)
    end
end

println("Valid trials (with enough interactions): ", length(interaction_matrices), "/", n_trials)

if isempty(interaction_matrices)
    println("Warning: No valid trials - try increasing scale")
    avg_matrix = zeros(n_guilds, n_guilds)
    edge_frequency = zeros(n_guilds, n_guilds)
else
    # Average interaction matrix
    avg_matrix = zeros(n_guilds, n_guilds)
    for mat in interaction_matrices
        avg_matrix .+= mat
    end
    avg_matrix ./= length(interaction_matrices)

    # Edge frequency: fraction of trials where each edge appeared
    edge_frequency = edge_presence ./ length(interaction_matrices)
end

# Filter edges by average count threshold (not frequency)
# This filters out weak edges that appear due to sampling noise
min_edge_count = 20.0  # Minimum average count to filter out low-probability noise
significant_edges = avg_matrix .>= min_edge_count
n_significant = count(significant_edges) - count(i -> significant_edges[i,i], 1:n_guilds)  # exclude diagonal
println("Edges with avg count ≥", min_edge_count, ": ", n_significant, " (excl. diagonal)")

println("\nAverage interaction counts (* = included in filtered graph, count ≥", min_edge_count, "):")
print("           ")
for j in 1:n_guilds
    print(rpad(guild_names[j][1:6], 10))
end
println()
for i in 1:n_guilds
    print(rpad(guild_names[i], 11))
    for j in 1:n_guilds
        count_str = string(round(avg_matrix[i, j], digits=1))
        # Mark significant edges (above threshold, excluding diagonal)
        if i != j && significant_edges[i, j]
            count_str = count_str * "*"
        end
        print(rpad(count_str, 10))
    end
    println()
end

# Create detailed figure with 2 panels
fig2 = Figure(size=(1100, 500))

# Panel 1: Average interaction counts (observed)
ax1 = Axis(fig2[1, 1],
    title = "Observed Interactions\n(avg over " * string(n_trials) * " trials)",
    xlabel = "Consumer guild",
    ylabel = "Resource guild",
    xticks = (1:n_guilds, guild_names),
    yticks = (1:n_guilds, guild_names),
    xticklabelrotation = π/4,
    yreversed = true)

hm = heatmap!(ax1, permutedims(avg_matrix), colormap=:YlOrRd)
Colorbar(fig2[1, 2], hm, label="Avg. interactions")

# Panel 2: Expected interaction strength
ax2 = Axis(fig2[1, 3],
    title = "Expected Interaction\n(g_resource · r_consumer)",
    xlabel = "Consumer guild",
    ylabel = "Resource guild",
    xticks = (1:n_guilds, guild_names),
    yticks = (1:n_guilds, guild_names),
    xticklabelrotation = π/4,
    yreversed = true)

expected_matrix = zeros(n_guilds, n_guilds)
for i in 1:n_guilds
    for j in 1:n_guilds
        expected_matrix[i, j] = dot(means_G[i], means_R[j])
    end
end

hm2 = heatmap!(ax2, permutedims(expected_matrix), colormap=:YlOrRd)
Colorbar(fig2[1, 4], hm2, label="g · r")

save("output/applications/ecological_4d_detailed.png", fig2)
println("\nSaved ecological_4d_detailed.png")

# Create filtered food web graph
println("\n--- Filtered Food Web (avg count ≥" * string(min_edge_count) * ") ---")
fig3 = Figure(size=(900, 700))
ax_graph = Axis(fig3[1, 1],
    title = "Filtered Food Web (edges with avg count ≥" * string(Int(min_edge_count)) * ")")
hidedecorations!(ax_graph)
hidespines!(ax_graph)

# Build filtered graph using count threshold
filtered_web = SimpleDiGraph(n_guilds)
filtered_edge_widths = Float64[]
filtered_arrow_sizes = Float64[]

for i in 1:n_guilds
    for j in 1:n_guilds
        if i != j && significant_edges[i, j]
            add_edge!(filtered_web, j, i)  # Arrow: consumer j → resource i
            w_normalized = avg_matrix[i, j] / maximum(avg_matrix)
            push!(filtered_edge_widths, 1.0 + 5.0 * w_normalized^0.6)
            push!(filtered_arrow_sizes, 8.0 + 18.0 * w_normalized^0.6)
        end
    end
end

println("Filtered edges: ", ne(filtered_web))

if ne(filtered_web) > 0
    node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_guilds]

    # Custom trophic layout: y = trophic level, x = spread within level
    # Producers (1) at bottom, Apex (5) at top
    trophic_levels = [0, 1, 1, 2, 3]  # Producers=0, Herbivores=1, SmallPred=2, Apex=3
    x_offsets = [0.0, -0.5, 0.5, 0.0, 0.0]  # Spread herbivores horizontally
    trophic_positions = [Point2f(x_offsets[i], trophic_levels[i]) for i in 1:n_guilds]

    graphplot!(ax_graph, filtered_web,
        layout=(_) -> trophic_positions,
        node_size=50,
        node_color=node_colors,
        edge_color=(:black, 0.7),
        edge_width=filtered_edge_widths,
        arrow_size=filtered_arrow_sizes,
        arrow_show=true)

    limits!(ax_graph, -1.5, 1.5, -0.5, 3.5)
    ax_graph.yticks = ([0, 1, 2, 3], ["Producers", "Herbivores", "Predators", "Apex"])
    hidespines!(ax_graph)
    ax_graph.yticklabelsvisible = true
    ax_graph.yticksvisible = false
    ax_graph.ygridvisible = false
end

save("output/applications/ecological_4d_filtered_web.png", fig3)
println("Saved ecological_4d_filtered_web.png")

println("\n" * "=" ^ 60)
println("Ecological Interpretation (4D):")
println("=" ^ 60)
println("""
In this 4D food web model, each dimension represents a trophic level:
  - Dimension 1: Producer niche (plants, phytoplankton)
  - Dimension 2: Herbivore niche (grazers, filter feeders)
  - Dimension 3: Predator niche (secondary consumers)
  - Dimension 4: Apex predator niche (top predators)

Key design principles:
  - Resource role (g): species peak in their own trophic dimension
  - Consumer role (r): species target the dimension BELOW them
  - Radius controls visibility: apex have low g-radius (rarely eaten)
  - Radius controls consumption: producers have near-zero r-radius

Expected trophic structure (matrix should be heavy ABOVE diagonal):
  - Row 1 (Producers) eaten by Cols 2-3 (Herbivores): HIGH
  - Rows 2-3 (Herbivores) eaten by Cols 4-5 (Predators): HIGH
  - Row 4 (Small Pred) eaten by Col 5 (Apex): HIGH
  - Row 5 (Apex) eaten by anyone: LOW (top of chain)
  - Column 1 (Producers as consumers): ~0 (autotrophs don't eat)
  - Diagonal (self-consumption): LOW (species don't eat themselves)
""")

println("\nDone!")
