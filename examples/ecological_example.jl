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

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Ecological Example: Food Web Modeling")
println("=" ^ 60)

# In this model:
# - The "green" coordinate represents an individual's role as a RESOURCE
#   (i.e., its propensity to be eaten, its nutritional profile)
# - The "red" coordinate represents an individual's role as a CONSUMER
#   (i.e., its foraging behavior, dietary preferences)
#
# Under the edge-centric interpretation, each sampled point (g, r) represents
# a trophic interaction: a consumer "at r" eating a resource "at g".

# We model a simple 2D ecosystem with two trophic niche axes:
# - Axis 1: Body size / energy content
# - Axis 2: Defensive capability / handling time

println("\n--- Setting Up Ecological Intensity ---")

# Resource intensity (green space): where individuals are positioned as food
# Primary producers have high nutritional value (high x1) and low defense (low x2)
ρ_G = BdPlusMixture(
    [0.7, 0.2, 0.1],  # Most biomass is at the base (producers)
    [
        [0.8, 0.2],  # Producers: high energy, low defense
        [0.5, 0.5],  # Herbivore resources
        [0.3, 0.7],  # Well-defended prey
    ],
    [15.0, 10.0, 10.0],  # Concentrations
    100.0  # Total resource availability
)

# Consumer intensity (red space): where individuals are positioned as eaters
# Consumers are positioned by their foraging preferences
ρ_R = BdPlusMixture(
    [0.3, 0.4, 0.3],  # More even distribution of consumers
    [
        [0.7, 0.3],  # Herbivores preferring high-energy resources
        [0.5, 0.5],  # Generalist consumers
        [0.4, 0.8],  # Specialists that can handle defended prey
    ],
    [10.0, 8.0, 12.0],
    50.0  # Total consumer activity
)

ρ = ProductIntensity(ρ_G, ρ_R)

# Compute statistics
stats = marginal_stats(ρ; rng=rng)

println("\nIntensity Statistics:")
println("  Resource intensity (c_G): ", round(stats.c_G, digits=2))
println("  Consumer intensity (c_R): ", round(stats.c_R, digits=2))
println("  Expected opportunities (E[N]): ", round(stats.E_N, digits=2))
println("  Mean resource position (μ̃_G): ", round.(stats.μ̃_G, digits=3))
println("  Mean consumer position (μ̃_R): ", round.(stats.μ̃_R, digits=3))
println("  Average connection probability: ", round(stats.avg_conn_prob, digits=4))
println("  Expected trophic interactions (E[L]): ", round(stats.E_edges_edge_centric, digits=1))

# Sample trophic interactions (edge-centric)
println("\n--- Sampling Trophic Interactions ---")

# First sample the opportunity structure
sites = sample_ppp_product(ρ; rng=rng)
println("Sampled ", length(sites), " potential interaction opportunities")

# Then determine which actually result in consumption (edge-centric)
trophic_interactions = generate_edge_centric(sites; rng=rng)
println("Realized trophic interactions: ", length(trophic_interactions))

# Analyze the interactions
println("\n--- Analyzing Trophic Structure ---")

# Discretize into "species" by clustering
n_species = 6  # Number of species-level clusters
graph, source_assign, target_assign = discretize_edge_centric(trophic_interactions, n_species; rng=rng)

println("\nDiscretized food web:")
println("  Species (clusters): ", n_species)
println("  Unique trophic links: ", ne(graph))

# Create visualization
println("\n--- Generating Plots ---")

fig = Figure(size=(1200, 800))

# Plot 1: Resource intensity
ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Resource Distribution (ρ_G)\n'Who gets eaten'")
plot_intensity_Bd_plus!(ax1, ρ_G; resolution=30)
text!(ax1, Point2f(1.05, 0.0), text="Energy", fontsize=10, align=(:left, :center))
text!(ax1, Point2f(0.0, 1.05), text="Defense", fontsize=10, align=(:center, :bottom))

# Plot 2: Consumer intensity
ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Consumer Distribution (ρ_R)\n'Who does the eating'")
plot_intensity_Bd_plus!(ax2, ρ_R; resolution=30)
text!(ax2, Point2f(1.05, 0.0), text="Energy pref.", fontsize=10, align=(:left, :center))
text!(ax2, Point2f(0.0, 1.05), text="Defense handling", fontsize=10, align=(:center, :bottom))

# Plot 3: Sampled resources (what gets eaten)
ax3 = Axis(fig[2, 1], aspect=DataAspect(), title="Resources Consumed\n(Source positions of interactions)")
draw_Bd_plus_boundary!(ax3)
if !isempty(trophic_interactions.sources)
    source_2d = [Bd_plus_to_2d(s) for s in trophic_interactions.sources]
    scatter!(ax3, source_2d, color=:green, markersize=5, alpha=0.5)
end

# Plot 4: Sampled consumers (who does the eating)
ax4 = Axis(fig[2, 2], aspect=DataAspect(), title="Active Consumers\n(Target positions of interactions)")
draw_Bd_plus_boundary!(ax4)
if !isempty(trophic_interactions.targets)
    target_2d = [Bd_plus_to_2d(t) for t in trophic_interactions.targets]
    scatter!(ax4, target_2d, color=:red, markersize=5, alpha=0.5)
end

save("output/applications/ecological_foodweb.png", fig)
println("Saved ecological_foodweb.png")

# Monte Carlo validation
println("\n" * "=" ^ 60)
println("Monte Carlo Validation")
println("=" ^ 60)

n_trials = 200
L_samples = Int[]

for _ in 1:n_trials
    sites_trial = sample_ppp_product(ρ; rng=rng)
    interactions_trial = generate_edge_centric(sites_trial; rng=rng)
    push!(L_samples, length(interactions_trial))
end

println("\nExpected trophic interactions E[L]:")
println("  Theory:    ", round(stats.E_edges_edge_centric, digits=1))
println("  Empirical: ", round(mean(L_samples), digits=1), " ± ", round(std(L_samples)/sqrt(n_trials), digits=1))

println("\n" * "=" ^ 60)
println("Ecological Interpretation:")
println("=" ^ 60)
println("""
In this food web model:

1. The intensity ρ(g, r) = ρ_G(g) · ρ_R(r) represents the 'opportunity structure'
   for trophic interactions - which resource-consumer pairs have the chance to meet.

2. The connection probability g · r represents whether an encounter results in
   predation, based on the alignment of resource availability and consumer preference.

3. The edge-centric interpretation is natural because each predation event involves
   ephemeral individuals - the consumed individual ceases to exist after the interaction.

4. Species-level food webs (as typically studied) are statistical summaries:
   they cluster individual interactions into species nodes and count link frequencies.

5. Evolutionary dynamics (e.g., adaptation, speciation) can be modeled as PDE
   evolution of the intensities ρ_G and ρ_R on B^d_+.
""")

println("\nDone!")
