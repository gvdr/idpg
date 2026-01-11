# Edge-Centric Model Comparison
# Demonstrates and validates the two edge-centric interpretations:
#
# 1. ENCOUNTER MODEL (generate_edge_centric):
#    - Sample N sites from PPP, each site (g, r) is an encounter opportunity
#    - Each encounter becomes an edge with probability g . r
#    - Expected edges: E[L] = E[N] * E[g.r]
#
# 2. PAIRING MODEL (sample_edge_centric):
#    - Sample L ~ Poisson(E[N]/2) opportunities
#    - Each opportunity pairs independent source and target from the marginal intensities
#    - Each pair connects with probability g_source . r_target
#    - Expected edges: E[L] = (E[N]/2) * E[g.r]
#
# Key insight: In the pairing model, each edge opportunity "consumes" 2 node-equivalents
# (1 source + 1 target), so from E[N] node-equivalents we get E[N]/2 opportunities.
#
# Comparison with node-centric:
#    - Node-centric: N sites become nodes, all N^2 pairs considered
#    - E[|E|] = E[N]^2 * E[g.r]
#    - Ratio: E[|E|]_node / E[L]_pairing = 2*E[N]

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs

rng = MersenneTwister(42)

println("=" ^ 70)
println("Edge-Centric Model Comparison: Encounter vs Pairing")
println("=" ^ 70)

# Create a product intensity
ρ_G = BdPlusMixture(
    [0.6, 0.4],
    [[0.7, 0.3], [0.3, 0.7]],
    [15.0, 15.0],
    40.0
)

ρ_R = BdPlusMixture(
    [0.5, 0.5],
    [[0.4, 0.6], [0.6, 0.4]],
    [12.0, 12.0],
    40.0
)

ρ = ProductIntensity(ρ_G, ρ_R)

# Compute theoretical statistics
stats = marginal_stats(ρ)

println("\n--- Theoretical Statistics ---")
println("c_G (green total intensity): " * string(round(stats.c_G, digits=2)))
println("c_R (red total intensity):   " * string(round(stats.c_R, digits=2)))
println("E[N] = c_G * c_R:             " * string(round(stats.E_N, digits=2)))
println("E[g.r] (avg connection prob): " * string(round(stats.avg_conn_prob, digits=4)))
println()

println("--- Expected Edges by Model ---")
println()
println("NODE-CENTRIC: N sites -> N^2 pairs")
println("  E[|E|] = E[N]^2 * E[g.r] = " * string(round(stats.E_edges_node_centric, digits=2)))
println()
println("EDGE-CENTRIC (Encounter): N sites, each accepts with prob g.r")
println("  E[L] = E[N] * E[g.r] = " * string(round(stats.E_N * stats.avg_conn_prob, digits=2)))
println()
println("EDGE-CENTRIC (Pairing): E[N]/2 opportunities, each accepts with prob g.r")
println("  E[L] = (E[N]/2) * E[g.r] = " * string(round(stats.E_edges_edge_centric, digits=2)))
println()
println("Ratio node/encounter = E[N] = " * string(round(stats.E_N, digits=2)))
println("Ratio node/pairing = 2*E[N] = " * string(round(2 * stats.E_N, digits=2)))

# Create symmetric edge intensity for pairing model
ei = SymmetricEdgeIntensity(ρ)
println("\nSymmetricEdgeIntensity: C_edge = E[N]/2 = " * string(round(ei.C_edge, digits=2)))

# Monte Carlo validation
println("\n" * "=" ^ 70)
println("Monte Carlo Validation")
println("=" ^ 70)

n_trials = 300

# Storage
N_samples = Int[]
E_node_samples = Int[]
L_encounter_samples = Int[]
L_pairing_samples = Int[]

for trial in 1:n_trials
    # Sample sites for node-centric and encounter model
    sites = sample_ppp_product(ρ; rng=MersenneTwister(trial))
    push!(N_samples, length(sites))

    # Node-centric: N sites -> N^2 pairs
    if length(sites) > 0
        graph_nc, _ = generate_node_centric(sites; rng=MersenneTwister(trial))
        push!(E_node_samples, ne(graph_nc))
    else
        push!(E_node_samples, 0)
    end

    # Encounter model: each site accepts with prob g.r
    if length(sites) > 0
        ec_encounter = generate_edge_centric(sites; rng=MersenneTwister(trial))
        push!(L_encounter_samples, length(ec_encounter))
    else
        push!(L_encounter_samples, 0)
    end

    # Pairing model: L ~ Poisson(E[N]/2) opportunities with independent source/target
    ec_pairing = sample_edge_centric(ei; rng=MersenneTwister(trial))
    push!(L_pairing_samples, length(ec_pairing))
end

# Results
println("\n--- Node Count ---")
println("  E[N] theory:    " * string(round(stats.E_N, digits=2)))
println("  E[N] empirical: " * string(round(mean(N_samples), digits=2)) * " +/- " * string(round(std(N_samples)/sqrt(n_trials), digits=2)))

println("\n--- Node-Centric Edges ---")
println("  E[|E|] theory:    " * string(round(stats.E_edges_node_centric, digits=2)))
println("  E[|E|] empirical: " * string(round(mean(E_node_samples), digits=2)) * " +/- " * string(round(std(E_node_samples)/sqrt(n_trials), digits=2)))

println("\n--- Encounter Model Edges ---")
E_encounter_theory = stats.E_N * stats.avg_conn_prob
println("  E[L] theory:    " * string(round(E_encounter_theory, digits=2)))
println("  E[L] empirical: " * string(round(mean(L_encounter_samples), digits=2)) * " +/- " * string(round(std(L_encounter_samples)/sqrt(n_trials), digits=2)))

println("\n--- Pairing Model Edges ---")
println("  E[L] theory:    " * string(round(stats.E_edges_edge_centric, digits=2)))
println("  E[L] empirical: " * string(round(mean(L_pairing_samples), digits=2)) * " +/- " * string(round(std(L_pairing_samples)/sqrt(n_trials), digits=2)))

println("\n--- Ratios ---")
if mean(L_encounter_samples) > 0
    println("  Node/Encounter: theory = " * string(round(stats.E_N, digits=2)) *
            ", empirical = " * string(round(mean(E_node_samples)/mean(L_encounter_samples), digits=2)))
end
if mean(L_pairing_samples) > 0
    println("  Node/Pairing:   theory = " * string(round(2*stats.E_N, digits=2)) *
            ", empirical = " * string(round(mean(E_node_samples)/mean(L_pairing_samples), digits=2)))
end
if mean(L_pairing_samples) > 0
    println("  Encounter/Pairing: theory = 2.0, empirical = " *
            string(round(mean(L_encounter_samples)/mean(L_pairing_samples), digits=2)))
end

# Create visualization
println("\n" * "=" ^ 70)
println("Creating Visualizations")
println("=" ^ 70)

fig = Figure(size=(1400, 900))

# Row 1: Intensity plots and distributions
ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Intensity rho_G (Green/Source)")
plot_intensity_Bd_plus!(ax1, ρ_G; resolution=30)

ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Intensity rho_R (Red/Target)")
plot_intensity_Bd_plus!(ax2, ρ_R; resolution=30)

# Node count distribution
ax3 = Axis(fig[1, 3], xlabel="Number of Sites (N)", ylabel="Count",
    title="Site Count Distribution")
hist!(ax3, N_samples, bins=25, color=:steelblue)
vlines!(ax3, [stats.E_N], color=:red, linestyle=:dash, linewidth=2, label="E[N]")
axislegend(ax3, position=:rt)

# Row 2: Edge count distributions for each model
ax4 = Axis(fig[2, 1], xlabel="Number of Edges", ylabel="Count",
    title="Node-Centric: E[|E|] = E[N]^2 * E[g.r]")
hist!(ax4, E_node_samples, bins=30, color=:forestgreen)
vlines!(ax4, [stats.E_edges_node_centric], color=:red, linestyle=:dash, linewidth=2, label="E[|E|]")
axislegend(ax4, position=:rt)

ax5 = Axis(fig[2, 2], xlabel="Number of Edges", ylabel="Count",
    title="Encounter: E[L] = E[N] * E[g.r]")
hist!(ax5, L_encounter_samples, bins=20, color=:purple)
vlines!(ax5, [E_encounter_theory], color=:red, linestyle=:dash, linewidth=2, label="E[L]")
axislegend(ax5, position=:rt)

ax6 = Axis(fig[2, 3], xlabel="Number of Edges", ylabel="Count",
    title="Pairing: E[L] = (E[N]/2) * E[g.r]")
hist!(ax6, L_pairing_samples, bins=20, color=:orange)
vlines!(ax6, [stats.E_edges_edge_centric], color=:red, linestyle=:dash, linewidth=2, label="E[L]")
axislegend(ax6, position=:rt)

# Row 3: Comparison plots
ax7 = Axis(fig[3, 1:2], xlabel="E[N] (site count)", ylabel="E[edges]",
    title="Scaling Comparison: Node-Centric vs Edge-Centric Models")

# Show scaling relationships
N_range = collect(10:10:200)
E_node_scaling = [n^2 * stats.avg_conn_prob for n in N_range]
L_encounter_scaling = [n * stats.avg_conn_prob for n in N_range]
L_pairing_scaling = [(n/2) * stats.avg_conn_prob for n in N_range]

lines!(ax7, N_range, E_node_scaling, color=:forestgreen, linewidth=3, label="Node: N^2 * p")
lines!(ax7, N_range, L_encounter_scaling, color=:purple, linewidth=3, label="Encounter: N * p")
lines!(ax7, N_range, L_pairing_scaling, color=:orange, linewidth=3, label="Pairing: (N/2) * p")

# Add actual data point
scatter!(ax7, [mean(N_samples)], [mean(E_node_samples)], color=:forestgreen, markersize=15)
scatter!(ax7, [mean(N_samples)], [mean(L_encounter_samples)], color=:purple, markersize=15)
scatter!(ax7, [mean(N_samples)], [mean(L_pairing_samples)], color=:orange, markersize=15)

axislegend(ax7, position=:lt)

# Summary text
ax8 = Axis(fig[3, 3], limits=(0, 1, 0, 1))
hidedecorations!(ax8)
hidespines!(ax8)

summary_text = """
Model Summary:

Node-Centric: Long-lived entities
  N sites become nodes
  All N^2 pairs considered
  E[|E|] = E[N]^2 * E[g.r]

Encounter: Sites are opportunities
  N sites, each accepts with p = g.r
  E[L] = E[N] * E[g.r]

Pairing: Edge opportunities
  L ~ Poisson(E[N]/2) opportunities
  Each pairs source + target
  E[L] = (E[N]/2) * E[g.r]

Ratio: Node/Pairing = 2*E[N]
"""

text!(ax8, 0.05, 0.95, text=summary_text, align=(:left, :top), fontsize=12)

save("output/validation/edge_centric_comparison.png", fig)
println("Saved edge_centric_comparison.png")

# Model interpretation
println("\n" * "=" ^ 70)
println("Interpretation")
println("=" ^ 70)
println("""
The three models represent different assumptions about entities and interactions:

1. NODE-CENTRIC (long-lived entities):
   - Sites represent persistent entities (nodes) in a network
   - Every pair of nodes has the opportunity to connect
   - Typical for social networks, citation networks, etc.
   - E[|E|] = E[N]^2 * E[g.r] (quadratic in N)

2. ENCOUNTER MODEL (ephemeral interactions):
   - Sites represent encounter opportunities, not persistent entities
   - Each encounter (g, r) becomes an edge with probability g.r
   - The same site provides both source role (g) and target role (r)
   - Typical for animal interactions, disease transmission
   - E[L] = E[N] * E[g.r] (linear in N)

3. PAIRING MODEL (independent source-target):
   - Edge opportunities are sampled directly
   - Source and target positions are independent draws
   - Each opportunity "consumes" 2 node-equivalents
   - From E[N] node-equivalents, get E[N]/2 opportunities
   - E[L] = (E[N]/2) * E[g.r] (half of encounter model)

The pairing model is appropriate when:
  - Edges are the primary objects of interest
  - Source and target roles are fundamentally different
  - The "node" doesn't exist as a persistent entity

The encounter model is appropriate when:
  - Each site represents an encounter opportunity
  - The same individual provides both roles in an interaction
  - E.g., predator-prey where each site is an encounter event
""")

println("\nDone!")
