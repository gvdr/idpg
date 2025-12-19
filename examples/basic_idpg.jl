# Basic IDPG Example
# Demonstrates sampling from a product intensity and generating graphs
# under both node-centric and edge-centric interpretations.

using IDPG
using Distributions
using Random
using CairoMakie
using Graphs

# Set random seed for reproducibility
rng = MersenneTwister(42)

# Set up 2-dimensional latent space (B^2_+)
d = 2

# Define product intensity using BdPlusMixture
# ρ_G: intensity on the "green" (source/proposing) space
# ρ_R: intensity on the "red" (target/accepting) space

# Means are positions in B^d_+ (non-negative unit ball)
# Concentrations control how peaked the distribution is around the mean
ρ_G = BdPlusMixture(
    [0.6, 0.4],  # Mixture weights
    [[0.8, 0.2], [0.2, 0.8]],  # Mean positions in B^2_+
    [10.0, 10.0],  # Concentrations
    50.0  # Total intensity (expected number of green coordinates)
)

ρ_R = BdPlusMixture(
    [0.5, 0.5],
    [[0.3, 0.7], [0.5, 0.5]],
    [10.0, 5.0],
    50.0
)

ρ = ProductIntensity(ρ_G, ρ_R)

# Compute theoretical statistics
stats = marginal_stats(ρ; rng=rng)
println("=== Theoretical Statistics ===")
println("c_G (total green intensity): ", round(stats.c_G, digits=2))
println("c_R (total red intensity): ", round(stats.c_R, digits=2))
println("E[N] (expected nodes): ", round(stats.E_N, digits=2))
println("μ̃_G (normalized green mean): ", stats.μ̃_G)
println("μ̃_R (normalized red mean): ", stats.μ̃_R)
println("Average connection probability: ", round(stats.avg_conn_prob, digits=4))
println("E[|E|] (node-centric edges): ", round(stats.E_edges_node_centric, digits=2))
println("E[L] (edge-centric edges): ", round(stats.E_edges_edge_centric, digits=2))
println()

# Sample interaction sites from PPP
println("=== Sampling ===")
sites = sample_ppp_product(ρ; rng=rng)
println("Sampled ", length(sites), " interaction sites")

# Generate node-centric graph
println("\n=== Node-Centric Graph ===")
graph_nc, _ = generate_node_centric(sites; rng=rng)
println("Nodes: ", nv(graph_nc))
println("Edges: ", ne(graph_nc))

# Generate edge-centric sample
println("\n=== Edge-Centric Sample ===")
sample_ec = generate_edge_centric(sites; rng=rng)
println("Edges: ", length(sample_ec))

# Visualize intensities
println("\n=== Generating Plots ===")

fig = Figure(size=(1200, 500))

# Plot ρ_G
ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Intensity ρ_G (Green/Source)")
plot_intensity_Bd_plus!(ax1, ρ_G; resolution=30)

# Plot ρ_R
ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Intensity ρ_R (Red/Target)")
plot_intensity_Bd_plus!(ax2, ρ_R; resolution=30)

# Plot sampled sites (green coordinates)
ax3 = Axis(fig[1, 3], aspect=DataAspect(), title="Sampled Sites (Green coords)")
draw_Bd_plus_boundary!(ax3)

g_points = [Bd_plus_to_2d(site.g) for site in sites]
scatter!(ax3, g_points, color=:green, markersize=6, alpha=0.7)

save("output/basics/basic_idpg.png", fig)
println("Saved basic_idpg.png")

# Run multiple trials to validate formulas
println("\n=== Formula Validation ===")
n_trials = 100

N_samples = Int[]
E_nc_samples = Int[]
L_samples = Int[]

for _ in 1:n_trials
    sites_trial = sample_ppp_product(ρ; rng=rng)
    push!(N_samples, length(sites_trial))

    if length(sites_trial) > 0
        graph_trial, _ = generate_node_centric(sites_trial; rng=rng)
        push!(E_nc_samples, ne(graph_trial))

        ec_trial = generate_edge_centric(sites_trial; rng=rng)
        push!(L_samples, length(ec_trial))
    else
        push!(E_nc_samples, 0)
        push!(L_samples, 0)
    end
end

using Statistics: mean, std

println("E[N]: Theory = ", round(stats.E_N, digits=2),
        ", Empirical = ", round(mean(N_samples), digits=2),
        " ± ", round(std(N_samples)/sqrt(n_trials), digits=2))

println("E[|E|] (node-centric): Theory = ", round(stats.E_edges_node_centric, digits=2),
        ", Empirical = ", round(mean(E_nc_samples), digits=2),
        " ± ", round(std(E_nc_samples)/sqrt(n_trials), digits=2))

println("E[L] (edge-centric): Theory = ", round(stats.E_edges_edge_centric, digits=2),
        ", Empirical = ", round(mean(L_samples), digits=2),
        " ± ", round(std(L_samples)/sqrt(n_trials), digits=2))

println("\nDone!")
