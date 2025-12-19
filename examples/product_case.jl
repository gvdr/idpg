# Product Case Analysis
# Demonstrates the key formulas for product intensities:
# - E[N] = c_G · c_R
# - E[|E|] = E[N]² · (μ̃_G · μ̃_R)  [node-centric]
# - E[L] = E[N] · (μ̃_G · μ̃_R)     [edge-centric]
# - E[|E|] / E[L] = E[N]

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Product Case: Formula Validation")
println("=" ^ 60)

# Test with different intensity scales, mean positions, and connectivity levels
# Note: scale parameter in BdPlusMixture controls total intensity of Gaussian kernels
# Connectivity is controlled by the dot product of normalized means (μ̃_G · μ̃_R)

test_cases = [
    # BASELINE: moderate connectivity
    (name="Baseline", scale_G=60.0, scale_R=60.0, mean_G=[0.5, 0.5], mean_R=[0.5, 0.5], conc=10.0),

    # SPARSE: interior orthogonal means (0.7,0.1)·(0.1,0.7)=0.14, ||mean||=0.71)
    # More interior means less truncation even with high κ
    (name="Sparse interior, κ=30", scale_G=80.0, scale_R=80.0, mean_G=[0.7, 0.1], mean_R=[0.1, 0.7], conc=30.0),

    # SPARSE: same means, higher concentration
    (name="Sparse interior, κ=80", scale_G=80.0, scale_R=80.0, mean_G=[0.7, 0.1], mean_R=[0.1, 0.7], conc=80.0),

    # VERY SPARSE: nearly orthogonal but more interior (||mean||=0.64)
    # (0.6, 0.2)·(0.2, 0.6) = 0.12 + 0.12 = 0.24
    (name="Very sparse, κ=50", scale_G=100.0, scale_R=100.0, mean_G=[0.6, 0.2], mean_R=[0.2, 0.6], conc=50.0),

    # EXTREMELY SPARSE: use small dot product interior means
    # (0.55, 0.1)·(0.1, 0.55) = 0.055 + 0.055 = 0.11
    (name="Extremely sparse, κ=100", scale_G=120.0, scale_R=120.0, mean_G=[0.55, 0.1], mean_R=[0.1, 0.55], conc=100.0),

    # LARGE + SPARSE: many nodes but low connectivity
    (name="Large sparse, κ=60", scale_G=150.0, scale_R=150.0, mean_G=[0.6, 0.15], mean_R=[0.15, 0.6], conc=60.0),
]

n_trials = 200

# Storage for plotting
theoretical_N = Float64[]
empirical_N = Float64[]
theoretical_E = Float64[]
empirical_E = Float64[]
theoretical_L = Float64[]
empirical_L = Float64[]

for (i, case) in enumerate(test_cases)
    println("\n--- Test Case $i: ", case.name, " ---")
    println("scale_G = ", case.scale_G, ", scale_R = ", case.scale_R)
    println("mean_G = ", case.mean_G, ", mean_R = ", case.mean_R)

    # Create intensities using BdPlusMixture
    ρ_G = BdPlusMixture([1.0], [case.mean_G], [case.conc], case.scale_G)
    ρ_R = BdPlusMixture([1.0], [case.mean_R], [case.conc], case.scale_R)
    ρ = ProductIntensity(ρ_G, ρ_R)

    # Theoretical values
    stats = marginal_stats(ρ; rng=rng)

    push!(theoretical_N, stats.E_N)
    push!(theoretical_E, stats.E_edges_node_centric)
    push!(theoretical_L, stats.E_edges_edge_centric)

    # Empirical sampling
    N_samples = Int[]
    E_samples = Int[]
    L_samples = Int[]

    for _ in 1:n_trials
        sites = sample_ppp_product(ρ; rng=rng)
        push!(N_samples, length(sites))

        if length(sites) > 0
            graph, _ = generate_node_centric(sites; rng=rng)
            push!(E_samples, ne(graph))

            ec = generate_edge_centric(sites; rng=rng)
            push!(L_samples, length(ec))
        else
            push!(E_samples, 0)
            push!(L_samples, 0)
        end
    end

    push!(empirical_N, mean(N_samples))
    push!(empirical_E, mean(E_samples))
    push!(empirical_L, mean(L_samples))

    # Report
    println("\nExpected Nodes E[N]:")
    println("  Theory:    ", round(stats.E_N, digits=2))
    println("  Empirical: ", round(mean(N_samples), digits=2), " ± ", round(std(N_samples)/sqrt(n_trials), digits=2))

    println("\nNode-centric E[|E|] = E[N]² · (μ̃_G · μ̃_R):")
    println("  Theory:    ", round(stats.E_edges_node_centric, digits=2))
    println("  Empirical: ", round(mean(E_samples), digits=2), " ± ", round(std(E_samples)/sqrt(n_trials), digits=2))

    println("\nEdge-centric E[L] = E[N] · (μ̃_G · μ̃_R):")
    println("  Theory:    ", round(stats.E_edges_edge_centric, digits=2))
    println("  Empirical: ", round(mean(L_samples), digits=2), " ± ", round(std(L_samples)/sqrt(n_trials), digits=2))

    println("\nRatio E[|E|]/E[L] ≈ E[N]:")
    if mean(L_samples) > 0
        println("  Theory:    ", round(stats.E_N, digits=2))
        println("  Empirical: ", round(mean(E_samples)/mean(L_samples), digits=2))
    end

    # Average degree (out-degree in directed graph)
    avg_degree_theory = stats.E_edges_node_centric / stats.E_N
    avg_degree_empirical = mean(E_samples) / mean(N_samples)
    println("\nAverage out-degree (E[|E|]/E[N]):")
    println("  Theory:    ", round(avg_degree_theory, digits=2), " (conn. prob = ", round(stats.avg_conn_prob, digits=3), ")")
    println("  Empirical: ", round(avg_degree_empirical, digits=2))
end

# Create validation plots with log-log scales
println("\n" * "=" ^ 60)
println("Creating validation plots...")

fig = Figure(size=(1200, 400))

# Helper to get log-scale limits
function log_lims(vals)
    mn = minimum(vals) * 0.5
    mx = maximum(vals) * 2.0
    return (mn, mx)
end

ax1 = Axis(fig[1, 1],
    xlabel="Theoretical E[N]",
    ylabel="Empirical E[N]",
    title="Node Count Validation",
    xscale=log10, yscale=log10)
scatter!(ax1, theoretical_N, empirical_N, color=:blue, markersize=15)
lims = log_lims(vcat(theoretical_N, empirical_N))
lines!(ax1, [lims[1], lims[2]], [lims[1], lims[2]], color=:red, linestyle=:dash, linewidth=2)

ax2 = Axis(fig[1, 2],
    xlabel="Theoretical E[|E|]",
    ylabel="Empirical E[|E|]",
    title="Node-Centric Edge Count",
    xscale=log10, yscale=log10)
scatter!(ax2, theoretical_E, empirical_E, color=:green, markersize=15)
lims = log_lims(vcat(theoretical_E, empirical_E))
lines!(ax2, [lims[1], lims[2]], [lims[1], lims[2]], color=:red, linestyle=:dash, linewidth=2)

ax3 = Axis(fig[1, 3],
    xlabel="Theoretical E[L]",
    ylabel="Empirical E[L]",
    title="Edge-Centric Edge Count",
    xscale=log10, yscale=log10)
scatter!(ax3, theoretical_L, empirical_L, color=:purple, markersize=15)
lims = log_lims(vcat(theoretical_L, empirical_L))
lines!(ax3, [lims[1], lims[2]], [lims[1], lims[2]], color=:red, linestyle=:dash, linewidth=2)

save("output/validation/product_case_validation.png", fig)
println("Saved product_case_validation.png")

println("\n" * "=" ^ 60)
println("Done!")
