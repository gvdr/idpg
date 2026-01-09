# IDPG vs Erdős-Rényi Comparison
# Tests whether IDPG graphs have structure beyond simple random graphs

using IDPG
using Random
using Graphs
using Statistics
using CairoMakie
using LinearAlgebra: dot

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG vs Erdős-Rényi: Is there additional structure?")
println("=" ^ 60)

# --- Metric functions (same as mesoscale_metrics.jl) ---

function compute_reciprocity(g::SimpleDiGraph)
    ne(g) == 0 && return 0.0
    reciprocal_count = sum(has_edge(g, dst(e), src(e)) for e in edges(g))
    return reciprocal_count / ne(g)
end

function degree_correlation(g::SimpleDiGraph)
    nv(g) < 2 && return NaN
    in_degs, out_degs = Float64.(indegree(g)), Float64.(outdegree(g))
    (std(in_degs) == 0 || std(out_degs) == 0) && return 0.0
    return cor(in_degs, out_degs)
end

function directed_clustering(g::SimpleDiGraph)
    nv(g) < 3 && return 0.0
    triangles, triplets = 0, 0
    for v in vertices(g)
        all_neigh = union(Set(outneighbors(g, v)), Set(inneighbors(g, v)))
        for u in all_neigh, w in all_neigh
            if u != w
                triplets += 1
                (has_edge(g, u, w) || has_edge(g, w, u)) && (triangles += 1)
            end
        end
    end
    triplets == 0 ? 0.0 : triangles / triplets
end

function directed_assortativity(g::SimpleDiGraph)
    ne(g) == 0 && return NaN
    src_out = Float64[outdegree(g, src(e)) for e in edges(g)]
    dst_in = Float64[indegree(g, dst(e)) for e in edges(g)]
    (std(src_out) == 0 || std(dst_in) == 0) && return 0.0
    return cor(src_out, dst_in)
end

function degree_variance_ratio(g::SimpleDiGraph)
    # Ratio of degree variance to mean (coefficient of dispersion)
    # ER expects variance ≈ mean (Poisson-like), heterogeneous models have higher variance
    out_degs = Float64.(outdegree(g))
    mean_d = mean(out_degs)
    var_d = var(out_degs)
    mean_d == 0 ? 0.0 : var_d / mean_d
end

function largest_scc_fraction(g::SimpleDiGraph)
    sccs = strongly_connected_components(g)
    maximum(length.(sccs)) / nv(g)
end

"""
Generate an Erdős-Rényi directed graph matching IDPG graph statistics.
"""
function generate_matched_er(n::Int, m::Int; rng=Random.default_rng())
    # Edge probability to get expected m edges in directed graph with n nodes
    # Expected edges = n*(n-1)*p (no self-loops)
    p = m / (n * (n - 1))
    p = clamp(p, 0.0, 1.0)

    g = SimpleDiGraph(n)
    for i in 1:n, j in 1:n
        i != j && rand(rng) < p && add_edge!(g, i, j)
    end
    return g
end

# --- Setup IDPG intensity ---
println("\n--- IDPG Setup ---")

# Test multiple intensity configurations
configs = [
    (
        name = "Concentrated (orthogonal)",
        ρ_G = BdPlusMixture([1.0], [[0.65, 0.15]], [40.0], 80.0),
        ρ_R = BdPlusMixture([1.0], [[0.15, 0.65]], [40.0], 80.0)
    ),
    (
        name = "Spread (aligned)",
        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 40.0),  # Reduced scale
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 40.0)
    ),
    (
        name = "Bimodal G",
        ρ_G = BdPlusMixture([0.5, 0.5], [[0.8, 0.1], [0.1, 0.8]], [30.0, 30.0], 60.0),  # Reduced
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [10.0], 60.0)
    ),
]

n_realizations = 30

# Storage for all results
all_results = []

for config in configs
    println("\n" * "=" ^ 50)
    println("Configuration: ", config.name)
    println("=" ^ 50)

    ρ = ProductIntensity(config.ρ_G, config.ρ_R)
    stats = marginal_stats(ρ; rng=rng)
    println("Expected N: ", round(stats.E_N, digits=1), ", conn prob: ", round(stats.avg_conn_prob, digits=3))

    # Metrics storage
    idpg_metrics = (
        reciprocity = Float64[], clustering = Float64[], assortativity = Float64[],
        in_out_corr = Float64[], degree_var_ratio = Float64[], scc_frac = Float64[],
        in_degree_std = Float64[], out_degree_std = Float64[]
    )
    er_metrics = (
        reciprocity = Float64[], clustering = Float64[], assortativity = Float64[],
        in_out_corr = Float64[], degree_var_ratio = Float64[], scc_frac = Float64[],
        in_degree_std = Float64[], out_degree_std = Float64[]
    )

    for i in 1:n_realizations
        # Generate IDPG graph
        sites = sample_ppp_product(ρ; rng=rng)
        length(sites) < 3 && continue

        g_idpg, _ = generate_node_centric(sites; rng=rng)
        n, m = nv(g_idpg), ne(g_idpg)

        # Generate matched ER graph
        g_er = generate_matched_er(n, m; rng=rng)

        # Compute metrics for both
        push!(idpg_metrics.reciprocity, compute_reciprocity(g_idpg))
        push!(idpg_metrics.clustering, directed_clustering(g_idpg))
        push!(idpg_metrics.assortativity, directed_assortativity(g_idpg))
        push!(idpg_metrics.in_out_corr, degree_correlation(g_idpg))
        push!(idpg_metrics.degree_var_ratio, degree_variance_ratio(g_idpg))
        push!(idpg_metrics.scc_frac, largest_scc_fraction(g_idpg))
        push!(idpg_metrics.in_degree_std, std(indegree(g_idpg)))
        push!(idpg_metrics.out_degree_std, std(outdegree(g_idpg)))

        push!(er_metrics.reciprocity, compute_reciprocity(g_er))
        push!(er_metrics.clustering, directed_clustering(g_er))
        push!(er_metrics.assortativity, directed_assortativity(g_er))
        push!(er_metrics.in_out_corr, degree_correlation(g_er))
        push!(er_metrics.degree_var_ratio, degree_variance_ratio(g_er))
        push!(er_metrics.scc_frac, largest_scc_fraction(g_er))
        push!(er_metrics.in_degree_std, std(indegree(g_er)))
        push!(er_metrics.out_degree_std, std(outdegree(g_er)))
    end

    # Report comparison
    println("\nMetric                  IDPG              ER              Diff")
    println("-" ^ 65)

    metrics_to_compare = [
        ("Reciprocity", idpg_metrics.reciprocity, er_metrics.reciprocity),
        ("Clustering", idpg_metrics.clustering, er_metrics.clustering),
        ("Assortativity", idpg_metrics.assortativity, er_metrics.assortativity),
        ("In/Out Corr", idpg_metrics.in_out_corr, er_metrics.in_out_corr),
        ("Degree Var/Mean", idpg_metrics.degree_var_ratio, er_metrics.degree_var_ratio),
        ("Giant SCC Frac", idpg_metrics.scc_frac, er_metrics.scc_frac),
        ("In-degree Std", idpg_metrics.in_degree_std, er_metrics.in_degree_std),
        ("Out-degree Std", idpg_metrics.out_degree_std, er_metrics.out_degree_std),
    ]

    for (name, idpg_vals, er_vals) in metrics_to_compare
        idpg_mean, idpg_std = mean(idpg_vals), std(idpg_vals)
        er_mean, er_std = mean(er_vals), std(er_vals)
        diff = idpg_mean - er_mean
        # Effect size (Cohen's d)
        pooled_std = sqrt((idpg_std^2 + er_std^2) / 2)
        effect = pooled_std > 0 ? abs(diff) / pooled_std : 0.0

        marker = effect > 0.5 ? " **" : (effect > 0.2 ? " *" : "")
        println(rpad(name, 18),
            rpad(string(round(idpg_mean, digits=3), " ± ", round(idpg_std, digits=3)), 18),
            rpad(string(round(er_mean, digits=3), " ± ", round(er_std, digits=3)), 16),
            round(diff, digits=3), marker)
    end
    println("\n* = small effect (d>0.2), ** = medium effect (d>0.5)")

    push!(all_results, (config=config, idpg=idpg_metrics, er=er_metrics))
end

# --- Visualization ---
println("\n" * "=" ^ 60)
println("Creating comparison visualizations...")
println("=" ^ 60)

fig = Figure(size=(1400, 1000))

# For each configuration, show key metric comparisons
metric_names = ["Clustering", "Degree Var/Mean", "In-degree Std", "Assortativity"]
metric_keys = [:clustering, :degree_var_ratio, :in_degree_std, :assortativity]

for (row, result) in enumerate(all_results)
    for (col, (mname, mkey)) in enumerate(zip(metric_names, metric_keys))
        ax = Axis(fig[row, col],
            title = row == 1 ? mname : "",
            xlabel = mkey == :assortativity ? "Value" : "",
            ylabel = col == 1 ? result.config.name : "")

        idpg_vals = getfield(result.idpg, mkey)
        er_vals = getfield(result.er, mkey)

        # Paired comparison: show distributions
        hist!(ax, idpg_vals, bins=12, color=(:steelblue, 0.6), label="IDPG")
        hist!(ax, er_vals, bins=12, color=(:orange, 0.6), label="ER")

        if row == 1 && col == 1
            axislegend(ax, position=:rt)
        end
    end
end

save("output/mesoscale/idpg_vs_er_comparison.png", fig)
println("Saved idpg_vs_er_comparison.png")

# --- Focused plot: Degree distributions for one realization ---
println("\n--- Degree Distribution Comparison (single realization) ---")

fig2 = Figure(size=(1200, 400))

# Use the bimodal config (should show most difference)
config = configs[3]
ρ = ProductIntensity(config.ρ_G, config.ρ_R)

sites = sample_ppp_product(ρ; rng=MersenneTwister(123))
g_idpg, _ = generate_node_centric(sites; rng=MersenneTwister(123))
g_er = generate_matched_er(nv(g_idpg), ne(g_idpg); rng=MersenneTwister(456))

ax1 = Axis(fig2[1, 1], xlabel="Out-Degree", ylabel="Count", title="Out-Degree Distribution")
hist!(ax1, outdegree(g_idpg), bins=15, color=(:steelblue, 0.7), label="IDPG")
hist!(ax1, outdegree(g_er), bins=15, color=(:orange, 0.7), label="ER")
axislegend(ax1)

ax2 = Axis(fig2[1, 2], xlabel="In-Degree", ylabel="Count", title="In-Degree Distribution")
hist!(ax2, indegree(g_idpg), bins=15, color=(:steelblue, 0.7), label="IDPG")
hist!(ax2, indegree(g_er), bins=15, color=(:orange, 0.7), label="ER")
axislegend(ax2)

ax3 = Axis(fig2[1, 3], xlabel="Out-Degree", ylabel="In-Degree",
    title="In vs Out Degree Scatter")
scatter!(ax3, outdegree(g_idpg), indegree(g_idpg), color=:steelblue, markersize=10, label="IDPG")
scatter!(ax3, outdegree(g_er), indegree(g_er), color=:orange, markersize=10, label="ER")
axislegend(ax3, position=:lt)

save("output/mesoscale/idpg_vs_er_degrees.png", fig2)
println("Saved idpg_vs_er_degrees.png")

# --- Key insight: heterogeneity in connection probabilities ---
println("\n" * "=" ^ 60)
println("Connection Probability Heterogeneity Analysis")
println("=" ^ 60)

# For IDPG, the connection probability varies: P(i→j) = g_i · r_j
# For ER, it's constant: p for all pairs
# This should create degree heterogeneity in IDPG

config = configs[3]  # Bimodal
ρ = ProductIntensity(config.ρ_G, config.ρ_R)

sites = sample_ppp_product(ρ; rng=MersenneTwister(999))
n = length(sites)

if n > 2
    # Compute all pairwise connection probabilities
    conn_probs = Float64[]
    out_propensities = Float64[]  # Sum of g_i · r_j over j for each i
    in_propensities = Float64[]   # Sum of g_j · r_i over j for each i

    for i in 1:n
        out_prop = 0.0
        in_prop = 0.0
        for j in 1:n
            if i != j
                p_ij = dot(sites[i].g, sites[j].r)
                p_ji = dot(sites[j].g, sites[i].r)
                push!(conn_probs, p_ij)
                out_prop += p_ij
                in_prop += p_ji
            end
        end
        push!(out_propensities, out_prop)
        push!(in_propensities, in_prop)
    end

    println("\nConnection probability distribution:")
    println("  Mean: ", round(mean(conn_probs), digits=3))
    println("  Std:  ", round(std(conn_probs), digits=3))
    println("  Min:  ", round(minimum(conn_probs), digits=3))
    println("  Max:  ", round(maximum(conn_probs), digits=3))
    println("  CV:   ", round(std(conn_probs)/mean(conn_probs), digits=3), " (coefficient of variation)")

    println("\nExpected out-degree heterogeneity:")
    println("  Std of node out-propensities: ", round(std(out_propensities), digits=2))
    println("  (ER would have std = 0)")

    # Visualize
    fig3 = Figure(size=(1000, 400))

    ax1 = Axis(fig3[1, 1], xlabel="Connection Probability P(i→j)", ylabel="Count",
        title="Heterogeneity in Connection Probabilities\n(ER would be a single spike)")
    hist!(ax1, conn_probs, bins=30, color=:purple)

    ax2 = Axis(fig3[1, 2], xlabel="Expected Out-Degree", ylabel="Expected In-Degree",
        title="Node Propensity Heterogeneity\n(ER: all nodes at same point)")
    scatter!(ax2, out_propensities, in_propensities, color=:teal, markersize=12)

    save("output/mesoscale/connection_probability_heterogeneity.png", fig3)
    println("\nSaved connection_probability_heterogeneity.png")
end

println("\n" * "=" ^ 60)
println("Done!")
println("\nConclusion: IDPG ≠ ER when there is heterogeneity in latent positions.")
println("The key differentiator is the VARIANCE in connection probabilities.")
