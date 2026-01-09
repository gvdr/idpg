# Mesoscale Network Metrics Example
# Explores structural properties beyond global statistics (nodes, edges)
# Computes metrics across multiple network realizations to assess consistency

using IDPG
using Random
using Graphs
using Statistics
using CairoMakie
using LinearAlgebra: dot

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Mesoscale Network Metrics")
println("=" ^ 60)

# --- Setup: Medium density case ---
println("\n--- Setting Up Intensity ---")

# G concentrated near (0.65, 0.15)
ρ_G = BdPlusMixture(
    [1.0],
    [[0.65, 0.15]],
    [40.0],
    80.0
)

# R concentrated near (0.15, 0.65) - nearly orthogonal to G
ρ_R = BdPlusMixture(
    [1.0],
    [[0.15, 0.65]],
    [40.0],
    80.0
)

ρ = ProductIntensity(ρ_G, ρ_R)
stats = marginal_stats(ρ; rng=rng)

println("Expected nodes E[N]: ", round(stats.E_N, digits=2))
println("Expected edges E[|E|]: ", round(stats.E_edges_node_centric, digits=1))
println("Average connection prob: ", round(stats.avg_conn_prob, digits=3))

# --- Metric computation functions ---

"""
Compute reciprocity: fraction of edges where both i→j and j→i exist.
"""
function compute_reciprocity(g::SimpleDiGraph)
    if ne(g) == 0
        return 0.0
    end
    reciprocal_count = 0
    for e in edges(g)
        if has_edge(g, dst(e), src(e))
            reciprocal_count += 1
        end
    end
    return reciprocal_count / ne(g)
end

"""
Compute correlation between in-degree and out-degree across nodes.
"""
function degree_correlation(g::SimpleDiGraph)
    if nv(g) < 2
        return NaN
    end
    in_degs = Float64.(indegree(g))
    out_degs = Float64.(outdegree(g))

    # Pearson correlation
    if std(in_degs) == 0 || std(out_degs) == 0
        return 0.0
    end
    return cor(in_degs, out_degs)
end

"""
Compute directed clustering coefficient (transitivity).
For directed graphs: fraction of potential triangles that are closed.
"""
function directed_clustering(g::SimpleDiGraph)
    if nv(g) < 3
        return 0.0
    end

    triangles = 0
    triplets = 0

    for v in vertices(g)
        # Get all neighbors (in or out)
        out_neigh = Set(outneighbors(g, v))
        in_neigh = Set(inneighbors(g, v))
        all_neigh = union(out_neigh, in_neigh)

        # Count directed 2-paths through v and closed triangles
        for u in all_neigh
            for w in all_neigh
                if u != w
                    # Check if u-v-w forms a 2-path
                    # and if there's an edge u-w (in either direction)
                    triplets += 1
                    if has_edge(g, u, w) || has_edge(g, w, u)
                        triangles += 1
                    end
                end
            end
        end
    end

    if triplets == 0
        return 0.0
    end
    return triangles / triplets
end

"""
Compute degree assortativity for directed graph.
Measures if high out-degree nodes connect to high in-degree nodes.
"""
function directed_assortativity(g::SimpleDiGraph)
    if ne(g) == 0
        return NaN
    end

    # For each edge, get out-degree of source and in-degree of target
    src_out = Float64[]
    dst_in = Float64[]

    for e in edges(g)
        push!(src_out, outdegree(g, src(e)))
        push!(dst_in, indegree(g, dst(e)))
    end

    if std(src_out) == 0 || std(dst_in) == 0
        return 0.0
    end
    return cor(src_out, dst_in)
end

"""
Compute statistics about strongly connected components.
"""
function scc_stats(g::SimpleDiGraph)
    sccs = strongly_connected_components(g)
    sizes = [length(scc) for scc in sccs]

    return (
        n_components = length(sccs),
        largest_scc = maximum(sizes),
        largest_scc_fraction = maximum(sizes) / nv(g),
        mean_scc_size = mean(sizes),
        n_singleton = count(s -> s == 1, sizes)
    )
end

"""
Compute PageRank and return statistics about its distribution.
"""
function pagerank_stats(g::SimpleDiGraph)
    if nv(g) == 0
        return (mean=NaN, std=NaN, max=NaN, gini=NaN)
    end

    pr = pagerank(g)

    # Gini coefficient for inequality
    sorted_pr = sort(pr)
    n = length(sorted_pr)
    gini = (2 * sum((1:n) .* sorted_pr) - (n + 1) * sum(sorted_pr)) / (n * sum(sorted_pr))

    return (
        mean = mean(pr),
        std = std(pr),
        max = maximum(pr),
        gini = gini
    )
end

# --- Run simulations ---
println("\n--- Running Simulations ---")

n_realizations = 50

# Storage for metrics
metrics = (
    n_nodes = Int[],
    n_edges = Int[],
    reciprocity = Float64[],
    in_out_correlation = Float64[],
    clustering = Float64[],
    assortativity = Float64[],
    n_scc = Int[],
    largest_scc_frac = Float64[],
    pagerank_gini = Float64[],
    avg_in_degree = Float64[],
    avg_out_degree = Float64[],
    in_degree_std = Float64[],
    out_degree_std = Float64[]
)

for i in 1:n_realizations
    # Sample sites and generate graph
    sites = sample_ppp_product(ρ; rng=rng)

    if length(sites) < 2
        continue
    end

    graph, _ = generate_node_centric(sites; rng=rng)

    # Basic counts
    push!(metrics.n_nodes, nv(graph))
    push!(metrics.n_edges, ne(graph))

    # Mesoscale metrics
    push!(metrics.reciprocity, compute_reciprocity(graph))
    push!(metrics.in_out_correlation, degree_correlation(graph))
    push!(metrics.clustering, directed_clustering(graph))
    push!(metrics.assortativity, directed_assortativity(graph))

    # SCC stats
    scc = scc_stats(graph)
    push!(metrics.n_scc, scc.n_components)
    push!(metrics.largest_scc_frac, scc.largest_scc_fraction)

    # PageRank stats
    pr = pagerank_stats(graph)
    push!(metrics.pagerank_gini, pr.gini)

    # Degree statistics
    in_degs = indegree(graph)
    out_degs = outdegree(graph)
    push!(metrics.avg_in_degree, mean(in_degs))
    push!(metrics.avg_out_degree, mean(out_degs))
    push!(metrics.in_degree_std, std(in_degs))
    push!(metrics.out_degree_std, std(out_degs))

    if i % 10 == 0
        println("  Completed ", i, "/", n_realizations, " realizations")
    end
end

# --- Report results ---
println("\n" * "=" ^ 60)
println("RESULTS SUMMARY (", length(metrics.n_nodes), " valid realizations)")
println("=" ^ 60)

println("\n--- Global Metrics ---")
println("Nodes:  ", round(mean(metrics.n_nodes), digits=1), " ± ", round(std(metrics.n_nodes), digits=1),
    "  (theory: ", round(stats.E_N, digits=1), ")")
println("Edges:  ", round(mean(metrics.n_edges), digits=1), " ± ", round(std(metrics.n_edges), digits=1),
    "  (theory: ", round(stats.E_edges_node_centric, digits=1), ")")

println("\n--- Reciprocity ---")
println("Observed:  ", round(mean(metrics.reciprocity), digits=3), " ± ", round(std(metrics.reciprocity), digits=3))
# Expected reciprocity: P(j→i | i→j) = P(j→i) = avg_conn_prob (independence)
println("Expected:  ", round(stats.avg_conn_prob, digits=3), "  (= avg connection probability)")

println("\n--- In/Out Degree Correlation ---")
println("Observed:  ", round(mean(metrics.in_out_correlation), digits=3), " ± ", round(std(metrics.in_out_correlation), digits=3))
println("Expected:  ~0  (g and r are independent in product intensity)")

println("\n--- Clustering Coefficient ---")
println("Observed:  ", round(mean(metrics.clustering), digits=3), " ± ", round(std(metrics.clustering), digits=3))

println("\n--- Degree Assortativity ---")
println("Observed:  ", round(mean(metrics.assortativity), digits=3), " ± ", round(std(metrics.assortativity), digits=3))

println("\n--- Strongly Connected Components ---")
println("Number of SCCs:      ", round(mean(metrics.n_scc), digits=1), " ± ", round(std(metrics.n_scc), digits=1))
println("Largest SCC fraction: ", round(mean(metrics.largest_scc_frac), digits=3), " ± ", round(std(metrics.largest_scc_frac), digits=3))

println("\n--- PageRank Inequality (Gini) ---")
println("Gini coefficient:  ", round(mean(metrics.pagerank_gini), digits=3), " ± ", round(std(metrics.pagerank_gini), digits=3))
println("  (0 = equal, 1 = maximally unequal)")

println("\n--- Degree Distributions ---")
println("In-degree:   mean=", round(mean(metrics.avg_in_degree), digits=2), ", std=", round(mean(metrics.in_degree_std), digits=2))
println("Out-degree:  mean=", round(mean(metrics.avg_out_degree), digits=2), ", std=", round(mean(metrics.out_degree_std), digits=2))

# --- Visualizations ---
println("\n--- Creating Visualizations ---")

fig = Figure(size=(1400, 1000))

# Panel 1: Reciprocity distribution
ax1 = Axis(fig[1, 1], xlabel="Reciprocity", ylabel="Count", title="Reciprocity Distribution")
hist!(ax1, metrics.reciprocity, bins=15, color=:steelblue)
vlines!(ax1, [stats.avg_conn_prob], color=:red, linestyle=:dash, linewidth=2, label="Expected")
axislegend(ax1, position=:rt)

# Panel 2: In/Out degree correlation
ax2 = Axis(fig[1, 2], xlabel="In-Out Degree Correlation", ylabel="Count", title="Degree Independence Test")
hist!(ax2, metrics.in_out_correlation, bins=15, color=:forestgreen)
vlines!(ax2, [0.0], color=:red, linestyle=:dash, linewidth=2, label="Expected (0)")
axislegend(ax2, position=:rt)

# Panel 3: Clustering coefficient
ax3 = Axis(fig[1, 3], xlabel="Clustering Coefficient", ylabel="Count", title="Clustering Distribution")
hist!(ax3, metrics.clustering, bins=15, color=:purple)

# Panel 4: Largest SCC fraction
ax4 = Axis(fig[2, 1], xlabel="Largest SCC / N", ylabel="Count", title="Giant SCC Fraction")
hist!(ax4, metrics.largest_scc_frac, bins=15, color=:orange)

# Panel 5: Assortativity
ax5 = Axis(fig[2, 2], xlabel="Assortativity", ylabel="Count", title="Degree Assortativity")
hist!(ax5, metrics.assortativity, bins=15, color=:teal)

# Panel 6: PageRank Gini
ax6 = Axis(fig[2, 3], xlabel="Gini Coefficient", ylabel="Count", title="PageRank Inequality")
hist!(ax6, metrics.pagerank_gini, bins=15, color=:crimson)

# Panel 7: Example degree distribution from one realization
# Generate one more graph for visualization
sites_example = sample_ppp_product(ρ; rng=rng)
graph_example, _ = generate_node_centric(sites_example; rng=rng)

ax7 = Axis(fig[3, 1], xlabel="Degree", ylabel="Count", title="Example In-Degree Distribution")
hist!(ax7, indegree(graph_example), bins=20, color=:steelblue, label="In-degree")

ax8 = Axis(fig[3, 2], xlabel="Degree", ylabel="Count", title="Example Out-Degree Distribution")
hist!(ax8, outdegree(graph_example), bins=20, color=:forestgreen, label="Out-degree")

# Panel 9: In vs Out degree scatter for one realization
ax9 = Axis(fig[3, 3], xlabel="Out-Degree", ylabel="In-Degree",
    title="In vs Out Degree (one realization)\nr = " * string(round(cor(Float64.(outdegree(graph_example)), Float64.(indegree(graph_example))), digits=3)))
scatter!(ax9, outdegree(graph_example), indegree(graph_example), color=:navy, markersize=8, alpha=0.6)
# Add diagonal line
lims = (0, max(maximum(outdegree(graph_example)), maximum(indegree(graph_example))) + 1)
lines!(ax9, [lims[1], lims[2]], [lims[1], lims[2]], color=:red, linestyle=:dash, linewidth=1)

save("output/mesoscale/mesoscale_metrics.png", fig)
println("Saved mesoscale_metrics.png")

# --- Additional: Compare across different intensity parameters ---
println("\n" * "=" ^ 60)
println("Comparing mesoscale structure across different intensities...")
println("=" ^ 60)

# Test cases with different connectivity patterns
comparison_cases = [
    (name="Aligned (high conn)", mean_G=[0.5, 0.5], mean_R=[0.5, 0.5], conc=20.0, scale=60.0),
    (name="Orthogonal (medium)", mean_G=[0.65, 0.15], mean_R=[0.15, 0.65], conc=40.0, scale=80.0),
    (name="Very orthogonal (low)", mean_G=[0.7, 0.1], mean_R=[0.1, 0.7], conc=60.0, scale=100.0),
]

n_compare = 30

fig2 = Figure(size=(1200, 400))

for (idx, case) in enumerate(comparison_cases)
    ρ_G_c = BdPlusMixture([1.0], [case.mean_G], [case.conc], case.scale)
    ρ_R_c = BdPlusMixture([1.0], [case.mean_R], [case.conc], case.scale)
    ρ_c = ProductIntensity(ρ_G_c, ρ_R_c)
    stats_c = marginal_stats(ρ_c; rng=rng)

    recip_vals = Float64[]
    cluster_vals = Float64[]

    for _ in 1:n_compare
        sites_c = sample_ppp_product(ρ_c; rng=rng)
        if length(sites_c) < 2
            continue
        end
        graph_c, _ = generate_node_centric(sites_c; rng=rng)
        push!(recip_vals, compute_reciprocity(graph_c))
        push!(cluster_vals, directed_clustering(graph_c))
    end

    println("\n", case.name, ":")
    println("  Avg conn prob: ", round(stats_c.avg_conn_prob, digits=3))
    println("  Reciprocity:   ", round(mean(recip_vals), digits=3), " (expected: ", round(stats_c.avg_conn_prob, digits=3), ")")
    println("  Clustering:    ", round(mean(cluster_vals), digits=3))

    ax = Axis(fig2[1, idx], xlabel="Reciprocity", ylabel="Count",
        title=case.name * "\nconn=" * string(round(stats_c.avg_conn_prob, digits=2)))
    hist!(ax, recip_vals, bins=12, color=[:steelblue, :forestgreen, :purple][idx])
    vlines!(ax, [stats_c.avg_conn_prob], color=:red, linestyle=:dash, linewidth=2)
end

save("output/mesoscale/mesoscale_comparison.png", fig2)
println("\nSaved mesoscale_comparison.png")

println("\n" * "=" ^ 60)
println("Done!")
