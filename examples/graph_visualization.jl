# Graph Visualization Example
# Shows the actual realized graphs from IDPG sampling
# Uses directed edges (arrows) to show the asymmetric nature of connections

using IDPG
using Random
using CairoMakie
using GraphMakie
using Graphs
using NetworkLayout
using LinearAlgebra: norm

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Graph Visualization Example")
println("=" ^ 60)

# Medium density case: ~30-50 nodes, average degree ~5-10
# Use interior orthogonal means for lower connectivity

println("\n--- Setting Up Intensity ---")

# G concentrated near (0.65, 0.15) - interior position
ρ_G = BdPlusMixture(
    [1.0],
    [[0.65, 0.15]],
    [40.0],  # High concentration for tighter clustering
    80.0
)

# R concentrated near (0.15, 0.65) - nearly orthogonal to G
# (0.65, 0.15) · (0.15, 0.65) = 0.0975 + 0.0975 = 0.195
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
avg_degree = stats.E_edges_node_centric / stats.E_N
println("Expected avg out-degree: ", round(avg_degree, digits=2))

# Sample sites
println("\n--- Sampling ---")
sites = sample_ppp_product(ρ; rng=rng)
println("Sampled ", length(sites), " interaction sites (nodes)")

# Generate node-centric graph (directed!)
graph_nc, _ = generate_node_centric(sites; rng=rng)
println("Node-centric graph: ", nv(graph_nc), " nodes, ", ne(graph_nc), " directed edges")
actual_avg_degree = ne(graph_nc) / nv(graph_nc)
println("Actual avg out-degree: ", round(actual_avg_degree, digits=2))

# Generate edge-centric sample
sample_ec = generate_edge_centric(sites; rng=rng)
println("Edge-centric edges: ", length(sample_ec))

# Create visualization
println("\n--- Generating Visualizations ---")

fig = Figure(size=(1600, 800))

# Panel 1: Node positions in latent space (G coordinates)
ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Node Positions in G-space (source propensity)")
draw_Bd_plus_boundary!(ax1)
g_points = [Bd_plus_to_2d(site.g) for site in sites]
scatter!(ax1, g_points, color=:forestgreen, markersize=12)
text!(ax1, Point2f(1.05, 0.0), text="g₁", fontsize=12, align=(:left, :center))
text!(ax1, Point2f(0.0, 1.05), text="g₂", fontsize=12, align=(:center, :bottom))

# Panel 2: Node positions in latent space (R coordinates)
ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Node Positions in R-space (target propensity)")
draw_Bd_plus_boundary!(ax2)
r_points = [Bd_plus_to_2d(site.r) for site in sites]
scatter!(ax2, r_points, color=:firebrick, markersize=12)
text!(ax2, Point2f(1.05, 0.0), text="r₁", fontsize=12, align=(:left, :center))
text!(ax2, Point2f(0.0, 1.05), text="r₂", fontsize=12, align=(:center, :bottom))

# Panel 3: Directed graph with spring layout and ARROWS
ax3 = Axis(fig[1, 3],
    title="Directed Graph (Spring Layout)\n" * string(nv(graph_nc)) * " nodes, " * string(ne(graph_nc)) * " edges, avg degree=" * string(round(actual_avg_degree, digits=1)))
hidedecorations!(ax3)
hidespines!(ax3)

if nv(graph_nc) > 0 && ne(graph_nc) > 0
    layout = Spring(; iterations=200)

    # Use the directed graph directly with arrows
    graphplot!(ax3, graph_nc,
        layout=layout,
        node_size=20,
        node_color=:steelblue,
        edge_color=(:gray60, 0.7),
        edge_width=1.5,
        arrow_size=15,
        arrow_show=true)
end

# Panel 4: Directed graph embedded in latent G-space with arrows
ax4 = Axis(fig[2, 1], aspect=DataAspect(),
    title="Directed Graph in G-space\nArrows: i → j if edge exists")
draw_Bd_plus_boundary!(ax4)

# Draw directed edges with arrows
if ne(graph_nc) > 0
    for e in edges(graph_nc)
        i, j = src(e), dst(e)
        p1, p2 = g_points[i], g_points[j]
        # Draw arrow from p1 to p2
        arrows!(ax4, [p1[1]], [p1[2]], [p2[1] - p1[1]], [p2[2] - p1[2]],
            color=(:purple, 0.3), linewidth=0.8, arrowsize=8)
    end
end
scatter!(ax4, g_points, color=:steelblue, markersize=10)

# Panel 5: Show a few example edges with their connection probabilities
ax5 = Axis(fig[2, 2], aspect=DataAspect(),
    title="Example: Why edges form\nP(i→j) = gᵢ · rⱼ")
draw_Bd_plus_boundary!(ax5)

# Show G and R positions for a few nodes
n_show = min(5, length(sites))
colors = [:red, :blue, :green, :orange, :purple]
for i in 1:n_show
    scatter!(ax5, [g_points[i]], color=colors[i], markersize=15, marker=:circle)
    scatter!(ax5, [r_points[i]], color=colors[i], markersize=15, marker=:diamond)
    # Draw line connecting g and r for same node
    lines!(ax5, [g_points[i], r_points[i]], color=(colors[i], 0.3), linestyle=:dash, linewidth=1)
end
text!(ax5, Point2f(0.5, -0.1), text="● = g position, ◆ = r position", fontsize=10, align=(:center, :top))

# Panel 6: Edge-centric visualization
ax6 = Axis(fig[2, 3], aspect=DataAspect(),
    title="Edge-Centric Sample\n" * string(length(sample_ec)) * " interaction events")
draw_Bd_plus_boundary!(ax6)
if length(sample_ec) > 0
    src_2d = [Bd_plus_to_2d(s) for s in sample_ec.sources]
    tgt_2d = [Bd_plus_to_2d(t) for t in sample_ec.targets]

    # Draw arrows from source to target positions
    for k in 1:length(sample_ec)
        arrows!(ax6, [src_2d[k][1]], [src_2d[k][2]],
            [tgt_2d[k][1] - src_2d[k][1]], [tgt_2d[k][2] - src_2d[k][2]],
            color=(:gray, 0.4), linewidth=0.5, arrowsize=6)
    end
    scatter!(ax6, src_2d, color=:forestgreen, markersize=6, marker=:circle)
    scatter!(ax6, tgt_2d, color=:firebrick, markersize=6, marker=:diamond)
end

save("output/graphs/graph_visualization.png", fig)
println("Saved graph_visualization.png")

# Also create a cleaner focused plot of just the directed graph
println("\n--- Creating Focused Directed Graph Plot ---")

fig2 = Figure(size=(900, 900))
ax = Axis(fig2[1, 1],
    title="IDPG Directed Graph\n" * string(nv(graph_nc)) * " nodes, " * string(ne(graph_nc)) * " directed edges")
hidedecorations!(ax)
hidespines!(ax)

if nv(graph_nc) > 0 && ne(graph_nc) > 0
    # Color nodes by their G coordinate (how much they "send")
    node_colors = [norm(sites[i].g) for i in 1:length(sites)]

    layout = Spring(; iterations=300, C=2.0)
    graphplot!(ax, graph_nc,
        layout=layout,
        node_size=25,
        node_color=node_colors,
        node_colormap=:viridis,
        edge_color=(:gray50, 0.5),
        edge_width=1.0,
        arrow_size=12,
        arrow_show=true)

    Colorbar(fig2[1, 2], limits=extrema(node_colors), colormap=:viridis, label="||g|| (source strength)")
end

save("output/graphs/graph_directed.png", fig2)
println("Saved graph_directed.png")

println("\n" * "=" ^ 60)
println("Done!")
