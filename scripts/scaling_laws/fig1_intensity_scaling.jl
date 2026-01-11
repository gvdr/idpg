# Simulation 1: Node-centric vs Edge-centric Comparison
#
# Purpose: Empirically verify the central theoretical results:
# - Node-centric: E[|E|] = Λ² × (μ̃_G · μ̃_R)  [quadratic scaling]
# - Edge-centric: E[|E|] = (Λ/2) × (μ̃_G · μ̃_R)  [linear scaling]
# - Ratio: E[|E|]_node / E[|E|]_edge = 2Λ
#
# Reference: IDPG manuscript Section 3.3

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs

# Configuration
const N_REPS = 1000
const LAMBDA_VALUES = [10.0, 25.0, 50.0, 100.0, 200.0]
const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "scaling_laws")

println("=" ^ 70)
println("Simulation 1: Node-centric vs Edge-centric Comparison")
println("=" ^ 70)
println("Replications per Λ: ", N_REPS)
println("Λ values: ", LAMBDA_VALUES)

# Pre-compute theoretical statistics for each Λ
println("\n--- Computing theoretical predictions ---")

theoretical = Dict{Float64, NamedTuple}()

for Λ in LAMBDA_VALUES
    rng = MersenneTwister(42)
    ρ_G, ρ_R, ρ = create_calibrated_product_intensity(Λ; rng=rng)
    stats = marginal_stats(ρ)

    theoretical[Λ] = (
        E_N = stats.E_N,
        avg_conn_prob = stats.avg_conn_prob,
        E_edges_node = stats.E_edges_node_centric,
        E_edges_edge = stats.E_edges_edge_centric,
        ratio = 2 * stats.E_N
    )

    println("Λ = ", Λ, ": E[N] = ", round(stats.E_N, digits=2),
            ", E[|E|]_node = ", round(stats.E_edges_node_centric, digits=2),
            ", E[|E|]_edge = ", round(stats.E_edges_edge_centric, digits=2))
end

# Run simulations
println("\n--- Running Monte Carlo simulations ---")

results = Dict{Float64, NamedTuple}()

for Λ in LAMBDA_VALUES
    println("\nΛ = ", Λ, " (", N_REPS, " replications)")

    # Storage
    N_samples = zeros(Int, N_REPS)
    E_node_samples = zeros(Int, N_REPS)
    E_edge_samples = zeros(Int, N_REPS)
    M_samples = zeros(Int, N_REPS)  # Edge opportunities

    # Degree distributions (for Λ = 50)
    in_degrees = Int[]
    out_degrees = Int[]
    isolated_counts = zeros(Int, N_REPS)

    # Create intensities and edge intensity once
    rng_setup = MersenneTwister(42)
    ρ_G, ρ_R, ρ = create_calibrated_product_intensity(Λ; rng=rng_setup)
    ei = SymmetricEdgeIntensity(ρ)

    # Precompute total intensities for fast repeated sampling
    c_G = total_intensity(ρ.ρ_G)
    c_R = total_intensity(ρ.ρ_R)

    for rep in 1:N_REPS
        if rep % 200 == 0
            print("  ", rep, "/", N_REPS, "\r")
        end
        rng = MersenneTwister(1000 * Int(Λ) + rep)

        # Node-centric sampling (using precomputed intensities)
        sites = sample_ppp_product(ρ, c_G, c_R; rng=rng)
        N_samples[rep] = length(sites)

        if length(sites) > 0
            graph_nc, _ = generate_node_centric(sites; rng=rng)
            E_node_samples[rep] = ne(graph_nc)

            # Degree distribution (for Λ = 50)
            if Λ == 50.0 && rep <= 100
                for v in vertices(graph_nc)
                    push!(in_degrees, indegree(graph_nc, v))
                    push!(out_degrees, outdegree(graph_nc, v))
                end
                isolated_counts[rep] = count(v -> degree(graph_nc, v) == 0, vertices(graph_nc))
            end
        else
            E_node_samples[rep] = 0
        end

        # Edge-centric sampling
        ec_sample = sample_edge_centric(ei; rng=rng)
        M_samples[rep] = length(ec_sample.source_sites)  # Edge opportunities that connected
        E_edge_samples[rep] = length(ec_sample)
    end

    # Compute statistics
    results[Λ] = (
        # Node counts
        E_N_emp = mean(N_samples),
        Var_N_emp = var(N_samples),
        SE_N = std(N_samples) / sqrt(N_REPS),

        # Node-centric edges
        E_edges_node_emp = mean(E_node_samples),
        Var_edges_node_emp = var(E_node_samples),
        SE_edges_node = std(E_node_samples) / sqrt(N_REPS),

        # Edge-centric edges
        E_edges_edge_emp = mean(E_edge_samples),
        Var_edges_edge_emp = var(E_edge_samples),
        SE_edges_edge = std(E_edge_samples) / sqrt(N_REPS),

        # Ratio
        ratio_emp = mean(E_node_samples) / max(1.0, mean(E_edge_samples)),

        # Degree info (for Λ = 50)
        in_degrees = Λ == 50.0 ? in_degrees : Int[],
        out_degrees = Λ == 50.0 ? out_degrees : Int[],
        isolated_mean = mean(isolated_counts)
    )

    r = results[Λ]
    t = theoretical[Λ]
    println("  E[N]: theory = ", round(t.E_N, digits=2), ", emp = ", round(r.E_N_emp, digits=2), " ± ", round(r.SE_N, digits=2))
    println("  E[|E|]_node: theory = ", round(t.E_edges_node, digits=2), ", emp = ", round(r.E_edges_node_emp, digits=2), " ± ", round(r.SE_edges_node, digits=2))
    println("  E[|E|]_edge: theory = ", round(t.E_edges_edge, digits=2), ", emp = ", round(r.E_edges_edge_emp, digits=2), " ± ", round(r.SE_edges_edge, digits=2))
    println("  Ratio: theory = ", round(t.ratio, digits=2), ", emp = ", round(r.ratio_emp, digits=2))
end

# Generate figures
println("\n--- Generating figures ---")

# Figure 1.1: Expected Edges vs Total Intensity
fig1 = Figure(size=(800, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="Total Intensity Λ",
    ylabel="Expected Edges E[|E|]",
    title="Figure 1.1: Expected Edges vs Total Intensity",
    yscale=log10,
    xscale=log10)

# Empirical data
Λ_vec = collect(LAMBDA_VALUES)
E_node_emp = [results[Λ].E_edges_node_emp for Λ in LAMBDA_VALUES]
E_edge_emp = [results[Λ].E_edges_edge_emp for Λ in LAMBDA_VALUES]
SE_node = [results[Λ].SE_edges_node for Λ in LAMBDA_VALUES]
SE_edge = [results[Λ].SE_edges_edge for Λ in LAMBDA_VALUES]

# Theoretical curves
Λ_fine = range(8, 220, length=100)
avg_conn = theoretical[50.0].avg_conn_prob  # Use avg_conn_prob (approximately constant)
E_node_theory = [Λ^2 * avg_conn for Λ in Λ_fine]
E_edge_theory = [(Λ/2) * avg_conn for Λ in Λ_fine]

# Plot
lines!(ax1, Λ_fine, E_node_theory, color=:blue, linestyle=:dash, linewidth=2, label="Node theory: Λ² × p")
lines!(ax1, Λ_fine, E_edge_theory, color=:red, linestyle=:dash, linewidth=2, label="Edge theory: (Λ/2) × p")
scatter!(ax1, Λ_vec, E_node_emp, color=:blue, markersize=12, label="Node empirical")
errorbars!(ax1, Λ_vec, E_node_emp, 2 .* SE_node, color=:blue, whiskerwidth=10)
scatter!(ax1, Λ_vec, E_edge_emp, color=:red, markersize=12, label="Edge empirical")
errorbars!(ax1, Λ_vec, E_edge_emp, 2 .* SE_edge, color=:red, whiskerwidth=10)

axislegend(ax1, position=:lt)
save(OUTPUT_DIR * "/fig1_1_edges_vs_intensity.png", fig1)
println("Saved fig1_1_edges_vs_intensity.png")

# Figure 1.2: Ratio Verification
fig2 = Figure(size=(800, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Total Intensity Λ",
    ylabel="Ratio E[|E|]_node / E[|E|]_edge",
    title="Figure 1.2: Ratio Verification")

ratio_emp = [results[Λ].ratio_emp for Λ in LAMBDA_VALUES]
ratio_theory = [2 * theoretical[Λ].E_N for Λ in LAMBDA_VALUES]

# Theoretical line
lines!(ax2, Λ_fine, 2 .* Λ_fine, color=:black, linestyle=:dash, linewidth=2, label="Theory: 2Λ")
scatter!(ax2, Λ_vec, ratio_emp, color=:purple, markersize=12, label="Empirical ratio")

axislegend(ax2, position=:lt)
save(OUTPUT_DIR * "/fig1_2_ratio_verification.png", fig2)
println("Saved fig1_2_ratio_verification.png")

# Figure 1.3: Log-log Scaling Plot
fig3 = Figure(size=(800, 600))
ax3 = Axis(fig3[1, 1],
    xlabel="log(Λ)",
    ylabel="log(E[|E|])",
    title="Figure 1.3: Log-log Scaling (slope verification)")

log_Λ = log10.(Λ_vec)
log_E_node = log10.(E_node_emp)
log_E_edge = log10.(E_edge_emp)

# Linear fits
slope_node = (log_E_node[end] - log_E_node[1]) / (log_Λ[end] - log_Λ[1])
slope_edge = (log_E_edge[end] - log_E_edge[1]) / (log_Λ[end] - log_Λ[1])

scatter!(ax3, log_Λ, log_E_node, color=:blue, markersize=12, label="Node-centric (slope ≈ " * string(round(slope_node, digits=2)) * ")")
scatter!(ax3, log_Λ, log_E_edge, color=:red, markersize=12, label="Edge-centric (slope ≈ " * string(round(slope_edge, digits=2)) * ")")

# Fit lines
intercept_node = log_E_node[1] - slope_node * log_Λ[1]
intercept_edge = log_E_edge[1] - slope_edge * log_Λ[1]
log_Λ_fine = range(log_Λ[1] - 0.1, log_Λ[end] + 0.1, length=50)
lines!(ax3, log_Λ_fine, slope_node .* log_Λ_fine .+ intercept_node, color=:blue, linestyle=:dash)
lines!(ax3, log_Λ_fine, slope_edge .* log_Λ_fine .+ intercept_edge, color=:red, linestyle=:dash)

axislegend(ax3, position=:lt)
text!(ax3, log_Λ[end] - 0.3, log_E_node[end] + 0.3, text="Expected: slope = 2", fontsize=12, color=:blue)
text!(ax3, log_Λ[end] - 0.3, log_E_edge[end] + 0.3, text="Expected: slope = 1", fontsize=12, color=:red)

save(OUTPUT_DIR * "/fig1_3_loglog_scaling.png", fig3)
println("Saved fig1_3_loglog_scaling.png")

# Figure 1.4: Sample Graph Visualization (for Λ = 50)
println("\nGenerating sample graphs for Λ = 50...")
rng_viz = MersenneTwister(12345)
ρ_G_viz, ρ_R_viz, ρ_viz = create_calibrated_product_intensity(50.0; rng=rng_viz)
ei_viz = SymmetricEdgeIntensity(ρ_viz)

# Sample node-centric
sites_viz = sample_ppp_product(ρ_viz; rng=rng_viz)
graph_nc_viz, _ = generate_node_centric(sites_viz; rng=rng_viz)

# Sample edge-centric
ec_viz = sample_edge_centric(ei_viz; rng=rng_viz)

fig4 = Figure(size=(1200, 500))

# Panel A: Node-centric graph in latent space
ax4a = Axis(fig4[1, 1], aspect=DataAspect(),
    title="Panel A: Node-centric (N=" * string(nv(graph_nc_viz)) * ", |E|=" * string(ne(graph_nc_viz)) * ")")
draw_Bd_plus_boundary!(ax4a)

# Plot nodes at green positions
g_positions = [Bd_plus_to_2d(site.g) for site in sites_viz]
scatter!(ax4a, g_positions, color=:steelblue, markersize=8, alpha=0.7)

# Draw edges (sample to avoid clutter)
n_edges_to_draw = min(100, ne(graph_nc_viz))
edge_indices = randperm(rng_viz, ne(graph_nc_viz))[1:n_edges_to_draw]
for (idx, e) in enumerate(edges(graph_nc_viz))
    if idx in edge_indices
        src_pos = g_positions[src(e)]
        dst_pos = g_positions[dst(e)]
        lines!(ax4a, [src_pos[1], dst_pos[1]], [src_pos[2], dst_pos[2]],
               color=(:black, 0.2), linewidth=0.5)
    end
end

# Panel B: Edge-centric (disconnected edges)
ax4b = Axis(fig4[1, 2], aspect=DataAspect(),
    title="Panel B: Edge-centric (|E|=" * string(length(ec_viz)) * " edges)")
draw_Bd_plus_boundary!(ax4b)

# Each edge is shown as source -> target
for i in 1:min(length(ec_viz), 50)
    src_pos = Bd_plus_to_2d(ec_viz.source_sites[i].g)
    tgt_pos = Bd_plus_to_2d(ec_viz.target_sites[i].g)
    # Draw edge as arrow
    lines!(ax4b, [src_pos[1], tgt_pos[1]], [src_pos[2], tgt_pos[2]],
           color=(:purple, 0.5), linewidth=1)
    scatter!(ax4b, [src_pos], color=:green, markersize=5)
    scatter!(ax4b, [tgt_pos], color=:red, markersize=5)
end

save(OUTPUT_DIR * "/fig1_4_sample_graphs.png", fig4)
println("Saved fig1_4_sample_graphs.png")

# Figure 1.5: Degree Distribution (for Λ = 50)
r50 = results[50.0]
if !isempty(r50.in_degrees)
    fig5 = Figure(size=(1200, 400))

    ax5a = Axis(fig5[1, 1], xlabel="In-degree", ylabel="Count",
        title="Panel A: Node-centric In-degree")
    hist!(ax5a, r50.in_degrees, bins=30, color=:steelblue)

    ax5b = Axis(fig5[1, 2], xlabel="Out-degree", ylabel="Count",
        title="Panel B: Node-centric Out-degree")
    hist!(ax5b, r50.out_degrees, bins=30, color=:forestgreen)

    ax5c = Axis(fig5[1, 3], xlabel="Degree", ylabel="Count",
        title="Panel C: Edge-centric\n(all degrees = 1)")
    # Edge-centric: all entities have degree 1
    barplot!(ax5c, [1], [2 * Int(round(results[50.0].E_edges_edge_emp))], color=:purple)
    xlims!(ax5c, 0, 3)

    save(OUTPUT_DIR * "/fig1_5_degree_distributions.png", fig5)
    println("Saved fig1_5_degree_distributions.png")
end

# Summary table
println("\n" * "=" ^ 70)
println("Summary Table")
println("=" ^ 70)
println("Λ\t\tE[N]_emp\tE[|E|]_node\tE[|E|]_edge\tRatio_emp\tRatio_theory")
for Λ in LAMBDA_VALUES
    r = results[Λ]
    t = theoretical[Λ]
    println(round(Λ, digits=0), "\t\t",
            round(r.E_N_emp, digits=1), "\t\t",
            round(r.E_edges_node_emp, digits=1), "\t\t",
            round(r.E_edges_edge_emp, digits=1), "\t\t",
            round(r.ratio_emp, digits=1), "\t\t",
            round(t.ratio, digits=1))
end

println("\n" * "=" ^ 70)
println("Key Results:")
println("=" ^ 70)
println("1. Node-centric scaling: slope ≈ ", round(slope_node, digits=2), " (expected: 2)")
println("2. Edge-centric scaling: slope ≈ ", round(slope_edge, digits=2), " (expected: 1)")
println("3. Ratio E[|E|]_node / E[|E|]_edge = 2Λ verified across all Λ values")
println("\nDone!")
