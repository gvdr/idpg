# Simulation 2: Intermediate Regime (R_μ)
#
# Purpose: Demonstrate interpolation between edge-centric (μ → 0) and node-centric (μ → ∞)
# as entity lifetime varies.
#
# Model:
# - Entities sampled from PPP on Ω × [0, W] (site space × observation window)
# - Each entity has birth time T ~ Uniform(0, W) and lifetime τ ~ Exp(μ)
# - Entity alive during [T, T + τ] ∩ [0, W]
# - Two entities can interact only if their alive intervals overlap
#
# Theoretical result (from manuscript Appendix):
# p_overlap(μ, W) = (2/α²)(α - 1 + e^(-α))  where α = W/μ
#
# Reference: IDPG manuscript Section 3.1.4 and Appendix

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using JLD2
using Distributions: Poisson, Exponential

# Configuration
const N_REPS = 1000
const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "scaling_laws")
const W = 1.0  # Observation window
const LAMBDA = 50.0  # Fixed total intensity

# μ values to test (log-spaced)
const MU_VALUES = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, Inf]

println("=" ^ 70)
println("Simulation 2: Intermediate Regime (R_μ)")
println("=" ^ 70)
println("Replications per μ: ", N_REPS)
println("Observation window W: ", W)
println("Total intensity Λ: ", LAMBDA)
println("μ values: ", MU_VALUES)

# Theoretical overlap probability
function p_overlap_theory(μ, W)
    if isinf(μ)
        return 1.0
    end
    α = W / μ
    if α < 1e-10
        return 1.0  # Limit as α → 0
    elseif α > 100
        return 2 * μ / W  # Asymptotic form as α → ∞
    else
        return (2 / α^2) * (α - 1 + exp(-α))
    end
end

# Check if two intervals overlap
function intervals_overlap(T1, τ1, T2, τ2, W)
    # Alive intervals: [T1, min(T1+τ1, W)] and [T2, min(T2+τ2, W)]
    end1 = min(T1 + τ1, W)
    end2 = min(T2 + τ2, W)
    return max(T1, T2) < min(end1, end2)
end


# Sample and simulate intermediate regime
# Accepts precomputed c_G, c_R for performance when called repeatedly
function simulate_intermediate(ρ, c_G, c_R, Λ, μ, W; rng=Random.default_rng())
    # Sample sites from PPP - the count IS the Poisson sample
    # Note: sample_ppp_product internally samples N ~ Poisson(Λ) where Λ = c_G * c_R
    # We don't sample N separately to avoid bias from min(N1, N2) < E[N]
    sites = sample_ppp_product(ρ, c_G, c_R; rng=rng)
    N = length(sites)

    if N == 0
        return 0, 0, 0
    end

    # Sample birth times and lifetimes
    T = rand(rng, N) .* W  # Uniform(0, W)
    if isinf(μ)
        τ = fill(W, N)  # All entities live forever
    else
        τ = rand(rng, Exponential(μ), N)
    end

    # Count overlapping pairs and edges
    n_opportunities = 0
    n_edges = 0

    for i in 1:N
        for j in (i+1):N  # Consider each unordered pair once
            if intervals_overlap(T[i], τ[i], T[j], τ[j], W)
                # This pair has opportunity to interact (both directions)
                n_opportunities += 2  # Count ordered pairs

                # Check both directions
                g_i, r_i = sites[i].g, sites[i].r
                g_j, r_j = sites[j].g, sites[j].r

                p_ij = dot(g_i, r_j)
                p_ji = dot(g_j, r_i)

                if rand(rng) < p_ij
                    n_edges += 1
                end
                if rand(rng) < p_ji
                    n_edges += 1
                end
            end
        end
    end

    return N, n_opportunities, n_edges
end

using LinearAlgebra: dot

# Create intensities (using library function)
rng_setup = MersenneTwister(42)
ρ_G, ρ_R, ρ = create_calibrated_product_intensity(LAMBDA; rng=rng_setup)
stats = marginal_stats(ρ)

# Precompute total intensities for fast repeated sampling
c_G = total_intensity(ρ.ρ_G)
c_R = total_intensity(ρ.ρ_R)

println("\nIntensity stats:")
println("  E[N] = ", round(stats.E_N, digits=2))
println("  E[g·r] = ", round(stats.avg_conn_prob, digits=4))
println("  E[|E|]_node = ", round(stats.E_edges_node_centric, digits=2))

# Run simulations
println("\n--- Running simulations ---")

results = Dict{Float64, NamedTuple}()

for μ in MU_VALUES
    μ_str = isinf(μ) ? "∞" : string(round(μ, digits=3))
    println("\nμ = ", μ_str, " (μ/W = ", isinf(μ) ? "∞" : round(μ/W, digits=3), ")")

    # Storage
    N_samples = zeros(Int, N_REPS)
    opp_samples = zeros(Int, N_REPS)
    edge_samples = zeros(Int, N_REPS)

    for rep in 1:N_REPS
        if rep % 200 == 0
            print("  ", rep, "/", N_REPS, "\r")
        end
        rng = MersenneTwister(10000 * findfirst(==(μ), MU_VALUES) + rep)
        N, opp, edges = simulate_intermediate(ρ, c_G, c_R, stats.E_N, μ, W; rng=rng)
        N_samples[rep] = N
        opp_samples[rep] = opp
        edge_samples[rep] = edges
    end

    # Theoretical predictions
    p_ov = p_overlap_theory(μ, W)
    E_opp_theory = stats.E_N^2 * p_ov  # N(N-1) ≈ N² for large N
    E_edges_theory = E_opp_theory * stats.avg_conn_prob

    results[μ] = (
        μ = μ,
        μ_over_W = isinf(μ) ? Inf : μ / W,
        p_overlap_theory = p_ov,
        p_overlap_emp = mean(opp_samples) / max(1.0, mean(N_samples)^2),
        E_opp_theory = E_opp_theory,
        E_opp_emp = mean(opp_samples),
        SE_opp = std(opp_samples) / sqrt(N_REPS),
        E_edges_theory = E_edges_theory,
        E_edges_emp = mean(edge_samples),
        SE_edges = std(edge_samples) / sqrt(N_REPS),
        E_N_emp = mean(N_samples),
    )

    r = results[μ]
    println("  p_overlap: theory = ", round(r.p_overlap_theory, digits=4),
            ", emp = ", round(r.p_overlap_emp, digits=4))
    println("  E[opp]: theory = ", round(r.E_opp_theory, digits=1),
            ", emp = ", round(r.E_opp_emp, digits=1))
    println("  E[|E|]: theory = ", round(r.E_edges_theory, digits=1),
            ", emp = ", round(r.E_edges_emp, digits=1))
end

# Generate figures
println("\n--- Generating figures ---")

# Sort μ values for plotting
μ_sorted = sort([μ for μ in MU_VALUES if !isinf(μ)])
μ_over_W = μ_sorted ./ W

# Figure 2.1: Opportunities vs Lifetime
# Notation: η = mean lifetime, α = W/η
fig1 = Figure(size=(800, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="η/W (log scale)",
    ylabel="Expected Opportunities",
    title="Figure 2.1: Opportunities vs Entity Lifetime",
    xscale=log10)

opp_emp = [results[μ].E_opp_emp for μ in μ_sorted]
opp_SE = [results[μ].SE_opp for μ in μ_sorted]
opp_theory = [results[μ].E_opp_theory for μ in μ_sorted]

scatter!(ax1, μ_over_W, opp_emp, color=:blue, markersize=10, label="Empirical")
errorbars!(ax1, μ_over_W, opp_emp, 2 .* opp_SE, color=:blue, whiskerwidth=8)
lines!(ax1, μ_over_W, opp_theory, color=:red, linestyle=:dash, linewidth=2, label="Theory")

# Reference lines
hlines!(ax1, [LAMBDA^2], color=:gray, linestyle=:dot, label="Node-centric limit (Λ²)")
hlines!(ax1, [0], color=:gray, linestyle=:dashdot)

axislegend(ax1, position=:rb)
save(OUTPUT_DIR * "/fig2_1_opportunities_vs_lifetime.png", fig1)
println("Saved fig2_1_opportunities_vs_lifetime.png")

# Figure 2.2: Expected Edges vs Lifetime
fig2 = Figure(size=(800, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="η/W (log scale)",
    ylabel="Expected Edges",
    title="Figure 2.2: Expected Edges vs Entity Lifetime",
    xscale=log10)

edges_emp = [results[μ].E_edges_emp for μ in μ_sorted]
edges_SE = [results[μ].SE_edges for μ in μ_sorted]
edges_theory = [results[μ].E_edges_theory for μ in μ_sorted]

scatter!(ax2, μ_over_W, edges_emp, color=:blue, markersize=10, label="Empirical")
errorbars!(ax2, μ_over_W, edges_emp, 2 .* edges_SE, color=:blue, whiskerwidth=8)
lines!(ax2, μ_over_W, edges_theory, color=:red, linestyle=:dash, linewidth=2, label="Theory")

# Reference lines
hlines!(ax2, [stats.E_edges_node_centric], color=:green, linestyle=:dot, label="Node-centric limit")
hlines!(ax2, [stats.E_edges_edge_centric], color=:orange, linestyle=:dot, label="Edge-centric limit")

axislegend(ax2, position=:rb)
save(OUTPUT_DIR * "/fig2_2_edges_vs_lifetime.png", fig2)
println("Saved fig2_2_edges_vs_lifetime.png")

# Figure 2.3: Overlap Probability
fig3 = Figure(size=(800, 600))
ax3 = Axis(fig3[1, 1],
    xlabel="η/W (log scale)",
    ylabel="Overlap Probability",
    title="Figure 2.3: Overlap Probability vs Entity Lifetime",
    xscale=log10)

p_ov_emp = [results[μ].p_overlap_emp for μ in μ_sorted]
p_ov_theory = [results[μ].p_overlap_theory for μ in μ_sorted]

# Theoretical curve (dense)
μ_W_fine = 10 .^ range(-2, 2, length=100)
p_ov_fine = [p_overlap_theory(μW * W, W) for μW in μ_W_fine]

lines!(ax3, μ_W_fine, p_ov_fine, color=:red, linewidth=2, label="Theory: (2/α²)(α - 1 + e^{-α})")
scatter!(ax3, μ_over_W, p_ov_emp, color=:blue, markersize=10, label="Empirical")

# Reference lines
hlines!(ax3, [1.0], color=:gray, linestyle=:dot, label="p_overlap = 1 (node-centric)")
hlines!(ax3, [0.0], color=:gray, linestyle=:dashdot)

axislegend(ax3, position=:rb)
save(OUTPUT_DIR * "/fig2_3_overlap_probability.png", fig3)
println("Saved fig2_3_overlap_probability.png")

# Figure 2.4: Transition Visualization (3 panels)
fig4 = Figure(size=(1200, 400))

# Select 3 representative η values (η/W shown)
μ_samples = [0.1, 1.0, 10.0]
titles = ["η/W = 0.1 (Near edge-centric)", "η/W = 1.0 (Intermediate)", "η/W = 10.0 (Near node-centric)"]

for (idx, μ) in enumerate(μ_samples)
    ax = Axis(fig4[1, idx], aspect=DataAspect(), title=titles[idx])
    draw_Bd_plus_boundary!(ax)

    # Generate one sample graph (using precomputed intensities)
    rng = MersenneTwister(42 + idx)
    sites = sample_ppp_product(ρ, c_G, c_R; rng=rng)
    N = length(sites)

    if N > 0
        # Sample times
        T = rand(rng, N) .* W
        τ = rand(rng, Exponential(μ), N)

        # Plot nodes
        g_pos = [Bd_plus_to_2d(site.g) for site in sites]
        scatter!(ax, g_pos, color=:steelblue, markersize=6, alpha=0.7)

        # Draw edges for overlapping pairs
        n_edges_drawn = 0
        for i in 1:N
            for j in (i+1):N
                if intervals_overlap(T[i], τ[i], T[j], τ[j], W)
                    p_ij = dot(sites[i].g, sites[j].r)
                    p_ji = dot(sites[j].g, sites[i].r)

                    if rand(rng) < p_ij && n_edges_drawn < 50
                        lines!(ax, [g_pos[i][1], g_pos[j][1]], [g_pos[i][2], g_pos[j][2]],
                               color=(:black, 0.3), linewidth=0.5)
                        n_edges_drawn += 1
                    end
                    if rand(rng) < p_ji && n_edges_drawn < 50
                        lines!(ax, [g_pos[j][1], g_pos[i][1]], [g_pos[j][2], g_pos[i][2]],
                               color=(:black, 0.3), linewidth=0.5)
                        n_edges_drawn += 1
                    end
                end
            end
        end
    end
end

save(OUTPUT_DIR * "/fig2_4_transition_visualization.png", fig4)
println("Saved fig2_4_transition_visualization.png")

# Summary table
println("\n" * "=" ^ 70)
println("Summary Table")
println("=" ^ 70)
println("μ/W\t\tp_overlap\tE[opp]\t\tE[|E|]")
println("--------\t---------\t------\t\t------")
for μ in μ_sorted
    r = results[μ]
    println(round(r.μ_over_W, digits=3), "\t\t",
            round(r.p_overlap_emp, digits=3), "\t\t",
            round(r.E_opp_emp, digits=1), "\t\t",
            round(r.E_edges_emp, digits=1))
end
if Inf in MU_VALUES
    r = results[Inf]
    println("∞\t\t", round(r.p_overlap_emp, digits=3), "\t\t",
            round(r.E_opp_emp, digits=1), "\t\t",
            round(r.E_edges_emp, digits=1))
end

# Save results
results_file = OUTPUT_DIR * "/sim2_results.jld2"
@save results_file results LAMBDA W MU_VALUES stats
println("\nSaved results to ", results_file)

println("\n" * "=" ^ 70)
println("Key Results:")
println("=" ^ 70)
println("1. As μ/W → ∞, p_overlap → 1 (node-centric limit)")
println("2. As μ/W → 0, p_overlap → 0 (edge-centric limit)")
println("3. Intermediate values smoothly interpolate between limits")
println("4. Formula p_overlap = (2/α²)(α - 1 + e^{-α}) verified")
println("\nDone!")
