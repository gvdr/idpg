# Figure: Desire Operator Spectral Consistency Validation
#
# Validates the theorem: σ_k(A_N)/N → σ_k(D̃) as N → ∞
# where σ_k(D̃) = √(λ_k(Σ_G Σ_R))
#
# Panels:
#   (a) Convergence of σ_k(A)/N vs N for k=1,2,3,4,5
#   (b) Bias scaling: shows bias = O(1/√Λ)
#   (c) Singular value spectrum comparison (theory vs empirical)
#
# Usage:
#   julia --project=. scripts/desire_validation/figure_desire_spectral.jl

using IDPG
using CairoMakie
using LaTeXStrings
using LinearAlgebra
using Statistics
using Random
using Graphs
using JLD2

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "desire_validation")
const DATA_FILE = joinpath(OUTPUT_DIR, "desire_spectral_data.jld2")

# 2-component mixture in d=4 for well-separated singular values
const D = 4
const G_MEANS = [[0.7, 0.5, 0.15, 0.1], [0.15, 0.1, 0.5, 0.7]]
const R_MEANS = [[0.5, 0.7, 0.1, 0.15], [0.1, 0.15, 0.7, 0.5]]
const G_WEIGHTS = [0.5, 0.5]
const R_WEIGHTS = [0.5, 0.5]
const CONCENTRATION = 10.0

# Test parameters
const LAMBDAS_CONVERGENCE = [100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0]
const LAMBDAS_BIAS = [100.0, 200.0, 500.0, 1000.0, 2000.0]
const N_REPS_CONVERGENCE = 10
const M_GRAPHS_BIAS = 100

# Multi-graph convergence parameters
const LAMBDA_FIXED = 300.0  # Fixed Λ for multi-graph test
const M_VALUES = [1, 2, 5, 10, 20, 50, 100]  # Number of graphs to average

# ============================================================================
# Helper Functions
# ============================================================================

function create_test_intensity()
    ρ_G = BdPlusMixture(G_WEIGHTS, G_MEANS, fill(CONCENTRATION, length(G_MEANS)), 1.0)
    ρ_R = BdPlusMixture(R_WEIGHTS, R_MEANS, fill(CONCENTRATION, length(R_MEANS)), 1.0)
    return ProductIntensity(ρ_G, ρ_R)
end

function sample_graph_singular_values(Λ; rng=Random.default_rng(), n_sv=5)
    ρ_G = BdPlusMixture(G_WEIGHTS, G_MEANS, fill(CONCENTRATION, length(G_MEANS)), 1.0)
    ρ_R = BdPlusMixture(R_WEIGHTS, R_MEANS, fill(CONCENTRATION, length(R_MEANS)), 1.0)

    # Scale to target Λ
    c_G = total_intensity(ρ_G)
    c_R = total_intensity(ρ_R)
    scale = sqrt(Λ / (c_G * c_R))

    ρ_G_scaled = BdPlusMixture(G_WEIGHTS, G_MEANS, fill(CONCENTRATION, length(G_MEANS)), scale)
    ρ_R_scaled = BdPlusMixture(R_WEIGHTS, R_MEANS, fill(CONCENTRATION, length(R_MEANS)), scale)
    ρ = ProductIntensity(ρ_G_scaled, ρ_R_scaled)

    sites = sample_ppp_product(ρ; rng=rng)
    N = length(sites)

    if N < n_sv
        return (N=N, σ_scaled=zeros(n_sv))
    end

    graph, _ = generate_node_centric(sites; rng=rng)
    A = Float64.(Matrix(adjacency_matrix(graph)))

    F = svd(A)
    σ_A = F.S

    σ_scaled = zeros(n_sv)
    for k in 1:min(n_sv, length(σ_A))
        σ_scaled[k] = σ_A[k] / N
    end

    return (N=N, σ_scaled=σ_scaled)
end

# ============================================================================
# Data Computation
# ============================================================================

function compute_data()
    println("="^70)
    println("Computing Desire Operator Spectral Validation Data")
    println("="^70)

    Random.seed!(42)

    # Theoretical values
    ρ = create_test_intensity()
    σ_theory = desire_operator_singular_values(ρ)
    stats = desire_stats(ρ)

    println("\nTheoretical singular values:")
    for k in 1:D
        println("  σ_", k, "(D̃) = ", round(σ_theory[k], digits=4))
    end

    # Test (a): Convergence with growing Λ
    println("\nComputing convergence data (growing Λ)...")
    convergence_data = Dict{Float64, Any}()

    for Λ in LAMBDAS_CONVERGENCE
        print("  Λ = ", Int(Λ), "...")
        all_σ = [zeros(5) for _ in 1:N_REPS_CONVERGENCE]
        all_N = zeros(Int, N_REPS_CONVERGENCE)

        for rep in 1:N_REPS_CONVERGENCE
            result = sample_graph_singular_values(Λ; rng=MersenneTwister(42 + rep + Int(Λ)), n_sv=5)
            all_N[rep] = result.N
            all_σ[rep] = result.σ_scaled
        end

        convergence_data[Λ] = (
            mean_N = mean(all_N),
            mean_σ = mean(all_σ),
            std_σ = std(all_σ),
            se_σ = std(all_σ) / sqrt(N_REPS_CONVERGENCE)
        )
        println(" done (mean N = ", round(convergence_data[Λ].mean_N, digits=0), ")")
    end

    # Test (b): Bias scaling
    println("\nComputing bias data (fixed Λ, m graphs)...")
    bias_data = Dict{Float64, Any}()

    for Λ in LAMBDAS_BIAS
        print("  Λ = ", Int(Λ), " (m=", M_GRAPHS_BIAS, ")...")
        all_σ = [zeros(5) for _ in 1:M_GRAPHS_BIAS]
        all_N = zeros(Int, M_GRAPHS_BIAS)

        for l in 1:M_GRAPHS_BIAS
            result = sample_graph_singular_values(Λ; rng=MersenneTwister(Int(Λ)*1000 + l), n_sv=5)
            all_N[l] = result.N
            all_σ[l] = result.σ_scaled
        end

        # Compute bias (only for k ≤ d where theory exists)
        mean_σ_vec = mean(all_σ)
        bias_vec = zeros(5)
        for k in 1:min(D, 5)
            bias_vec[k] = mean_σ_vec[k] - σ_theory[k]
        end
        # For k > d, bias is just the empirical value (theory is 0)
        for k in (D+1):5
            bias_vec[k] = mean_σ_vec[k]
        end

        bias_data[Λ] = (
            mean_N = mean(all_N),
            mean_σ = mean_σ_vec,
            std_σ = std(all_σ),
            bias = bias_vec
        )
        println(" done")
    end

    # Test (c): Multi-graph convergence (fixed Λ, growing m)
    println("\nComputing multi-graph convergence data (fixed Λ=", Int(LAMBDA_FIXED), ", growing m)...")
    multigraph_data = Dict{Int, Any}()

    for m in M_VALUES
        print("  m = ", m, "...")

        # For each m, we compute the averaged estimator multiple times to get its distribution
        n_trials = 50  # Number of times to compute the m-graph average
        trial_means = [zeros(5) for _ in 1:n_trials]

        for trial in 1:n_trials
            # Sample m graphs and average their scaled singular values
            graph_σs = [zeros(5) for _ in 1:m]
            for l in 1:m
                result = sample_graph_singular_values(LAMBDA_FIXED;
                    rng=MersenneTwister(trial * 10000 + m * 1000 + l), n_sv=5)
                graph_σs[l] = result.σ_scaled
            end
            trial_means[trial] = mean(graph_σs)
        end

        # Statistics across trials
        multigraph_data[m] = (
            mean_σ = mean(trial_means),
            std_σ = std(trial_means),
            se_σ = std(trial_means) / sqrt(n_trials)
        )
        println(" done")
    end

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        σ_theory = σ_theory,
        Σ_G = stats.Σ_G,
        Σ_R = stats.Σ_R,
        convergence_data = convergence_data,
        bias_data = bias_data,
        multigraph_data = multigraph_data,
        lambdas_convergence = LAMBDAS_CONVERGENCE,
        lambdas_bias = LAMBDAS_BIAS,
        m_values = M_VALUES,
        lambda_fixed = LAMBDA_FIXED,
    )
    @save DATA_FILE data
    println("\nSaved data: ", DATA_FILE)

    return data
end

function load_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

function plot_figure(data)
    println("="^70)
    println("Generating Figure")
    println("="^70)

    σ_theory = data.σ_theory
    conv_data = data.convergence_data
    bias_data = data.bias_data
    mg_data = data.multigraph_data
    λ_fixed = data.lambda_fixed

    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=16))) do
        fig = Figure(size=(1400, 400))

        # Colors for singular values
        colors = [:blue, :red, :green, :orange, :purple]

        # ===== Panel (a): Single graph convergence (growing Λ) =====
        ax_a = Axis(fig[1, 1],
            xlabel = L"Total intensity $\Lambda$",
            ylabel = L"Scaled singular value $\sigma_k(A)/N$",
            title = L"(a) Single graph: $\Lambda \to \infty$",
            xscale = log10,
            yscale = log10)

        lambdas = sort(collect(keys(conv_data)))

        for k in 1:5
            means = [conv_data[Λ].mean_σ[k] for Λ in lambdas]
            ses = [conv_data[Λ].se_σ[k] for Λ in lambdas]

            # Empirical points with error bars
            scatter!(ax_a, lambdas, means, color=colors[k], markersize=8,
                    label=L"\sigma_%$k")
            errorbars!(ax_a, lambdas, means, ses, color=colors[k], whiskerwidth=6)

            # Theoretical line
            if k <= 4
                hlines!(ax_a, [σ_theory[k]], color=colors[k], linestyle=:dash, alpha=0.7)
            end
        end

        # Noise floor
        Λ_dense = 10 .^ range(log10(minimum(lambdas)), log10(maximum(lambdas)), length=50)
        lines!(ax_a, Λ_dense, 1 ./ sqrt.(Λ_dense), color=:gray, linestyle=:dot,
               label=L"1/\sqrt{\Lambda}")

        Legend(fig[1, 2], ax_a, nbanks=2)

        # ===== Panel (b): Bias scaling (single graph property) =====
        ax_b = Axis(fig[1, 3],
            xlabel = L"Total intensity $\Lambda$",
            ylabel = L"Bias: $|\bar{\sigma}_k - \sigma_k(\tilde{D})|$",
            title = L"(b) Single graph bias: $O(1/\sqrt{\Lambda})$",
            xscale = log10,
            yscale = log10)

        lambdas_b = sort(collect(keys(bias_data)))

        for k in 1:4  # All four singular values
            biases = [abs(bias_data[Λ].bias[k]) for Λ in lambdas_b]
            scatter!(ax_b, lambdas_b, biases, color=colors[k], markersize=8,
                    label=L"\sigma_%$k")
        end

        # Reference line 1/√Λ
        Λ_dense = 10 .^ range(log10(minimum(lambdas_b)), log10(maximum(lambdas_b)), length=50)
        lines!(ax_b, Λ_dense, 0.5 ./ sqrt.(Λ_dense), color=:gray, linestyle=:dash,
               label=L"c/\sqrt{\Lambda}")

        Legend(fig[1, 4], ax_b)

        # ===== Panel (c): Variance vs m (multiple graphs) =====
        ms = sort(collect(keys(mg_data)))
        ax_c = Axis(fig[1, 5],
            xlabel = L"Number of graphs $m$",
            ylabel = L"Std of $\bar{\sigma}_k$",
            title = L"(c) Multi-graph variance ($\Lambda = %$(Int(λ_fixed))$): $O(1/\sqrt{m})$",
            xscale = log10,
            yscale = log10,
            xticks = (ms, string.(ms)))

        for k in 1:2  # Only σ₁ and σ₂ (others are noise-dominated)
            stds = [mg_data[m].std_σ[k] for m in ms]
            scatter!(ax_c, ms, stds, color=colors[k], markersize=8,
                    label=L"\sigma_%$k")
        end

        # Reference line 1/√m
        ms_dense = 10 .^ range(log10(minimum(ms)), log10(maximum(ms)), length=50)
        # Scale to match σ₁ at m=1
        c_ref = mg_data[minimum(ms)].std_σ[1] * sqrt(minimum(ms))
        lines!(ax_c, ms_dense, c_ref ./ sqrt.(ms_dense), color=:gray, linestyle=:dash,
               label=L"c/\sqrt{m}")

        Legend(fig[1, 6], ax_c)

        # Title
        Label(fig[0, 1:6],
              L"Desire Operator Spectral Consistency ($d=4$, 2-component mixture)",
              fontsize=20)

        fig
    end

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_desire_spectral.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: ", output_file)

    # Also save PDF
    output_pdf = joinpath(OUTPUT_DIR, "figure_desire_spectral.pdf")
    save(output_pdf, fig)
    println("Saved: ", output_pdf)

    return fig
end

# ============================================================================
# Summary Statistics
# ============================================================================

function print_summary(data)
    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)

    σ_theory = data.σ_theory
    println("\nTheoretical σ(D̃):")
    for k in 1:length(σ_theory)
        println("  σ_", k, " = ", round(σ_theory[k], digits=4))
    end

    println("\nConvergence (largest Λ = ", maximum(data.lambdas_convergence), "):")
    Λ_max = maximum(data.lambdas_convergence)
    conv = data.convergence_data[Λ_max]
    for k in 1:4
        err = abs(conv.mean_σ[k] - σ_theory[k])
        println("  σ_", k, ": theory=", round(σ_theory[k], digits=4),
                ", empirical=", round(conv.mean_σ[k], digits=4),
                ", error=", round(err, digits=4))
    end
    println("  σ_5: empirical=", round(conv.mean_σ[5], digits=4),
            " (should be ≈ 1/√N = ", round(1/sqrt(conv.mean_N), digits=4), ")")

    println("\nBias scaling verification:")
    for Λ in sort(collect(keys(data.bias_data)))
        bias2 = data.bias_data[Λ].bias[2]
        expected = 1/sqrt(Λ)
        println("  Λ=", Int(Λ), ": bias₂=", round(bias2, digits=4),
                ", 1/√Λ=", round(expected, digits=4),
                ", ratio=", round(bias2/expected, digits=2))
    end
end

# ============================================================================
# Main
# ============================================================================

function main(; recompute=false)
    if recompute || !isfile(DATA_FILE)
        data = compute_data()
    else
        println("Loading cached data from ", DATA_FILE)
        data = load_data()
    end

    plot_figure(data)
    print_summary(data)

    return data
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main(recompute=true)
end
