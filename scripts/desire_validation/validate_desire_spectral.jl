# Validation Script: Desire Operator Spectral Consistency
#
# Tests the theorem: σ_k(A_N) / N → σ_k(D̃) as N → ∞
# where σ_k(D̃) = √(λ_k(Σ_G Σ_R))
#
# Two test regimes:
# (a) Growing Λ (single large graph): verify σ_k(A)/N → σ_k(D̃)
# (b) Fixed Λ, growing m: verify bias = O(1/√Λ), variance = O(1/√(mΛ))
#
# Usage:
#   julia --project=. scripts/desire_validation/validate_desire_spectral.jl

using IDPG
using LinearAlgebra
using Statistics
using Random
using Graphs

# ============================================================================
# Helper Functions
# ============================================================================

"""
Sample a graph and compute σ_k(A)/N.
"""
function sample_graph_singular_values(mean_G, mean_R, conc, Λ; rng=Random.default_rng())
    ρ_G_cal, ρ_R_cal, ρ_cal = create_calibrated_product_intensity(
        Λ; mean_G=mean_G, mean_R=mean_R, conc=conc, rng=rng
    )

    sites = sample_ppp_product(ρ_cal; rng=rng)
    N = length(sites)

    if N < 2
        return (N=N, σ_scaled=zeros(2))
    end

    graph, _ = generate_node_centric(sites; rng=rng)
    A = Float64.(Matrix(adjacency_matrix(graph)))

    F = svd(A)
    σ_A = F.S

    # Return top 2 singular values scaled by N
    σ_scaled = zeros(2)
    σ_scaled[1] = σ_A[1] / N
    σ_scaled[2] = length(σ_A) >= 2 ? σ_A[2] / N : 0.0

    return (N=N, σ_scaled=σ_scaled)
end

# ============================================================================
# Test (a): Growing Λ
# ============================================================================

function test_growing_lambda(; lambdas=[100.0, 500.0, 1000.0, 2000.0, 5000.0],
                              mean_G=[0.5, 0.5], mean_R=[0.5, 0.5], conc=15.0,
                              n_reps=5)
    println("\n" * "="^80)
    println("TEST (a): Growing Λ — Single Large Graph Convergence")
    println("="^80)

    # Compute theoretical values
    ρ_G = BdPlusMixture([1.0], [mean_G], [conc], 1.0)
    ρ_R = BdPlusMixture([1.0], [mean_R], [conc], 1.0)
    ρ = ProductIntensity(ρ_G, ρ_R)

    σ_theory = desire_operator_singular_values(ρ)
    println("\nGeometry: mean_G=", mean_G, ", mean_R=", mean_R, ", κ=", conc)
    println("Theoretical σ(D̃) = ", round.(σ_theory, digits=4))
    println("  σ₁ = ", round(σ_theory[1], digits=4))
    println("  σ₂ = ", round(σ_theory[2], digits=4))

    println("\nΛ\t\tE[N]\t\tσ₁(A)/N ± SE\t\tσ₂(A)/N ± SE\t\t1/√N")
    println("-"^80)

    for Λ in lambdas
        all_σ = [zeros(2) for _ in 1:n_reps]
        all_N = zeros(Int, n_reps)

        for rep in 1:n_reps
            result = sample_graph_singular_values(mean_G, mean_R, conc, Λ;
                                                   rng=MersenneTwister(42 + rep + Int(Λ)))
            all_N[rep] = result.N
            all_σ[rep] = result.σ_scaled
        end

        mean_N = mean(all_N)
        mean_σ = mean(all_σ)
        se_σ = std(all_σ) / sqrt(n_reps)
        noise_floor = 1.0 / sqrt(mean_N)

        println(Int(Λ), "\t\t", round(mean_N, digits=0), "\t\t",
                round(mean_σ[1], digits=4), " ± ", round(se_σ[1], digits=4), "\t\t",
                round(mean_σ[2], digits=4), " ± ", round(se_σ[2], digits=4), "\t\t",
                round(noise_floor, digits=4))
    end

    println("\nTheory: σ₁ = ", round(σ_theory[1], digits=4),
            ", σ₂ = ", round(σ_theory[2], digits=4))
end

# ============================================================================
# Test (b): Fixed Λ, growing m — Bias Verification
# ============================================================================

function test_bias_vs_lambda(; lambdas=[100.0, 250.0, 500.0, 1000.0, 2000.0],
                              mean_G=[0.5, 0.5], mean_R=[0.5, 0.5], conc=15.0,
                              m=100)
    println("\n" * "="^80)
    println("TEST (b): Bias vs Λ — Verifying O(1/√Λ) Scaling")
    println("="^80)

    # Compute theoretical values
    ρ_G = BdPlusMixture([1.0], [mean_G], [conc], 1.0)
    ρ_R = BdPlusMixture([1.0], [mean_R], [conc], 1.0)
    ρ = ProductIntensity(ρ_G, ρ_R)

    σ_theory = desire_operator_singular_values(ρ)
    println("\nTheory: σ₁ = ", round(σ_theory[1], digits=4),
            ", σ₂ = ", round(σ_theory[2], digits=4))
    println("Averaging m=", m, " graphs for each Λ")

    println("\nΛ\t\t1/√Λ\t\tmean(σ₁/N)\tmean(σ₂/N)\tbias₂\t\tbias₂/theory")
    println("-"^90)

    for Λ in lambdas
        all_σ1 = Float64[]
        all_σ2 = Float64[]

        for l in 1:m
            result = sample_graph_singular_values(mean_G, mean_R, conc, Λ;
                                                   rng=MersenneTwister(Int(Λ)*1000 + l))
            if result.N >= 2
                push!(all_σ1, result.σ_scaled[1])
                push!(all_σ2, result.σ_scaled[2])
            end
        end

        mean_σ1 = mean(all_σ1)
        mean_σ2 = mean(all_σ2)
        bias2 = mean_σ2 - σ_theory[2]
        expected_bias = 1/sqrt(Λ)
        bias_ratio = bias2 / expected_bias

        println(Int(Λ), "\t\t", round(expected_bias, digits=4), "\t\t",
                round(mean_σ1, digits=4), "\t\t",
                round(mean_σ2, digits=4), "\t\t",
                round(bias2, digits=4), "\t\t",
                round(bias_ratio, digits=2))
    end

    println("\nNote: bias₂ should scale as O(1/√Λ)")
    println("The ratio bias₂/(1/√Λ) should be roughly constant")
end

# ============================================================================
# Test (c): Multiple Test Cases
# ============================================================================

function test_multiple_geometries()
    println("\n" * "="^80)
    println("TEST (c): Multiple Geometries")
    println("="^80)

    test_cases = [
        (name="Symmetric", mean_G=[0.5, 0.5], mean_R=[0.5, 0.5], conc=15.0),
        (name="Asymmetric", mean_G=[0.7, 0.3], mean_R=[0.3, 0.7], conc=15.0),
        (name="Concentrated", mean_G=[0.6, 0.4], mean_R=[0.4, 0.6], conc=50.0),
    ]

    Λ = 1000.0
    n_reps = 10

    for tc in test_cases
        ρ_G = BdPlusMixture([1.0], [tc.mean_G], [tc.conc], 1.0)
        ρ_R = BdPlusMixture([1.0], [tc.mean_R], [tc.conc], 1.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        σ_theory = desire_operator_singular_values(ρ)

        all_σ = [zeros(2) for _ in 1:n_reps]
        for rep in 1:n_reps
            result = sample_graph_singular_values(tc.mean_G, tc.mean_R, tc.conc, Λ;
                                                   rng=MersenneTwister(hash(tc.name) + rep))
            all_σ[rep] = result.σ_scaled
        end

        mean_σ = mean(all_σ)
        se_σ = std(all_σ) / sqrt(n_reps)

        println("\n", tc.name, ":")
        println("  Theory:   σ₁ = ", round(σ_theory[1], digits=4),
                ", σ₂ = ", round(σ_theory[2], digits=4))
        println("  Empirical: σ₁ = ", round(mean_σ[1], digits=4), " ± ", round(se_σ[1], digits=4),
                ", σ₂ = ", round(mean_σ[2], digits=4), " ± ", round(se_σ[2], digits=4))
        println("  Error:    |Δσ₁| = ", round(abs(mean_σ[1] - σ_theory[1]), digits=4),
                ", |Δσ₂| = ", round(abs(mean_σ[2] - σ_theory[2]), digits=4))
    end
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("="^80)
    println("DESIRE OPERATOR SPECTRAL CONSISTENCY VALIDATION")
    println("="^80)
    println("\nTheorem: σ_k(A_N)/N → σ_k(D̃) as N → ∞")
    println("where σ_k(D̃) = √(λ_k(Σ_G Σ_R))")
    println("\nKey predictions:")
    println("  (a) Single large graph: convergence as N → ∞")
    println("  (b) Multiple small graphs: bias = O(1/√Λ), variance = O(1/√(mΛ))")

    test_growing_lambda()
    test_bias_vs_lambda()
    test_multiple_geometries()

    println("\n" * "="^80)
    println("VALIDATION COMPLETE")
    println("="^80)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
