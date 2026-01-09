# Benchmark: Integrals.jl (HCubature) vs Monte Carlo
# Compares accuracy and performance for total_intensity and intensity_weighted_mean

using IDPG
using Random
using Statistics
using Printf

println("=" ^ 60)
println("BENCHMARK: Integration Methods (Quadrature vs Monte Carlo)")
println("=" ^ 60)

# Create test intensity functions
function create_test_intensities()
    # Simple single-component
    ρ_simple = BdPlusMixture([1.0], [[0.5, 0.5]], [10.0], 50.0)

    # Multi-component mixture
    ρ_mixture = BdPlusMixture(
        [0.4, 0.3, 0.3],
        [[0.3, 0.2], [0.7, 0.5], [0.4, 0.8]],
        [15.0, 20.0, 12.0],
        100.0
    )

    return ρ_simple, ρ_mixture
end

ρ_simple, ρ_mixture = create_test_intensities()

# =============================================================================
# Accuracy Comparison
# =============================================================================
println("\n--- ACCURACY COMPARISON ---\n")

# Use quadrature with good precision as "ground truth"
println("Computing ground truth with quadrature (reltol=1e-8)...")
c_simple_truth = total_intensity(ρ_simple; method=:quadrature, reltol=1e-8)
c_mixture_truth = total_intensity(ρ_mixture; method=:quadrature, reltol=1e-8)
μ_simple_truth = intensity_weighted_mean(ρ_simple; method=:quadrature, reltol=1e-8)
μ_mixture_truth = intensity_weighted_mean(ρ_mixture; method=:quadrature, reltol=1e-8)

println("\nGround truth values:")
@printf("  Simple:  c = %.8f, μ = [%.6f, %.6f]\n", c_simple_truth, μ_simple_truth[1], μ_simple_truth[2])
@printf("  Mixture: c = %.8f, μ = [%.6f, %.6f]\n", c_mixture_truth, μ_mixture_truth[1], μ_mixture_truth[2])

# Compare quadrature at different tolerances
println("\n[Quadrature at various tolerances]")
for reltol in [1e-4, 1e-6, 1e-8]
    c = total_intensity(ρ_mixture; method=:quadrature, reltol=reltol)
    err = abs(c - c_mixture_truth) / c_mixture_truth
    @printf("  reltol=%.0e: c = %.8f, rel_error = %.2e\n", reltol, c, err)
end

# Compare Monte Carlo at different sample sizes
println("\n[Monte Carlo at various sample sizes]")
for n_samples in [1000, 10000, 50000]
    # Run multiple times to get mean and std
    errors = Float64[]
    for seed in 1:10
        rng = Random.MersenneTwister(seed)
        c = total_intensity(ρ_mixture; method=:montecarlo, n_samples=n_samples, rng=rng)
        push!(errors, abs(c - c_mixture_truth) / c_mixture_truth)
    end
    @printf("  n_samples=%6d: mean_rel_error = %.2e ± %.2e\n", n_samples, mean(errors), std(errors))
end

# =============================================================================
# Performance Comparison
# =============================================================================
println("\n--- PERFORMANCE COMPARISON ---\n")

function benchmark_method(f, n_runs=10)
    # Warmup
    f()
    # Time
    times = [(@elapsed f()) for _ in 1:n_runs]
    return median(times), minimum(times), maximum(times)
end

println("[total_intensity on mixture]")

med, min_t, max_t = benchmark_method(() -> total_intensity(ρ_mixture; method=:quadrature, reltol=1e-6))
@printf("  Quadrature (reltol=1e-6): %.4f ms (min: %.4f, max: %.4f)\n", med*1000, min_t*1000, max_t*1000)

med, min_t, max_t = benchmark_method(() -> total_intensity(ρ_mixture; method=:montecarlo, n_samples=10000))
@printf("  Monte Carlo (n=10000):    %.4f ms (min: %.4f, max: %.4f)\n", med*1000, min_t*1000, max_t*1000)

med, min_t, max_t = benchmark_method(() -> total_intensity(ρ_mixture; method=:montecarlo, n_samples=100000))
@printf("  Monte Carlo (n=100000):   %.4f ms (min: %.4f, max: %.4f)\n", med*1000, min_t*1000, max_t*1000)

println("\n[intensity_weighted_mean on mixture]")

med, min_t, max_t = benchmark_method(() -> intensity_weighted_mean(ρ_mixture; method=:quadrature, reltol=1e-6))
@printf("  Quadrature (reltol=1e-6): %.4f ms (min: %.4f, max: %.4f)\n", med*1000, min_t*1000, max_t*1000)

med, min_t, max_t = benchmark_method(() -> intensity_weighted_mean(ρ_mixture; method=:montecarlo, n_samples=10000))
@printf("  Monte Carlo (n=10000):    %.4f ms (min: %.4f, max: %.4f)\n", med*1000, min_t*1000, max_t*1000)

# =============================================================================
# Summary
# =============================================================================
println("\n" * "=" ^ 60)
println("SUMMARY")
println("=" ^ 60)
println("""
Quadrature (HCubatureJL):
  + Deterministic (reproducible results)
  + Controllable accuracy via reltol
  + Generally faster for smooth integrands in low dimensions
  - May struggle with discontinuities at B^d_+ boundary

Monte Carlo:
  + Dimension-independent convergence rate
  + Robust to integrand discontinuities
  - Stochastic (results vary with RNG seed)
  - Slower convergence (O(1/√n))
""")
