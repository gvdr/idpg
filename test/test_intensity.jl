using Test
using IDPG
using Distributions
using StaticArrays
using LinearAlgebra
using Random

@testset "Intensity Functions" begin

    @testset "BdPlusMixture single component" begin
        # Single component using mean/concentration constructor
        means = [[0.5, 0.5]]
        concentrations = [5.0]
        scale = 10.0

        bm = BdPlusMixture([1.0], means, concentrations, scale)

        # Should evaluate to positive value at interior points
        point = SVector{2, Float64}(0.5, 0.5)
        val = bm(point)
        @test val > 0
        @test isfinite(val)
    end

    @testset "BdPlusMixture with multiple components" begin
        means = [[0.8, 0.2], [0.2, 0.8]]
        concentrations = [10.0, 10.0]
        scale = 50.0

        bm = BdPlusMixture([0.3, 0.7], means, concentrations, scale)

        # Should evaluate correctly
        point = SVector{2, Float64}(0.5, 0.5)
        val = bm(point)
        @test val > 0
        @test isfinite(val)
    end

    @testset "BdPlusMixture total intensity" begin
        rng = MersenneTwister(42)

        means = [[0.5, 0.5]]
        concentrations = [5.0]
        scale = 10.0

        bm = BdPlusMixture([1.0], means, concentrations, scale)

        # Monte Carlo estimate of total intensity
        c = total_intensity(bm; n_samples=5000, rng=rng)

        # Should be positive and of similar order to scale
        @test c > 0
        @test isfinite(c)
    end

    @testset "ProductIntensity" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)

        ρ = ProductIntensity(ρ_G, ρ_R)

        # Product of marginals
        g = SVector{2, Float64}(0.5, 0.5)
        r = SVector{2, Float64}(0.3, 0.4)

        @test isapprox(ρ(g, r), ρ_G(g) * ρ_R(r), rtol=1e-10)
    end

    @testset "marginal_stats" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        stats = marginal_stats(ρ; n_samples=5000, rng=rng)

        # E[N] should be c_G * c_R
        @test stats.E_N > 0
        @test isfinite(stats.E_N)

        # Normalized means should be in B^d_+
        @test in_Bd_plus(stats.μ̃_G; tol=0.1)
        @test in_Bd_plus(stats.μ̃_R; tol=0.1)

        # Average connection probability should be in [0, 1]
        @test 0 <= stats.avg_conn_prob <= 1
    end

    @testset "intensity_weighted_mean" begin
        rng = MersenneTwister(42)

        # Concentrated distribution - mean should be near the mode
        bm = BdPlusMixture([1.0], [[0.7, 0.3]], [20.0], 1.0)
        μ = intensity_weighted_mean(bm; n_samples=5000, rng=rng)

        # Mean should be positive and finite
        @test all(x -> x >= 0, μ)
        @test all(isfinite, μ)
    end

    @testset "sample_from_mixture" begin
        rng = MersenneTwister(42)

        bm = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)

        for _ in 1:100
            point = sample_from_mixture(bm; rng=rng)
            @test length(point) == 2
            @test in_Bd_plus(point; tol=0.1)  # Allow some tolerance
            @test all(x -> x >= -0.01, point)  # Nearly non-negative
        end
    end

    @testset "3D BdPlusMixture" begin
        rng = MersenneTwister(42)

        means = [[0.5, 0.3, 0.3]]
        concentrations = [5.0]
        scale = 10.0

        bm = BdPlusMixture([1.0], means, concentrations, scale)

        # Should evaluate correctly
        point = SVector{3, Float64}(0.4, 0.3, 0.3)
        val = bm(point)
        @test val > 0
        @test isfinite(val)

        # Sampling should work
        for _ in 1:50
            p = sample_from_mixture(bm; rng=rng)
            @test length(p) == 3
            @test all(x -> x >= -0.01, p)
        end
    end

end
