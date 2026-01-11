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
        means = [[0.5, 0.5]]
        concentrations = [5.0]
        scale = 10.0

        bm = BdPlusMixture([1.0], means, concentrations, scale)

        # Quadrature estimate of total intensity
        c = total_intensity(bm)

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
        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        stats = marginal_stats(ρ)

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
        # Concentrated distribution - mean should be near the mode
        bm = BdPlusMixture([1.0], [[0.7, 0.3]], [20.0], 1.0)
        μ = intensity_weighted_mean(bm)

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

    @testset "ScaledProductEdgeIntensity" begin
        # Create product intensity
        ρ_G = BdPlusMixture([1.0], [[0.6, 0.4]], [10.0], 20.0)
        ρ_R = BdPlusMixture([1.0], [[0.4, 0.6]], [10.0], 20.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        # Create symmetric edge intensity
        ei = SymmetricEdgeIntensity(ρ)

        # C_edge should be E[N]/2
        E_N = total_intensity(ρ)
        @test isapprox(ei.C_edge, E_N / 2, rtol=0.2)

        # edge_intensity should return C_edge
        @test edge_intensity(ei) == ei.C_edge
    end

    @testset "SymmetricEdgeIntensity C_edge formula" begin
        # Create product intensity with known total intensity
        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [10.0], 30.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [10.0], 30.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        # Create edge intensity
        ei = SymmetricEdgeIntensity(ρ)

        # Verify C_edge = E[N]/2
        stats = marginal_stats(ρ)
        @test isapprox(ei.C_edge, stats.E_N / 2, rtol=0.15)
    end

    @testset "ScaledProductEdgeIntensity asymmetric" begin
        # Create separate source and target intensities
        ρ_source = ProductIntensity(
            BdPlusMixture([1.0], [[0.7, 0.3]], [10.0], 25.0),
            BdPlusMixture([1.0], [[0.3, 0.7]], [10.0], 25.0)
        )
        ρ_target = ProductIntensity(
            BdPlusMixture([1.0], [[0.4, 0.6]], [10.0], 35.0),
            BdPlusMixture([1.0], [[0.6, 0.4]], [10.0], 35.0)
        )

        # Create edge intensity
        ei = ScaledProductEdgeIntensity(ρ_source, ρ_target)

        # C_edge should be (C_source + C_target) / 2 for asymmetric case
        C_S = total_intensity(ρ_source)
        C_T = total_intensity(ρ_target)
        @test isapprox(ei.C_edge, (C_S + C_T) / 2, rtol=0.2)
    end

    @testset "marginal_stats for edge intensity" begin
        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [10.0], 25.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [10.0], 25.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        ei = SymmetricEdgeIntensity(ρ)
        stats = marginal_stats(ei)

        # Check all expected fields exist
        @test haskey(stats, :C_edge)
        @test haskey(stats, :E_edges)
        @test haskey(stats, :avg_conn_prob)

        # C_edge should match
        @test stats.C_edge == ei.C_edge

        # E_edges = C_edge * avg_conn_prob
        @test isapprox(stats.E_edges, stats.C_edge * stats.avg_conn_prob, rtol=0.01)

        # avg_conn_prob should be in [0, 1]
        @test 0 <= stats.avg_conn_prob <= 1
    end

end
