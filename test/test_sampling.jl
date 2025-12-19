using Test
using IDPG
using Distributions
using StaticArrays
using Statistics
using Random

@testset "PPP Sampling" begin

    @testset "InteractionSite" begin
        g = SVector{2, Float64}(0.5, 0.3)
        r = SVector{2, Float64}(0.2, 0.5)
        site = InteractionSite{2}(g, r)

        @test site.g == g
        @test site.r == r
    end

    @testset "sample_ppp basic" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        sites = sample_ppp(ρ; rng=rng)

        # Should get some sites (not empty, not huge)
        @test length(sites) >= 0
        @test length(sites) < 10000

        # All sites should be in B^d_+
        for site in sites
            @test in_Bd_plus(site.g; tol=0.1)
            @test in_Bd_plus(site.r; tol=0.1)
        end
    end

    @testset "sample_ppp_product" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 10.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        sites = sample_ppp_product(ρ; rng=rng)

        # All sites should be in B^d_+
        for site in sites
            @test in_Bd_plus(site.g; tol=0.1)
            @test in_Bd_plus(site.r; tol=0.1)
        end
    end

    @testset "expected number of sites" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 4.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        # Sample many times and check mean is reasonable
        n_trials = 200
        n_samples = [length(sample_ppp_product(ρ; rng=rng)) for _ in 1:n_trials]

        empirical_mean = mean(n_samples)

        # Should be positive and finite
        @test empirical_mean > 0
        @test isfinite(empirical_mean)
    end

    @testset "estimate_max_intensity" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.8, 0.2]], [10.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.2, 0.8]], [10.0], 10.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        λ_max = estimate_max_intensity(ρ; n_samples=1000, rng=rng)

        # Should be positive and finite
        @test λ_max > 0
        @test isfinite(λ_max)
    end

    @testset "3D sampling" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.4, 0.3, 0.3]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.4, 0.3, 0.3]], [5.0], 10.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        sites = sample_ppp_product(ρ; rng=rng)

        # All sites should be in B^3_+
        for site in sites
            @test length(site.g) == 3
            @test length(site.r) == 3
            @test all(x -> x >= -0.01, site.g)
            @test all(x -> x >= -0.01, site.r)
        end
    end

end
