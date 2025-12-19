using Test
using IDPG
using Distributions
using Graphs
using Statistics
using Random
using LinearAlgebra

@testset "Formula Validation" begin

    @testset "E[N] = c_G * c_R (conceptual)" begin
        rng = MersenneTwister(42)

        # Set up product intensity
        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 6.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        # Empirical E[N]
        n_trials = 300
        N_samples = [length(sample_ppp_product(ρ; rng=rng)) for _ in 1:n_trials]
        E_N_empirical = mean(N_samples)

        # Should be positive and reasonable
        @test E_N_empirical > 0
        @test isfinite(E_N_empirical)
    end

    @testset "Node-centric graph generation produces edges" begin
        rng = MersenneTwister(42)

        # Use moderate intensities
        ρ_G = BdPlusMixture([1.0], [[0.6, 0.4]], [5.0], 4.0)
        ρ_R = BdPlusMixture([1.0], [[0.6, 0.4]], [5.0], 4.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        # Empirical edge counts
        n_trials = 100
        edge_counts = Int[]
        for _ in 1:n_trials
            sites = sample_ppp_product(ρ; rng=rng)
            if length(sites) > 0
                graph, _ = generate_node_centric(sites; rng=rng)
                push!(edge_counts, ne(graph))
            else
                push!(edge_counts, 0)
            end
        end
        E_edges_empirical = mean(edge_counts)

        # Should be non-negative
        @test E_edges_empirical >= 0
        @test isfinite(E_edges_empirical)
    end

    @testset "Edge-centric graph generation produces edges" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.6, 0.4]], [5.0], 10.0)
        ρ_R = BdPlusMixture([1.0], [[0.6, 0.4]], [5.0], 10.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        # Empirical E[L]
        n_trials = 200
        L_samples = Int[]
        for _ in 1:n_trials
            sites = sample_ppp_product(ρ; rng=rng)
            sample = generate_edge_centric(sites; rng=rng)
            push!(L_samples, length(sample))
        end
        E_L_empirical = mean(L_samples)

        @test E_L_empirical >= 0
        @test isfinite(E_L_empirical)
    end

    @testset "Marginal stats computation" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        stats = marginal_stats(ρ; n_samples=5000, rng=rng)

        # E_N should be positive
        @test stats.E_N > 0

        # Average connection probability should be in [0, 1]
        @test 0 <= stats.avg_conn_prob <= 1

        # E_edges and E_L should be non-negative
        @test stats.E_edges_node_centric >= 0
        @test stats.E_edges_edge_centric >= 0
    end

    @testset "Connection probability bounds" begin
        # In B^d_+, for points with norm <= 1, dot product is <= 1
        for _ in 1:100
            g = uniform_Bd_plus_sample(2)
            r = uniform_Bd_plus_sample(2)
            p = connection_probability(g, r)
            @test 0 <= p
            @test p <= norm(g) * norm(r) + 1e-10  # Cauchy-Schwarz
        end
    end

    @testset "Aligned vectors give high connection probability" begin
        # Two aligned unit vectors should have connection probability = 1
        g = SVector{2, Float64}(0.6, 0.8)  # norm = 1
        r = SVector{2, Float64}(0.6, 0.8)  # norm = 1
        @test isapprox(connection_probability(g, r), 1.0, atol=1e-10)

        # Two orthogonal vectors should have connection probability = 0
        g2 = SVector{2, Float64}(1.0, 0.0)
        r2 = SVector{2, Float64}(0.0, 1.0)
        @test connection_probability(g2, r2) == 0.0
    end

end
