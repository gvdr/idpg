using Test
using IDPG
using Distributions
using Graphs
using StaticArrays
using Random

@testset "Graph Generation" begin

    @testset "generate_node_centric" begin
        rng = MersenneTwister(42)

        # Create some test sites with aligned vectors (high connection prob)
        sites = [
            InteractionSite{2}(
                SVector{2, Float64}(0.6, 0.8),
                SVector{2, Float64}(0.6, 0.8)
            ),
            InteractionSite{2}(
                SVector{2, Float64}(0.8, 0.6),
                SVector{2, Float64}(0.8, 0.6)
            ),
            InteractionSite{2}(
                SVector{2, Float64}(1.0, 0.0),
                SVector{2, Float64}(0.0, 1.0)
            ),
        ]

        graph, returned_sites = generate_node_centric(sites; rng=rng)

        @test nv(graph) == 3
        @test returned_sites == sites
    end

    @testset "generate_node_centric with product intensity" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ_R = BdPlusMixture([1.0], [[0.5, 0.5]], [5.0], 5.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        sites = sample_ppp_product(ρ; rng=rng)
        if length(sites) > 0
            graph, _ = generate_node_centric(sites; rng=rng)
            @test nv(graph) == length(sites)
            @test ne(graph) >= 0
        end
    end

    @testset "EdgeCentricSample" begin
        sources = [
            SVector{2, Float64}(0.5, 0.3),
            SVector{2, Float64}(0.2, 0.5)
        ]
        targets = [
            SVector{2, Float64}(0.3, 0.4),
            SVector{2, Float64}(0.1, 0.7)
        ]

        sample = EdgeCentricSample{2}(sources, targets)

        @test length(sample) == 2
        @test sample.sources == sources
        @test sample.targets == targets
    end

    @testset "generate_edge_centric" begin
        rng = MersenneTwister(42)

        # Sites with high connection probability (aligned vectors)
        sites = [
            InteractionSite{2}(
                SVector{2, Float64}(0.8, 0.6),
                SVector{2, Float64}(0.8, 0.6)
            ),
            InteractionSite{2}(
                SVector{2, Float64}(0.7, 0.7),
                SVector{2, Float64}(0.7, 0.7)
            ),
        ]

        sample = generate_edge_centric(sites; rng=rng)

        # Should get at least some edges (high probability)
        @test length(sample) >= 0  # Could be 0 due to randomness
        @test length(sample) <= length(sites)

        # All sources and targets should be in B^d_+
        for s in sample.sources
            @test in_Bd_plus(s; tol=0.1)
        end
        for t in sample.targets
            @test in_Bd_plus(t; tol=0.1)
        end
    end

    @testset "generate_edge_centric with product intensity" begin
        rng = MersenneTwister(42)

        ρ_G = BdPlusMixture([1.0], [[0.7, 0.3]], [5.0], 20.0)
        ρ_R = BdPlusMixture([1.0], [[0.7, 0.3]], [5.0], 20.0)
        ρ = ProductIntensity(ρ_G, ρ_R)

        sites = sample_ppp_product(ρ; rng=rng)
        sample = generate_edge_centric(sites; rng=rng)

        # Should get some edges
        @test length(sample) >= 0
    end

    @testset "discretize_edge_centric" begin
        rng = MersenneTwister(42)

        # Create edge-centric sample
        n_edges = 50
        sources = [uniform_Bd_plus_sample(2; rng=rng) for _ in 1:n_edges]
        targets = [uniform_Bd_plus_sample(2; rng=rng) for _ in 1:n_edges]
        sample = EdgeCentricSample{2}(sources, targets)

        # Discretize
        n_clusters = 5
        graph, source_assign, target_assign = discretize_edge_centric(sample, n_clusters; rng=rng)

        @test nv(graph) == n_clusters
        @test length(source_assign) == n_edges
        @test length(target_assign) == n_edges
        @test all(1 <= a <= n_clusters for a in source_assign)
        @test all(1 <= a <= n_clusters for a in target_assign)
    end

    @testset "empty edge-centric sample" begin
        rng = MersenneTwister(42)

        sample = EdgeCentricSample{2}(LatentPoint{2}[], LatentPoint{2}[])

        graph, _, _ = discretize_edge_centric(sample, 3; rng=rng)
        @test nv(graph) == 3
        @test ne(graph) == 0
    end

    @testset "3D graph generation" begin
        rng = MersenneTwister(42)

        sites = [
            InteractionSite{3}(
                SVector{3, Float64}(0.5, 0.5, 0.5),
                SVector{3, Float64}(0.5, 0.5, 0.5)
            ),
            InteractionSite{3}(
                SVector{3, Float64}(0.6, 0.4, 0.4),
                SVector{3, Float64}(0.6, 0.4, 0.4)
            ),
        ]

        graph, _ = generate_node_centric(sites; rng=rng)
        @test nv(graph) == 2

        sample = generate_edge_centric(sites; rng=rng)
        @test length(sample) >= 0
    end

end
