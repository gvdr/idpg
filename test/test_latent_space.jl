using Test
using IDPG
using StaticArrays
using LinearAlgebra
using Random

@testset "Latent Space (B^d_+)" begin

    @testset "in_Bd_plus" begin
        # Valid B^d_+ points
        @test in_Bd_plus([0.0, 0.0])           # Origin
        @test in_Bd_plus([1.0, 0.0])           # On x-axis boundary
        @test in_Bd_plus([0.0, 1.0])           # On y-axis boundary
        @test in_Bd_plus([0.5, 0.5])           # Interior (norm = sqrt(0.5) < 1)
        @test in_Bd_plus([0.6, 0.8])           # On sphere boundary (norm = 1)
        @test in_Bd_plus([0.3, 0.4])           # Interior

        # Invalid points
        @test !in_Bd_plus([0.8, 0.8])          # norm > 1
        @test !in_Bd_plus([1.5, 0.0])          # norm > 1
        @test !in_Bd_plus([-0.1, 0.5])         # Negative component
        @test !in_Bd_plus([0.5, -0.1])         # Negative component

        # 3D cases
        @test in_Bd_plus([0.5, 0.5, 0.5])      # norm = sqrt(0.75) < 1
        @test !in_Bd_plus([0.6, 0.6, 0.6])     # norm = sqrt(1.08) > 1
    end

    @testset "on_Bd_plus_boundary" begin
        # On sphere boundary (norm = 1)
        @test on_Bd_plus_boundary([1.0, 0.0])
        @test on_Bd_plus_boundary([0.0, 1.0])
        @test on_Bd_plus_boundary([0.6, 0.8])
        @test on_Bd_plus_boundary([1/sqrt(2), 1/sqrt(2)])

        # On coordinate boundary (one coord = 0)
        @test on_Bd_plus_boundary([0.0, 0.5])
        @test on_Bd_plus_boundary([0.5, 0.0])
        @test on_Bd_plus_boundary([0.0, 0.0])  # Origin is on boundary

        # Interior points
        @test !on_Bd_plus_boundary([0.3, 0.3])
        @test !on_Bd_plus_boundary([0.2, 0.2])
    end

    @testset "uniform_Bd_plus_sample" begin
        rng = MersenneTwister(42)

        for d in [2, 3, 4]
            for _ in 1:100
                point = uniform_Bd_plus_sample(d; rng=rng)
                @test length(point) == d
                @test in_Bd_plus(point)
                @test all(x -> x >= 0, point)
                @test norm(point) <= 1.0 + 1e-10
            end
        end
    end

    @testset "project_to_Bd_plus" begin
        # Already in B^d_+
        p = [0.3, 0.4]
        proj = project_to_Bd_plus(p)
        @test isapprox(proj, p, atol=1e-10)

        # Outside sphere - needs scaling
        p_out = [1.5, 0.0]
        proj_out = project_to_Bd_plus(p_out)
        @test in_Bd_plus(proj_out)
        @test isapprox(norm(proj_out), 1.0, atol=1e-10)

        # Negative values - needs clamping
        p_neg = [-0.5, 0.8]
        proj_neg = project_to_Bd_plus(p_neg)
        @test in_Bd_plus(proj_neg)
        @test all(x -> x >= 0, proj_neg)

        # Both negative and outside
        p_both = [-1.0, 2.0]
        proj_both = project_to_Bd_plus(p_both)
        @test in_Bd_plus(proj_both)
    end

    @testset "Bd_plus_volume" begin
        # Volume of B^d_+ = volume of d-ball / 2^d
        # For d=2: pi * r^2 / 4 = pi/4 (quarter circle)
        @test isapprox(Bd_plus_volume(2), pi/4, atol=1e-10)

        # For d=3: (4/3)*pi*r^3 / 8 = pi/6 (eighth of sphere)
        @test isapprox(Bd_plus_volume(3), pi/6, atol=1e-10)
    end

    @testset "connection_probability" begin
        # Aligned vectors
        g = SVector{2, Float64}(0.6, 0.8)
        r = SVector{2, Float64}(0.6, 0.8)
        @test isapprox(connection_probability(g, r), 0.36 + 0.64, atol=1e-10)

        # Orthogonal vectors
        g2 = SVector{2, Float64}(1.0, 0.0)
        r2 = SVector{2, Float64}(0.0, 1.0)
        @test connection_probability(g2, r2) == 0.0

        # Partial overlap
        g3 = SVector{2, Float64}(0.5, 0.5)
        r3 = SVector{2, Float64}(0.5, 0.0)
        @test isapprox(connection_probability(g3, r3), 0.25, atol=1e-10)

        # In B^d_+, dot products are always <= 1 when both vectors have norm <= 1
        g4 = SVector{2, Float64}(0.6, 0.8)  # norm = 1
        r4 = SVector{2, Float64}(0.8, 0.6)  # norm = 1
        @test connection_probability(g4, r4) <= 1.0
    end

    @testset "radial_coordinate" begin
        @test radial_coordinate([0.0, 0.0]) == 0.0
        @test radial_coordinate([1.0, 0.0]) == 1.0
        @test radial_coordinate([0.0, 1.0]) == 1.0
        @test isapprox(radial_coordinate([0.6, 0.8]), 1.0, atol=1e-10)
        @test isapprox(radial_coordinate([0.3, 0.4]), 0.5, atol=1e-10)
    end

    @testset "angular_coordinates" begin
        # On axes
        @test isapprox(angular_coordinates([1.0, 0.0]), [1.0, 0.0], atol=1e-10)
        @test isapprox(angular_coordinates([0.0, 1.0]), [0.0, 1.0], atol=1e-10)

        # Diagonal
        diag = angular_coordinates([1.0, 1.0])
        @test isapprox(diag, [1/sqrt(2), 1/sqrt(2)], atol=1e-10)

        # General point
        ang = angular_coordinates([0.6, 0.8])
        @test isapprox(ang, [0.6, 0.8], atol=1e-10)  # Already unit norm
    end

end
