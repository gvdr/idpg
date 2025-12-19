using Test
using IDPG
using StaticArrays
using Random

@testset "PDE Evolution" begin

    @testset "BdPlusGrid creation" begin
        grid = create_Bd_plus_grid(2, 10)

        @test grid.resolution == 10
        @test length(grid.points) > 0

        # All points should be in B^d_+
        for p in grid.points
            @test in_Bd_plus(p; tol=0.1)  # Larger tolerance for grid points
        end
    end

    @testset "Grid points coverage" begin
        grid = create_Bd_plus_grid(2, 20)

        # Should have reasonable number of points
        # For B^2_+, expect roughly (pi/4) * n^2 points
        @test length(grid.points) > 10
        @test length(grid.points) < 1000
    end

    @testset "evolve_diffusion! basic" begin
        grid = create_Bd_plus_grid(2, 15)
        n_points = length(grid.points)

        # Initial condition: concentrated at one point
        ρ_values = zeros(n_points)
        ρ_values[1] = 1.0

        # Evolve
        D = 0.1
        dt = 0.001
        n_steps = 100

        ρ_evolved = evolve_diffusion!(copy(ρ_values), grid, D, dt, n_steps)

        # Should still be non-negative
        @test all(ρ_evolved .>= 0)

        # Total mass should decrease or stay same (absorbing boundary)
        # With numerical errors, allow small increase
        @test sum(ρ_evolved) <= sum(ρ_values) * 1.1
    end

    @testset "evolve_diffusion! spreads mass" begin
        grid = create_Bd_plus_grid(2, 15)
        n_points = length(grid.points)

        # Initial condition: all mass at one interior point
        ρ_values = zeros(n_points)
        center_idx = div(n_points, 2)
        ρ_values[center_idx] = 10.0

        initial_max = maximum(ρ_values)

        # Evolve for longer time
        D = 0.5
        dt = 0.0001
        n_steps = 1000

        ρ_evolved = evolve_diffusion!(copy(ρ_values), grid, D, dt, n_steps)

        # Peak should decrease (mass spreads)
        @test maximum(ρ_evolved) < initial_max
    end

    @testset "evolve_advection! basic" begin
        grid = create_Bd_plus_grid(2, 15)
        n_points = length(grid.points)

        # Uniform initial condition
        ρ_values = ones(n_points)

        v = SVector{2, Float64}(0.1, 0.0)
        dt = 0.001
        n_steps = 100

        ρ_evolved = evolve_advection!(copy(ρ_values), grid, v, dt, n_steps)

        # Should still be non-negative
        @test all(ρ_evolved .>= 0)
    end

    @testset "evolve_reaction_diffusion! with growth" begin
        grid = create_Bd_plus_grid(2, 10)
        n_points = length(grid.points)

        # Initial condition
        ρ_values = ones(n_points) * 0.5

        # Logistic growth: f(ρ) = r*ρ*(1 - ρ/K)
        r = 0.5
        K = 1.0
        f(ρ) = r * ρ * (1 - ρ / K)

        D = 0.01
        dt = 0.001
        n_steps = 100

        ρ_evolved = evolve_reaction_diffusion!(copy(ρ_values), grid, D, f, dt, n_steps)

        # Should still be non-negative
        @test all(ρ_evolved .>= 0)

        # With growth, mean should increase (or stay same if at carrying capacity)
        # Since we start below K, expect some growth
        @test mean(ρ_evolved) >= mean(ρ_values) * 0.9  # Allow for boundary losses
    end

    @testset "evolve_and_track" begin
        grid = create_Bd_plus_grid(2, 10)

        # Initial intensity function
        ρ_initial(p) = 1.0

        results = evolve_and_track(
            ρ_initial, grid;
            pde_type=:diffusion,
            D=0.1,
            dt=0.01,
            t_final=0.5,
            sample_interval=0.1
        )

        @test length(results.times) > 1
        @test results.times[1] == 0.0
        @test results.times[end] >= 0.4
        @test length(results.ρ_history) == length(results.times)
        @test length(results.total_intensity) == length(results.times)
    end

    @testset "3D grid creation" begin
        grid = create_Bd_plus_grid(3, 8)

        @test grid.resolution == 8
        @test length(grid.points) > 0

        # All points should be in B^3_+
        for p in grid.points
            @test length(p) == 3
            @test all(x -> x >= -0.01, p)
            @test norm(p) <= 1.0 + 0.1
        end
    end

end

# Helper function for testing
using Statistics: mean
using LinearAlgebra: norm
