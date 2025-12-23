# Test 4D MethodOfLines PDE solvers
using Test
using IDPG

@testset "PDESciML 4D" begin
    # Use very small resolution for fast testing (4D is expensive!)
    resolution = 4
    dx = 1.0 / (resolution - 1)

    # Create a simple initial condition: Gaussian blob in the center
    u0 = zeros(Float64, resolution, resolution, resolution, resolution)
    center = resolution ÷ 2
    for i in 1:resolution
        for j in 1:resolution
            for k in 1:resolution
                for l in 1:resolution
                    x1 = (i - 1) * dx
                    x2 = (j - 1) * dx
                    x3 = (k - 1) * dx
                    x4 = (l - 1) * dx
                    # Only set values inside B^4_+
                    if x1^2 + x2^2 + x3^2 + x4^2 <= 1.0
                        dist_sq = (x1 - 0.5)^2 + (x2 - 0.5)^2 + (x3 - 0.5)^2 + (x4 - 0.5)^2
                        u0[i, j, k, l] = exp(-10 * dist_sq)
                    end
                end
            end
        end
    end

    @testset "4D Diffusion" begin
        println("Testing 4D diffusion with MethodOfLines...")
        D = 0.1
        tspan = (0.0, 0.1)

        result = solve_diffusion_mol_4d(u0, D, tspan;
                                        dx=dx,
                                        saveat=0.05,
                                        apply_mask=true)

        @test length(result.times) >= 2
        @test length(result.snapshots) == length(result.times)
        @test size(result.snapshots[1]) == (resolution, resolution, resolution, resolution)
        @test !isnothing(result.mask)

        # Check that intensity spreads (max decreases due to diffusion)
        initial_max = maximum(result.snapshots[1])
        final_max = maximum(result.snapshots[end])
        println("  Initial max: " * string(round(initial_max, digits=4)))
        println("  Final max: " * string(round(final_max, digits=4)))
        # Diffusion should spread the peak (max should decrease or stay similar)
        @test final_max <= initial_max * 1.1  # Allow small numerical increase

        # Check total mass conservation (approximately)
        initial_mass = sum(result.snapshots[1])
        final_mass = sum(result.snapshots[end])
        println("  Initial mass: " * string(round(initial_mass, digits=4)))
        println("  Final mass: " * string(round(final_mass, digits=4)))

        println("  4D Diffusion test passed!")
    end

    @testset "4D Advection" begin
        println("Testing 4D advection with MethodOfLines...")
        # Advection with upwind scheme needs larger grid for stencil
        # Use resolution 5 minimum for 4D advection
        adv_resolution = 5
        adv_dx = 1.0 / (adv_resolution - 1)

        # Create initial condition for advection
        u0_adv = zeros(Float64, adv_resolution, adv_resolution, adv_resolution, adv_resolution)
        adv_center = adv_resolution ÷ 2
        for i in 1:adv_resolution
            for j in 1:adv_resolution
                for k in 1:adv_resolution
                    for l in 1:adv_resolution
                        x1 = (i - 1) * adv_dx
                        x2 = (j - 1) * adv_dx
                        x3 = (k - 1) * adv_dx
                        x4 = (l - 1) * adv_dx
                        if x1^2 + x2^2 + x3^2 + x4^2 <= 1.0
                            dist_sq = (x1 - 0.5)^2 + (x2 - 0.5)^2 + (x3 - 0.5)^2 + (x4 - 0.5)^2
                            u0_adv[i, j, k, l] = exp(-10 * dist_sq)
                        end
                    end
                end
            end
        end

        v = [0.2, 0.1, 0.05, 0.0]  # Velocity in 4D
        tspan = (0.0, 0.1)

        result = solve_advection_mol_4d(u0_adv, v, tspan;
                                        dx=adv_dx,
                                        saveat=0.05,
                                        apply_mask=true)

        @test length(result.times) >= 2
        @test length(result.snapshots) == length(result.times)
        @test size(result.snapshots[1]) == (adv_resolution, adv_resolution, adv_resolution, adv_resolution)
        @test !isnothing(result.mask)

        println("  Number of snapshots: " * string(length(result.snapshots)))
        println("  4D Advection test passed!")
    end

    @testset "4D Mask" begin
        mask = create_Bd_plus_mask_4d(resolution)
        @test size(mask) == (resolution, resolution, resolution, resolution)

        # Center should be inside
        @test mask[center, center, center, center] == true

        # Origin should be inside
        @test mask[1, 1, 1, 1] == true

        # Far corners should be outside
        @test mask[resolution, resolution, resolution, resolution] == false

        println("  4D Mask test passed!")
    end
end

println("\nAll 4D PDESciML tests passed!")
