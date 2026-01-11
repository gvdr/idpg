# Verification: Run simulations to check behaviors are as expected

using IDPG
using CairoMakie
using LinearAlgebra
using Printf

output_dir = joinpath(@__DIR__, "figures")
mkpath(output_dir)

println("=" ^ 60)
println("SIMULATION VERIFICATION")
println("=" ^ 60)
println("Output directory: ", output_dir)

# =============================================================================
# Test 1: Diffusion spreads a localized intensity
# =============================================================================
println("\n--- TEST 1: Diffusion Spreading ---")

d = 2
resolution = 30
grid = create_Bd_plus_grid(d, resolution)

# Localized initial condition
μ₀ = [0.6, 0.4]
σ = 0.1
ρ₀ = [exp(-sum((p .- μ₀).^2) / (2σ^2)) for p in grid.points]

D = 0.02
tspan = (0.0, 1.0)
sol = evolve_diffusion(ρ₀, grid, D, tspan; saveat=[0.0, 0.25, 0.5, 1.0])

fig = Figure(size=(1000, 300))
for (i, t) in enumerate(sol.t)
    ax = Axis(fig[1, i], title="t = " * string(t), aspect=1)
    ρ = max.(sol.u[i], 0.0)

    # Plot as scatter with intensity as color
    xs = [p[1] for p in grid.points]
    ys = [p[2] for p in grid.points]
    scatter!(ax, xs, ys, color=ρ, colormap=:viridis, markersize=8)

    # Draw B²₊ boundary
    θ = range(0, π/2, 50)
    lines!(ax, cos.(θ), sin.(θ), color=:black, linewidth=2)

    xlims!(ax, -0.05, 1.05)
    ylims!(ax, -0.05, 1.05)
end
save(joinpath(output_dir, "diffusion_spreading.png"), fig)
println("  Saved: diffusion_spreading.png")
println("  Expected: Intensity should spread out and flatten over time")

# Check: variance should increase
variances = Float64[]
for ρ in sol.u
    ρ_pos = max.(ρ, 0.0)
    total = sum(ρ_pos)
    if total > 1e-10
        μ = sum(ρᵢ * xᵢ for (ρᵢ, xᵢ) in zip(ρ_pos, grid.points)) ./ total
        var = sum(ρᵢ * sum((xᵢ .- μ).^2) for (ρᵢ, xᵢ) in zip(ρ_pos, grid.points)) / total
        push!(variances, var)
    end
end
println("  Variances over time: ", round.(variances, digits=4))
@assert all(diff(variances) .>= -1e-6) "Variance should not decrease under diffusion"
println("  ✓ Variance increases (diffusion working correctly)")

# =============================================================================
# Test 2: Advection moves intensity in velocity direction
# =============================================================================
println("\n--- TEST 2: Advection Transport ---")

v⃗ = [0.3, 0.2]
ρ₀_adv = [exp(-sum((p .- [0.3, 0.3]).^2) / (2*0.1^2)) for p in grid.points]
sol_adv = evolve_advection(ρ₀_adv, grid, v⃗, (0.0, 0.5); saveat=[0.0, 0.25, 0.5])

fig2 = Figure(size=(800, 300))
centers = []
for (i, t) in enumerate(sol_adv.t)
    ax = Axis(fig2[1, i], title="t = " * string(t), aspect=1)
    ρ = max.(sol_adv.u[i], 0.0)

    xs = [p[1] for p in grid.points]
    ys = [p[2] for p in grid.points]
    scatter!(ax, xs, ys, color=ρ, colormap=:plasma, markersize=8)

    # Compute and mark center of mass
    total = sum(ρ)
    if total > 1e-10
        μ = sum(ρᵢ * xᵢ for (ρᵢ, xᵢ) in zip(ρ, grid.points)) ./ total
        push!(centers, μ)
        scatter!(ax, [μ[1]], [μ[2]], color=:red, markersize=15, marker=:x)
    end

    θ = range(0, π/2, 50)
    lines!(ax, cos.(θ), sin.(θ), color=:black, linewidth=2)
    xlims!(ax, -0.05, 1.05)
    ylims!(ax, -0.05, 1.05)
end
save(joinpath(output_dir, "advection_transport.png"), fig2)
println("  Saved: advection_transport.png")
println("  Centers of mass: ", [round.(c, digits=3) for c in centers])
println("  Velocity: ", v⃗)
println("  Expected displacement in 0.5s: ", 0.5 .* v⃗)

# Check direction of movement
if length(centers) >= 2
    displacement = centers[end] .- centers[1]
    v_normalized = v⃗ ./ norm(v⃗)
    disp_normalized = displacement ./ norm(displacement)
    alignment = dot(v_normalized, disp_normalized)
    @printf("  Movement alignment with velocity: %.3f (should be ~1.0)\n", alignment)
    @assert alignment > 0.9 "Mass should move in velocity direction"
    println("  ✓ Advection moves mass in correct direction")
end

# =============================================================================
# Test 3: Reaction-diffusion with growth
# =============================================================================
println("\n--- TEST 3: Reaction-Diffusion (Logistic Growth) ---")

# Logistic growth: f(ρ) = r*ρ*(1 - ρ/K)
r = 0.5  # growth rate
K = 5.0  # carrying capacity
f_logistic(ρ) = r * ρ * (1 - ρ / K)

ρ₀_rd = [0.1 * exp(-sum((p .- [0.5, 0.5]).^2) / (2*0.15^2)) for p in grid.points]
sol_rd = evolve_reaction_diffusion(ρ₀_rd, grid, 0.005, f_logistic, (0.0, 5.0); saveat=0.0:1.0:5.0)

fig3 = Figure(size=(1200, 250))
total_intensities = Float64[]
for (i, t) in enumerate(sol_rd.t)
    ax = Axis(fig3[1, i], title="t = " * string(Int(t)), aspect=1)
    ρ = max.(sol_rd.u[i], 0.0)
    push!(total_intensities, sum(ρ) * grid.h^d)

    xs = [p[1] for p in grid.points]
    ys = [p[2] for p in grid.points]
    scatter!(ax, xs, ys, color=ρ, colormap=:inferno, markersize=8,
             colorrange=(0, K))

    θ = range(0, π/2, 50)
    lines!(ax, cos.(θ), sin.(θ), color=:white, linewidth=2)
    xlims!(ax, -0.05, 1.05)
    ylims!(ax, -0.05, 1.05)
end
save(joinpath(output_dir, "reaction_diffusion.png"), fig3)
println("  Saved: reaction_diffusion.png")
println("  Total intensities over time: ", round.(total_intensities, digits=3))
println("  Expected: Intensity should grow and spread, approaching carrying capacity")

# Check growth occurred
@assert total_intensities[end] > total_intensities[1] "Total intensity should increase with logistic growth"
println("  ✓ Reaction-diffusion shows growth")

# =============================================================================
# Test 4: Integration accuracy for marginal_stats
# =============================================================================
println("\n--- TEST 4: Integration for marginal_stats ---")

ρ_G = BdPlusMixture([1.0], [[0.6, 0.4]], [15.0], 10.0)
ρ_R = BdPlusMixture([1.0], [[0.5, 0.6]], [12.0], 8.0)
ρ_product = ProductIntensity(ρ_G, ρ_R)

stats = marginal_stats(ρ_product)
println("  c_G = ", round(stats.c_G, digits=4))
println("  c_R = ", round(stats.c_R, digits=4))
println("  E[N] = ", round(stats.E_N, digits=4))
println("  μ̃_G = ", round.(stats.μ̃_G, digits=4))
println("  μ̃_R = ", round.(stats.μ̃_R, digits=4))
println("  avg_conn_prob = ", round(stats.avg_conn_prob, digits=4))
println("  E[edges] (node-centric) = ", round(stats.E_edges_node_centric, digits=2))
println("  E[edges] (edge-centric) = ", round(stats.E_edges_edge_centric, digits=4))

# Sanity checks
@assert stats.c_G > 0 && stats.c_R > 0 "Intensities should be positive"
@assert all(0 .<= stats.μ̃_G .<= 1) "Normalized mean should be in B^d_+"
@assert 0 <= stats.avg_conn_prob <= 1 "Connection probability should be in [0,1]"
println("  ✓ marginal_stats produces valid values")

# =============================================================================
# Summary
# =============================================================================
println("\n" * "=" ^ 60)
println("ALL VERIFICATION TESTS PASSED")
println("=" ^ 60)
println("\nFigures saved to: ", output_dir)
println("  - diffusion_spreading.png")
println("  - advection_transport.png")
println("  - reaction_diffusion.png")
