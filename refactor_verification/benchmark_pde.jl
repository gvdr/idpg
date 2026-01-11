# Benchmark: OrdinaryDiffEq vs Explicit Euler (archived implementation)
# Compares accuracy and performance for PDE evolution

using IDPG
using LinearAlgebra
using Statistics
using Printf

println("=" ^ 60)
println("BENCHMARK: PDE Evolution (OrdinaryDiffEq vs Explicit Euler)")
println("=" ^ 60)

# =============================================================================
# Setup: Create grid and initial condition
# =============================================================================

d = 2
resolution = 20
grid = create_Bd_plus_grid(d, resolution)

println("\nGrid info:")
println("  Dimension: ", d)
println("  Resolution: ", resolution)
println("  Points in B^d_+: ", length(grid.points))
println("  Grid spacing h: ", grid.h)

# Gaussian initial condition centered at (0.5, 0.5)
μ₀ = [0.5, 0.5]
σ = 0.15
ρ₀ = [exp(-sum((p .- μ₀).^2) / (2σ^2)) for p in grid.points]
ρ₀ ./= sum(ρ₀) * grid.h^d  # Normalize to integrate to ~1

println("  Initial total intensity: ", sum(ρ₀) * grid.h^d)

# =============================================================================
# Test 1: Diffusion Evolution
# =============================================================================
println("\n--- DIFFUSION TEST (D=0.01, t=0→0.5) ---\n")

D = 0.01
tspan = (0.0, 0.5)

# New implementation (OrdinaryDiffEq)
println("Running OrdinaryDiffEq solver...")
t_ode = @elapsed begin
    sol_ode = evolve_diffusion(copy(ρ₀), grid, D, tspan; saveat=0.1)
end
ρ_final_ode = max.(sol_ode[end], 0.0)

# Explicit Euler (manual implementation for comparison)
println("Running explicit Euler (for comparison)...")
function explicit_euler_diffusion(ρ₀, grid, D, dt, n_steps)
    ρ = copy(ρ₀)
    ρ_new = similar(ρ)
    for _ in 1:n_steps
        for i in eachindex(ρ)
            ∇²ρ = laplacian_stencil(grid, ρ, i)
            ρ_new[i] = ρ[i] + dt * D * ∇²ρ
        end
        @. ρ_new = max(ρ_new, 0.0)
        copyto!(ρ, ρ_new)
    end
    return ρ
end

# For explicit Euler, we need small dt for stability
# CFL condition: dt < h² / (2*d*D)
dt_stable = 0.4 * grid.h^2 / (2 * d * D)
n_steps = ceil(Int, 0.5 / dt_stable)
t_euler = @elapsed begin
    ρ_final_euler = explicit_euler_diffusion(copy(ρ₀), grid, D, dt_stable, n_steps)
end

# Compare results
diff_norm = norm(ρ_final_ode - ρ_final_euler)
rel_diff = diff_norm / norm(ρ_final_euler)

println("\nResults:")
@printf("  OrdinaryDiffEq time:     %.4f ms\n", t_ode * 1000)
@printf("  Explicit Euler time:     %.4f ms (dt=%.2e, %d steps)\n", t_euler * 1000, dt_stable, n_steps)
@printf("  Speedup:                 %.2fx\n", t_euler / t_ode)
println()
@printf("  Final total intensity (ODE):   %.6f\n", sum(ρ_final_ode) * grid.h^d)
@printf("  Final total intensity (Euler): %.6f\n", sum(ρ_final_euler) * grid.h^d)
@printf("  Relative difference:           %.2e\n", rel_diff)

# =============================================================================
# Test 2: Advection Evolution
# =============================================================================
println("\n--- ADVECTION TEST (v⃗=[0.5, 0.3], t=0→0.3) ---\n")

v⃗ = [0.5, 0.3]
tspan_adv = (0.0, 0.3)

# OrdinaryDiffEq
println("Running OrdinaryDiffEq solver...")
t_ode_adv = @elapsed begin
    sol_ode_adv = evolve_advection(copy(ρ₀), grid, v⃗, tspan_adv; saveat=0.1)
end
ρ_final_ode_adv = max.(sol_ode_adv[end], 0.0)

# Explicit Euler for advection
function explicit_euler_advection(ρ₀, grid, v⃗, dt, n_steps)
    ρ = copy(ρ₀)
    ρ_new = similar(ρ)
    d = length(grid.points[1])
    for _ in 1:n_steps
        for i in eachindex(ρ)
            advection = 0.0
            for dim in 1:d
                ∂ρ = gradient_component(grid, ρ, i, dim)
                advection -= v⃗[dim] * ∂ρ
            end
            ρ_new[i] = ρ[i] + dt * advection
        end
        @. ρ_new = max(ρ_new, 0.0)
        copyto!(ρ, ρ_new)
    end
    return ρ
end

# CFL for advection: dt < h / max(|v|)
dt_stable_adv = 0.4 * grid.h / maximum(abs.(v⃗))
n_steps_adv = ceil(Int, 0.3 / dt_stable_adv)
t_euler_adv = @elapsed begin
    ρ_final_euler_adv = explicit_euler_advection(copy(ρ₀), grid, v⃗, dt_stable_adv, n_steps_adv)
end

diff_norm_adv = norm(ρ_final_ode_adv - ρ_final_euler_adv)
rel_diff_adv = diff_norm_adv / norm(ρ_final_euler_adv)

println("\nResults:")
@printf("  OrdinaryDiffEq time:     %.4f ms\n", t_ode_adv * 1000)
@printf("  Explicit Euler time:     %.4f ms (dt=%.2e, %d steps)\n", t_euler_adv * 1000, dt_stable_adv, n_steps_adv)
@printf("  Speedup:                 %.2fx\n", t_euler_adv / t_ode_adv)
println()
@printf("  Final total intensity (ODE):   %.6f\n", sum(ρ_final_ode_adv) * grid.h^d)
@printf("  Final total intensity (Euler): %.6f\n", sum(ρ_final_euler_adv) * grid.h^d)
@printf("  Relative difference:           %.2e\n", rel_diff_adv)

# =============================================================================
# Test 3: Different solvers comparison
# =============================================================================
println("\n--- SOLVER COMPARISON (diffusion, various ODE solvers) ---\n")

using OrdinaryDiffEq: Tsit5, ROCK4, Rodas5P, RK4

solvers = [
    ("Tsit5 (default)", Tsit5()),
    ("RK4 (fixed step)", RK4()),
    ("ROCK4 (stiff)", ROCK4()),
]

for (name, solver) in solvers
    t = @elapsed begin
        sol = evolve_diffusion(copy(ρ₀), grid, D, tspan; solver=solver)
    end
    ρ_final = max.(sol[end], 0.0)
    total_int = sum(ρ_final) * grid.h^d
    @printf("  %-20s: %.4f ms, final_intensity = %.6f\n", name, t * 1000, total_int)
end

# =============================================================================
# Summary
# =============================================================================
println("\n" * "=" ^ 60)
println("SUMMARY")
println("=" ^ 60)
println("""
OrdinaryDiffEq advantages:
  + Adaptive time-stepping (no manual CFL calculation)
  + Error control built-in
  + Choice of solvers for different problem types
  + Stiff solvers available (ROCK4, Rodas5P)
  + Generally faster due to larger adaptive steps

Explicit Euler:
  + Simple implementation
  - Requires manual CFL stability calculation
  - Many small steps needed for stability
  - No error control
""")
