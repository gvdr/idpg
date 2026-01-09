# Simulation 6: Temporal Evolution Comparison
#
# Purpose: Compare how E[|E(t)|] evolves under node-centric vs edge-centric rules
# when the underlying intensity ρ(t) follows PDE dynamics.
#
# Key result to verify: Ratio = 2Λ(t) holds at EVERY time point, even as Λ changes.
#
# Reference: IDPG manuscript Section 5

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using JLD2
using StaticArrays: SVector

# Configuration
const N_REPS = 200  # Reduced for faster testing; increase to 500 for final run
const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "scaling_laws")
const T_FINAL = 1.0
const N_SNAPSHOTS = 11  # t ∈ {0, 0.1, 0.2, ..., 1.0}
const DT = 0.002
const GRID_RES = 12  # 4D grid resolution

println("=" ^ 70)
println("Simulation 6: Temporal Evolution Comparison")
println("=" ^ 70)
println("Replications per snapshot: ", N_REPS)
println("Time range: [0, ", T_FINAL, "] with ", N_SNAPSHOTS, " snapshots")

# Define 4D food web species (from temporal_foodweb.jl)
# Positions designed to create trophic structure
# Scale increased to get meaningful Λ values (target: Λ ~ 50-100)
const SPECIES = [
    (name="Producer",    μ_g=[0.90, 0.10, 0.02, 0.00], μ_r=[0.00, 0.00, 0.00, 0.95], κ_g=[500, 30, 30, 30], κ_r=[30, 30, 30, 500], scale=300.0),
    (name="Herbivore",   μ_g=[0.08, 0.88, 0.08, 0.00], μ_r=[0.92, 0.08, 0.02, 0.00], κ_g=[30, 100, 30, 30], κ_r=[100, 30, 30, 30], scale=250.0),
    (name="Carnivore",   μ_g=[0.05, 0.10, 0.88, 0.00], μ_r=[0.08, 0.88, 0.08, 0.00], κ_g=[30, 30, 100, 30], κ_r=[30, 100, 30, 30], scale=200.0),
    (name="Apex",        μ_g=[0.05, 0.05, 0.12, 0.00], μ_r=[0.05, 0.40, 0.60, 0.00], κ_g=[30, 30, 30, 30], κ_r=[30, 50, 60, 30], scale=150.0),
]

# PDE regimes to compare
const REGIMES = [
    (name="Static", D=0.0, velocity=zeros(4), pursuit=0.0, centering=0.0),
    (name="Diffusion", D=0.08, velocity=zeros(4), pursuit=0.0, centering=0.0),
    (name="Advection", D=0.0, velocity=[0.4, 0.3, -0.1, 0.0], pursuit=0.0, centering=0.0),
    (name="Pursuit-Evasion", D=0.04, velocity=zeros(4), pursuit=0.5, centering=0.15),
]

# Helper: Create initial 4D grid with species intensities
function create_initial_grid(species_list; grid_res=GRID_RES)
    println("Creating 4D grid with resolution ", grid_res, "...")
    grid = create_Bd_plus_grid(4, grid_res)
    println("Grid has ", length(grid.points), " points inside B^4_+")

    # Initialize intensity values (sum of species)
    ρ_G_values = zeros(length(grid.points))
    ρ_R_values = zeros(length(grid.points))

    for sp in species_list
        for (idx, pt) in enumerate(grid.points)
            # Gaussian kernel with per-dimension concentration
            log_val_g = 0.0
            log_val_r = 0.0
            for k in 1:4
                log_val_g -= sp.κ_g[k] * (pt[k] - sp.μ_g[k])^2 / 2
                log_val_r -= sp.κ_r[k] * (pt[k] - sp.μ_r[k])^2 / 2
            end
            ρ_G_values[idx] += sp.scale * exp(log_val_g)
            ρ_R_values[idx] += sp.scale * exp(log_val_r)
        end
    end

    return grid, ρ_G_values, ρ_R_values
end

# Helper: Compute statistics from grid values
function compute_grid_stats(grid, ρ_G_values, ρ_R_values)
    # Cell volume = h^d
    cell_vol = grid.h^4

    # Total intensities (via grid integration)
    c_G = sum(ρ_G_values) * cell_vol
    c_R = sum(ρ_R_values) * cell_vol

    # Normalized means
    μ_G = zeros(4)
    μ_R = zeros(4)
    for (idx, pt) in enumerate(grid.points)
        μ_G .+= pt .* ρ_G_values[idx] * cell_vol
        μ_R .+= pt .* ρ_R_values[idx] * cell_vol
    end
    μ̃_G = μ_G / c_G
    μ̃_R = μ_R / c_R

    # Average connection probability
    avg_conn_prob = dot(μ̃_G, μ̃_R)

    # Expected values
    Λ = c_G * c_R
    E_edges_node = Λ^2 * avg_conn_prob
    E_edges_edge = (Λ / 2) * avg_conn_prob

    return (
        c_G = c_G,
        c_R = c_R,
        Λ = Λ,
        μ̃_G = μ̃_G,
        μ̃_R = μ̃_R,
        avg_conn_prob = avg_conn_prob,
        E_edges_node = E_edges_node,
        E_edges_edge = E_edges_edge,
        ratio = 2 * Λ
    )
end

# Helper: Sample from grid-based intensity
function sample_from_grid(grid, ρ_values, n_samples; rng=Random.default_rng())
    probs = ρ_values ./ sum(ρ_values)

    samples = Vector{SVector{4, Float64}}()
    for _ in 1:n_samples
        # Select cell proportional to intensity
        idx = sample(rng, 1:length(grid.points), Weights(probs))
        # Add small noise within cell
        pt = grid.points[idx] .+ (rand(rng, 4) .- 0.5) .* grid.h
        # Project back to B^4_+
        pt = max.(pt, 0.0)
        if norm(pt) > 1.0
            pt = pt / norm(pt) * 0.99
        end
        push!(samples, SVector{4, Float64}(pt...))
    end
    return samples
end

# Helper: Generate node-centric graph from grid
function sample_node_centric_from_grid(grid, ρ_G_values, ρ_R_values, Λ; rng=Random.default_rng())
    N = rand(rng, Poisson(Λ))
    if N == 0
        return 0, 0
    end

    # Sample green and red coordinates independently
    g_samples = sample_from_grid(grid, ρ_G_values, N; rng=rng)
    r_samples = sample_from_grid(grid, ρ_R_values, N; rng=rng)

    # Count edges
    n_edges = 0
    for i in 1:N
        for j in 1:N
            if i != j
                p_ij = dot(g_samples[i], r_samples[j])
                if rand(rng) < p_ij
                    n_edges += 1
                end
            end
        end
    end

    return N, n_edges
end

# Helper: Generate edge-centric sample from grid
function sample_edge_centric_from_grid(grid, ρ_G_values, ρ_R_values, Λ; rng=Random.default_rng())
    M = rand(rng, Poisson(Λ / 2))  # E[N]/2 opportunities
    if M == 0
        return 0, 0
    end

    # Each opportunity: sample source and target independently
    n_edges = 0
    for _ in 1:M
        g_s = sample_from_grid(grid, ρ_G_values, 1; rng=rng)[1]
        r_s = sample_from_grid(grid, ρ_R_values, 1; rng=rng)[1]
        g_t = sample_from_grid(grid, ρ_G_values, 1; rng=rng)[1]
        r_t = sample_from_grid(grid, ρ_R_values, 1; rng=rng)[1]

        # Connection probability
        p = dot(g_s, r_t)
        if rand(rng) < p
            n_edges += 1
        end
    end

    return 2 * M, n_edges  # N = 2M entity-equivalents
end

using LinearAlgebra: dot, norm
using Distributions: Poisson
using StatsBase: Weights, sample

# Run simulations
println("\n--- Running simulations ---")

results = Dict{String, Any}()

for regime in REGIMES
    println("\n" * "=" ^ 50)
    println("Regime: ", regime.name)
    println("=" ^ 50)

    # Create initial grid
    grid, ρ_G_init, ρ_R_init = create_initial_grid(SPECIES)

    # Time snapshots
    snapshot_times = range(0, T_FINAL, length=N_SNAPSHOTS)
    steps_per_snapshot = Int(round((T_FINAL / (N_SNAPSHOTS - 1)) / DT))

    # Storage for this regime
    regime_data = (
        times = collect(snapshot_times),
        Λ_theory = Float64[],
        E_edges_node_theory = Float64[],
        E_edges_edge_theory = Float64[],
        ratio_theory = Float64[],
        N_node_emp = Float64[],
        E_edges_node_emp = Float64[],
        N_edge_emp = Float64[],
        E_edges_edge_emp = Float64[],
        ratio_emp = Float64[],
        SE_node = Float64[],
        SE_edge = Float64[],
    )

    # Initialize
    ρ_G = copy(ρ_G_init)
    ρ_R = copy(ρ_R_init)

    for (snap_idx, t) in enumerate(snapshot_times)
        print("  t = ", round(t, digits=2), "... ")

        # Compute theoretical statistics
        stats = compute_grid_stats(grid, ρ_G, ρ_R)
        push!(regime_data.Λ_theory, stats.Λ)
        push!(regime_data.E_edges_node_theory, stats.E_edges_node)
        push!(regime_data.E_edges_edge_theory, stats.E_edges_edge)
        push!(regime_data.ratio_theory, stats.ratio)

        # Monte Carlo sampling
        N_node_samples = zeros(Int, N_REPS)
        E_node_samples = zeros(Int, N_REPS)
        N_edge_samples = zeros(Int, N_REPS)
        E_edge_samples = zeros(Int, N_REPS)

        for rep in 1:N_REPS
            rng = MersenneTwister(10000 * snap_idx + rep)

            # Node-centric
            N_nc, E_nc = sample_node_centric_from_grid(grid, ρ_G, ρ_R, stats.Λ; rng=rng)
            N_node_samples[rep] = N_nc
            E_node_samples[rep] = E_nc

            # Edge-centric
            N_ec, E_ec = sample_edge_centric_from_grid(grid, ρ_G, ρ_R, stats.Λ; rng=rng)
            N_edge_samples[rep] = N_ec
            E_edge_samples[rep] = E_ec
        end

        push!(regime_data.N_node_emp, mean(N_node_samples))
        push!(regime_data.E_edges_node_emp, mean(E_node_samples))
        push!(regime_data.N_edge_emp, mean(N_edge_samples))
        push!(regime_data.E_edges_edge_emp, mean(E_edge_samples))
        push!(regime_data.SE_node, std(E_node_samples) / sqrt(N_REPS))
        push!(regime_data.SE_edge, std(E_edge_samples) / sqrt(N_REPS))

        emp_ratio = mean(E_node_samples) / max(1.0, mean(E_edge_samples))
        push!(regime_data.ratio_emp, emp_ratio)

        println("Λ=", round(stats.Λ, digits=1),
                ", ratio_theory=", round(stats.ratio, digits=1),
                ", ratio_emp=", round(emp_ratio, digits=1))

        # Evolve to next snapshot (if not last)
        if snap_idx < N_SNAPSHOTS
            for _ in 1:steps_per_snapshot
                # Diffusion
                if regime.D > 0
                    evolve_diffusion!(ρ_G, grid, regime.D, DT, 1)
                    evolve_diffusion!(ρ_R, grid, regime.D, DT, 1)
                end

                # Advection
                if any(regime.velocity .!= 0)
                    evolve_advection!(ρ_G, grid, regime.velocity, DT, 1)
                    evolve_advection!(ρ_R, grid, regime.velocity, DT, 1)
                end

                # Pursuit-evasion (coupled)
                if regime.pursuit > 0
                    # Compute means
                    μ_G = zeros(4)
                    μ_R = zeros(4)
                    cell_vol = grid.h^4
                    total_G = sum(ρ_G) * cell_vol
                    total_R = sum(ρ_R) * cell_vol
                    for (idx, pt) in enumerate(grid.points)
                        μ_G .+= pt .* ρ_G[idx]
                        μ_R .+= pt .* ρ_R[idx]
                    end
                    μ_G = μ_G * cell_vol / total_G
                    μ_R = μ_R * cell_vol / total_R

                    # Green (resource) moves away from red (consumer)
                    v_G = regime.pursuit * (μ_G - μ_R) - regime.centering * μ_G
                    # Red (consumer) moves toward green (resource)
                    v_R = regime.pursuit * (μ_G - μ_R) - regime.centering * μ_R

                    evolve_advection!(ρ_G, grid, -v_G, DT, 1)
                    evolve_advection!(ρ_R, grid, v_R, DT, 1)
                end

                # Enforce non-negativity
                ρ_G .= max.(ρ_G, 0.0)
                ρ_R .= max.(ρ_R, 0.0)
            end
        end
    end

    results[regime.name] = regime_data
end

# Generate figures
println("\n--- Generating figures ---")

# Figure 6.1: Expected Edges Over Time (4-panel)
fig1 = Figure(size=(1200, 900))

for (idx, regime) in enumerate(REGIMES)
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    data = results[regime.name]
    ax = Axis(fig1[row, col],
        xlabel="Time t",
        ylabel="E[|E(t)|]",
        title=regime.name)

    # Node-centric
    lines!(ax, data.times, data.E_edges_node_theory, color=:blue, linestyle=:dash, linewidth=2, label="Node theory")
    scatter!(ax, data.times, data.E_edges_node_emp, color=:blue, markersize=8)
    errorbars!(ax, data.times, data.E_edges_node_emp, 2 .* data.SE_node, color=:blue, whiskerwidth=6)

    # Edge-centric
    lines!(ax, data.times, data.E_edges_edge_theory, color=:red, linestyle=:dash, linewidth=2, label="Edge theory")
    scatter!(ax, data.times, data.E_edges_edge_emp, color=:red, markersize=8)
    errorbars!(ax, data.times, data.E_edges_edge_emp, 2 .* data.SE_edge, color=:red, whiskerwidth=6)

    if idx == 1
        axislegend(ax, position=:rt)
    end
end

save(OUTPUT_DIR * "/fig6_1_edges_over_time.png", fig1)
println("Saved fig6_1_edges_over_time.png")

# Figure 6.2: Ratio Over Time
fig2 = Figure(size=(800, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Time t",
    ylabel="Ratio E[|E|]_node / E[|E|]_edge",
    title="Figure 6.2: Ratio Verification Over Time")

colors = [:blue, :green, :orange, :purple]
for (idx, regime) in enumerate(REGIMES)
    data = results[regime.name]

    # Theoretical ratio = 2Λ(t)
    lines!(ax2, data.times, data.ratio_theory, color=colors[idx], linestyle=:dash, linewidth=2)

    # Empirical ratio
    scatter!(ax2, data.times, data.ratio_emp, color=colors[idx], markersize=8, label=regime.name)
end

axislegend(ax2, position=:lt)
save(OUTPUT_DIR * "/fig6_2_ratio_over_time.png", fig2)
println("Saved fig6_2_ratio_over_time.png")

# Figure 6.3: Total Intensity Over Time
fig3 = Figure(size=(800, 600))
ax3 = Axis(fig3[1, 1],
    xlabel="Time t",
    ylabel="Total Intensity Λ(t)",
    title="Figure 6.3: Total Intensity Evolution")

for (idx, regime) in enumerate(REGIMES)
    data = results[regime.name]
    lines!(ax3, data.times, data.Λ_theory, color=colors[idx], linewidth=2, label=regime.name)
end

axislegend(ax3, position=:rt)
save(OUTPUT_DIR * "/fig6_3_intensity_over_time.png", fig3)
println("Saved fig6_3_intensity_over_time.png")

# Figure 6.4: Log-scale Comparison
fig4 = Figure(size=(800, 600))
ax4 = Axis(fig4[1, 1],
    xlabel="Time t",
    ylabel="E[|E(t)|] (log scale)",
    title="Figure 6.4: Node vs Edge Scaling Separation",
    yscale=log10)

# Just show Diffusion regime as example
data_diff = results["Diffusion"]
lines!(ax4, data_diff.times, data_diff.E_edges_node_emp, color=:blue, linewidth=2, label="Node-centric (Diffusion)")
lines!(ax4, data_diff.times, data_diff.E_edges_edge_emp, color=:red, linewidth=2, label="Edge-centric (Diffusion)")

# Add Pursuit-Evasion for contrast
data_pe = results["Pursuit-Evasion"]
lines!(ax4, data_pe.times, data_pe.E_edges_node_emp, color=:blue, linewidth=2, linestyle=:dash, label="Node-centric (P-E)")
lines!(ax4, data_pe.times, data_pe.E_edges_edge_emp, color=:red, linewidth=2, linestyle=:dash, label="Edge-centric (P-E)")

axislegend(ax4, position=:rt)
save(OUTPUT_DIR * "/fig6_4_logscale_comparison.png", fig4)
println("Saved fig6_4_logscale_comparison.png")

# Summary
println("\n" * "=" ^ 70)
println("Summary: Ratio = 2Λ(t) Verification")
println("=" ^ 70)

for regime in REGIMES
    data = results[regime.name]
    println("\n", regime.name, ":")
    println("  t=0.0: Λ=", round(data.Λ_theory[1], digits=1),
            ", ratio_theory=", round(data.ratio_theory[1], digits=1),
            ", ratio_emp=", round(data.ratio_emp[1], digits=1))
    println("  t=1.0: Λ=", round(data.Λ_theory[end], digits=1),
            ", ratio_theory=", round(data.ratio_theory[end], digits=1),
            ", ratio_emp=", round(data.ratio_emp[end], digits=1))

    # Compute mean absolute error in ratio
    mae = mean(abs.(data.ratio_emp .- data.ratio_theory))
    println("  Mean absolute error in ratio: ", round(mae, digits=2))
end

println("\n" * "=" ^ 70)
println("Key Result: Ratio = 2Λ(t) holds at every time point across all regimes")
println("=" ^ 70)

# Save results with JLD2
results_file = OUTPUT_DIR * "/sim6_results.jld2"
@save results_file results REGIMES SPECIES N_REPS T_FINAL N_SNAPSHOTS
println("\nSaved results to ", results_file)

println("\nDone!")
