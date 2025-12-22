# Temporal Food Web Example: PDE Evolution of Ecological Niches
# Demonstrates how food web structure changes as species niches evolve
# under different dynamical regimes (diffusion, advection, pursuit-evasion).
#
# Key insight: The intensity functions ρ_G and ρ_R evolve according to PDEs,
# and the resulting food web structure emerges from the evolved distributions.

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using GraphMakie
using NetworkLayout
using LinearAlgebra

rng = MersenneTwister(42)

println("=" ^ 70)
println("IDPG Temporal Food Web: PDE Evolution of Ecological Niches")
println("=" ^ 70)

# =============================================================================
# Guild Configuration (same as ecological_4d_example.jl)
# =============================================================================

guild_names = ["Producers", "Small Herb.", "Large Herb.", "Small Pred.", "Apex Pred."]
n_guilds = 5

# Angle constants for orthogonal trophic levels
const small_angle = π/60
const large_angle = 29π/60

# Resource roles (g) - where species sit when being consumed
g_hyperspherical = [
    (0.95, [small_angle, small_angle, small_angle]),  # Producers: peak in dim 1
    (0.90, [large_angle, small_angle, small_angle]),  # Small herbivores: peak in dim 2
    (0.85, [large_angle, small_angle, small_angle]),  # Large herbivores: peak in dim 2
    (0.75, [large_angle, large_angle, small_angle]),  # Small predators: peak in dim 3
    (0.20, [large_angle, large_angle, large_angle]),  # Apex predators: peak in dim 4
]

# Consumer roles (r) - each level targets the level below
r_hyperspherical = [
    (0.02, [small_angle, small_angle, small_angle]),   # Producers: nearly zero consumption
    (0.95, [small_angle, small_angle, small_angle]),   # Small herbivores: target dim 1
    (0.90, [small_angle, small_angle, small_angle]),   # Large herbivores: target dim 1
    (0.90, [large_angle, small_angle, small_angle]),   # Small predators: target dim 2
    (0.95, [large_angle, large_angle, small_angle]),   # Apex: target dim 3
]

# Convert to Cartesian
means_G = [Vector(Bd_plus_from_hyperspherical(r, angles)) for (r, angles) in g_hyperspherical]
means_R = [Vector(Bd_plus_from_hyperspherical(r, angles)) for (r, angles) in r_hyperspherical]

# Weights (abundance)
weights_G = [0.35, 0.25, 0.20, 0.15, 0.05]
weights_R = [0.05, 0.25, 0.20, 0.30, 0.20]

println("\nGuild means (4D Cartesian):")
for i in 1:n_guilds
    println("  ", guild_names[i], ": g=", round.(means_G[i], digits=2), ", r=", round.(means_R[i], digits=2))
end

# =============================================================================
# Grid Setup for 4D
# =============================================================================

println("\n" * "=" ^ 70)
println("Creating 4D Grid")
println("=" ^ 70)

# Resolution for 4D grid
resolution_4d = 12

println("Creating 4D grid with resolution ", resolution_4d, "...")
grid_4d = create_Bd_plus_grid(4, resolution_4d)
n_grid_points = length(grid_4d.points)
println("Grid has ", n_grid_points, " points inside B^4_+")

# =============================================================================
# Helper Functions
# =============================================================================

# Assign point to nearest guild mean
function assign_to_nearest_guild(point, guild_means)
    min_dist = Inf
    best_guild = 1
    for (i, mean) in enumerate(guild_means)
        dist = norm(Vector(point) .- mean)
        if dist < min_dist
            min_dist = dist
            best_guild = i
        end
    end
    return best_guild
end

# Compute food web matrix from EdgeCentricSample
function compute_foodweb_matrix(sample::EdgeCentricSample, means_G, means_R)
    n_guilds = length(means_G)
    edge_weights = zeros(n_guilds, n_guilds)

    for k in 1:length(sample)
        src_guild = assign_to_nearest_guild(sample.sources[k], means_G)
        tgt_guild = assign_to_nearest_guild(sample.targets[k], means_R)
        edge_weights[src_guild, tgt_guild] += 1
    end

    return edge_weights
end

# =============================================================================
# Velocity Fields for Pursuit-Evasion
# =============================================================================

# Center of B^4_+ (approximate barycenter)
center_4d = [0.25, 0.25, 0.25, 0.25]

# Velocity field: move TOWARD target with centering force
function make_velocity_toward(μ_target, μ_center, speed, centering_strength)
    function v_toward(x)
        # Direction toward target
        dir_target = μ_target .- x
        dist_target = norm(dir_target)
        if dist_target > 0.01
            dir_target = dir_target ./ dist_target
        end

        # Centering force
        dir_center = μ_center .- x

        return speed .* dir_target .+ centering_strength .* dir_center
    end
    return v_toward
end

# Velocity field: move AWAY from threat with centering force
function make_velocity_away(μ_threat, μ_center, speed, centering_strength)
    function v_away(x)
        # Direction away from threat
        dir_away = x .- μ_threat
        dist_threat = norm(dir_away)
        if dist_threat > 0.01
            dir_away = dir_away ./ dist_threat
        end

        # Centering force
        dir_center = μ_center .- x

        return speed .* dir_away .+ centering_strength .* dir_center
    end
    return v_away
end

# =============================================================================
# PDE Regimes
# =============================================================================

println("\n" * "=" ^ 70)
println("Defining PDE Regimes")
println("=" ^ 70)

# Parameters
κ_base = 30.0
scale_base = 2000.0  # Scale for grid-based sampling (higher for more interactions)
D_diffusion = 0.08  # Diffusion coefficient (increased for visible spreading)
v_advection = [0.4, 0.3, -0.1, 0.0]  # Advection velocity (stronger, asymmetric shift)
pursuit_speed = 0.5  # Speed for pursuit-evasion (increased)
centering_strength = 0.15  # Pull toward center (reduced to allow more movement)

# Time parameters
dt = 0.002  # Smaller dt for stability with larger coefficients
t_final = 1.0  # Longer evolution time
n_snapshots = 5
snapshot_times = collect(range(0, t_final, length=n_snapshots))
steps_between_snapshots = Int(round((t_final / (n_snapshots - 1)) / dt))

regimes = [
    (name = "Static", type = :static),
    (name = "Diffusion", type = :diffusion),
    (name = "Advection", type = :advection),
    (name = "Pursuit-Evasion", type = :pursuit_evasion),
]

println("Regimes: ", [r.name for r in regimes])
println("Time: 0 to ", t_final, " with ", n_snapshots, " snapshots")
println("Steps between snapshots: ", steps_between_snapshots)

# =============================================================================
# Run Evolution for Each Regime
# =============================================================================

println("\n" * "=" ^ 70)
println("Running PDE Evolution")
println("=" ^ 70)

# Storage for results
all_results = Dict()
κ_vals = fill(κ_base, n_guilds)

for regime in regimes
    println("\n--- Regime: ", regime.name, " ---")

    # Initialize intensities from mixture using IDPG function
    ρ_G = initialize_grid_from_mixture(grid_4d, weights_G, means_G, κ_vals, scale_base)
    ρ_R = initialize_grid_from_mixture(grid_4d, weights_R, means_R, κ_vals, scale_base)

    # Storage for this regime
    foodweb_history = []
    interaction_counts = Float64[]
    mean_G_history = []
    mean_R_history = []

    for (snap_idx, t) in enumerate(snapshot_times)
        # Compute current means using IDPG function
        μ_G = compute_mean_position(ρ_G, grid_4d)
        μ_R = compute_mean_position(ρ_R, grid_4d)
        push!(mean_G_history, Vector(μ_G))
        push!(mean_R_history, Vector(μ_R))

        # Sample food web at this time point (multiple samples for averaging)
        n_samples = 15
        avg_foodweb = zeros(n_guilds, n_guilds)
        total_interactions = 0

        for s in 1:n_samples
            sample = sample_from_grid(ρ_G, ρ_R, grid_4d; rng=MersenneTwister(1000*snap_idx + s))
            if length(sample) > 0
                fw = compute_foodweb_matrix(sample, means_G, means_R)
                avg_foodweb .+= fw
                total_interactions += length(sample)
            end
        end
        avg_foodweb ./= n_samples

        push!(foodweb_history, avg_foodweb)
        push!(interaction_counts, total_interactions / n_samples)

        println("  t=", round(t, digits=2), ": ", round(total_interactions / n_samples, digits=1), " avg interactions")

        # Evolve to next snapshot (if not last)
        if snap_idx < n_snapshots
            if regime.type == :static
                # No evolution
            elseif regime.type == :diffusion
                evolve_diffusion!(ρ_G, grid_4d, D_diffusion, dt, steps_between_snapshots)
                evolve_diffusion!(ρ_R, grid_4d, D_diffusion, dt, steps_between_snapshots)
            elseif regime.type == :advection
                evolve_advection!(ρ_G, grid_4d, v_advection, dt, steps_between_snapshots)
                evolve_advection!(ρ_R, grid_4d, v_advection, dt, steps_between_snapshots)
            elseif regime.type == :pursuit_evasion
                # Coupled dynamics: prey flee consumers, consumers chase prey
                # ρ_G (resources) move AWAY from consumer mean
                v_flee = make_velocity_away(Vector(μ_R), center_4d, pursuit_speed, centering_strength)
                evolve_advection_field!(ρ_G, grid_4d, v_flee, dt, steps_between_snapshots)

                # ρ_R (consumers) move TOWARD resource mean
                v_chase = make_velocity_toward(Vector(μ_G), center_4d, pursuit_speed, centering_strength)
                evolve_advection_field!(ρ_R, grid_4d, v_chase, dt, steps_between_snapshots)
            end
        end
    end

    all_results[regime.name] = (
        foodweb_history = foodweb_history,
        interaction_counts = interaction_counts,
        mean_G_history = mean_G_history,
        mean_R_history = mean_R_history,
    )
end

# =============================================================================
# Visualization
# =============================================================================

println("\n" * "=" ^ 70)
println("Creating Visualizations")
println("=" ^ 70)

# Create output directory
mkpath("output/temporal_foodweb")

# Figure 1: Food web heatmaps at start and end for each regime
fig1 = Figure(size=(1400, 700))
Label(fig1[0, :], "Food Web Evolution: Start vs End", fontsize=20, font=:bold)

for (col, regime) in enumerate(regimes)
    results = all_results[regime.name]

    # Start (t=0)
    ax1 = Axis(fig1[1, col],
        title = regime.name * "\nt = 0",
        xlabel = col == 1 ? "Consumer" : "",
        ylabel = col == 1 ? "Resource" : "",
        xticks = 1:n_guilds,
        yticks = 1:n_guilds,
        yreversed = true)

    start_matrix = copy(results.foodweb_history[1])
    for i in 1:n_guilds
        start_matrix[i, i] = 0
    end
    heatmap!(ax1, permutedims(start_matrix), colormap=:YlOrRd)

    # End (t=t_final)
    ax2 = Axis(fig1[2, col],
        title = "t = " * string(t_final),
        xlabel = "Consumer",
        ylabel = col == 1 ? "Resource" : "",
        xticks = 1:n_guilds,
        yticks = 1:n_guilds,
        yreversed = true)

    end_matrix = copy(results.foodweb_history[end])
    for i in 1:n_guilds
        end_matrix[i, i] = 0
    end
    heatmap!(ax2, permutedims(end_matrix), colormap=:YlOrRd)
end

save("output/temporal_foodweb/regime_comparison_heatmaps.png", fig1)
println("Saved regime_comparison_heatmaps.png")

# Figure 2: Interaction counts over time
fig2 = Figure(size=(800, 500))
ax = Axis(fig2[1, 1],
    title = "Total Interactions Over Time",
    xlabel = "Time",
    ylabel = "Average Interactions per Sample")

colors = [:gray, :blue, :orange, :red]
for (idx, regime) in enumerate(regimes)
    results = all_results[regime.name]
    lines!(ax, snapshot_times, results.interaction_counts,
           color=colors[idx], linewidth=2, label=regime.name)
    scatter!(ax, snapshot_times, results.interaction_counts,
             color=colors[idx], markersize=8)
end
axislegend(ax, position=:rt)

save("output/temporal_foodweb/interaction_counts.png", fig2)
println("Saved interaction_counts.png")

# Figure 3: Food web graphs at final time
fig3 = Figure(size=(1200, 350))
Label(fig3[0, :], "Food Web Structure at t = " * string(t_final), fontsize=18, font=:bold)

# Trophic layout
trophic_levels = [0, 1, 1, 2, 3]
x_offsets = [0.0, -0.5, 0.5, 0.0, 0.0]
trophic_positions = [Point2f(x_offsets[i], trophic_levels[i]) for i in 1:n_guilds]
node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_guilds]

min_edge_threshold = 0.5

for (col, regime) in enumerate(regimes)
    results = all_results[regime.name]
    end_matrix = results.foodweb_history[end]

    ax_fw = Axis(fig3[1, col], aspect=DataAspect(), title=regime.name)
    hidedecorations!(ax_fw)
    hidespines!(ax_fw)
    limits!(ax_fw, -1.5, 1.5, -0.5, 3.5)

    # Build graph
    food_web = SimpleDiGraph(n_guilds)
    edge_widths = Float64[]
    arrow_sizes = Float64[]
    max_w = max(1.0, maximum(end_matrix))

    for i in 1:n_guilds
        for j in 1:n_guilds
            if i != j && end_matrix[i, j] >= min_edge_threshold
                add_edge!(food_web, j, i)  # consumer j → resource i
                w_norm = end_matrix[i, j] / max_w
                push!(edge_widths, 0.5 + 4.0 * w_norm^0.6)
                push!(arrow_sizes, 6.0 + 14.0 * w_norm^0.6)
            end
        end
    end

    if ne(food_web) > 0
        graphplot!(ax_fw, food_web,
            layout=(_) -> trophic_positions,
            node_size=35,
            node_color=node_colors,
            edge_color=(:black, 0.7),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    else
        for (i, pos) in enumerate(trophic_positions)
            scatter!(ax_fw, [pos], color=node_colors[i], markersize=35)
        end
    end

    if col == 1
        ax_fw.yticks = ([0, 1, 2, 3], ["Prod.", "Herb.", "Pred.", "Apex"])
        ax_fw.yticklabelsvisible = true
        ax_fw.yticksvisible = false
    end
end

save("output/temporal_foodweb/final_foodwebs.png", fig3)
println("Saved final_foodwebs.png")

# Figure 4: Mean position trajectories (first 2 dimensions projected)
fig4 = Figure(size=(1000, 500))
Label(fig4[0, :], "Niche Mean Trajectories (dims 1-2 projection)", fontsize=18, font=:bold)

ax_G = Axis(fig4[1, 1], title="Resource Means (ρ_G)", xlabel="Dim 1", ylabel="Dim 2")
ax_R = Axis(fig4[1, 2], title="Consumer Means (ρ_R)", xlabel="Dim 1", ylabel="Dim 2")

for (idx, regime) in enumerate(regimes)
    results = all_results[regime.name]

    xs_G = [μ[1] for μ in results.mean_G_history]
    ys_G = [μ[2] for μ in results.mean_G_history]
    xs_R = [μ[1] for μ in results.mean_R_history]
    ys_R = [μ[2] for μ in results.mean_R_history]

    lines!(ax_G, xs_G, ys_G, color=colors[idx], linewidth=2, label=regime.name)
    scatter!(ax_G, [xs_G[1]], [ys_G[1]], color=colors[idx], markersize=12, marker=:circle)
    scatter!(ax_G, [xs_G[end]], [ys_G[end]], color=colors[idx], markersize=12, marker=:star5)

    lines!(ax_R, xs_R, ys_R, color=colors[idx], linewidth=2, label=regime.name)
    scatter!(ax_R, [xs_R[1]], [ys_R[1]], color=colors[idx], markersize=12, marker=:circle)
    scatter!(ax_R, [xs_R[end]], [ys_R[end]], color=colors[idx], markersize=12, marker=:star5)
end

axislegend(ax_G, position=:rt)
text!(ax_G, Point2f(0.02, 0.98), text="● = start, ★ = end", fontsize=10,
      align=(:left, :top), space=:relative)

save("output/temporal_foodweb/mean_trajectories.png", fig4)
println("Saved mean_trajectories.png")

# =============================================================================
# Figure 5: Temporal evolution filmstrip for pursuit-evasion
# =============================================================================

println("\n--- Creating Pursuit-Evasion Filmstrip ---")

fig5 = Figure(size=(1400, 400))
Label(fig5[0, :], "Pursuit-Evasion: Food Web Over Time", fontsize=18, font=:bold)

results_pe = all_results["Pursuit-Evasion"]

for (col, snap_idx) in enumerate(1:n_snapshots)
    t = snapshot_times[snap_idx]
    fw_matrix = results_pe.foodweb_history[snap_idx]

    ax_pe = Axis(fig5[1, col], aspect=DataAspect(),
        title = "t = " * string(round(t, digits=2)))
    hidedecorations!(ax_pe)
    hidespines!(ax_pe)
    limits!(ax_pe, -1.5, 1.5, -0.5, 3.5)

    # Build graph
    food_web = SimpleDiGraph(n_guilds)
    edge_widths = Float64[]
    arrow_sizes = Float64[]
    max_w = max(1.0, maximum(fw_matrix))

    for i in 1:n_guilds
        for j in 1:n_guilds
            if i != j && fw_matrix[i, j] >= min_edge_threshold
                add_edge!(food_web, j, i)
                w_norm = fw_matrix[i, j] / max_w
                push!(edge_widths, 0.5 + 4.0 * w_norm^0.6)
                push!(arrow_sizes, 6.0 + 14.0 * w_norm^0.6)
            end
        end
    end

    if ne(food_web) > 0
        graphplot!(ax_pe, food_web,
            layout=(_) -> trophic_positions,
            node_size=30,
            node_color=node_colors,
            edge_color=(:black, 0.7),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    else
        for (i, pos) in enumerate(trophic_positions)
            scatter!(ax_pe, [pos], color=node_colors[i], markersize=30)
        end
    end

    if col == 1
        ax_pe.yticks = ([0, 1, 2, 3], ["P", "H", "Pr", "A"])
        ax_pe.yticklabelsvisible = true
        ax_pe.yticksvisible = false
    end
end

save("output/temporal_foodweb/pursuit_evasion_filmstrip.png", fig5)
println("Saved pursuit_evasion_filmstrip.png")

# =============================================================================
# Ecological Interpretation
# =============================================================================

println("\n" * "=" ^ 70)
println("Ecological Interpretation")
println("=" ^ 70)
println("""
This example demonstrates how food web structure changes under different
dynamical regimes for species niches:

1. STATIC: Baseline case - no niche evolution, food web is stable.

2. DIFFUSION: Niches spread over time (becoming more generalist).
   - Trophic interactions become less specific
   - Cross-trophic "noise" interactions may increase
   - Represents niche expansion / reduced specialization

3. ADVECTION: Niches shift uniformly in a direction.
   - All species shift together in latent space
   - Food web structure may be preserved if shift is parallel
   - Represents climate-driven niche displacement

4. PURSUIT-EVASION: Resources flee consumers, consumers chase resources.
   - Creates ecological "arms race" dynamics
   - Centering force prevents escape to corners of B^d_+
   - Food web structure may destabilize as niches diverge
   - Represents co-evolutionary dynamics

The grid-based PDE approach allows exact evolution of any dynamical regime,
with food web structure emerging from the evolved intensity distributions.
""")

println("\nDone!")
