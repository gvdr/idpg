# Temporal Food Web Example: PDE Evolution of Ecological Niches
# Demonstrates how food web structure changes as species niches evolve
# under different dynamical regimes (diffusion, advection, pursuit-evasion).
#
# Key insight: The intensity functions ρ_G and ρ_R evolve according to PDEs,
# and the resulting food web structure emerges from the evolved distributions.
#
# Structure:
#   Phase 1: Simulation - run PDE evolution and save results
#   Phase 2: Visualization - load results and create figures
#
# Run with:
#   julia --project=. examples/temporal_foodweb.jl           # Full run
#   julia --project=. examples/temporal_foodweb.jl --viz     # Viz only

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using GraphMakie
using NetworkLayout
using LinearAlgebra
using Serialization
using Distributions
using StatsBase

# Check for --viz flag
VIZ_ONLY = "--viz" in ARGS

# Output directories
DATA_DIR = "output/temporal_foodweb/data"
mkpath(DATA_DIR)
mkpath("output/temporal_foodweb")

rng = MersenneTwister(42)

println("=" ^ 70)
println("IDPG Temporal Food Web: PDE Evolution of Ecological Niches")
if VIZ_ONLY
    println("  Mode: Visualization only (loading saved data)")
else
    println("  Mode: Full run (simulation + visualization)")
end
println("=" ^ 70)

# =============================================================================
# Guild Configuration (MixtureOfProductIntensities-consistent design)
# =============================================================================
# Uses orthogonal trophic design from ecological_4d_example.jl where:
# - Dimension 4 is "null consumption" - producers point r here, no species has g here
# - This ensures producers have zero consumption (g · r_producer = 0 for all species)
# - Each species has a single proportion (not separate G and R weights)

guild_names = ["Producers", "Small Herb.", "Large Herb.", "Small Pred.", "Apex Pred."]
n_guilds = 5

# Resource roles (g) - where species sit when being consumed
# Each trophic level peaks in its designated dimension
# g[4] ≈ 0 for all species so producers (who target dim 4) consume nothing
means_G = [
    [0.90, 0.10, 0.02, 0.00],   # Producers: strong in dim 1
    [0.08, 0.88, 0.08, 0.00],   # Small herbivores: strong in dim 2
    [0.08, 0.82, 0.15, 0.00],   # Large herbivores: dim 2, some dim 3
    [0.05, 0.10, 0.88, 0.00],   # Small predators: strong in dim 3
    [0.05, 0.05, 0.12, 0.00],   # Apex predators: small g overall (rarely eaten)
]

# Consumer roles (r) - each level targets the dimension below
# Producers point to dim 4 which has NO resources → zero consumption
means_R = [
    [0.00, 0.00, 0.00, 0.95],   # Producers: ONLY dim 4 (orthogonal to all g!)
    [0.92, 0.08, 0.02, 0.00],   # Small herbivores: target dim 1 (producers)
    [0.88, 0.12, 0.02, 0.00],   # Large herbivores: target dim 1, some dim 2
    [0.08, 0.88, 0.08, 0.00],   # Small predators: target dim 2 (herbivores)
    [0.05, 0.40, 0.60, 0.00],   # Apex: target dims 2-3 (herbivores AND predators)
]

# Normalize to ensure they're in B^4_+ (norm ≤ 1)
for i in 1:n_guilds
    g_norm = sqrt(sum(means_G[i].^2))
    r_norm = sqrt(sum(means_R[i].^2))
    if g_norm > 1.0
        means_G[i] = means_G[i] ./ g_norm
    end
    if r_norm > 1.0
        means_R[i] = means_R[i] ./ r_norm
    end
end

# Species proportions (same for G and R - species have coupled niches)
# This is consistent with MixtureOfProductIntensities: γ_m = c_{G,m} · c_{R,m}
species_proportions = [0.30, 0.25, 0.20, 0.15, 0.10]  # Ecological pyramid

println("\nGuild means (4D Cartesian):")
for i in 1:n_guilds
    println("  ", guild_names[i], ": g=", round.(means_G[i], digits=2), ", r=", round.(means_R[i], digits=2))
end

# =============================================================================
# Grid Setup for 4D (only needed for simulation)
# =============================================================================

if !VIZ_ONLY
    println("\n" * "=" ^ 70)
    println("Creating 4D Grid")
    println("=" ^ 70)

    # Resolution for 4D grid
    resolution_4d = 12

    println("Creating 4D grid with resolution ", resolution_4d, "...")
    grid_4d = create_Bd_plus_grid(4, resolution_4d)
    n_grid_points = length(grid_4d.points)
    println("Grid has ", n_grid_points, " points inside B^4_+")
end

# =============================================================================
# Helper Functions
# =============================================================================

# Assign site to nearest guild using FULL (g, r) signature
# guild_full_means should be a vector of 2d-dimensional vectors: [means_G[i]; means_R[i]]
function assign_site_to_guild(site::InteractionSite, guild_full_means)
    site_full = vcat(Vector(site.g), Vector(site.r))  # 2d-dimensional
    min_dist = Inf
    best_guild = 1
    for (i, full_mean) in enumerate(guild_full_means)
        dist = norm(site_full .- full_mean)
        if dist < min_dist
            min_dist = dist
            best_guild = i
        end
    end
    return best_guild
end

# Build full (g, r) guild means from separate G and R means
function build_full_guild_means(means_G, means_R)
    return [vcat(means_G[i], means_R[i]) for i in 1:length(means_G)]
end

# Compute food web matrix from FullEdgeCentricSample using full (g,r) centroids
# guild_full_means is a vector of 2d-dimensional centroids (one per guild)
function compute_foodweb_matrix(sample::FullEdgeCentricSample, guild_full_means::Vector)
    n_guilds = length(guild_full_means)
    edge_weights = zeros(n_guilds, n_guilds)

    for k in 1:length(sample)
        # Assign source and target sites using their FULL (g, r) signatures
        src_guild = assign_site_to_guild(sample.source_sites[k], guild_full_means)
        tgt_guild = assign_site_to_guild(sample.target_sites[k], guild_full_means)
        edge_weights[src_guild, tgt_guild] += 1
    end

    return edge_weights
end

# Convenience method that builds full means from separate G and R
function compute_foodweb_matrix(sample::FullEdgeCentricSample, means_G, means_R)
    guild_full_means = build_full_guild_means(means_G, means_R)
    return compute_foodweb_matrix(sample, guild_full_means)
end

"""
Sample from per-species grids with COUPLED (g, r) per species.

Algorithm:
1. Compute species intensities γ_m = c_{G,m} · c_{R,m}
2. Sample N ~ Poisson(Σ γ_m)
3. For each site: sample species m, then sample g from ρ_{G,m} and r from ρ_{R,m}
4. For each pair (source, target): accept with probability g_source · r_target
"""
function sample_from_species_grids(ρ_G_species::Vector{Vector{Float64}},
                                    ρ_R_species::Vector{Vector{Float64}},
                                    grid::BdPlusGrid{d};
                                    rng::AbstractRNG=Random.default_rng()) where d
    n_species = length(ρ_G_species)
    h_d = grid.h^d

    # Compute species intensities γ_m = c_{G,m} · c_{R,m}
    c_G = [sum(ρ_G_species[m]) * h_d for m in 1:n_species]
    c_R = [sum(ρ_R_species[m]) * h_d for m in 1:n_species]
    γ = c_G .* c_R
    C = sum(γ)

    if C < 1e-10
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    # Sample number of site pairs (interaction opportunities)
    N = rand(rng, Poisson(C))
    if N == 0
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    # Normalize for sampling
    probs = γ ./ C
    p_G = [ρ_G_species[m] ./ max(1e-10, sum(ρ_G_species[m])) for m in 1:n_species]
    p_R = [ρ_R_species[m] ./ max(1e-10, sum(ρ_R_species[m])) for m in 1:n_species]

    source_sites = InteractionSite{d}[]
    target_sites = InteractionSite{d}[]

    for _ in 1:N
        # Sample source species and position
        m_src = sample(rng, 1:n_species, Weights(probs))
        g_idx = sample(rng, 1:length(grid.points), Weights(p_G[m_src]))
        r_idx = sample(rng, 1:length(grid.points), Weights(p_R[m_src]))
        source_site = InteractionSite{d}(grid.points[g_idx], grid.points[r_idx])

        # Sample target species and position
        m_tgt = sample(rng, 1:n_species, Weights(probs))
        g_idx = sample(rng, 1:length(grid.points), Weights(p_G[m_tgt]))
        r_idx = sample(rng, 1:length(grid.points), Weights(p_R[m_tgt]))
        target_site = InteractionSite{d}(grid.points[g_idx], grid.points[r_idx])

        # Accept with probability g_source · r_target
        p_connect = connection_probability(source_site.g, target_site.r)
        if rand(rng) < p_connect
            push!(source_sites, source_site)
            push!(target_sites, target_site)
        end
    end

    return FullEdgeCentricSample{d}(source_sites, target_sites)
end

# Legacy support: Compute food web matrix from EdgeCentricSample (partial info only)
function compute_foodweb_matrix(sample::EdgeCentricSample, means_G, means_R)
    n_guilds = length(means_G)
    edge_weights = zeros(n_guilds, n_guilds)

    # Fallback: assign sources by g, targets by r (loses half the info)
    for k in 1:length(sample)
        # Find nearest mean for source (using g only)
        src_point = Vector(sample.sources[k])
        src_guild = 1
        min_dist = Inf
        for (i, mean) in enumerate(means_G)
            dist = norm(src_point .- mean)
            if dist < min_dist
                min_dist = dist
                src_guild = i
            end
        end

        # Find nearest mean for target (using r only)
        tgt_point = Vector(sample.targets[k])
        tgt_guild = 1
        min_dist = Inf
        for (i, mean) in enumerate(means_R)
            dist = norm(tgt_point .- mean)
            if dist < min_dist
                min_dist = dist
                tgt_guild = i
            end
        end

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

RESULTS_FILE = DATA_DIR * "/temporal_results.jls"

# =============================================================================
# Run Evolution for Each Regime (or load saved data)
# =============================================================================

if !VIZ_ONLY
    println("\n" * "=" ^ 70)
    println("Running PDE Evolution")
    println("=" ^ 70)

    # Storage for results
    all_results = Dict()
    κ_vals = fill(κ_base, n_guilds)

    for regime in regimes
        println("\n--- Regime: ", regime.name, " ---")

        # Initialize PER-SPECIES intensity grids (preserves species coupling)
        # Each species has its own ρ_G and ρ_R that evolve together
        ρ_G_species = Vector{Vector{Float64}}(undef, n_guilds)
        ρ_R_species = Vector{Vector{Float64}}(undef, n_guilds)
        for m in 1:n_guilds
            # Single-component mixture for each species
            ρ_G_species[m] = initialize_grid_from_mixture(grid_4d, [1.0], [means_G[m]], [κ_base], scale_base * species_proportions[m])
            ρ_R_species[m] = initialize_grid_from_mixture(grid_4d, [1.0], [means_R[m]], [κ_base], scale_base * species_proportions[m])
        end

        # Track per-species centroids as they evolve
        current_means_G = [copy(m) for m in means_G]
        current_means_R = [copy(m) for m in means_R]

        # Storage for this regime
        foodweb_history = []
        interaction_counts = Float64[]
        mean_G_history = []
        mean_R_history = []
        evolved_centroids_G_history = []
        evolved_centroids_R_history = []

        for (snap_idx, t) in enumerate(snapshot_times)
            # Compute overall means from per-species grids (for pursuit-evasion dynamics)
            # Weighted average across species
            ρ_G_total = sum(ρ_G_species)
            ρ_R_total = sum(ρ_R_species)
            μ_G = compute_mean_position(ρ_G_total, grid_4d)
            μ_R = compute_mean_position(ρ_R_total, grid_4d)
            push!(mean_G_history, Vector(μ_G))
            push!(mean_R_history, Vector(μ_R))

            # Store current species centroids
            push!(evolved_centroids_G_history, [copy(m) for m in current_means_G])
            push!(evolved_centroids_R_history, [copy(m) for m in current_means_R])

            # Build full (g, r) centroids for guild assignment
            current_full_means = build_full_guild_means(current_means_G, current_means_R)

            # Sample food web using COUPLED per-species sampling
            n_samples_per_trial = 15
            avg_foodweb = zeros(n_guilds, n_guilds)
            total_interactions = 0

            for s in 1:n_samples_per_trial
                rng_s = MersenneTwister(1000*snap_idx + s)

                # Sample with species coupling: for each site, pick species first,
                # then sample (g, r) from that species's distributions
                sample = sample_from_species_grids(ρ_G_species, ρ_R_species, grid_4d; rng=rng_s)

                if length(sample) > 0
                    fw = compute_foodweb_matrix(sample, current_full_means)
                    avg_foodweb .+= fw
                    total_interactions += length(sample)
                end
            end
            avg_foodweb ./= n_samples_per_trial

            push!(foodweb_history, avg_foodweb)
            push!(interaction_counts, total_interactions / n_samples_per_trial)

            println("  t=", round(t, digits=2), ": ", round(total_interactions / n_samples_per_trial, digits=1), " avg interactions")

            # Evolve per-species grids and centroids to next snapshot
            if snap_idx < n_snapshots
                Δt_total = dt * steps_between_snapshots

                if regime.type == :static
                    # No evolution
                elseif regime.type == :diffusion
                    for m in 1:n_guilds
                        evolve_diffusion!(ρ_G_species[m], grid_4d, D_diffusion, dt, steps_between_snapshots)
                        evolve_diffusion!(ρ_R_species[m], grid_4d, D_diffusion, dt, steps_between_snapshots)
                    end
                    # Diffusion: centroids don't move
                elseif regime.type == :advection
                    for m in 1:n_guilds
                        evolve_advection!(ρ_G_species[m], grid_4d, v_advection, dt, steps_between_snapshots)
                        evolve_advection!(ρ_R_species[m], grid_4d, v_advection, dt, steps_between_snapshots)
                        # Centroids move with velocity
                        current_means_G[m] = project_to_Bd_plus(current_means_G[m] .+ v_advection .* Δt_total)
                        current_means_R[m] = project_to_Bd_plus(current_means_R[m] .+ v_advection .* Δt_total)
                    end
                elseif regime.type == :pursuit_evasion
                    # All species flee/chase based on overall mean
                    v_flee = make_velocity_away(Vector(μ_R), center_4d, pursuit_speed, centering_strength)
                    v_chase = make_velocity_toward(Vector(μ_G), center_4d, pursuit_speed, centering_strength)

                    for m in 1:n_guilds
                        evolve_advection_field!(ρ_G_species[m], grid_4d, v_flee, dt, steps_between_snapshots)
                        evolve_advection_field!(ρ_R_species[m], grid_4d, v_chase, dt, steps_between_snapshots)

                        # Evolve centroids
                        v_G = v_flee(current_means_G[m])
                        v_R = v_chase(current_means_R[m])
                        current_means_G[m] = project_to_Bd_plus(current_means_G[m] .+ v_G .* Δt_total)
                        current_means_R[m] = project_to_Bd_plus(current_means_R[m] .+ v_R .* Δt_total)
                    end
                end
            end
        end

        all_results[regime.name] = (
            foodweb_history = foodweb_history,
            interaction_counts = interaction_counts,
            mean_G_history = mean_G_history,
            mean_R_history = mean_R_history,
            evolved_centroids_G = evolved_centroids_G_history,
            evolved_centroids_R = evolved_centroids_R_history,
        )
    end

    # Save results
    serialize(RESULTS_FILE, all_results)
    println("\nSaved simulation results to ", RESULTS_FILE)

else
    # Load saved results
    println("\n" * "=" ^ 70)
    println("Loading saved results")
    println("=" ^ 70)
    all_results = deserialize(RESULTS_FILE)
    println("Loaded results from ", RESULTS_FILE)
end

# =============================================================================
# Visualization
# =============================================================================

println("\n" * "=" ^ 70)
println("Creating Visualizations")
println("=" ^ 70)

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

# Figure 3: Food web graphs at final time (using per-capita rates)
fig3 = Figure(size=(1200, 350))
Label(fig3[0, :], "Food Web Structure at t = " * string(t_final), fontsize=18, font=:bold)

# Trophic layout
trophic_levels = [0, 1, 1, 2, 3]
x_offsets = [0.0, -0.5, 0.5, 0.0, 0.0]
trophic_positions = [Point2f(x_offsets[i], trophic_levels[i]) for i in 1:n_guilds]
node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_guilds]

# Per-capita rate threshold (same as ecological example)
min_percapita_rate = 0.25

for (col, regime) in enumerate(regimes)
    results = all_results[regime.name]
    end_matrix = results.foodweb_history[end]

    # Estimate guild counts from species proportions
    # For grid-based sampling, total interactions ≈ (c_G × c_R) × E[g·r]
    # We approximate guild counts as proportional to species_proportions
    total_counts = sum(end_matrix)
    est_guild_counts = species_proportions .* sqrt(total_counts)  # Approx sqrt(N) per guild

    # Compute per-capita matrix
    percapita_matrix = zeros(n_guilds, n_guilds)
    for i in 1:n_guilds
        for j in 1:n_guilds
            denom = est_guild_counts[i] * est_guild_counts[j]
            if denom > 0
                percapita_matrix[i, j] = end_matrix[i, j] / denom
            end
        end
    end

    ax_fw = Axis(fig3[1, col], aspect=DataAspect(), title=regime.name)
    hidedecorations!(ax_fw)
    hidespines!(ax_fw)
    limits!(ax_fw, -1.5, 1.5, -0.5, 3.5)

    # Build graph using per-capita rates
    food_web = SimpleDiGraph(n_guilds)
    edge_widths = Float64[]
    arrow_sizes = Float64[]

    # Filter by per-capita rate (same as ecological example)
    for i in 1:n_guilds
        for j in 1:n_guilds
            if i != j && percapita_matrix[i, j] >= min_percapita_rate
                # Arrow: resource i → consumer j (energy flow direction)
                add_edge!(food_web, i, j)
                rate = percapita_matrix[i, j]
                push!(edge_widths, 0.5 + 4.0 * rate^0.6)
                push!(arrow_sizes, 6.0 + 14.0 * rate^0.6)
            end
        end
    end

    if ne(food_web) > 0
        graphplot!(ax_fw, food_web,
            layout=(_) -> trophic_positions,
            node_size=35,
            node_color=node_colors,
            edge_color=(:black, 0.5),
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

    # Compute per-capita matrix for this snapshot
    total_counts = sum(fw_matrix)
    est_guild_counts = species_proportions .* sqrt(max(1.0, total_counts))
    percapita_matrix = zeros(n_guilds, n_guilds)
    for i in 1:n_guilds
        for j in 1:n_guilds
            denom = est_guild_counts[i] * est_guild_counts[j]
            if denom > 0
                percapita_matrix[i, j] = fw_matrix[i, j] / denom
            end
        end
    end

    ax_pe = Axis(fig5[1, col], aspect=DataAspect(),
        title = "t = " * string(round(t, digits=2)))
    hidedecorations!(ax_pe)
    hidespines!(ax_pe)
    limits!(ax_pe, -1.5, 1.5, -0.5, 3.5)

    # Build graph using per-capita rates
    food_web = SimpleDiGraph(n_guilds)
    edge_widths = Float64[]
    arrow_sizes = Float64[]

    for i in 1:n_guilds
        for j in 1:n_guilds
            if i != j && percapita_matrix[i, j] >= min_percapita_rate
                # Arrow: resource i → consumer j (energy flow direction)
                add_edge!(food_web, i, j)
                rate = percapita_matrix[i, j]
                push!(edge_widths, 0.5 + 4.0 * rate^0.6)
                push!(arrow_sizes, 6.0 + 14.0 * rate^0.6)
            end
        end
    end

    if ne(food_web) > 0
        graphplot!(ax_pe, food_web,
            layout=(_) -> trophic_positions,
            node_size=30,
            node_color=node_colors,
            edge_color=(:black, 0.5),
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
