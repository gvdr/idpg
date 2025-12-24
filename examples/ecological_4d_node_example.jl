# Ecological Example: NODE-CENTRIC Food Web in 4D Latent Space
#
# This example uses the NODE-CENTRIC interpretation of IDPG:
#   1. Sample N sites from MixtureOfProductIntensities (each site has coupled (g, r))
#   2. Sites become "nodes" - ALL pairs (i, j) are considered for edges
#   3. Each pair forms an edge with probability g_i · r_j
#   4. Full (g, r) preserved for both source and target for clustering
#
# This is O(N²) complexity: N sites produce up to N² potential edges.
#
# Key model structure:
#   ρ(g,r) = Σ_m ρ_{G,m}(g) · ρ_{R,m}(r)
#
# Each species m has:
#   - ρ_{G,m}: niche in resource space (where it sits when being consumed)
#   - ρ_{R,m}: niche in consumer space (what it targets when consuming)
#   - γ_m = c_{G,m} · c_{R,m}: species-specific intensity (abundance)
#
# Structure:
#   Phase 1: Simulation - run models and save results to files
#   Phase 2: Visualization - load results and create figures
#
# Run with:
#   julia --project=. examples/ecological_4d_node_example.jl           # Run both phases
#   julia --project=. examples/ecological_4d_node_example.jl --viz     # Viz only (uses saved data)

using IDPG
using Random
using CairoMakie
using Statistics
using Graphs
using GraphMakie
using NetworkLayout
using LinearAlgebra
using Serialization

# Check for --viz flag (visualization only, skip simulation)
VIZ_ONLY = "--viz" in ARGS

# Output directory for saved data
DATA_DIR = "output/applications/data"
mkpath(DATA_DIR)

rng = MersenneTwister(161)

println("=" ^ 60)
println("IDPG Ecological Example: 4D NODE-CENTRIC Food Web")
println("  (Sites become nodes, all pairs considered for edges)")
if VIZ_ONLY
    println("  Mode: Visualization only (loading saved data)")
else
    println("  Mode: Full run (simulation + visualization)")
end
println("=" ^ 60)

# In 4D hyperspherical coordinates, we specify:
#   r ∈ [0, 1]: radius (network activity level)
#   φ₁, φ₂, φ₃ ∈ [0, π/2]: angles defining niche position
#
# Key insight for trophic structure:
#   - Small angle (≈0) concentrates weight in current dimension
#   - Large angle (≈π/2) passes weight to later dimensions
#
# To point to dimension k, use (k-1) large angles followed by small angles:
#   dim 1: [small, small, small]  → producers (food source)
#   dim 2: [large, small, small]  → herbivores
#   dim 3: [large, large, small]  → small predators
#   dim 4: [large, large, large]  → apex predators

guild_names = ["Producers", "Small Herb.", "Large Herb.", "Small Pred.", "Apex Pred."]

# =============================================================================
# Orthogonal Trophic Design
# =============================================================================
# Key insight: use orthogonality in 4D to enforce trophic structure
#
# Dimension usage:
#   Dim 1: Producer resources (g)
#   Dim 2: Herbivore resources (g)
#   Dim 3: Predator resources (g)
#   Dim 4: "Null consumption" - producers point r here, minimal g here
#
# This ensures:
#   - Producers don't consume (r orthogonal to all g)
#   - Each trophic level targets the one below
#   - Some noise/overlap for realism

# Resource roles (g) - where species sit when being consumed
# Each trophic level peaks in its designated dimension
# g[4] ≈ 0 for all species so producers (who target dim 4) consume almost nothing
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
for i in 1:5
    g_norm = sqrt(sum(means_G[i].^2))
    r_norm = sqrt(sum(means_R[i].^2))
    if g_norm > 1.0
        means_G[i] = means_G[i] ./ g_norm
    end
    if r_norm > 1.0
        means_R[i] = means_R[i] ./ r_norm
    end
end

# Note: Species scales will be computed from target E[N] using compute_scales_for_target_EN()
# This ensures we get meaningful numbers of interactions regardless of dimension

println("\nGuild structure (4D positions):")
for i in 1:5
    println("  ", guild_names[i], ":")
    println("    g = ", round.(means_G[i], digits=2))
    println("    r = ", round.(means_R[i], digits=2))
    g_dot_r = dot(means_G[i], means_R[i])
    println("    self-interaction (g·r) = ", round(g_dot_r, digits=3))
end

# Print expected interactions between guilds
println("\nExpected interaction strengths (g_i · r_j):")
print("           ")
for j in 1:5
    print(rpad(guild_names[j][1:6], 10))
end
println()
for i in 1:5
    print(rpad(guild_names[i], 11))
    for j in 1:5
        interaction = dot(means_G[i], means_R[j])
        print(rpad(string(round(interaction, digits=2)), 10))
    end
    println()
end

# Helper function to create MixtureOfProductIntensities from means and scales
function create_species_mixture(means_G, means_R, κ::Float64, scales)
    d = length(means_G[1])
    n_species = length(means_G)

    # Create per-species BdPlusMixture (single component each)
    ρ_Gs = BdPlusMixture{d}[]
    ρ_Rs = BdPlusMixture{d}[]

    for m in 1:n_species
        push!(ρ_Gs, BdPlusMixture([1.0], [means_G[m]], [κ], scales[m]))
        push!(ρ_Rs, BdPlusMixture([1.0], [means_R[m]], [κ], scales[m]))
    end

    return MixtureOfProductIntensities(ρ_Gs, ρ_Rs)
end

"""
Compute scales to achieve target E[N] with given relative abundances.

Theory:
  γ_m = c_{G,m} · c_{R,m} = s_m² · I_{G,m} · I_{R,m}
  E[N] = Σ_m γ_m

Given target E[N] and proportions p_m (Σp_m = 1):
  s_m = √(p_m · E[N] / (I_{G,m} · I_{R,m}))

where I_{G,m} = ∫ exp(-κ||g-μ_{G,m}||²) dg is computed via MC.
"""
function compute_scales_for_target_EN(means_G, means_R, κ::Float64,
                                       target_EN::Float64, proportions;
                                       n_samples::Int=50000,
                                       rng::AbstractRNG=Random.default_rng())
    d = length(means_G[1])
    M = length(means_G)
    vol = Bd_plus_volume(d)

    # Compute geometry-dependent integrals I_{G,m} and I_{R,m} via MC
    # Using unit scale (scale=1) to get the "shape" integral
    I_G = zeros(M)
    I_R = zeros(M)

    for _ in 1:n_samples
        x = uniform_Bd_plus_sample(d; rng=rng)
        for m in 1:M
            # Evaluate unnormalized Gaussian kernel (scale=1)
            dist_sq_G = sum((x[i] - means_G[m][i])^2 for i in 1:d)
            dist_sq_R = sum((x[i] - means_R[m][i])^2 for i in 1:d)
            I_G[m] += exp(-κ * dist_sq_G)
            I_R[m] += exp(-κ * dist_sq_R)
        end
    end

    # Convert to integrals
    I_G .*= vol / n_samples
    I_R .*= vol / n_samples

    # Compute scales: s_m = √(p_m · E[N] / (I_G[m] · I_R[m]))
    scales = zeros(M)
    for m in 1:M
        product = I_G[m] * I_R[m]
        if product > 1e-10
            scales[m] = sqrt(proportions[m] * target_EN / product)
        else
            scales[m] = 1.0  # Fallback for edge cases
            @warn "Species " * string(m) * " has near-zero integral product"
        end
    end

    return scales, I_G, I_R
end

# Compare 2D vs 4D using MixtureOfProductIntensities
# Scales computed from target E[N] and species proportions

# 2D scenario setup
means_G_2d = [Vector(Bd_plus_from_hyperspherical(0.85, [π/10])),
              Vector(Bd_plus_from_hyperspherical(0.85, [π/4])),
              Vector(Bd_plus_from_hyperspherical(0.85, [π/3]))]
means_R_2d = [Vector(Bd_plus_from_hyperspherical(0.85, [π/6])),
              Vector(Bd_plus_from_hyperspherical(0.85, [π/5])),
              Vector(Bd_plus_from_hyperspherical(0.85, [π/4]))]
κ_2d = 20.0
target_EN_2d = 50.0
proportions_2d = [0.5, 0.3, 0.2]

println("\nComputing scales for 2D scenario (target E[N] = ", target_EN_2d, ")...")
scales_2d, I_G_2d, I_R_2d = compute_scales_for_target_EN(
    means_G_2d, means_R_2d, κ_2d, target_EN_2d, proportions_2d; rng=rng)
println("  Geometry integrals I_G: ", round.(I_G_2d, digits=4))
println("  Geometry integrals I_R: ", round.(I_R_2d, digits=4))
println("  Computed scales: ", round.(scales_2d, digits=2))

# 4D scenario setup - use the guild means defined above
# With orthogonal design, moderate κ is sufficient
κ_4d = 50.0  # Moderate concentration (σ ≈ 0.1)
target_EN_4d = 100.0  # Target expected sites
proportions_4d = [0.30, 0.25, 0.20, 0.15, 0.10]  # Ecological pyramid, apex at 10% for visibility

println("\nComputing scales for 4D scenario (target E[N] = ", target_EN_4d, ")...")
scales_4d, I_G_4d, I_R_4d = compute_scales_for_target_EN(
    means_G, means_R, κ_4d, target_EN_4d, proportions_4d; rng=rng)
println("  Geometry integrals I_G: ", round.(I_G_4d, digits=4))
println("  Geometry integrals I_R: ", round.(I_R_4d, digits=4))
println("  Computed scales: ", round.(scales_4d, digits=2))

scenarios = [
    (
        name = "2D Latent Space",
        dim = 2,
        means_G = means_G_2d,
        means_R = means_R_2d,
        species_scales = scales_2d,
        κ = κ_2d,
    ),
    (
        name = "4D Latent Space",
        dim = 4,
        means_G = means_G,
        means_R = means_R,
        species_scales = scales_4d,
        κ = κ_4d,
    ),
]

# Helper function: assign site to nearest guild using FULL (g, r) signature
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

# Legacy: assign point to nearest guild mean (for backward compatibility)
function assign_to_nearest_guild(point, guild_means)
    min_dist = Inf
    best_guild = 1
    for (i, mean) in enumerate(guild_means)
        dist = norm(point .- mean)
        if dist < min_dist
            min_dist = dist
            best_guild = i
        end
    end
    return best_guild
end

fig = Figure(size=(1000, 900))

println("\n" * "=" ^ 60)
println("Comparing 2D vs 4D Latent Space (Node-Centric Model)")
println("=" ^ 60)

# Minimum edge count for graph visualization (filters noise)
comparison_min_count = 3.0

for (col, scenario) in enumerate(scenarios)
    println("\n--- ", scenario.name, " ---")

    n_guilds = length(scenario.species_scales)

    # Create MixtureOfProductIntensities (species have coupled G, R niches)
    ρ = create_species_mixture(scenario.means_G, scenario.means_R, scenario.κ, scenario.species_scales)

    # Use more MC samples for 4D case
    n_mc = scenario.dim == 4 ? 50000 : 10000
    stats = marginal_stats(ρ; n_samples=n_mc, rng=rng)
    println("  Dimension: ", scenario.dim)
    println("  Number of species: ", n_species(ρ))
    println("  Species intensities γ: ", round.(stats.γ, digits=1))
    println("  Total intensity C: ", round(stats.C, digits=1))
    println("  Expected interactions E[L]: ", round(stats.E_edges_edge_centric, digits=1))

    # Sample sites from MixtureOfProductIntensities (species labels + sites)
    labeled_sites = sample_ppp_mixture(ρ; n_samples=n_mc, rng=MersenneTwister(161 + col))
    sites = [site for (_, site) in labeled_sites]  # Extract just the InteractionSites
    # NODE-CENTRIC: sites become nodes, ALL N² pairs considered for edges
    interactions = generate_edge_centric_full(sites; rng=MersenneTwister(161 + col))
    println("  Sampled sites (nodes): ", length(sites))
    println("  Realized interactions: ", length(interactions))

    # Build full guild means for clustering
    guild_full_means = build_full_guild_means(scenario.means_G, scenario.means_R)

    # Assign to guilds using FULL (g, r) signature
    n_clusters = n_guilds
    edge_weights = zeros(Int, n_clusters, n_clusters)
    for k in 1:length(interactions)
        src_guild = assign_site_to_guild(interactions.source_sites[k], guild_full_means)
        tgt_guild = assign_site_to_guild(interactions.target_sites[k], guild_full_means)
        edge_weights[src_guild, tgt_guild] += 1
    end

    # Count filtered edges (above threshold, excluding diagonal)
    n_links = 0
    for i in 1:n_clusters
        for j in 1:n_clusters
            if i != j && edge_weights[i, j] >= comparison_min_count
                n_links += 1
            end
        end
    end
    println("  Guilds: ", n_clusters, ", Links (count≥", comparison_min_count, "): ", n_links)

    # Row 1: Food web structure (filtered by count threshold)
    ax1 = Axis(fig[1, col], aspect=DataAspect(),
        title = scenario.name * "\nFood web (" * string(n_links) * " links, count≥" * string(Int(comparison_min_count)) * ")")
    hidedecorations!(ax1)
    hidespines!(ax1)
    limits!(ax1, -1.5, 1.5, -1.5, 1.5)

    if n_clusters > 0 && n_links > 0
        food_web = SimpleDiGraph(n_clusters)
        edge_widths = Float64[]
        arrow_sizes = Float64[]
        max_weight = maximum(edge_weights)

        for i in 1:n_clusters
            for j in 1:n_clusters
                if i != j && edge_weights[i, j] >= comparison_min_count
                    add_edge!(food_web, i, j)  # Arrow: resource i → consumer j (energy flow)
                    w_normalized = edge_weights[i, j] / max_weight
                    width = 0.5 + 6.0 * w_normalized^0.6
                    push!(edge_widths, width)
                    push!(arrow_sizes, 5.0 + 18.0 * w_normalized^0.6)
                end
            end
        end

        node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_clusters]

        graphplot!(ax1, food_web,
            layout=Shell(),
            node_size=35,
            node_color=node_colors,
            nlabels=string.(1:n_clusters),
            nlabels_align=(:center, :center),
            nlabels_color=:white,
            nlabels_fontsize=14,
            edge_color=(:black, 0.8),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    end

    # Row 2: Edge weight matrix as heatmap (with colorbar below)
    ax2 = Axis(fig[2, col],
        title = "Interaction counts",
        xlabel = "Consumer",
        ylabel = "Resource",
        xticks = 1:n_clusters,
        yticks = 1:n_clusters,
        yreversed = true)

    # Exclude diagonal for visualization
    weights_no_diag = copy(edge_weights)
    for i in 1:n_clusters
        weights_no_diag[i, i] = 0
    end

    hm = heatmap!(ax2, permutedims(weights_no_diag), colormap=:viridis)
    Colorbar(fig[3, col], hm, label="Count", vertical=false)
end

# Add title
Label(fig[0, :], "Node-Centric Food Web: 2D vs 4D Latent Space", fontsize=18, font=:bold)

save("output/applications/ecological_4d_node_comparison.png", fig)
println("\nSaved ecological_4d_node_comparison.png")

# --- 2x2 Grid: Variability across realizations ---
println("\n" * "=" ^ 60)
println("4D Node-Centric Food Web Variability (2x2 grid of samples)")
println("=" ^ 60)

# Use the 4D scenario settings with MixtureOfProductIntensities
scenario_4d = scenarios[2]
n_guilds_4d = length(scenario_4d.species_scales)

ρ_4d = create_species_mixture(scenario_4d.means_G, scenario_4d.means_R, scenario_4d.κ, scenario_4d.species_scales)

fig_var = Figure(size=(1000, 1000))
Label(fig_var[0, :], "4D Node-Centric Food Web: Variability Across Samples", fontsize=18, font=:bold)

# Trophic layout for consistent visualization
trophic_levels_4d = [0, 1, 1, 2, 3]
x_offsets_4d = [0.0, -0.5, 0.5, 0.0, 0.0]
trophic_positions_4d = [Point2f(x_offsets_4d[i], trophic_levels_4d[i]) for i in 1:n_guilds_4d]
node_colors_4d = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_guilds_4d]

sample_seeds = [1001, 1002, 1003, 1004]  # Different seeds for variety

# Build full guild means for 4D scenario
guild_full_means_4d = build_full_guild_means(scenario_4d.means_G, scenario_4d.means_R)

for (idx, seed) in enumerate(sample_seeds)
    row = div(idx - 1, 2) + 1
    col = mod(idx - 1, 2) + 1

    # Sample sites from MixtureOfProductIntensities
    labeled_sites = sample_ppp_mixture(ρ_4d; n_samples=50000, rng=MersenneTwister(seed))
    sites = [site for (_, site) in labeled_sites]
    # NODE-CENTRIC: sites become nodes, ALL N² pairs considered for edges
    interactions = generate_edge_centric_full(sites; rng=MersenneTwister(seed))

    # Assign to guilds using full (g, r) signature
    edge_weights = zeros(Int, n_guilds_4d, n_guilds_4d)
    for k in 1:length(interactions)
        src_guild = assign_site_to_guild(interactions.source_sites[k], guild_full_means_4d)
        tgt_guild = assign_site_to_guild(interactions.target_sites[k], guild_full_means_4d)
        edge_weights[src_guild, tgt_guild] += 1
    end

    # Count edges (no threshold - show raw sample)
    n_edges = 0
    for i in 1:n_guilds_4d
        for j in 1:n_guilds_4d
            if i != j && edge_weights[i, j] > 0
                n_edges += 1
            end
        end
    end

    ax = Axis(fig_var[row, col], aspect=DataAspect(),
        title = "Sample " * string(idx) * " (" * string(length(interactions)) * " interactions, " * string(n_edges) * " edges)")
    hidedecorations!(ax)
    hidespines!(ax)
    limits!(ax, -1.5, 1.5, -0.5, 3.5)

    # Build and plot graph
    sample_web = SimpleDiGraph(n_guilds_4d)
    edge_widths = Float64[]
    arrow_sizes = Float64[]
    max_w = max(1, maximum(edge_weights))

    for i in 1:n_guilds_4d
        for j in 1:n_guilds_4d
            if i != j && edge_weights[i, j] > 0
                add_edge!(sample_web, i, j)  # Energy flow direction
                w_norm = edge_weights[i, j] / max_w
                push!(edge_widths, 0.5 + 4.0 * w_norm^0.6)
                push!(arrow_sizes, 6.0 + 14.0 * w_norm^0.6)
            end
        end
    end

    if ne(sample_web) > 0
        graphplot!(ax, sample_web,
            layout=(_) -> trophic_positions_4d,
            node_size=40,
            node_color=node_colors_4d,
            edge_color=(:black, 0.7),
            edge_width=edge_widths,
            arrow_size=arrow_sizes,
            arrow_show=true)
    else
        # Just show nodes if no edges
        for (i, pos) in enumerate(trophic_positions_4d)
            scatter!(ax, [pos], color=node_colors_4d[i], markersize=40)
        end
    end

    # Add trophic level labels on leftmost column
    if col == 1
        ax.yticks = ([0, 1, 2, 3], ["Prod.", "Herb.", "Pred.", "Apex"])
        ax.yticklabelsvisible = true
        ax.yticksvisible = false
    end

    println("  Sample ", idx, ": ", length(interactions), " interactions, ", n_edges, " edges")
end

save("output/applications/ecological_4d_node_variability.png", fig_var)
println("Saved ecological_4d_node_variability.png")

# --- Detailed 4D analysis ---
println("\n" * "=" ^ 60)
println("Detailed 4D Node-Centric Food Web Analysis")
println("=" ^ 60)

n_guilds = 5
n_trials = 50
DETAILED_DATA_FILE = DATA_DIR * "/ecological_4d_node_detailed.jls"

# Expected interaction matrix (computed from means - always available)
expected_matrix = zeros(n_guilds, n_guilds)
for i in 1:n_guilds
    for j in 1:n_guilds
        expected_matrix[i, j] = dot(means_G[i], means_R[j])
    end
end

if !VIZ_ONLY
    # =========== SIMULATION PHASE ===========
    κ_detailed = 50.0
    target_EN_detailed = 500.0
    proportions_detailed = [0.30, 0.25, 0.20, 0.15, 0.10]

    println("Computing scales for detailed analysis (target E[N] = ", target_EN_detailed, ")...")
    detailed_scales, _, _ = compute_scales_for_target_EN(
        means_G, means_R, κ_detailed, target_EN_detailed, proportions_detailed; rng=rng)
    println("  Computed scales: ", round.(detailed_scales, digits=2))

    ρ = create_species_mixture(means_G, means_R, κ_detailed, detailed_scales)

    stats = marginal_stats(ρ; n_samples=50000, rng=rng)
    println("Number of species: ", n_species(ρ))
    println("Species intensities γ: ", round.(stats.γ, digits=1))
    println("Total intensity C: ", round(stats.C, digits=1))
    println("Expected interactions E[L]: ", round(stats.E_edges_edge_centric, digits=1))

    # Build full guild means for detailed analysis
    guild_full_means_detailed = build_full_guild_means(means_G, means_R)

    interaction_matrices = []
    guild_counts_list = []
    edge_presence = zeros(n_guilds, n_guilds)

    for t in 1:n_trials
        labeled_sites = sample_ppp_mixture(ρ; n_samples=50000, rng=MersenneTwister(1610 + t))
        sites = [site for (_, site) in labeled_sites]
        # NODE-CENTRIC: sites become nodes, ALL N² pairs considered for edges
        interactions = generate_edge_centric_full(sites; rng=MersenneTwister(1610 + t))

        if length(interactions) >= 5
            guild_counts = zeros(Int, n_guilds)
            for site in sites
                g = assign_site_to_guild(site, guild_full_means_detailed)
                guild_counts[g] += 1
            end
            push!(guild_counts_list, guild_counts)

            edge_weights = zeros(n_guilds, n_guilds)
            for k in 1:length(interactions)
                src_guild = assign_site_to_guild(interactions.source_sites[k], guild_full_means_detailed)
                tgt_guild = assign_site_to_guild(interactions.target_sites[k], guild_full_means_detailed)
                edge_weights[src_guild, tgt_guild] += 1
            end
            push!(interaction_matrices, edge_weights)
            edge_presence .+= (edge_weights .> 0)
        end
    end

    println("Valid trials (with enough interactions): ", length(interaction_matrices), "/", n_trials)

    if isempty(interaction_matrices)
        println("Warning: No valid trials - try increasing scale")
        avg_matrix = zeros(n_guilds, n_guilds)
        avg_guild_counts = ones(n_guilds)
        edge_frequency = zeros(n_guilds, n_guilds)
    else
        avg_matrix = zeros(n_guilds, n_guilds)
        for mat in interaction_matrices
            avg_matrix .+= mat
        end
        avg_matrix ./= length(interaction_matrices)

        avg_guild_counts = zeros(n_guilds)
        for counts in guild_counts_list
            avg_guild_counts .+= counts
        end
        avg_guild_counts ./= length(guild_counts_list)
        println("Average guild counts: ", round.(avg_guild_counts, digits=1))

        edge_frequency = edge_presence ./ length(interaction_matrices)
    end

    # Per-capita normalization
    percapita_matrix = zeros(n_guilds, n_guilds)
    for i in 1:n_guilds
        for j in 1:n_guilds
            denom = avg_guild_counts[i] * avg_guild_counts[j]
            if denom > 0
                percapita_matrix[i, j] = avg_matrix[i, j] / denom
            end
        end
    end

    # Save results
    results = Dict(
        "avg_matrix" => avg_matrix,
        "avg_guild_counts" => avg_guild_counts,
        "percapita_matrix" => percapita_matrix,
        "edge_frequency" => edge_frequency,
        "n_trials" => n_trials,
    )
    serialize(DETAILED_DATA_FILE, results)
    println("Saved simulation results to ", DETAILED_DATA_FILE)

else
    # =========== LOAD SAVED DATA ===========
    println("Loading saved data from ", DETAILED_DATA_FILE)
    results = deserialize(DETAILED_DATA_FILE)
    avg_matrix = results["avg_matrix"]
    avg_guild_counts = results["avg_guild_counts"]
    percapita_matrix = results["percapita_matrix"]
    edge_frequency = results["edge_frequency"]
    n_trials = results["n_trials"]
    println("Loaded results from ", n_trials, " trials")
end

# =========== VISUALIZATION PHASE ===========
# (This runs in both modes)

# Visualization parameters (tweak these without re-running simulation!)
min_edge_count = 50.0  # Minimum average count for filtered graph (was 20.0)
min_percapita_rate = 0.25  # Alternative: filter by per-capita rate

significant_edges = avg_matrix .>= min_edge_count
n_significant = count(significant_edges) - count(i -> significant_edges[i,i], 1:n_guilds)
println("\nEdges with avg count ≥", min_edge_count, ": ", n_significant, " (excl. diagonal)")

println("\nAverage interaction counts (* = included in filtered graph, count ≥", min_edge_count, "):")
print("           ")
for j in 1:n_guilds
    print(rpad(guild_names[j][1:6], 10))
end
println()
for i in 1:n_guilds
    print(rpad(guild_names[i], 11))
    for j in 1:n_guilds
        count_str = string(round(avg_matrix[i, j], digits=1))
        # Mark significant edges (above threshold, excluding diagonal)
        if i != j && significant_edges[i, j]
            count_str = count_str * "*"
        end
        print(rpad(count_str, 10))
    end
    println()
end

println("\nPer-capita interaction rates (observed / (n_source × n_target)):")
print("           ")
for j in 1:n_guilds
    print(rpad(guild_names[j][1:6], 10))
end
println()
for i in 1:n_guilds
    print(rpad(guild_names[i], 11))
    for j in 1:n_guilds
        print(rpad(string(round(percapita_matrix[i, j], digits=3)), 10))
    end
    println()
end

# Create detailed figure with 3 panels
fig2 = Figure(size=(1400, 500))

# Panel 1: Average interaction counts (observed - raw)
ax1 = Axis(fig2[1, 1],
    title = "Observed Counts\n(avg over " * string(n_trials) * " trials)",
    xlabel = "Consumer guild",
    ylabel = "Resource guild",
    xticks = (1:n_guilds, guild_names),
    yticks = (1:n_guilds, guild_names),
    xticklabelrotation = π/4,
    yreversed = true)

hm = heatmap!(ax1, permutedims(avg_matrix), colormap=:YlOrRd)
Colorbar(fig2[1, 2], hm, label="Avg. count")

# Panel 2: Per-capita normalized (makes it comparable to g·r)
ax2 = Axis(fig2[1, 3],
    title = "Per-Capita Rate\n(count / (n_src × n_tgt))",
    xlabel = "Consumer guild",
    ylabel = "Resource guild",
    xticks = (1:n_guilds, guild_names),
    yticks = (1:n_guilds, guild_names),
    xticklabelrotation = π/4,
    yreversed = true)

hm2 = heatmap!(ax2, permutedims(percapita_matrix), colormap=:YlOrRd)
Colorbar(fig2[1, 4], hm2, label="Rate")

# Panel 3: Expected interaction strength (computed earlier from means_G, means_R)
ax3 = Axis(fig2[1, 5],
    title = "Expected g · r\n(theoretical)",
    xlabel = "Consumer guild",
    ylabel = "Resource guild",
    xticks = (1:n_guilds, guild_names),
    yticks = (1:n_guilds, guild_names),
    xticklabelrotation = π/4,
    yreversed = true)

hm3 = heatmap!(ax3, permutedims(expected_matrix), colormap=:YlOrRd)
Colorbar(fig2[1, 6], hm3, label="g · r")

save("output/applications/ecological_4d_node_detailed.png", fig2)
println("\nSaved ecological_4d_node_detailed.png")

# Create filtered food web graph using per-capita rate threshold
# This filters by interaction probability, not raw counts (better for trophic structure)
significant_by_rate = percapita_matrix .>= min_percapita_rate
n_significant_rate = count(significant_by_rate) - count(i -> significant_by_rate[i,i], 1:n_guilds)
println("\n--- Filtered Food Web (per-capita rate ≥", min_percapita_rate, ") ---")
println("Edges passing rate threshold: ", n_significant_rate, " (excl. diagonal)")

fig3 = Figure(size=(900, 700))
ax_graph = Axis(fig3[1, 1],
    title = "Filtered Food Web (per-capita rate ≥" * string(min_percapita_rate) * ")")
hidedecorations!(ax_graph)
hidespines!(ax_graph)

# Build filtered graph using per-capita rate threshold
filtered_web = SimpleDiGraph(n_guilds)
filtered_edge_widths = Float64[]
filtered_arrow_sizes = Float64[]

for i in 1:n_guilds
    for j in 1:n_guilds
        # Filter by per-capita rate (better reflects true interaction probability)
        if i != j && percapita_matrix[i, j] >= min_percapita_rate
            add_edge!(filtered_web, i, j)  # Arrow: resource i → consumer j (energy flow)
            # Scale by per-capita rate (not raw count)
            rate = percapita_matrix[i, j]
            push!(filtered_edge_widths, 1.0 + 8.0 * rate)
            push!(filtered_arrow_sizes, 10.0 + 20.0 * rate)
        end
    end
end

println("Filtered edges: ", ne(filtered_web))

if ne(filtered_web) > 0
    node_colors = [Makie.wong_colors()[mod1(i, 7)] for i in 1:n_guilds]

    # Custom trophic layout: y = trophic level, x = spread within level
    # Producers (1) at bottom, Apex (5) at top
    trophic_levels = [0, 1, 1, 2, 3]  # Producers=0, Herbivores=1, SmallPred=2, Apex=3
    x_offsets = [0.0, -0.5, 0.5, 0.0, 0.0]  # Spread herbivores horizontally
    trophic_positions = [Point2f(x_offsets[i], trophic_levels[i]) for i in 1:n_guilds]

    graphplot!(ax_graph, filtered_web,
        layout=(_) -> trophic_positions,
        node_size=50,
        node_color=node_colors,
        edge_color=(:black, 0.7),
        edge_width=filtered_edge_widths,
        arrow_size=filtered_arrow_sizes,
        arrow_show=true)

    limits!(ax_graph, -1.5, 1.5, -0.5, 3.5)
    ax_graph.yticks = ([0, 1, 2, 3], ["Producers", "Herbivores", "Predators", "Apex"])
    hidespines!(ax_graph)
    ax_graph.yticklabelsvisible = true
    ax_graph.yticksvisible = false
    ax_graph.ygridvisible = false
end

save("output/applications/ecological_4d_node_filtered_web.png", fig3)
println("Saved ecological_4d_node_filtered_web.png")

println("\n" * "=" ^ 60)
println("Ecological Interpretation (4D):")
println("=" ^ 60)
println("""
In this 4D food web model, each dimension represents a trophic level:
  - Dimension 1: Producer niche (plants, phytoplankton)
  - Dimension 2: Herbivore niche (grazers, filter feeders)
  - Dimension 3: Predator niche (secondary consumers)
  - Dimension 4: Apex predator niche (top predators)

Key design principles:
  - Resource role (g): species peak in their own trophic dimension
  - Consumer role (r): species target the dimension BELOW them
  - Radius controls visibility: apex have low g-radius (rarely eaten)
  - Radius controls consumption: producers have near-zero r-radius

Expected trophic structure (matrix should be heavy ABOVE diagonal):
  - Row 1 (Producers) eaten by Cols 2-3 (Herbivores): HIGH
  - Rows 2-3 (Herbivores) eaten by Cols 4-5 (Predators): HIGH
  - Row 4 (Small Pred) eaten by Col 5 (Apex): HIGH
  - Row 5 (Apex) eaten by anyone: LOW (top of chain)
  - Column 1 (Producers as consumers): ~0 (autotrophs don't eat)
  - Diagonal (self-consumption): LOW (species don't eat themselves)
""")

println("\nDone!")
