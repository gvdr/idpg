# Poisson Point Process sampling on Ω = G × R ⊆ B^d_+ × B^d_+

"""
    InteractionSite{d}

An interaction site sampled from the PPP on Ω = G × R.

# Fields
- `g::LatentPoint{d}`: Green coordinate (propensity to propose connections)
- `r::LatentPoint{d}`: Red coordinate (propensity to accept connections)
"""
struct InteractionSite{d}
    g::LatentPoint{d}
    r::LatentPoint{d}
end

function Base.show(io::IO, site::InteractionSite{d}) where d
    print(io, "InteractionSite{" * string(d) * "}(g=" * string(site.g) * ", r=" * string(site.r) * ")")
end

"""
    estimate_max_intensity(ρ::AbstractIntensity{d}; n_samples=1000, rng=Random.default_rng()) -> Float64

Estimate the supremum of intensity ρ over Ω via sampling.
Used for thinning algorithm.
"""
function estimate_max_intensity(ρ::AbstractIntensity{d};
                                n_samples::Int=1000,
                                rng::AbstractRNG=Random.default_rng()) where d
    max_val = 0.0
    for _ in 1:n_samples
        g = uniform_Bd_plus_sample(d; rng=rng)
        r = uniform_Bd_plus_sample(d; rng=rng)
        val = ρ(g, r)
        max_val = max(max_val, val)
    end
    # Add safety margin
    return 1.5 * max_val
end

"""
    estimate_max_intensity(ρ::ProductIntensity{d}; n_samples=1000, rng=Random.default_rng()) -> Float64

Estimate max intensity for product case more efficiently.
"""
function estimate_max_intensity(ρ::ProductIntensity{d};
                                n_samples::Int=1000,
                                rng::AbstractRNG=Random.default_rng()) where d
    # For product intensity, max is max_G * max_R
    max_G = 0.0
    max_R = 0.0
    for _ in 1:n_samples
        g = uniform_Bd_plus_sample(d; rng=rng)
        r = uniform_Bd_plus_sample(d; rng=rng)
        max_G = max(max_G, ρ.ρ_G(g))
        max_R = max(max_R, ρ.ρ_R(r))
    end
    return 1.5 * max_G * max_R
end

"""
    sample_ppp(ρ::AbstractIntensity{d}; λ_max=nothing, rng=Random.default_rng()) -> Vector{InteractionSite{d}}

Sample from inhomogeneous PPP on Ω = B^d_+ × B^d_+ using thinning.

# Algorithm
1. Find λ_max = sup ρ(g,r) over Ω
2. Sample N ~ Poisson(λ_max · |Ω|) candidate points uniformly on Ω
3. Accept each point (g,r) with probability ρ(g,r) / λ_max

# Arguments
- `ρ`: Intensity function
- `λ_max`: Upper bound on intensity (estimated if not provided)
- `rng`: Random number generator

# Returns
Vector of accepted interaction sites.
"""
function sample_ppp(ρ::AbstractIntensity{d};
                    λ_max::Union{Nothing, Float64}=nothing,
                    rng::AbstractRNG=Random.default_rng()) where d

    # Estimate λ_max if not provided
    if isnothing(λ_max)
        λ_max = estimate_max_intensity(ρ; rng=rng)
    end

    # Volume of Ω = B^d_+ × B^d_+
    vol_Bd_plus = Bd_plus_volume(d)
    vol_Ω = vol_Bd_plus^2

    # Sample number of candidates from Poisson
    n_candidates = rand(rng, Poisson(λ_max * vol_Ω))

    # Sample candidates uniformly and thin
    accepted = InteractionSite{d}[]

    for _ in 1:n_candidates
        g = uniform_Bd_plus_sample(d; rng=rng)
        r = uniform_Bd_plus_sample(d; rng=rng)

        # Accept with probability ρ(g,r) / λ_max
        if rand(rng) < ρ(g, r) / λ_max
            push!(accepted, InteractionSite{d}(g, r))
        end
    end

    return accepted
end

"""
    sample_ppp_product(ρ::ProductIntensity{d}; rng=Random.default_rng()) -> Vector{InteractionSite{d}}

Sample from product intensity more efficiently by sampling from each marginal independently.

For ProductIntensity, we can sample more efficiently:
1. Sample N ~ Poisson(c_G · c_R) total sites
2. For each site, sample g from ρ_G (normalized) and r from ρ_R (normalized) independently
"""
function sample_ppp_product(ρ::ProductIntensity{d};
                            rng::AbstractRNG=Random.default_rng()) where d

    # Get total intensities
    c_G = total_intensity(ρ.ρ_G; rng=rng)
    c_R = total_intensity(ρ.ρ_R; rng=rng)

    # Expected number of sites
    E_N = c_G * c_R

    # Sample actual number from Poisson
    n_sites = rand(rng, Poisson(E_N))

    # Sample each site
    sites = Vector{InteractionSite{d}}(undef, n_sites)
    for i in 1:n_sites
        g = sample_from_mixture(ρ.ρ_G; rng=rng)
        r = sample_from_mixture(ρ.ρ_R; rng=rng)
        sites[i] = InteractionSite{d}(g, r)
    end

    return sites
end

"""
    sample_ppp_temporal(ρ, t_start, t_end; dt=0.01, rng=Random.default_rng())

Sample from time-varying intensity over a time interval.
Returns list of (site, time) tuples.

NOTE: This is a simple implementation that discretizes time.
For more accurate temporal point process simulation, consider using JumpProcesses.jl.
"""
function sample_ppp_temporal(ρ, t_start::Real, t_end::Real;
                             dt::Float64=0.01,
                             rng::AbstractRNG=Random.default_rng())
    results = Tuple{InteractionSite, Float64}[]

    t = t_start
    while t < t_end
        # Sample from frozen intensity at time t
        # This is an approximation - proper implementation would use thinning over space-time
        sites = sample_ppp(ρ.ρ; λ_max=nothing, rng=rng)
        for site in sites
            if rand(rng) < dt  # Thin by dt to account for time slice
                push!(results, (site, t))
            end
        end
        t += dt
    end

    return results
end

"""
    sample_from_grid(ρ_G_values::Vector{Float64}, ρ_R_values::Vector{Float64},
                     grid::BdPlusGrid{d}; rng=Random.default_rng()) -> EdgeCentricSample{d}

Sample edge-centric interactions from grid-discretized intensity functions.

This is useful when intensity evolves via PDE and is stored as grid values
rather than as a callable function.

# Arguments
- `ρ_G_values`: Intensity values for source (resource) distribution at each grid point
- `ρ_R_values`: Intensity values for target (consumer) distribution at each grid point
- `grid`: The B^d_+ grid structure
- `rng`: Random number generator

# Returns
EdgeCentricSample containing source and target positions of realized interactions.

# Algorithm
1. Compute c_G = ∫ρ_G dx ≈ Σ ρ_G[i] * h^d (total source intensity)
2. Compute c_R = ∫ρ_R dx ≈ Σ ρ_R[i] * h^d (total target intensity)
3. Sample N ~ Poisson(c_G * c_R) interaction opportunities
4. For each opportunity:
   - Sample source position g from ρ_G (normalized)
   - Sample target position r from ρ_R (normalized)
   - Accept edge with probability g · r
"""
function sample_from_grid(ρ_G_values::Vector{Float64}, ρ_R_values::Vector{Float64},
                          grid::BdPlusGrid{d}; rng::AbstractRNG=Random.default_rng()) where d
    # Compute volume element
    h_d = grid.h^d

    # Compute total intensities (approximate integrals)
    c_G = sum(ρ_G_values) * h_d
    c_R = sum(ρ_R_values) * h_d

    # Expected number of interaction opportunities
    E_N = c_G * c_R

    # Handle edge case of zero intensity
    if E_N < 1e-10
        return EdgeCentricSample{d}(LatentPoint{d}[], LatentPoint{d}[])
    end

    # Sample number of opportunities from Poisson
    N = rand(rng, Poisson(E_N))

    if N == 0
        return EdgeCentricSample{d}(LatentPoint{d}[], LatentPoint{d}[])
    end

    # Normalize to get probability distributions over grid points
    p_G = ρ_G_values ./ sum(ρ_G_values)
    p_R = ρ_R_values ./ sum(ρ_R_values)

    # Sample interactions
    sources = LatentPoint{d}[]
    targets = LatentPoint{d}[]

    for _ in 1:N
        # Sample source position from ρ_G
        g_idx = sample(rng, 1:length(grid.points), Weights(p_G))
        g = grid.points[g_idx]

        # Sample target position from ρ_R
        r_idx = sample(rng, 1:length(grid.points), Weights(p_R))
        r = grid.points[r_idx]

        # Accept with connection probability g · r
        p_connect = connection_probability(g, r)
        if rand(rng) < p_connect
            push!(sources, g)
            push!(targets, r)
        end
    end

    return EdgeCentricSample{d}(sources, targets)
end

"""
    sample_from_grid_full(ρ_G_values, ρ_R_values, grid; rng=Random.default_rng()) -> FullEdgeCentricSample{d}

Sample interactions from grid-discretized intensities with full site information.

Unlike `sample_from_grid`, this samples FULL sites for both source and target:
- Source site: (g_s, r_s) where g_s ~ ρ_G, r_s ~ ρ_R
- Target site: (g_t, r_t) where g_t ~ ρ_G, r_t ~ ρ_R
- Connection probability: g_s · r_t

This preserves all 4 coordinates per edge, allowing clustering based on full
site characteristics rather than just the connection-relevant coordinates.

# Arguments
- `ρ_G_values`: Intensity values for source (resource) distribution
- `ρ_R_values`: Intensity values for target (consumer) distribution
- `grid`: The B^d_+ grid structure
- `rng`: Random number generator

# Returns
FullEdgeCentricSample containing full (g, r) for both source and target sites.
"""
function sample_from_grid_full(ρ_G_values::Vector{Float64}, ρ_R_values::Vector{Float64},
                                grid::BdPlusGrid{d}; rng::AbstractRNG=Random.default_rng()) where d
    # Compute volume element
    h_d = grid.h^d

    # Compute total intensities
    c_G = sum(ρ_G_values) * h_d
    c_R = sum(ρ_R_values) * h_d

    # Expected number of interaction opportunities
    E_N = c_G * c_R

    if E_N < 1e-10
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    # Sample number of opportunities
    N = rand(rng, Poisson(E_N))

    if N == 0
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    # Normalize to probability distributions
    p_G = ρ_G_values ./ sum(ρ_G_values)
    p_R = ρ_R_values ./ sum(ρ_R_values)

    source_sites = InteractionSite{d}[]
    target_sites = InteractionSite{d}[]

    for _ in 1:N
        # Sample FULL source site: g_s from ρ_G, r_s from ρ_R
        g_s_idx = sample(rng, 1:length(grid.points), Weights(p_G))
        r_s_idx = sample(rng, 1:length(grid.points), Weights(p_R))
        g_s = grid.points[g_s_idx]
        r_s = grid.points[r_s_idx]
        source_site = InteractionSite{d}(g_s, r_s)

        # Sample FULL target site: g_t from ρ_G, r_t from ρ_R
        g_t_idx = sample(rng, 1:length(grid.points), Weights(p_G))
        r_t_idx = sample(rng, 1:length(grid.points), Weights(p_R))
        g_t = grid.points[g_t_idx]
        r_t = grid.points[r_t_idx]
        target_site = InteractionSite{d}(g_t, r_t)

        # Connection probability uses g_source · r_target
        p_connect = connection_probability(g_s, r_t)
        if rand(rng) < p_connect
            push!(source_sites, source_site)
            push!(target_sites, target_site)
        end
    end

    return FullEdgeCentricSample{d}(source_sites, target_sites)
end

"""
    initialize_grid_from_mixture(grid::BdPlusGrid{d}, weights, means, κ_vals, scale) -> Vector{Float64}

Initialize intensity values on a grid from a mixture distribution.

Creates intensity values at each grid point by evaluating a Gaussian-like mixture:
    ρ(x) = scale * Σ_k weights[k] * exp(-κ[k] * ||x - means[k]||² / 2)

# Arguments
- `grid`: The B^d_+ grid
- `weights`: Mixture weights (should sum to 1)
- `means`: Vector of mean positions for each component
- `κ_vals`: Concentration parameters (like inverse variance)
- `scale`: Overall scaling factor

# Returns
Vector of intensity values at each grid point.
"""
function initialize_grid_from_mixture(grid::BdPlusGrid{d}, weights, means, κ_vals, scale) where d
    n_components = length(weights)
    ρ_values = zeros(length(grid.points))

    for (idx, point) in enumerate(grid.points)
        val = 0.0
        for k in 1:n_components
            diff = Vector(point) .- means[k]
            dist_sq = sum(diff.^2)
            val += weights[k] * exp(-κ_vals[k] * dist_sq / 2)
        end
        ρ_values[idx] = scale * val
    end

    return ρ_values
end

# ============================================================================
# Sampling from MixtureOfProductIntensities
# ============================================================================

"""
    sample_site_from_mixture(mop::MixtureOfProductIntensities; rng=Random.default_rng()) -> InteractionSite

Sample an InteractionSite from MixtureOfProductIntensities.

Algorithm:
1. Compute γ_m = c_{G,m} · c_{R,m} for each species
2. Sample species m with probability γ_m / C
3. Sample g from ρ_{G,m} and r from ρ_{R,m}
"""
function sample_site_from_mixture(mop::MixtureOfProductIntensities{d};
                                   n_samples::Int=10000,
                                   rng::AbstractRNG=Random.default_rng()) where d
    _, g, r = sample_from_mixture(mop; n_samples=n_samples, rng=rng)
    return InteractionSite{d}(g, r)
end

"""
    sample_ppp_mixture(mop::MixtureOfProductIntensities; rng=Random.default_rng()) -> Vector{Tuple{Int, InteractionSite}}

Sample from a Poisson Point Process with MixtureOfProductIntensities.

Returns a vector of (species_index, site) tuples.

Algorithm:
1. Compute γ_m = c_{G,m} · c_{R,m} and C = Σ_m γ_m
2. Sample N ~ Poisson(C)
3. For each of N sites:
   - Sample species m with probability γ_m / C
   - Sample g from ρ_{G,m}, r from ρ_{R,m}
"""
function sample_ppp_mixture(mop::MixtureOfProductIntensities{d};
                             n_samples::Int=10000,
                             rng::AbstractRNG=Random.default_rng()) where d
    # Compute species intensities γ_m
    γ = species_intensities(mop; n_samples=n_samples, rng=rng)
    C = sum(γ)
    probs = γ ./ C

    # Sample total number of sites
    N = rand(rng, Poisson(C))

    # Sample each site
    sites = Vector{Tuple{Int, InteractionSite{d}}}(undef, N)
    for i in 1:N
        # Sample species
        m = sample(rng, 1:length(mop.species), Weights(probs))
        species = mop.species[m]

        # Sample (g, r) from the chosen species
        g = sample_from_mixture(species.ρ_G; rng=rng)
        r = sample_from_mixture(species.ρ_R; rng=rng)

        sites[i] = (m, InteractionSite{d}(g, r))
    end

    return sites
end

"""
    sample_ppp_mixture_sites_only(mop::MixtureOfProductIntensities; rng=Random.default_rng()) -> Vector{InteractionSite}

Sample sites from MixtureOfProductIntensities PPP, discarding species labels.
Convenience wrapper when species identity is not needed.
"""
function sample_ppp_mixture_sites_only(mop::MixtureOfProductIntensities{d};
                                        n_samples::Int=10000,
                                        rng::AbstractRNG=Random.default_rng()) where d
    labeled = sample_ppp_mixture(mop; n_samples=n_samples, rng=rng)
    return [site for (_, site) in labeled]
end
