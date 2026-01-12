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
    valdim = Val{d}()
    max_val = maximum(
        ρ(uniform_Bd_plus_sample(valdim; rng=rng), uniform_Bd_plus_sample(valdim; rng=rng))
        for _ in 1:n_samples
    )
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
    valdim = Val{d}()
    max_G = 0.0
    max_R = 0.0
    for _ in 1:n_samples
        g = uniform_Bd_plus_sample(valdim; rng=rng)
        r = uniform_Bd_plus_sample(valdim; rng=rng)
        max_G = max(max_G, ρ.ρ_G(g))
        max_R = max(max_R, ρ.ρ_R(r))
    end
    return 1.5 * max_G * max_R
end

"""
    sample_ppp(ρ::AbstractIntensity{d}, λ_max::Float64; rng=Random.default_rng()) -> Vector{InteractionSite{d}}

Sample from inhomogeneous PPP on Ω = B^d_+ × B^d_+ using thinning.

# Algorithm
1. Use provided λ_max = sup ρ(g,r) over Ω
2. Sample N ~ Poisson(λ_max · |Ω|) candidate points uniformly on Ω
3. Accept each point (g,r) with probability ρ(g,r) / λ_max

# Arguments
- `ρ`: Intensity function
- `λ_max`: Upper bound on intensity
- `rng`: Random number generator

# Returns
Vector of accepted interaction sites.

For repeated sampling, precompute λ_max once:
```julia
λ_max = estimate_max_intensity(ρ)
for _ in 1:n_reps
    sites = sample_ppp(ρ, λ_max; rng=rng)
end
```
"""
function sample_ppp(ρ::AbstractIntensity{d}, λ_max::Float64;
                    rng::AbstractRNG=Random.default_rng()) where d
    # Volume of Ω = B^d_+ × B^d_+
    vol_Bd_plus = Bd_plus_volume(d)
    vol_Ω = vol_Bd_plus^2

    # Sample number of candidates from Poisson
    n_candidates = rand(rng, Poisson(λ_max * vol_Ω))

    # Sample candidates uniformly and thin
    valdim = Val{d}()
    accepted = InteractionSite{d}[]

    for _ in 1:n_candidates
        g = uniform_Bd_plus_sample(valdim; rng=rng)
        r = uniform_Bd_plus_sample(valdim; rng=rng)

        # Accept with probability ρ(g,r) / λ_max
        if rand(rng) < ρ(g, r) / λ_max
            push!(accepted, InteractionSite{d}(g, r))
        end
    end

    return accepted
end

"""
    sample_ppp(ρ::AbstractIntensity{d}; rng=Random.default_rng()) -> Vector{InteractionSite{d}}

Convenience method that estimates λ_max automatically.

See `sample_ppp(ρ, λ_max; rng)` for the core method.
"""
function sample_ppp(ρ::AbstractIntensity{d};
                    rng::AbstractRNG=Random.default_rng()) where d
    λ_max = estimate_max_intensity(ρ; rng=rng)
    return sample_ppp(ρ, λ_max; rng=rng)
end

"""
    sample_ppp_product(ρ::ProductIntensity{d}; rng=Random.default_rng()) -> Vector{InteractionSite{d}}

Sample from product intensity more efficiently by sampling from each marginal independently.

For ProductIntensity, we can sample more efficiently:
1. Sample N ~ Poisson(c_G · c_R) total sites
2. For each site, sample g from ρ_G (normalized) and r from ρ_R (normalized) independently

Two methods are provided:
- `sample_ppp_product(ρ; rng)` computes total intensities via quadrature (convenient but slow if called repeatedly)
- `sample_ppp_product(ρ, c_G, c_R; rng)` uses precomputed intensities (fast for repeated sampling)

For performance-critical code with many samples from the same intensity, precompute once:
```julia
c_G = total_intensity(ρ.ρ_G)
c_R = total_intensity(ρ.ρ_R)
for _ in 1:n_reps
    sites = sample_ppp_product(ρ, c_G, c_R; rng=rng)
end
```
"""
function sample_ppp_product(ρ::ProductIntensity{d}, c_G::Float64, c_R::Float64;
                            rng::AbstractRNG=Random.default_rng()) where d
    # Expected number of sites
    E_N = c_G * c_R

    # Sample actual number from Poisson
    n_sites = rand(rng, Poisson(E_N))

    # Sample each site
    sites = [InteractionSite{d}(sample_from_mixture(ρ.ρ_G; rng=rng),
                                 sample_from_mixture(ρ.ρ_R; rng=rng))
             for _ in 1:n_sites]

    return sites
end

# Convenience method that computes total intensities via quadrature
function sample_ppp_product(ρ::ProductIntensity{d};
                            rng::AbstractRNG=Random.default_rng()) where d
    c_G = total_intensity(ρ.ρ_G)
    c_R = total_intensity(ρ.ρ_R)
    return sample_ppp_product(ρ, c_G, c_R; rng=rng)
end

"""
Wrapper struct to freeze a time-varying intensity at a specific time t.
Allows passing to sample_ppp as an AbstractIntensity.
"""
struct FrozenTimeIntensity{d, F} <: AbstractIntensity{d}
    ρ_func::F
    t::Float64
end

function (fti::FrozenTimeIntensity{d})(g::AbstractVector, r::AbstractVector) where d
    return fti.ρ_func(g, r, fti.t)
end

"""
    sample_ppp_temporal(ρ::TimeVaryingIntensity, t_start, t_end; dt=0.01, rng=Random.default_rng())

Sample from time-varying intensity over a time interval.
Returns list of (site, time) tuples.

NOTE: This is a simple implementation that discretizes time.
For more accurate temporal point process simulation, consider using JumpProcesses.jl.
"""
function sample_ppp_temporal(ρ::TimeVaryingIntensity{d}, t_start::Real, t_end::Real;
                             dt::Float64=0.01,
                             rng::AbstractRNG=Random.default_rng()) where d
    results = Tuple{InteractionSite{d}, Float64}[]

    t = t_start
    while t < t_end
        # Create a frozen intensity wrapper for this time slice
        frozen_ρ = FrozenTimeIntensity{d, typeof(ρ.ρ)}(ρ.ρ, t)

        # Sample from the frozen intensity
        sites = sample_ppp(frozen_ρ; rng=rng)
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

# Edge-Centric Semantics
Each edge opportunity "consumes" two node-equivalents (one source, one target).
From E[N] = c_G · c_R node-equivalents, we get E[N]/2 edge opportunities.

Expected edges: E[L] = (E[N]/2) · E[g·r]

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
3. Compute E[N] = c_G · c_R (node-equivalent intensity)
4. Sample L ~ Poisson(E[N]/2) edge opportunities
5. For each opportunity:
   - Sample source position g from ρ_G (normalized)
   - Sample target position r from ρ_R (normalized)
   - Accept edge with probability g · r
"""
function sample_from_grid(ρ_G_values::Vector{Float64}, ρ_R_values::Vector{Float64},
                          grid::BdPlusGrid{d}; rng::AbstractRNG=Random.default_rng()) where d
    # Compute volume element
    h_d = grid.h^d

    # Compute sums once, reuse for both total intensity and normalization
    sum_G = sum(ρ_G_values)
    sum_R = sum(ρ_R_values)

    # Compute total intensities (approximate integrals)
    c_G = sum_G * h_d
    c_R = sum_R * h_d

    # Expected number of node-equivalents
    E_N = c_G * c_R

    # Handle edge case of zero intensity
    if E_N < 1e-10
        return EdgeCentricSample{d}(LatentPoint{d}[], LatentPoint{d}[])
    end

    # EDGE-CENTRIC: Sample E[N]/2 edge opportunities (not E[N])
    # Each opportunity consumes 2 node-equivalents (1 source + 1 target)
    L = rand(rng, Poisson(E_N / 2))

    if L == 0
        return EdgeCentricSample{d}(LatentPoint{d}[], LatentPoint{d}[])
    end

    # Normalize to get probability distributions over grid points
    p_G = ρ_G_values ./ sum_G
    p_R = ρ_R_values ./ sum_R

    # Sample interactions
    sources = LatentPoint{d}[]
    targets = LatentPoint{d}[]

    for _ in 1:L
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

Sample edge-centric interactions from grid-discretized intensities with full site information.

# Edge-Centric Semantics
Each edge opportunity "consumes" two node-equivalents (one source, one target).
From E[N] = c_G · c_R node-equivalents, we get E[N]/2 edge opportunities.

Expected edges: E[L] = (E[N]/2) · E[g·r]

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

    # Compute sums once, reuse for both total intensity and normalization
    sum_G = sum(ρ_G_values)
    sum_R = sum(ρ_R_values)

    # Compute total intensities
    c_G = sum_G * h_d
    c_R = sum_R * h_d

    # Expected number of node-equivalents
    E_N = c_G * c_R

    if E_N < 1e-10
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    # EDGE-CENTRIC: Sample E[N]/2 edge opportunities (not E[N])
    # Each opportunity consumes 2 node-equivalents (1 source + 1 target)
    L = rand(rng, Poisson(E_N / 2))

    if L == 0
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    # Normalize to probability distributions
    p_G = ρ_G_values ./ sum_G
    p_R = ρ_R_values ./ sum_R

    source_sites = InteractionSite{d}[]
    target_sites = InteractionSite{d}[]

    for _ in 1:L
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

    ρ_values = [
        scale * sum(
            weights[k] * exp(-κ_vals[k] * sum((point .- means[k]).^2) / 2)
            for k in 1:n_components
        )
        for point in grid.points
    ]

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

Two methods are provided:
- `sample_ppp_mixture(mop; rng)` computes species intensities via quadrature (convenient but slow if called repeatedly)
- `sample_ppp_mixture(mop, γ; rng)` uses precomputed species intensities (fast for repeated sampling)

For performance-critical code, precompute species intensities once:
```julia
γ = species_intensities(mop)
for _ in 1:n_reps
    sites = sample_ppp_mixture(mop, γ; rng=rng)
end
```
"""
function sample_ppp_mixture(mop::MixtureOfProductIntensities{d}, γ::Vector{Float64};
                             rng::AbstractRNG=Random.default_rng()) where d
    C = sum(γ)
    probs = γ ./ C

    # Sample total number of sites
    N = rand(rng, Poisson(C))

    # Sample each site: choose species, then sample (g, r) from that species
    sites = [
        let m = sample(rng, 1:length(mop.species), Weights(probs))
            species = mop.species[m]
            (m, InteractionSite{d}(sample_from_mixture(species.ρ_G; rng=rng),
                                   sample_from_mixture(species.ρ_R; rng=rng)))
        end
        for _ in 1:N
    ]

    return sites
end

# Convenience method that computes species intensities via quadrature
function sample_ppp_mixture(mop::MixtureOfProductIntensities{d};
                             rng::AbstractRNG=Random.default_rng()) where d
    γ = species_intensities(mop)
    return sample_ppp_mixture(mop, γ; rng=rng)
end

"""
    sample_ppp_mixture_sites_only(mop::MixtureOfProductIntensities; rng=Random.default_rng()) -> Vector{InteractionSite}

Sample sites from MixtureOfProductIntensities PPP, discarding species labels.
Convenience wrapper when species identity is not needed.

Two methods are provided:
- `sample_ppp_mixture_sites_only(mop; rng)` computes species intensities via quadrature
- `sample_ppp_mixture_sites_only(mop, γ; rng)` uses precomputed species intensities
"""
function sample_ppp_mixture_sites_only(mop::MixtureOfProductIntensities{d}, γ::Vector{Float64};
                                        rng::AbstractRNG=Random.default_rng()) where d
    labeled = sample_ppp_mixture(mop, γ; rng=rng)
    return [site for (_, site) in labeled]
end

# Convenience method
function sample_ppp_mixture_sites_only(mop::MixtureOfProductIntensities{d};
                                        rng::AbstractRNG=Random.default_rng()) where d
    labeled = sample_ppp_mixture(mop; rng=rng)
    return [site for (_, site) in labeled]
end

# ============================================================================
# True Edge-Centric Sampling from EdgeIntensity
# ============================================================================

"""
    sample_edge_centric(ei::ScaledProductEdgeIntensity; rng=Random.default_rng()) -> FullEdgeCentricSample

Sample edges from a ScaledProductEdgeIntensity using true edge-centric semantics.

# Edge-Centric Semantics
- C_edge = E[N]/2 edge opportunities (not E[N] or E[N]²)
- Each opportunity "consumes" 2 node-equivalents (1 source + 1 target)
- Source and target are sampled independently from their respective distributions
- Edge is accepted with probability g_source · r_target

# Algorithm
1. Sample L ~ Poisson(C_edge) edge opportunities
2. For each opportunity:
   a. Sample source site s = (g_s, r_s) from ρ_source (normalized)
   b. Sample target site t = (g_t, r_t) from ρ_target (normalized)
   c. Accept edge with probability g_s · r_t
3. Return all accepted edges as FullEdgeCentricSample

# Returns
FullEdgeCentricSample containing full (g, r) for both source and target of each edge.

# Example
```julia
# Create edge intensity from product intensities
ρ_G = BdPlusMixture([1.0], [[0.7, 0.3]], [30.0], 50.0)
ρ_R = BdPlusMixture([1.0], [[0.3, 0.7]], [30.0], 50.0)
ρ_source = ProductIntensity(ρ_G, ρ_R)
ρ_target = ProductIntensity(ρ_G, ρ_R)

ei = ScaledProductEdgeIntensity(ρ_source, ρ_target)
edges = sample_edge_centric(ei)
```
"""
function sample_edge_centric(ei::ScaledProductEdgeIntensity{d};
                              rng::AbstractRNG=Random.default_rng()) where d
    # Sample number of edge opportunities
    L = rand(rng, Poisson(ei.C_edge))

    if L == 0
        return FullEdgeCentricSample{d}(InteractionSite{d}[], InteractionSite{d}[])
    end

    source_sites = InteractionSite{d}[]
    target_sites = InteractionSite{d}[]

    for _ in 1:L
        # Sample source site from ρ_source
        source_site = _sample_site_from_intensity(ei.ρ_source; rng=rng)

        # Sample target site from ρ_target
        target_site = _sample_site_from_intensity(ei.ρ_target; rng=rng)

        # Accept with probability g_source · r_target
        p_connect = connection_probability(source_site.g, target_site.r)
        if rand(rng) < p_connect
            push!(source_sites, source_site)
            push!(target_sites, target_site)
        end
    end

    return FullEdgeCentricSample{d}(source_sites, target_sites)
end

"""
Internal helper to sample an InteractionSite from an AbstractIntensity.
"""
function _sample_site_from_intensity(ρ::ProductIntensity{d};
                                      rng::AbstractRNG=Random.default_rng()) where d
    g = sample_from_mixture(ρ.ρ_G; rng=rng)
    r = sample_from_mixture(ρ.ρ_R; rng=rng)
    return InteractionSite{d}(g, r)
end

function _sample_site_from_intensity(ρ::MixtureOfProductIntensities{d};
                                      rng::AbstractRNG=Random.default_rng()) where d
    _, g, r = sample_from_mixture(ρ; rng=rng)
    return InteractionSite{d}(g, r)
end

"""
    sample_edge_centric_ppp(ρ_source, ρ_target, C_edge; rng=Random.default_rng()) -> FullEdgeCentricSample

Sample edge-centric edges from source/target intensities with precomputed C_edge.

For repeated sampling, precompute C_edge once:
```julia
C_S = total_intensity(ρ_source)
C_T = total_intensity(ρ_target)
C_edge = (C_S + C_T) / 2
for _ in 1:n_reps
    edges = sample_edge_centric_ppp(ρ_source, ρ_target, C_edge; rng=rng)
end
```
"""
function sample_edge_centric_ppp(ρ_source::AbstractIntensity{d}, ρ_target::AbstractIntensity{d},
                                  C_edge::Float64;
                                  rng::AbstractRNG=Random.default_rng()) where d
    ei = ScaledProductEdgeIntensity{d}(ρ_source, ρ_target, C_edge)
    return sample_edge_centric(ei; rng=rng)
end

"""
    sample_edge_centric_ppp(ρ_source, ρ_target; rng=Random.default_rng()) -> FullEdgeCentricSample

Convenience method that computes C_edge automatically.

See `sample_edge_centric_ppp(ρ_source, ρ_target, C_edge; rng)` for the core method.
"""
function sample_edge_centric_ppp(ρ_source::AbstractIntensity{d}, ρ_target::AbstractIntensity{d};
                                  rng::AbstractRNG=Random.default_rng()) where d
    ei = ScaledProductEdgeIntensity(ρ_source, ρ_target)
    return sample_edge_centric(ei; rng=rng)
end

"""
    sample_edge_centric_symmetric(ρ, C_edge; rng=Random.default_rng()) -> FullEdgeCentricSample

Sample edge-centric edges with symmetric source/target distribution and precomputed C_edge.

For repeated sampling, precompute C_edge once:
```julia
C = total_intensity(ρ)
C_edge = C / 2
for _ in 1:n_reps
    edges = sample_edge_centric_symmetric(ρ, C_edge; rng=rng)
end
```
"""
function sample_edge_centric_symmetric(ρ::AbstractIntensity{d}, C_edge::Float64;
                                        rng::AbstractRNG=Random.default_rng()) where d
    ei = ScaledProductEdgeIntensity{d}(ρ, ρ, C_edge)
    return sample_edge_centric(ei; rng=rng)
end

"""
    sample_edge_centric_symmetric(ρ; rng=Random.default_rng()) -> FullEdgeCentricSample

Convenience method that computes C_edge automatically.

See `sample_edge_centric_symmetric(ρ, C_edge; rng)` for the core method.
"""
function sample_edge_centric_symmetric(ρ::AbstractIntensity{d};
                                        rng::AbstractRNG=Random.default_rng()) where d
    ei = SymmetricEdgeIntensity(ρ)
    return sample_edge_centric(ei; rng=rng)
end
