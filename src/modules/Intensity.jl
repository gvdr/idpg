# Intensity function representations on Ω = G × R ⊆ B^d_+ × B^d_+
# where B^d_+ is the non-negative unit ball

"""
Abstract type for intensity functions on Ω = G × R ⊆ B^d_+ × B^d_+.

Intensity functions are NOT probability densities - they need not integrate to 1.
The total intensity gives the expected number of sampled points.
"""
abstract type AbstractIntensity{d} end

"""
    BdPlusMixture{d}

A mixture of Gaussian-like kernels for intensity on B^d_+ (non-negative unit ball).

Each component is a truncated Gaussian centered at a mean position in B^d_+.
The intensity at x is: scale * Σᵢ wᵢ * exp(-κᵢ * ||x - μᵢ||²)

This is NOT a probability distribution - `scale` controls the total intensity
(expected number of points), not normalization.

# Fields
- `weights::Vector{Float64}`: Mixture weights (must sum to 1)
- `means::Vector{SVector{d, Float64}}`: Mean positions in B^d_+
- `concentrations::Vector{Float64}`: Concentration parameters (higher = more peaked)
- `scale::Float64`: Total intensity scaling
"""
struct BdPlusMixture{d}
    weights::Vector{Float64}
    means::Vector{SVector{d, Float64}}
    concentrations::Vector{Float64}
    scale::Float64

    function BdPlusMixture{d}(weights, means, concentrations, scale) where d
        @assert length(weights) == length(means) == length(concentrations)
        @assert abs(sum(weights) - 1.0) < 1e-10 "Mixture weights must sum to 1"
        @assert all(length(μ) == d for μ in means) "All means must have dimension d"
        @assert all(in_Bd_plus(μ) for μ in means) "All means must be in B^d_+"
        @assert all(κ > 0 for κ in concentrations) "All concentrations must be positive"
        @assert scale > 0 "Scale must be positive"
        new{d}(weights, means, concentrations, scale)
    end
end

"""
    BdPlusMixture(weights, means, concentrations, scale)

Construct a mixture of Gaussian kernels on B^d_+.

# Arguments
- `weights`: Mixture weights (must sum to 1)
- `means`: Vector of mean positions in B^d_+ (each a d-vector)
- `concentrations`: Concentration parameters (higher = more peaked around mean)
- `scale`: Total intensity scaling (controls expected number of points)

# Example
```julia
# Single component centered at (0.7, 0.3) with concentration 10
ρ = BdPlusMixture([1.0], [[0.7, 0.3]], [10.0], 50.0)
```
"""
function BdPlusMixture(weights::Vector{Float64},
                       means::Vector{<:AbstractVector},
                       concentrations::Vector{Float64},
                       scale::Float64)
    d = length(first(means))
    sv_means = [SVector{d, Float64}(μ) for μ in means]
    return BdPlusMixture{d}(weights, sv_means, concentrations, scale)
end

"""
Evaluate the BdPlusMixture intensity at point x in B^d_+.
Uses Gaussian kernel: intensity = scale * Σᵢ wᵢ * exp(-κᵢ * ||x - μᵢ||²)
"""
function (bm::BdPlusMixture{d})(x::AbstractVector) where d
    # Check if x is in B^d_+
    if !in_Bd_plus(x; tol=1e-6)
        return 0.0
    end

    # Compute mixture of Gaussian kernels
    p = sum(
        w * exp(-κ * sum((x[i] - μ[i])^2 for i in 1:d))
        for (w, μ, κ) in zip(bm.weights, bm.means, bm.concentrations)
    )

    return bm.scale * p
end

"""
    ProductIntensity{d, FG, FR}

Product intensity: ρ(g,r) = ρ_G(g) · ρ_R(r)

Under product intensity, a node's propensity to propose connections (position in G)
is independent of its propensity to accept connections (position in R).

# Fields
- `ρ_G::FG`: Intensity function on G (green/source space)
- `ρ_R::FR`: Intensity function on R (red/target space)
"""
struct ProductIntensity{d, FG, FR} <: AbstractIntensity{d}
    ρ_G::FG
    ρ_R::FR
end

function ProductIntensity(ρ_G::BdPlusMixture{d}, ρ_R::BdPlusMixture{d}) where d
    return ProductIntensity{d, BdPlusMixture{d}, BdPlusMixture{d}}(ρ_G, ρ_R)
end

"""
Evaluate product intensity at a point (g, r).
"""
function (ρ::ProductIntensity{d})(g::AbstractVector, r::AbstractVector) where d
    return ρ.ρ_G(g) * ρ.ρ_R(r)
end

"""
Evaluate product intensity for a tuple (g, r).
"""
function (ρ::ProductIntensity{d})(site::Tuple{<:AbstractVector, <:AbstractVector}) where d
    return ρ(site[1], site[2])
end

"""
    TimeVaryingIntensity{d, F}

Time-varying intensity for PDE evolution: ρ(g, r, t).

# Fields
- `ρ::F`: Function (g, r, t) -> intensity value
"""
struct TimeVaryingIntensity{d, F} <: AbstractIntensity{d}
    ρ::F
end

function (tvi::TimeVaryingIntensity{d})(g::AbstractVector, r::AbstractVector, t::Real) where d
    return tvi.ρ(g, r, t)
end

"""
    _select_quadrature_solver(d::Int)

Select the best quadrature solver based on dimension.
Based on benchmarks (January 2026):
- d ≤ 4: HCubatureJL (fast adaptive cubature)
- d = 5-6: CubaDivonne (good balance of speed and accuracy)
- d ≥ 7: CubaVegas (most robust for high dimensions)
"""
function _select_quadrature_solver(d::Int)
    if d <= 4
        return HCubatureJL()
    elseif d <= 6
        return CubaDivonne()
    else
        return CubaVegas()
    end
end

"""
    total_intensity(ρ::BdPlusMixture; reltol=1e-4) -> Float64

Compute the total intensity c = ∫ρ(x)dx over B^d_+.

# Arguments
- `ρ`: Intensity function (BdPlusMixture)
- `reltol`: Relative tolerance for quadrature (default: 1e-4)

# Implementation Note
Uses hyperspherical coordinates to avoid the discontinuity at the B^d_+ boundary.
Solver is automatically selected based on dimension:
- d ≤ 4: HCubatureJL (adaptive cubature)
- d = 5-6: CubaDivonne (Cuba library)
- d ≥ 7: CubaVegas (Vegas Monte Carlo from Cuba library)
"""
function total_intensity(ρ::BdPlusMixture{d}; reltol::Float64=1e-4) where d
    # 1D case: direct integration with scalar bounds
    if d == 1
        f1d(x, p) = ρ([x])
        prob = IntegralProblem(f1d, (0.0, 1.0))
        sol = solve(prob, HCubatureJL(); reltol=reltol)
        return sol.u
    end

    # d ≥ 2: use hyperspherical coordinates to avoid boundary discontinuity
    function integrand_c(u, p)
        r = u[1]
        r < 1e-12 && return 0.0  # avoid singularity at origin
        θ = @view u[2:end]
        x = hyperspherical_to_cartesian(r, θ)
        J = hyperspherical_jacobian(r, θ)
        return ρ(x) * J
    end

    # Bounds: r ∈ [0,1], θᵢ ∈ [0,π/2]
    lower = zeros(d)
    upper = vcat(1.0, fill(π/2, d-1))

    prob = IntegralProblem(integrand_c, (lower, upper))
    solver = _select_quadrature_solver(d)

    # HCubature uses reltol, Cuba solvers use defaults
    if d <= 4
        sol = solve(prob, solver; reltol=reltol)
    else
        sol = solve(prob, solver)
    end

    return sol.u
end

"""
    total_intensity(ρ::ProductIntensity; reltol=1e-4) -> Float64

Compute the total intensity C = c_G × c_R = E[N] for a product intensity.
"""
function total_intensity(ρ::ProductIntensity{d}; reltol::Float64=1e-4) where d
    c_G, c_R = marginal_total_intensity(ρ; reltol=reltol)
    return c_G * c_R
end

"""
    marginal_total_intensity(ρ::ProductIntensity; reltol=1e-4) -> Tuple{Float64, Float64}

Compute the marginal total intensities c_G = ∫ρ_G(g)dg and c_R = ∫ρ_R(r)dr.
Returns (c_G, c_R).
"""
function marginal_total_intensity(ρ::ProductIntensity{d}; reltol::Float64=1e-4) where d
    c_G = total_intensity(ρ.ρ_G; reltol=reltol)
    c_R = total_intensity(ρ.ρ_R; reltol=reltol)
    return (c_G, c_R)
end

"""
    intensity_weighted_mean(ρ::BdPlusMixture; reltol=1e-4) -> SVector

Compute the intensity-weighted mean position μ = ∫x·ρ(x)dx over B^d_+.

# Arguments
- `ρ`: Intensity function (BdPlusMixture)
- `reltol`: Relative tolerance for quadrature (used for d ≤ 4)

# Implementation Note
Uses hyperspherical coordinates to avoid the discontinuity at the B^d_+ boundary.
Solver is automatically selected based on dimension:
- d ≤ 4: HCubatureJL (adaptive cubature)
- d = 5-6: CubaDivonne (Cuba library)
- d ≥ 7: CubaVegas (Vegas Monte Carlo from Cuba library)
"""
function intensity_weighted_mean(ρ::BdPlusMixture{d}; reltol::Float64=1e-4) where d
    # 1D case: direct integration with scalar bounds
    if d == 1
        f1d(x, p) = x * ρ([x])
        prob = IntegralProblem(f1d, (0.0, 1.0))
        sol = solve(prob, HCubatureJL(); reltol=reltol)
        return SVector{1, Float64}(sol.u)
    end

    # d ≥ 2: integrate each component using hyperspherical coords
    μ = zeros(d)
    solver = _select_quadrature_solver(d)

    for k in 1:d
        function integrand_μk(u, p)
            r = u[1]
            r < 1e-12 && return 0.0  # avoid singularity at origin
            θ = @view u[2:end]
            x = hyperspherical_to_cartesian(r, θ)
            J = hyperspherical_jacobian(r, θ)
            return x[k] * ρ(x) * J
        end

        # Bounds: r ∈ [0,1], θᵢ ∈ [0,π/2]
        lower = zeros(d)
        upper = vcat(1.0, fill(π/2, d-1))
        prob = IntegralProblem(integrand_μk, (lower, upper))

        # HCubature uses reltol, Cuba solvers use defaults
        if d <= 4
            sol = solve(prob, solver; reltol=reltol)
        else
            sol = solve(prob, solver)
        end
        μ[k] = sol.u
    end

    return SVector{d, Float64}(μ)
end

"""
    normalized_mean(ρ::BdPlusMixture; reltol=1e-4) -> SVector

Compute the normalized mean μ̃ = μ/c where μ is the intensity-weighted mean
and c is the total intensity.
"""
function normalized_mean(ρ::BdPlusMixture{d}; reltol::Float64=1e-4) where d
    c = total_intensity(ρ; reltol=reltol)
    μ = intensity_weighted_mean(ρ; reltol=reltol)
    return μ ./ c
end

"""
    marginal_stats(ρ::ProductIntensity; reltol=1e-4)

Compute all marginal statistics for a product intensity:
- c_G, c_R: marginal total intensities
- μ_G, μ_R: intensity-weighted mean positions
- μ̃_G, μ̃_R: normalized means

# Edge-Centric vs Node-Centric
- Node-centric: N sites → N² pairs → E[|E|] = E[N]² · E[g·r]
- Edge-centric: E[N]/2 opportunities → E[L] = (E[N]/2) · E[g·r]

Each edge opportunity "consumes" two node-equivalents (one source, one target).

Returns a named tuple with all quantities.
"""
function marginal_stats(ρ::ProductIntensity{d}; reltol::Float64=1e-4) where d
    c_G = total_intensity(ρ.ρ_G; reltol=reltol)
    c_R = total_intensity(ρ.ρ_R; reltol=reltol)

    μ_G = intensity_weighted_mean(ρ.ρ_G; reltol=reltol)
    μ_R = intensity_weighted_mean(ρ.ρ_R; reltol=reltol)

    μ̃_G = μ_G ./ c_G
    μ̃_R = μ_R ./ c_R

    E_N = c_G * c_R
    avg_conn_prob = dot(μ̃_G, μ̃_R)

    return (
        c_G = c_G,
        c_R = c_R,
        μ_G = μ_G,
        μ_R = μ_R,
        μ̃_G = μ̃_G,
        μ̃_R = μ̃_R,
        E_N = E_N,
        avg_conn_prob = avg_conn_prob,
        E_edges_node_centric = E_N^2 * avg_conn_prob,
        # Edge-centric: E[N]/2 opportunities, each accepts with prob E[g·r]
        E_edges_edge_centric = (E_N / 2) * avg_conn_prob
    )
end

"""
    MixtureOfProductIntensities{d}

Mixture of Product intensities: ρ(g,r) = Σ_m ρ_{G,m}(g) · ρ_{R,m}(r)

This models M species/groups, each with its own (G, R) niche distribution.
Unlike ProductIntensity where G and R are independent, here each species
has a COUPLED (G, R) niche - a site's G and R coordinates come from the
same species.

Key property: No cross-species mixing. Species m contributes ρ_{G,m}(g)·ρ_{R,m}(r),
never ρ_{G,m}(g)·ρ_{R,n}(r) for n ≠ m.

# Fields
- `species::Vector{ProductIntensity{d}}`: Vector of species-specific (ρ_G, ρ_R) pairs
- `_cache`: Cached species intensities (computed lazily)

# Key quantities (computed via marginal_stats):
- γ_m = c_{G,m} · c_{R,m}: Species-specific intensity (abundance emerges from niche)
- C = Σ_m γ_m: Total intensity
- P(species=m) = γ_m / C: Probability of sampling from species m
"""
struct MixtureOfProductIntensities{d} <: AbstractIntensity{d}
    species::Vector{ProductIntensity{d, BdPlusMixture{d}, BdPlusMixture{d}}}

    function MixtureOfProductIntensities{d}(species) where d
        @assert length(species) > 0 "Must have at least one species"
        new{d}(species)
    end
end

"""
    MixtureOfProductIntensities(species::Vector{ProductIntensity{d}})

Construct a MixtureOfProductIntensities from a vector of species ProductIntensities.
"""
function MixtureOfProductIntensities(species::Vector{ProductIntensity{d, BdPlusMixture{d}, BdPlusMixture{d}}}) where d
    return MixtureOfProductIntensities{d}(species)
end

"""
    MixtureOfProductIntensities(ρ_Gs::Vector{BdPlusMixture{d}}, ρ_Rs::Vector{BdPlusMixture{d}})

Construct from paired vectors of G and R marginals for each species.
"""
function MixtureOfProductIntensities(ρ_Gs::Vector{BdPlusMixture{d}}, ρ_Rs::Vector{BdPlusMixture{d}}) where d
    @assert length(ρ_Gs) == length(ρ_Rs) "Must have same number of G and R distributions"
    species = [ProductIntensity(ρ_Gs[m], ρ_Rs[m]) for m in eachindex(ρ_Gs)]
    return MixtureOfProductIntensities{d}(species)
end

"""
Evaluate MixtureOfProductIntensities at (g, r).
ρ(g,r) = Σ_m ρ_{G,m}(g) · ρ_{R,m}(r)
"""
function (mop::MixtureOfProductIntensities{d})(g::AbstractVector, r::AbstractVector) where d
    return sum(species.ρ_G(g) * species.ρ_R(r) for species in mop.species)
end

"""
Evaluate MixtureOfProductIntensities for a tuple (g, r).
"""
function (mop::MixtureOfProductIntensities{d})(site::Tuple{<:AbstractVector, <:AbstractVector}) where d
    return mop(site[1], site[2])
end

"""
    n_species(mop::MixtureOfProductIntensities) -> Int

Return the number of species in the mixture.
"""
n_species(mop::MixtureOfProductIntensities) = length(mop.species)

"""
    species_intensities(mop::MixtureOfProductIntensities; reltol=1e-4)

Compute species-specific total intensities γ_m = c_{G,m} · c_{R,m}.
Returns a vector of γ values.
"""
function species_intensities(mop::MixtureOfProductIntensities{d}; reltol::Float64=1e-4) where d
    return [
        total_intensity(species.ρ_G; reltol=reltol) *
        total_intensity(species.ρ_R; reltol=reltol)
        for species in mop.species
    ]
end

"""
    total_intensity(mop::MixtureOfProductIntensities; reltol=1e-4)

Compute total intensity C = Σ_m γ_m = Σ_m c_{G,m} · c_{R,m}.
"""
function total_intensity(mop::MixtureOfProductIntensities{d}; reltol::Float64=1e-4) where d
    return sum(species_intensities(mop; reltol=reltol))
end

"""
    species_probabilities(mop::MixtureOfProductIntensities; reltol=1e-4)

Compute sampling probabilities for each species: P(species=m) = γ_m / C.
Returns a vector of probabilities that sum to 1.
"""
function species_probabilities(mop::MixtureOfProductIntensities{d}; reltol::Float64=1e-4) where d
    γ = species_intensities(mop; reltol=reltol)
    C = sum(γ)
    return γ ./ C
end

"""
    marginal_stats(mop::MixtureOfProductIntensities; reltol=1e-4)

Compute statistics for MixtureOfProductIntensities.

# Edge-Centric vs Node-Centric
- Node-centric: N sites → N² pairs → complex cross-species interactions
- Edge-centric: E[N]/2 opportunities → E[L] = (E[N]/2) · E[g·r]

Each edge opportunity "consumes" two node-equivalents (one source, one target).

Returns a named tuple with:
- γ: Vector of species-specific intensities γ_m = c_{G,m} · c_{R,m}
- C: Total intensity (sum of γ) = E[N]
- species_probs: Sampling probabilities P(m) = γ_m / C
- per_species: Vector of per-species stats (c_G, c_R, μ̃_G, μ̃_R, avg_conn_prob)
- E_N: Expected total sites = C
- E_edges_edge_centric: Expected edges = (E[N]/2) · E[g·r]
"""
function marginal_stats(mop::MixtureOfProductIntensities{d}; reltol::Float64=1e-4) where d
    M = length(mop.species)

    # Compute per-species statistics
    per_species = Vector{NamedTuple}(undef, M)
    γ = Vector{Float64}(undef, M)

    for (m, species) in enumerate(mop.species)
        c_G = total_intensity(species.ρ_G; reltol=reltol)
        c_R = total_intensity(species.ρ_R; reltol=reltol)
        μ_G = intensity_weighted_mean(species.ρ_G; reltol=reltol)
        μ_R = intensity_weighted_mean(species.ρ_R; reltol=reltol)
        μ̃_G = μ_G ./ c_G
        μ̃_R = μ_R ./ c_R
        avg_conn_prob = dot(μ̃_G, μ̃_R)

        γ[m] = c_G * c_R
        per_species[m] = (
            c_G = c_G,
            c_R = c_R,
            μ̃_G = μ̃_G,
            μ̃_R = μ̃_R,
            avg_conn_prob = avg_conn_prob
        )
    end

    C = sum(γ)
    species_probs = γ ./ C

    # Weighted average connection probability
    avg_conn_prob = sum(species_probs[m] * per_species[m].avg_conn_prob for m in 1:M)

    # Edge-centric: E[N]/2 opportunities, each accepts with prob E[g·r]
    E_edges = (C / 2) * avg_conn_prob

    return (
        γ = γ,
        C = C,
        species_probs = species_probs,
        per_species = per_species,
        E_N = C,
        avg_conn_prob = avg_conn_prob,
        E_edges_edge_centric = E_edges
    )
end

"""
    sample_from_mixture(mop::MixtureOfProductIntensities; rng=Random.default_rng()) -> (species_idx, g, r)

Sample a site from MixtureOfProductIntensities.

Algorithm:
1. Compute γ_m = c_{G,m} · c_{R,m} for each species (cached if possible)
2. Sample species m with probability γ_m / Σγ
3. Sample g from ρ_{G,m} and r from ρ_{R,m}

Returns a tuple (species_index, g, r) where g and r are LatentPoints.
"""
function sample_from_mixture(mop::MixtureOfProductIntensities{d};
                              reltol::Float64=1e-4,
                              rng::AbstractRNG=Random.default_rng()) where d
    # Compute species probabilities
    γ = species_intensities(mop; reltol=reltol)
    probs = γ ./ sum(γ)

    # Sample species
    m = sample(rng, 1:length(mop.species), Weights(probs))
    species = mop.species[m]

    # Sample (g, r) from the chosen species
    g = sample_from_mixture(species.ρ_G; rng=rng)
    r = sample_from_mixture(species.ρ_R; rng=rng)

    return (m, g, r)
end

"""
    sample_from_mixture(bm::BdPlusMixture; rng=Random.default_rng()) -> LatentPoint

Sample a single point from the BdPlusMixture using rejection sampling.
Samples from a Gaussian centered at a randomly chosen component mean,
then rejects if outside B^d_+.
"""
function sample_from_mixture(bm::BdPlusMixture{d};
                             rng::AbstractRNG=Random.default_rng()) where d
    # Choose component according to weights
    k = sample(rng, 1:length(bm.weights), Weights(bm.weights))

    μ = bm.means[k]
    κ = bm.concentrations[k]
    σ = 1.0 / sqrt(2 * κ)  # Standard deviation from concentration

    # Rejection sampling: sample from Gaussian, reject if outside B^d_+
    max_attempts = 1000
    for _ in 1:max_attempts
        # Sample from Gaussian centered at μ
        x = μ .+ σ .* randn(rng, d)

        # Accept if in B^d_+
        if in_Bd_plus(x)
            return SVector{d, Float64}(x)
        end
    end

    # Fallback: project to B^d_+ (shouldn't happen often)
    x = μ .+ σ .* randn(rng, d)
    return project_to_Bd_plus(SVector{d, Float64}(x))
end

# ============================================================================
# Edge-Centric Intensity: Joint distribution over source-target pairs
# ============================================================================

"""
    AbstractEdgeIntensity{d}

Abstract type for edge intensities on edge space E = (B^d_+)^4.

An edge intensity ρ_edge(s, t) defines a joint distribution over source-target pairs,
where s = (g_s, r_s) and t = (g_t, r_t) are full interaction sites.

Key property: C_edge = E[N]/2, where E[N] is the expected number of node-equivalents.
Each edge opportunity "consumes" 2 node-equivalents (1 source + 1 target).
"""
abstract type AbstractEdgeIntensity{d} end

"""
    ScaledProductEdgeIntensity{d}

Edge intensity where source and target are conditionally independent given an opportunity.

    ρ_edge(s, t) = (C_edge / (C_S · C_T)) · ρ_S(s) · ρ_T(t)

The scaling factor ensures the total intensity equals C_edge = E[N]/2, not C_S · C_T.

# Fields
- `ρ_source`: Intensity for source sites (can be ProductIntensity or MixtureOfProductIntensities)
- `ρ_target`: Intensity for target sites (can be ProductIntensity or MixtureOfProductIntensities)
- `C_edge`: Total edge opportunity intensity = E[N]/2

# Sampling
1. Sample L ~ Poisson(C_edge)
2. For each of L opportunities:
   - Sample source site s ~ ρ_S / C_S (normalized)
   - Sample target site t ~ ρ_T / C_T (normalized)
   - Accept edge with probability g_s · r_t
"""
struct ScaledProductEdgeIntensity{d, FS<:AbstractIntensity{d}, FT<:AbstractIntensity{d}} <: AbstractEdgeIntensity{d}
    ρ_source::FS
    ρ_target::FT
    C_edge::Float64

    function ScaledProductEdgeIntensity{d}(ρ_source::FS, ρ_target::FT, C_edge::Float64) where {d, FS<:AbstractIntensity{d}, FT<:AbstractIntensity{d}}
        @assert C_edge > 0 "C_edge must be positive"
        new{d, FS, FT}(ρ_source, ρ_target, C_edge)
    end
end

"""
    ScaledProductEdgeIntensity(ρ_source, ρ_target; reltol=1e-4)

Construct a ScaledProductEdgeIntensity from source and target intensities.

Computes C_edge = (C_source + C_target)/2, interpreting the total node-equivalents
as coming from two separate populations.

For symmetric case (same population), use SymmetricEdgeIntensity instead.

For repeated construction with known C_edge, use the direct constructor:
`ScaledProductEdgeIntensity{d}(ρ_source, ρ_target, C_edge)`
"""
function ScaledProductEdgeIntensity(ρ_source::AbstractIntensity{d}, ρ_target::AbstractIntensity{d};
                                     reltol::Float64=1e-4) where d
    # For separate source/target populations:
    # Total node-equivalents = C_source + C_target
    # Each edge consumes 2 (one from each), so C_edge = (C_source + C_target)/2
    C_S = total_intensity(ρ_source; reltol=reltol)
    C_T = total_intensity(ρ_target; reltol=reltol)
    C_edge = (C_S + C_T) / 2
    return ScaledProductEdgeIntensity{d}(ρ_source, ρ_target, C_edge)
end

"""
    SymmetricEdgeIntensity(ρ; reltol=1e-4)

Convenience constructor for symmetric case where ρ_source = ρ_target = ρ.

For a single population with intensity ρ:
- E[N] = total_intensity(ρ)
- C_edge = E[N]/2 (each edge consumes 2 node-equivalents from the same pool)
"""
function SymmetricEdgeIntensity(ρ::AbstractIntensity{d}; reltol::Float64=1e-4) where d
    # Single population: E[N] = C, so C_edge = C/2
    C = total_intensity(ρ; reltol=reltol)
    C_edge = C / 2
    return ScaledProductEdgeIntensity{d}(ρ, ρ, C_edge)
end

"""
    edge_intensity(ei::ScaledProductEdgeIntensity) -> Float64

Return the total edge opportunity intensity C_edge = E[N]/2.
"""
edge_intensity(ei::ScaledProductEdgeIntensity) = ei.C_edge

"""
    marginal_stats(ei::ScaledProductEdgeIntensity; reltol=1e-4)

Compute statistics for a ScaledProductEdgeIntensity.

Returns a named tuple with:
- C_edge: Edge opportunity intensity = E[N]/2
- C_source, C_target: Source and target total intensities
- E_N: Expected node-equivalents = C_source · C_target
- avg_conn_prob: Expected connection probability E[g_s · r_t]
- E_edges: Expected realized edges = C_edge · avg_conn_prob
"""
function marginal_stats(ei::ScaledProductEdgeIntensity{d}; reltol::Float64=1e-4) where d
    # Get stats from source and target
    stats_S = marginal_stats(ei.ρ_source; reltol=reltol)
    stats_T = marginal_stats(ei.ρ_target; reltol=reltol)

    # Compute cross-distribution connection probability using dispatch
    avg_conn_prob = cross_avg_conn_prob(stats_S, stats_T, ei.ρ_source, ei.ρ_target)

    # Extract total intensities (E_N present in both ProductIntensity and MixtureOfProductIntensities stats)
    C_source = stats_S.E_N
    C_target = stats_T.E_N

    return (
        C_edge = ei.C_edge,
        C_source = C_source,
        C_target = C_target,
        E_N = C_source * C_target,
        avg_conn_prob = avg_conn_prob,
        E_edges = ei.C_edge * avg_conn_prob
    )
end

"""
Compute expected connection probability E[g_source · r_target] for cross-distribution pairs.
Uses multiple dispatch to handle different intensity type combinations.
"""
function cross_avg_conn_prob(stats_S, stats_T, ::ProductIntensity, ::ProductIntensity)
    return dot(stats_S.μ̃_G, stats_T.μ̃_R)
end

function cross_avg_conn_prob(stats_S, stats_T, ::MixtureOfProductIntensities, ::MixtureOfProductIntensities)
    # Approximate: use geometric mean of within-population connection probs
    return sqrt(stats_S.avg_conn_prob * stats_T.avg_conn_prob)
end

function cross_avg_conn_prob(stats_S, stats_T, ::ProductIntensity, ::MixtureOfProductIntensities)
    return sqrt(stats_S.avg_conn_prob * stats_T.avg_conn_prob)
end

function cross_avg_conn_prob(stats_S, stats_T, ::MixtureOfProductIntensities, ::ProductIntensity)
    return sqrt(stats_S.avg_conn_prob * stats_T.avg_conn_prob)
end

# ============================================================================
# Utility Functions for 1D Intensities
# ============================================================================

"""
    truncated_gaussian_sigma(x, μ, σ)

Compute the normalized truncated Gaussian density on [0,1].

# Arguments
- `x`: Point to evaluate (returns 0 if outside [0,1])
- `μ`: Mean of the underlying Gaussian
- `σ`: Standard deviation of the underlying Gaussian

# Returns
The density value pdf(x)/Z where Z = CDF(1) - CDF(0) ensures normalization on [0,1].
"""
function truncated_gaussian_sigma(x::Real, μ::Real, σ::Real)
    if x < 0.0 || x > 1.0
        return 0.0
    end
    d = Normal(μ, σ)
    Z = cdf(d, 1.0) - cdf(d, 0.0)
    return pdf(d, x) / Z
end

"""
    truncated_gaussian_kappa(x, μ, κ)

Compute the normalized truncated Gaussian density on [0,1] using concentration parameter.

# Arguments
- `x`: Point to evaluate (returns 0 if outside [0,1])
- `μ`: Mean of the underlying Gaussian
- `κ`: Concentration parameter (higher = more peaked); σ = 1/√κ

# Returns
The density value pdf(x)/Z where Z = CDF(1) - CDF(0) ensures normalization on [0,1].
"""
function truncated_gaussian_kappa(x::Real, μ::Real, κ::Real)
    σ = 1.0 / sqrt(κ)
    return truncated_gaussian_sigma(x, μ, σ)
end

"""
    compute_1d_marginals(grid, μ_G, μ_R, σ)

Compute normalized marginal intensities ρ_G(g) and ρ_R(r) on a 1D grid.

# Arguments
- `grid`: Vector of grid points in [0,1]
- `μ_G`: Mean for the green (source) intensity
- `μ_R`: Mean for the red (target) intensity
- `σ`: Standard deviation for both

# Returns
Tuple (ρ_G, ρ_R) of vectors with truncated Gaussian values at grid points.
"""
function compute_1d_marginals(grid::AbstractVector, μ_G::Real, μ_R::Real, σ::Real)
    ρ_G = [truncated_gaussian_sigma(g, μ_G, σ) for g in grid]
    ρ_R = [truncated_gaussian_sigma(r, μ_R, σ) for r in grid]
    return ρ_G, ρ_R
end

"""
    compute_bound_heat_matrix(grid, ρ_G, ρ_R; normalize=false)

Compute bound heat matrix h̄(g,r) = K(g,r) · ρ_G(g) · ρ_R(r) on a 1D grid.

The kernel K(g,r) = g · r is the standard IDPG connection probability kernel.

# Arguments
- `grid`: Vector of grid points in [0,1]
- `ρ_G`: Green (source) marginal intensity values at grid points
- `ρ_R`: Red (target) marginal intensity values at grid points
- `normalize`: If true, normalize ρ_G and ρ_R to integrate to 1

# Returns
Matrix of bound heat values h̄[i,j] = grid[i] * grid[j] * ρ_G[i] * ρ_R[j].
"""
function compute_bound_heat_matrix(grid::AbstractVector, ρ_G::AbstractVector, ρ_R::AbstractVector; normalize::Bool=false)
    dx = length(grid) > 1 ? grid[2] - grid[1] : 1.0

    if normalize
        c_G = sum(ρ_G) * dx
        c_R = sum(ρ_R) * dx
        ρ_G_use = ρ_G ./ max(c_G, 1e-10)
        ρ_R_use = ρ_R ./ max(c_R, 1e-10)
    else
        ρ_G_use = ρ_G
        ρ_R_use = ρ_R
    end

    # Outer product: h̄[i,j] = grid[i] * ρ_G[i] * grid[j] * ρ_R[j]
    return (grid .* ρ_G_use) * (grid .* ρ_R_use)'
end

# ============================================================================
# Calibrated Product Intensity Construction
# ============================================================================

"""
    create_calibrated_product_intensity(target_lambda; mean_G, mean_R, conc, rng)

Create a ProductIntensity with BdPlusMixture marginals calibrated to achieve a target Λ = c_G × c_R.

This is useful for simulations where you want to control the expected number of sites
while specifying the shape of the marginal distributions.

# Arguments
- `target_lambda`: Target total intensity (c_G × c_R), which equals E[N] for node-centric
- `mean_G`: Mean position for G (source) intensity in B^d_+ (default: [0.6, 0.4])
- `mean_R`: Mean position for R (target) intensity in B^d_+ (default: [0.5, 0.5])
- `conc`: Concentration parameter κ for both marginals (default: 15.0)
- `rng`: Random number generator for calibration Monte Carlo

# Returns
Tuple (ρ_G, ρ_R, ρ) where:
- `ρ_G`: BdPlusMixture for source intensity
- `ρ_R`: BdPlusMixture for target intensity
- `ρ`: ProductIntensity(ρ_G, ρ_R)

# Example
```julia
ρ_G, ρ_R, ρ = create_calibrated_product_intensity(50.0)
stats = marginal_stats(ρ)
println("E[N] = ", stats.E_N)  # ≈ 50
```
"""
function create_calibrated_product_intensity(target_lambda::Real;
                                              mean_G::AbstractVector=[0.6, 0.4],
                                              mean_R::AbstractVector=[0.5, 0.5],
                                              conc::Real=15.0,
                                              rng::AbstractRNG=Random.default_rng())
    # Create unit-scale intensities to measure actual total intensity
    ρ_G_unit = BdPlusMixture([1.0], [mean_G], [conc], 1.0)
    ρ_R_unit = BdPlusMixture([1.0], [mean_R], [conc], 1.0)

    # Compute total intensities via quadrature
    c_G_unit = total_intensity(ρ_G_unit)
    c_R_unit = total_intensity(ρ_R_unit)

    # Compute scale factor to achieve target Λ
    # With equal scaling: (s × c_G_unit) × (s × c_R_unit) = target_lambda
    # So s² × c_G_unit × c_R_unit = target_lambda
    s = sqrt(target_lambda / (c_G_unit * c_R_unit))

    # Create scaled intensities
    ρ_G = BdPlusMixture([1.0], [mean_G], [conc], s)
    ρ_R = BdPlusMixture([1.0], [mean_R], [conc], s)
    ρ = ProductIntensity(ρ_G, ρ_R)

    return ρ_G, ρ_R, ρ
end
