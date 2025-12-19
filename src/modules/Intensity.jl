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
    p = 0.0
    for (w, μ, κ) in zip(bm.weights, bm.means, bm.concentrations)
        dist_sq = sum((x[i] - μ[i])^2 for i in 1:d)
        p += w * exp(-κ * dist_sq)
    end

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
    total_intensity(ρ::BdPlusMixture; n_samples=10000, rng=Random.default_rng()) -> Float64

Compute the total intensity c = ∫ρ(x)dx via Monte Carlo integration.
"""
function total_intensity(ρ::BdPlusMixture{d};
                         n_samples::Int=10000,
                         rng::AbstractRNG=Random.default_rng()) where d
    # Monte Carlo integration over B^d_+
    total = 0.0
    vol = Bd_plus_volume(d)

    for _ in 1:n_samples
        x = uniform_Bd_plus_sample(d; rng=rng)
        total += ρ(x)
    end

    return vol * total / n_samples
end

"""
    marginal_total_intensity(ρ::ProductIntensity; n_samples=10000) -> Tuple{Float64, Float64}

Compute the marginal total intensities c_G = ∫ρ_G(g)dg and c_R = ∫ρ_R(r)dr.
Returns (c_G, c_R).
"""
function marginal_total_intensity(ρ::ProductIntensity{d};
                                  n_samples::Int=10000,
                                  rng::AbstractRNG=Random.default_rng()) where d
    c_G = total_intensity(ρ.ρ_G; n_samples=n_samples, rng=rng)
    c_R = total_intensity(ρ.ρ_R; n_samples=n_samples, rng=rng)
    return (c_G, c_R)
end

"""
    intensity_weighted_mean(ρ::BdPlusMixture; n_samples=10000, rng=Random.default_rng()) -> SVector

Compute the intensity-weighted mean position μ = ∫x·ρ(x)dx via Monte Carlo.
"""
function intensity_weighted_mean(ρ::BdPlusMixture{d};
                                 n_samples::Int=10000,
                                 rng::AbstractRNG=Random.default_rng()) where d
    μ = zeros(d)
    vol = Bd_plus_volume(d)

    for _ in 1:n_samples
        x = uniform_Bd_plus_sample(d; rng=rng)
        μ .+= x .* ρ(x)
    end

    return SVector{d, Float64}(vol .* μ ./ n_samples)
end

"""
    normalized_mean(ρ::BdPlusMixture; n_samples=10000, rng=Random.default_rng()) -> SVector

Compute the normalized mean μ̃ = μ/c where μ is the intensity-weighted mean
and c is the total intensity.
"""
function normalized_mean(ρ::BdPlusMixture{d};
                         n_samples::Int=10000,
                         rng::AbstractRNG=Random.default_rng()) where d
    c = total_intensity(ρ; n_samples=n_samples, rng=rng)
    μ = intensity_weighted_mean(ρ; n_samples=n_samples, rng=rng)
    return μ ./ c
end

"""
    marginal_stats(ρ::ProductIntensity; n_samples=10000, rng=Random.default_rng())

Compute all marginal statistics for a product intensity:
- c_G, c_R: marginal total intensities
- μ_G, μ_R: intensity-weighted mean positions
- μ̃_G, μ̃_R: normalized means

Returns a named tuple with all quantities.
"""
function marginal_stats(ρ::ProductIntensity{d};
                        n_samples::Int=10000,
                        rng::AbstractRNG=Random.default_rng()) where d
    c_G = total_intensity(ρ.ρ_G; n_samples=n_samples, rng=rng)
    c_R = total_intensity(ρ.ρ_R; n_samples=n_samples, rng=rng)

    μ_G = intensity_weighted_mean(ρ.ρ_G; n_samples=n_samples, rng=rng)
    μ_R = intensity_weighted_mean(ρ.ρ_R; n_samples=n_samples, rng=rng)

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
        E_edges_edge_centric = E_N * avg_conn_prob
    )
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
