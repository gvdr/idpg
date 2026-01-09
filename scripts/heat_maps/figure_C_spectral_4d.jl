# Figure C': Spectral Structure (d=4, Food Web)
#
# Dimension: d = 4
# Purpose: Show spectral decomposition for ecologically realistic case
#
# Uses guild centroids from Phase 1
#
# Panels:
#   (a) Singular value spectrum {σ_n}
#   (b) Cumulative explained variance: Σ_{n≤k} σ_n² / Σ_n σ_n² vs k
#   (c) Guild loadings on first 2 singular modes
#   (d) Comparison: SVD spectrum of h̄ vs spectrum of target K*
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_C_spectral_4d.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_C_spectral_4d.jl --plot   # plot from saved data

using IDPG
using CairoMakie
using LaTeXStrings
using Distributions
using Statistics
using LinearAlgebra
using JLD2
using Random

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "heat_maps")
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_C_prime_data.jld2")
const CENTROIDS_FILE = joinpath(pkgdir(IDPG), "output", "heat_maps", "foodweb_centroids.jld2")

# Nyström parameters
const N_SAMPLES = 500
const KAPPA = 30.0  # Concentration for truncated Gaussian around centroids

# Random seed for reproducibility
Random.seed!(42)

# ============================================================================
# Load Centroids
# ============================================================================

"""
    load_centroids()

Load guild centroids from Phase 1 output.
"""
function load_centroids()
    @load CENTROIDS_FILE result
    return result
end

# ============================================================================
# Sampling from Mixture of Truncated Gaussians
# ============================================================================

# Note: Uses sample_guild_position and project_to_Bd_plus from IDPG library

"""
    sample_from_guild(μ, κ, n_samples)

Sample n points from truncated Gaussian centered at μ with concentration κ.
Uses library's sample_guild_position internally.
"""
function sample_from_guild(μ::Vector{Float64}, κ::Float64, n_samples::Int)
    samples = Vector{Vector{Float64}}(undef, n_samples)
    for i in 1:n_samples
        samples[i] = sample_guild_position(μ, κ)
    end
    return samples
end

"""
    sample_mixture(M, weights, κ, n_total)

Sample from mixture of guild distributions.
M: n_guilds × d matrix of centroids
weights: n_guilds vector of mixture weights
"""
function sample_mixture(M::Matrix{Float64}, weights::Vector{Float64}, κ::Float64, n_total::Int)
    n_guilds, d = size(M)

    # Normalize weights
    w = weights ./ sum(weights)

    # Determine samples per guild
    n_per_guild = round.(Int, w .* n_total)
    # Adjust to get exactly n_total
    n_per_guild[end] += n_total - sum(n_per_guild)

    # Sample from each guild
    samples = Vector{Vector{Float64}}()
    guild_labels = Int[]

    for g in 1:n_guilds
        μ = M[g, :]
        guild_samples = sample_from_guild(μ, κ, n_per_guild[g])
        append!(samples, guild_samples)
        append!(guild_labels, fill(g, n_per_guild[g]))
    end

    return samples, guild_labels
end

# ============================================================================
# Bound Heat Computation via Nyström
# ============================================================================

"""
    compute_bound_heat_nystrom(samples_G, samples_R, weights_G, weights_R)

Compute discretized bound heat matrix using Nyström approximation.
H[i,j] = K(g_i, r_j) * w_G[i] * w_R[j]
where K(g,r) = g · r (dot product kernel).
"""
function compute_bound_heat_nystrom(samples_G::Vector{Vector{Float64}},
                                      samples_R::Vector{Vector{Float64}},
                                      weights_G::Vector{Float64},
                                      weights_R::Vector{Float64})
    n_G = length(samples_G)
    n_R = length(samples_R)

    H = zeros(n_G, n_R)

    for i in 1:n_G
        for j in 1:n_R
            K = dot(samples_G[i], samples_R[j])
            H[i, j] = K * sqrt(weights_G[i]) * sqrt(weights_R[j])
        end
    end

    return H
end

"""
    compute_guild_loadings(U, guild_labels, n_guilds, n_modes)

Compute how strongly each guild loads onto each singular mode.
Returns n_guilds × n_modes matrix.
"""
function compute_guild_loadings(U::Matrix{Float64}, guild_labels::Vector{Int},
                                 n_guilds::Int, n_modes::Int)
    loadings = zeros(n_guilds, n_modes)

    for g in 1:n_guilds
        idx = findall(guild_labels .== g)
        if length(idx) > 0
            for k in 1:n_modes
                # RMS loading for guild g on mode k
                loadings[g, k] = sqrt(mean(U[idx, k].^2))
            end
        end
    end

    return loadings
end

# ============================================================================
# Data Computation
# ============================================================================

"""
    compute_figure_C_prime_data()

Compute all data for Figure C' and save to JLD2.
"""
function compute_figure_C_prime_data()
    println("=" ^ 60)
    println("Figure C': Computing Spectral Data (d=4, Food Web)")
    println("=" ^ 60)

    # Load centroids
    println("Loading centroids from Phase 1...")
    result = load_centroids()
    M_G = result.M_G  # 5 × 4 matrix
    M_R = result.M_R
    K_star = result.K_star
    K_approx = result.K_approx
    guild_names = result.guild_names

    n_guilds = size(M_G, 1)
    d = size(M_G, 2)

    println("  Guilds: " * join(guild_names, ", "))
    println("  Latent dimension: d = " * string(d))

    # Equal weights for all guilds
    weights = fill(1.0 / n_guilds, n_guilds)

    # Sample from guild distributions
    println("\nSampling from guild distributions (n=" * string(N_SAMPLES) * ", κ=" * string(KAPPA) * ")...")
    samples_G, guild_labels_G = sample_mixture(M_G, weights, KAPPA, N_SAMPLES)
    samples_R, guild_labels_R = sample_mixture(M_R, weights, KAPPA, N_SAMPLES)

    # Uniform sample weights for simplicity
    w_G = fill(1.0 / N_SAMPLES, N_SAMPLES)
    w_R = fill(1.0 / N_SAMPLES, N_SAMPLES)

    # Compute bound heat matrix
    println("Computing bound heat matrix...")
    H = compute_bound_heat_nystrom(samples_G, samples_R, w_G, w_R)

    # SVD of bound heat
    println("Computing SVD...")
    F = svd(H)
    σ = F.S

    println("  Top 10 singular values: " * string(round.(σ[1:min(10, length(σ))], digits=6)))
    println("  Rank (σ > 1e-10): " * string(sum(σ .> 1e-10)))

    # SVD of target affinity matrix
    F_K = svd(K_star)
    σ_K = F_K.S

    println("  K* singular values: " * string(round.(σ_K, digits=4)))

    # Cumulative explained variance
    total_var = sum(σ.^2)
    cum_var = cumsum(σ.^2) ./ total_var

    # Guild loadings on first modes
    loadings_G = compute_guild_loadings(Matrix(F.U), guild_labels_G, n_guilds, 3)
    loadings_R = compute_guild_loadings(Matrix(F.V), guild_labels_R, n_guilds, 3)

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        σ = σ,
        σ_K = σ_K,
        cum_var = cum_var,
        loadings_G = loadings_G,
        loadings_R = loadings_R,
        guild_names = guild_names,
        n_guilds = n_guilds,
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

"""
    load_figure_C_prime_data()

Load precomputed data for Figure C'.
"""
function load_figure_C_prime_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

"""
    plot_figure_C_prime(data)

Generate Figure C' visualization from precomputed data.
"""
function plot_figure_C_prime(data)
    println("=" ^ 60)
    println("Figure C': Generating Plot")
    println("=" ^ 60)

    σ = data.σ
    σ_K = data.σ_K
    cum_var = data.cum_var
    loadings_G = data.loadings_G
    guild_names = data.guild_names
    n_guilds = data.n_guilds

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
        fig = Figure(size=(1000, 900))

        # Panel (a): Singular value spectrum
        ax_a = Axis(fig[1, 1],
            xlabel = L"Index $n$",
            ylabel = L"Singular value $\sigma_n$",
            title = "(a) Singular Value Spectrum",
            yscale = log10)

        n_show = min(20, length(σ))
        scatter!(ax_a, 1:n_show, σ[1:n_show] .+ 1e-15, color=:blue, markersize=8, label=L"Bound heat $\bar{h}$")
        lines!(ax_a, 1:n_show, σ[1:n_show] .+ 1e-15, color=:blue, linewidth=1)

        # Panel (b): Cumulative explained variance
        ax_b = Axis(fig[1, 2],
            xlabel = L"Rank $k$",
            ylabel = "Cumulative variance explained",
            title = "(b) Explained Variance")

        n_cum = min(15, length(cum_var))
        lines!(ax_b, 1:n_cum, cum_var[1:n_cum], color=:blue, linewidth=2)
        scatter!(ax_b, 1:n_cum, cum_var[1:n_cum], color=:blue, markersize=8)
        hlines!(ax_b, [0.9, 0.95, 0.99], color=:red, linestyle=:dash, alpha=0.5)

        # Find rank for 95% variance
        k_95 = findfirst(cum_var .>= 0.95)
        if !isnothing(k_95)
            text!(ax_b, k_95, 0.95, text="95% at k=" * string(k_95), fontsize=12, align=(:left, :bottom))
        end

        # Panel (c): Guild loadings on first 2 modes
        ax_c = Axis(fig[2, 1],
            xlabel = "Mode 1 loading",
            ylabel = "Mode 2 loading",
            title = "(c) Guild Loadings (G coordinates)")

        colors = [:green, :orange, :purple, :red, :black]
        for g in 1:n_guilds
            scatter!(ax_c, [loadings_G[g, 1]], [loadings_G[g, 2]],
                     color=colors[g], markersize=15, label=guild_names[g])
        end
        axislegend(ax_c, position=:lt)

        # Panel (d): Compare spectra
        ax_d = Axis(fig[2, 2],
            xlabel = L"Index $n$",
            ylabel = "Singular value",
            title = "(d) Spectrum Comparison")

        # Normalize for comparison
        σ_norm = σ[1:5] ./ σ[1]
        σ_K_norm = σ_K ./ σ_K[1]

        scatter!(ax_d, 1:5, σ_norm, color=:blue, markersize=12, label=L"Bound heat $\bar{h}$ (scaled)")
        lines!(ax_d, 1:5, σ_norm, color=:blue, linewidth=2)

        scatter!(ax_d, 1:length(σ_K_norm), σ_K_norm, color=:red, markersize=12, marker=:diamond, label=L"Target $K^*$ (scaled)")
        lines!(ax_d, 1:length(σ_K_norm), σ_K_norm, color=:red, linewidth=2, linestyle=:dash)

        axislegend(ax_d, position=:rt)

        # Add overall title
        Label(fig[0, 1:2], L"Spectral Decomposition of Bound Heat $\bar{h}$ ($d=4$, Food Web)", fontsize=22)

        fig
    end

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_C_spectral_4d.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_C_prime_data()
        plot_figure_C_prime(data)
    else
        data = compute_figure_C_prime_data()
        plot_figure_C_prime(data)
    end
    return 0
end
