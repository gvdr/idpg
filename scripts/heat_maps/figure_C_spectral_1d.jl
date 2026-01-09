# Figure C: Spectral Structure (d=1)
#
# Dimension: d = 1
# Purpose: Illustrate the spectral decomposition of the bound heat operator
#
# Panels (2×2 grid):
#   (a) Singular value spectrum: σ_n vs n (log scale on y-axis)
#   (b) First 3 left singular functions: u_1(g), u_2(g), u_3(g)
#   (c) First 3 right singular functions: v_1(r), v_2(r), v_3(r)
#   (d) Reconstruction error: ||h - h_k||_L² vs rank k
#
# Usage:
#   julia --project=. scripts/heat_maps/figure_C_spectral_1d.jl          # compute + plot
#   julia --project=. scripts/heat_maps/figure_C_spectral_1d.jl --plot   # plot from saved data

using IDPG
using CairoMakie
using LaTeXStrings
using Distributions
using Statistics
using LinearAlgebra
using JLD2

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(pkgdir(IDPG), "output", "heat_maps")
const DATA_FILE = joinpath(OUTPUT_DIR, "figure_C_data.jld2")
const GRID_RES = 100

# Intensity parameters (d=1)
const MU_G = 0.3
const MU_R = 0.7
const SIGMA = 0.12

# Number of singular values to show
const N_SINGULAR = 20
const N_FUNCTIONS = 1  # Only show first singular function (u2,u3,v2,v3 are numerical noise)

# ============================================================================
# Intensity Functions
# ============================================================================

# Note: Uses truncated_gaussian_sigma from IDPG library
const truncated_gaussian = truncated_gaussian_sigma

"""
    compute_intensities(grid)

Compute normalized intensities on grid using library functions.
"""
function compute_intensities(grid::Vector{Float64})
    ρ_G, ρ_R = compute_1d_marginals(grid, MU_G, MU_R, SIGMA)

    # Normalize
    dx = grid[2] - grid[1]
    c_G = sum(ρ_G) * dx
    c_R = sum(ρ_R) * dx

    return ρ_G ./ c_G, ρ_R ./ c_R
end

# ============================================================================
# Bound Heat and SVD
# ============================================================================

"""
    compute_bound_heat_matrix(grid, ρ_G_norm, ρ_R_norm)

Compute bound heat matrix h(g,r) = (g·r) · ρ̃_G(g) · ρ̃_R(r)
Scaled by dx for proper L² discretization.
"""
function compute_bound_heat_matrix(grid::Vector{Float64}, ρ_G_norm::Vector{Float64}, ρ_R_norm::Vector{Float64})
    n = length(grid)
    dx = grid[2] - grid[1]

    H = zeros(n, n)
    for (i, g) in enumerate(grid)
        for (j, r) in enumerate(grid)
            K = g * r
            # Scale by sqrt(dx) on each side for proper L² inner product
            H[i, j] = K * sqrt(ρ_G_norm[i]) * sqrt(ρ_R_norm[j]) * dx
        end
    end

    return H
end

"""
    compute_spectral_decomposition(grid, ρ_G_norm, ρ_R_norm)

Compute SVD of bound heat operator.
Returns singular values and functions in physical (weighted) coordinates.
"""
function compute_spectral_decomposition(grid::Vector{Float64}, ρ_G_norm::Vector{Float64}, ρ_R_norm::Vector{Float64})
    n = length(grid)
    dx = grid[2] - grid[1]

    # Build bound heat matrix
    H = compute_bound_heat_matrix(grid, ρ_G_norm, ρ_R_norm)

    # SVD
    F = svd(H)

    # Singular values
    σ = F.S

    # Convert singular functions back to physical coordinates
    # U contains left singular vectors, V contains right singular vectors
    # The functions are U[:, n] and V[:, n]
    # We need to weight by sqrt(ρ) to get the L² functions
    U_phys = zeros(n, length(σ))
    V_phys = zeros(n, length(σ))

    for k in 1:min(N_SINGULAR, length(σ))
        # u_k(g) = U[i,k] / sqrt(ρ_G_norm[i] * dx) (approximately)
        # For visualization, we just use the raw singular vectors scaled
        U_phys[:, k] = F.U[:, k] ./ sqrt(dx)
        V_phys[:, k] = F.V[:, k] ./ sqrt(dx)
    end

    return σ, U_phys, V_phys, H
end

"""
    compute_reconstruction_errors(H, σ)

Compute ||H - H_k||_F for various ranks k.
"""
function compute_reconstruction_errors(H::Matrix{Float64}, σ::Vector{Float64})
    total_norm_sq = sum(σ.^2)
    errors = Float64[]

    for k in 0:min(N_SINGULAR, length(σ))
        if k == 0
            err_sq = total_norm_sq
        else
            err_sq = sum(σ[(k+1):end].^2)
        end
        push!(errors, sqrt(err_sq))
    end

    return errors
end

# ============================================================================
# Data Computation
# ============================================================================

"""
    compute_figure_C_data()

Compute all data for Figure C and save to JLD2.
"""
function compute_figure_C_data()
    println("=" ^ 60)
    println("Figure C: Computing Spectral Data (d=1)")
    println("=" ^ 60)

    # Setup grid
    grid = collect(range(0, 1, length=GRID_RES))
    dx = grid[2] - grid[1]

    # Compute intensities
    ρ_G_norm, ρ_R_norm = compute_intensities(grid)

    println("Intensity parameters:")
    println("  μ_G = " * string(MU_G) * ", μ_R = " * string(MU_R) * ", σ = " * string(SIGMA))

    # Compute spectral decomposition
    println("\nComputing SVD...")
    σ, U, V, H = compute_spectral_decomposition(grid, ρ_G_norm, ρ_R_norm)

    println("  Top 5 singular values: " * string(round.(σ[1:5], digits=6)))
    println("  Hilbert-Schmidt norm: " * string(round(sqrt(sum(σ.^2)), digits=6)))

    # Reconstruction errors
    errors = compute_reconstruction_errors(H, σ)
    println("  Reconstruction errors at k=1,3,5: " * string(round.(errors[[2, 4, 6]], digits=6)))

    # Save data
    mkpath(OUTPUT_DIR)
    data = (
        grid = grid,
        σ = σ,
        U = U,
        V = V,
        errors = errors,
        params = (μ_G = MU_G, μ_R = MU_R, σ = SIGMA),
    )
    @save DATA_FILE data
    println("\nSaved data: " * DATA_FILE)

    return data
end

"""
    load_figure_C_data()

Load precomputed data for Figure C.
"""
function load_figure_C_data()
    @load DATA_FILE data
    return data
end

# ============================================================================
# Visualization
# ============================================================================

"""
    plot_figure_C(data)

Generate Figure C visualization from precomputed data.
"""
function plot_figure_C(data)
    println("=" ^ 60)
    println("Figure C: Generating Plot")
    println("=" ^ 60)

    grid = data.grid
    σ = data.σ
    U = data.U
    V = data.V
    errors = data.errors

    # Use LaTeX fonts with larger sizes
    fig = with_theme(merge(theme_latexfonts(), Theme(fontsize=18))) do
        fig = Figure(size=(1000, 900))

        # Panel (a): Singular value spectrum
        ax_a = Axis(fig[1, 1],
            xlabel = L"Index $n$",
            ylabel = L"Singular value $\sigma_n$",
            title = "(a) Singular Value Spectrum",
            yscale = log10)

        n_show = min(N_SINGULAR, length(σ))
        scatter!(ax_a, 1:n_show, σ[1:n_show], color=:blue, markersize=8)
        lines!(ax_a, 1:n_show, σ[1:n_show], color=:blue, linewidth=1)

        # Panel (b): Left singular functions
        ax_b = Axis(fig[1, 2],
            xlabel = L"$g$",
            ylabel = L"$u_n(g)$",
            title = "(b) Left Singular Function")

        colors_b = [:blue, :red, :green]
        for k in 1:N_FUNCTIONS
            σ_str = string(round(σ[k], digits=4))
            lines!(ax_b, grid, U[:, k], color=colors_b[k], linewidth=2,
                   label="u" * string(k) * " (σ=" * σ_str * ")")
        end
        axislegend(ax_b, position=:rt)

        # Panel (c): Right singular functions
        ax_c = Axis(fig[2, 1],
            xlabel = L"$r$",
            ylabel = L"$v_n(r)$",
            title = "(c) Right Singular Function")

        colors_c = [:blue, :red, :green]
        for k in 1:N_FUNCTIONS
            lines!(ax_c, grid, V[:, k], color=colors_c[k], linewidth=2,
                   label="v" * string(k))
        end
        axislegend(ax_c, position=:rt)

        # Panel (d): Reconstruction error
        ax_d = Axis(fig[2, 2],
            xlabel = L"Rank $k$",
            ylabel = "Reconstruction error",
            title = "(d) Reconstruction Error")

        scatter!(ax_d, 0:(length(errors)-1), errors, color=:blue, markersize=8)
        lines!(ax_d, 0:(length(errors)-1), errors, color=:blue, linewidth=1)

        # Mark 90% and 99% reconstruction
        total_norm = errors[1]
        for thresh in [0.1, 0.01]
            for (k, err) in enumerate(errors)
                if err / total_norm < thresh
                    hlines!(ax_d, [err], color=:red, linestyle=:dash, alpha=0.5)
                    pct_str = string(round(100 * (1 - thresh), digits=0))
                    text!(ax_d, k-1, err * 1.1,
                          text = pct_str * "% at k=" * string(k-1),
                          fontsize = 12)
                    break
                end
            end
        end

        # Add overall title
        Label(fig[0, 1:2], "Spectral Decomposition of Bound Heat Operator (d=1)", fontsize=22)

        fig
    end

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "figure_C_spectral_1d.png")
    save(output_file, fig, px_per_unit=2)
    println("Saved: " * output_file)

    return fig
end

# ============================================================================
# Entry Point
# ============================================================================

function (@main)(args)
    if "--plot" in args
        data = load_figure_C_data()
        plot_figure_C(data)
    else
        data = compute_figure_C_data()
        plot_figure_C(data)
    end
    return 0
end
