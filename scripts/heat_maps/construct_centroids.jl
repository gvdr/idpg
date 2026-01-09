# Phase 1: Food Web Centroid Construction
#
# Goal: Find guild centroids μ_i^G and μ_j^R in B^d_+ such that
#       μ_i^G · μ_j^R ≈ K*_ij (target trophic affinity matrix)
#
# This script uses the construct_food_web_centroids function from the IDPG library.
#
# Reference: docs/simulation_for_heat_maps.md
#
# Usage: julia --project=. scripts/heat_maps/construct_centroids.jl

using IDPG
using JLD2

# ============================================================================
# Entry Point
# ============================================================================

"""
Main entry point for command-line execution.
"""
function (@main)(args)
    # Use library function to construct centroids
    result = construct_food_web_centroids()

    # Save results
    output_file = joinpath(pkgdir(IDPG), "output", "heat_maps", "foodweb_centroids.jld2")
    mkpath(dirname(output_file))
    @save output_file result
    println("\nSaved centroids to: " * output_file)

    return 0
end
