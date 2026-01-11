module IDPG

using LinearAlgebra
using StaticArrays
using Distributions
using Random
using Graphs
using StatsBase
using Clustering
using CairoMakie
using GraphMakie
using SpecialFunctions: gamma
using OrdinaryDiffEq
using Integrals
using Cuba

# Include component files (not submodules, just code)
# Note: Order matters - PDEEvolution before Sampling because Sampling uses BdPlusGrid
include("modules/LatentSpace.jl")
include("modules/Intensity.jl")
include("modules/PDEEvolution.jl")
include("modules/Sampling.jl")
include("modules/GraphGeneration.jl")
include("modules/Visualization.jl")
include("modules/EcologicalUtils.jl")

# Exports

# Latent Space (B^d_+ = non-negative unit ball)
export LatentPoint
export in_Bd_plus, on_Bd_plus_boundary, Bd_plus_outward_normal
export project_to_Bd_plus, uniform_Bd_plus_sample, Bd_plus_volume
export connection_probability, radial_coordinate, angular_coordinates
export hyperspherical_to_cartesian, cartesian_to_hyperspherical, hyperspherical_jacobian
export Bd_plus_from_hyperspherical, Bd_plus_to_hyperspherical

# Intensity
export AbstractIntensity, ProductIntensity, BdPlusMixture, TimeVaryingIntensity
export MixtureOfProductIntensities
export marginal_stats, total_intensity, sample_from_mixture
export marginal_total_intensity, intensity_weighted_mean, normalized_mean
export n_species, species_intensities, species_probabilities
export create_calibrated_product_intensity

# 1D intensity utilities
export truncated_gaussian_sigma, truncated_gaussian_kappa
export compute_1d_marginals, compute_bound_heat_matrix

# Edge-centric intensity (joint source-target distributions)
export AbstractEdgeIntensity, ScaledProductEdgeIntensity, SymmetricEdgeIntensity
export edge_intensity

# Sampling
export InteractionSite, sample_ppp, sample_ppp_product, estimate_max_intensity
export sample_from_grid, sample_from_grid_full, initialize_grid_from_mixture
export sample_site_from_mixture, sample_ppp_mixture, sample_ppp_mixture_sites_only

# Edge-centric sampling (true E[N]/2 semantics)
export sample_edge_centric, sample_edge_centric_ppp, sample_edge_centric_symmetric

# Graph generation
export EdgeCentricSample, FullEdgeCentricSample
export generate_node_centric, generate_node_centric_full
export generate_edge_centric
export discretize_edge_centric, discretize_edge_centric_joint, discretize_with_weights
export source_g, source_r, target_g, target_r, to_edge_centric

# PDE evolution (OrdinaryDiffEq-based)
export BdPlusGrid, create_Bd_plus_grid
export evolve_diffusion, evolve_diffusion!
export evolve_advection, evolve_advection!
export evolve_advection_field, evolve_advection_field!
export evolve_reaction_diffusion, evolve_reaction_diffusion!
export evolve_and_track
export compute_mean_position, gradient_component, get_neighbors, laplacian_stencil

# Visualization
export plot_intensity_Bd_plus, plot_intensity_Bd_plus!
export plot_node_centric_graph, plot_node_centric_graph!
export plot_edge_centric, plot_edge_centric!
export animate_evolution, plot_formula_validation
export Bd_plus_to_2d, draw_Bd_plus_boundary!, plot_sites_Bd_plus
export latex_figure_theme, with_latex_theme

# Ecological utilities
export assign_site_to_guild, assign_point_to_guild
export build_full_guild_means
export compute_foodweb_matrix, normalize_foodweb_matrix
export sample_guild_position
export compute_expected_guild_edges, compute_guild_affinity
export trophic_layout

# Food web centroid construction
export default_trophic_affinity, construct_food_web_centroids
export svd_initialize_centroids, verify_centroids
export project_rows_to_Bd_plus!

end # module IDPG
