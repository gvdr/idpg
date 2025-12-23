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
using ModelingToolkit
using MethodOfLines
using OrdinaryDiffEq
using DomainSets

# Include component files (not submodules, just code)
# Note: Order matters - PDEEvolution before Sampling because Sampling uses BdPlusGrid
include("modules/LatentSpace.jl")
include("modules/Intensity.jl")
include("modules/PDEEvolution.jl")
include("modules/PDESciML.jl")
include("modules/Sampling.jl")
include("modules/GraphGeneration.jl")
include("modules/Visualization.jl")

# Exports

# Latent Space (B^d_+ = non-negative unit ball)
export LatentPoint
export in_Bd_plus, on_Bd_plus_boundary, Bd_plus_outward_normal
export project_to_Bd_plus, uniform_Bd_plus_sample, Bd_plus_volume
export connection_probability, radial_coordinate, angular_coordinates
export hyperspherical_to_cartesian, cartesian_to_hyperspherical
export Bd_plus_from_hyperspherical, Bd_plus_to_hyperspherical

# Intensity
export AbstractIntensity, ProductIntensity, BdPlusMixture, TimeVaryingIntensity
export MixtureOfProductIntensities
export marginal_stats, total_intensity, sample_from_mixture
export marginal_total_intensity, intensity_weighted_mean, normalized_mean
export n_species, species_intensities, species_probabilities

# Sampling
export InteractionSite, sample_ppp, sample_ppp_product, estimate_max_intensity
export sample_from_grid, sample_from_grid_full, initialize_grid_from_mixture
export sample_site_from_mixture, sample_ppp_mixture, sample_ppp_mixture_sites_only

# Graph generation
export EdgeCentricSample, FullEdgeCentricSample
export generate_node_centric, generate_edge_centric, generate_edge_centric_full
export discretize_edge_centric, discretize_edge_centric_joint, discretize_with_weights
export source_g, source_r, target_g, target_r, to_edge_centric

# PDE evolution (hand-coded finite differences)
export BdPlusGrid, create_Bd_plus_grid
export evolve_diffusion!, evolve_advection!, evolve_reaction_diffusion!
export evolve_advection_field!
export evolve_and_track
export compute_mean_position, gradient_component, get_neighbors

# PDE evolution (SciML: MethodOfLines + OrdinaryDiffEq)
export create_Bd_plus_mask_2d, create_Bd_plus_mask_4d, apply_Bd_plus_mask!
export solve_diffusion_mol_2d, solve_advection_mol_2d
export solve_diffusion_mol_4d, solve_advection_mol_4d

# Visualization
export plot_intensity_Bd_plus, plot_intensity_Bd_plus!
export plot_node_centric_graph, plot_node_centric_graph!
export plot_edge_centric, plot_edge_centric!
export animate_evolution, plot_formula_validation
export Bd_plus_to_2d, draw_Bd_plus_boundary!, plot_sites_Bd_plus

end # module IDPG
