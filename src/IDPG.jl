module IDPG

using LinearAlgebra
using StaticArrays
using Distributions
using Random
using Graphs
using StatsBase
using CairoMakie
using GraphMakie
using SpecialFunctions: gamma

# Include component files (not submodules, just code)
include("modules/LatentSpace.jl")
include("modules/Intensity.jl")
include("modules/Sampling.jl")
include("modules/GraphGeneration.jl")
include("modules/PDEEvolution.jl")
include("modules/Visualization.jl")

# Exports

# Latent Space (B^d_+ = non-negative unit ball)
export LatentPoint
export in_Bd_plus, on_Bd_plus_boundary, Bd_plus_outward_normal
export project_to_Bd_plus, uniform_Bd_plus_sample, Bd_plus_volume
export connection_probability, radial_coordinate, angular_coordinates

# Intensity
export AbstractIntensity, ProductIntensity, BdPlusMixture, TimeVaryingIntensity
export marginal_stats, total_intensity, sample_from_mixture
export marginal_total_intensity, intensity_weighted_mean, normalized_mean

# Sampling
export InteractionSite, sample_ppp, sample_ppp_product, estimate_max_intensity

# Graph generation
export EdgeCentricSample
export generate_node_centric, generate_edge_centric, discretize_edge_centric

# PDE evolution
export BdPlusGrid, create_Bd_plus_grid
export evolve_diffusion!, evolve_advection!, evolve_reaction_diffusion!
export evolve_and_track

# Visualization
export plot_intensity_Bd_plus, plot_intensity_Bd_plus!
export plot_node_centric_graph, plot_node_centric_graph!
export plot_edge_centric, plot_edge_centric!
export animate_evolution, plot_formula_validation
export Bd_plus_to_2d, draw_Bd_plus_boundary!, plot_sites_Bd_plus

end # module IDPG
