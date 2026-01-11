#!/bin/bash
# Run all figure-producing scripts for the manuscript
# Usage: ./scripts/run_all_figures.sh

set -e  # Exit on first error

cd "$(dirname "$0")/.."  # Go to project root

echo "============================================================"
echo "IDPG Manuscript Figure Generation"
echo "Started: $(date)"
echo "============================================================"

# Heat maps
echo -e "\n--- Heat Maps ---"
julia --project=. scripts/heat_maps/figure_A_anatomy.jl
julia --project=. scripts/heat_maps/figure_A_prime_asymmetric.jl
julia --project=. scripts/heat_maps/figure_A_double_prime_nonproduct.jl
julia --project=. scripts/heat_maps/figure_B_pde_absorbing.jl
julia --project=. scripts/heat_maps/figure_B_pde_reflecting.jl
julia --project=. scripts/heat_maps/figure_C_spectral_1d.jl
julia --project=. scripts/heat_maps/figure_C_spectral_4d.jl
julia --project=. scripts/heat_maps/figure_D_foodweb_static.jl
julia --project=. scripts/heat_maps/figure_E_foodweb_dynamic.jl

# Scaling laws
echo -e "\n--- Scaling Laws ---"
julia --project=. scripts/scaling_laws/fig1_intensity_scaling.jl
julia --project=. scripts/scaling_laws/fig2_entity_lifetime.jl
julia --project=. scripts/scaling_laws/fig6_temporal_evolution.jl

echo -e "\n============================================================"
echo "All figures generated successfully!"
echo "Finished: $(date)"
echo "============================================================"
