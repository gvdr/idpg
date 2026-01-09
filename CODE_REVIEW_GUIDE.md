# IDPG Code Review Guide

This document provides a structured guide for reviewing the IDPG (Intensity Dot Product Graph) codebase. It explains what each component does, what to check during review, and how the pieces fit together.

---

## Table of Contents

1. [Repository Overview](#repository-overview)
2. [Library Modules (src/)](#library-modules)
3. [Scaling Laws Scripts](#scaling-laws-scripts)
4. [Heat Map Scripts](#heat-map-scripts)
5. [Tests](#tests)
6. [Review Checklist](#review-checklist)

---

## Repository Overview

```
idpg/
├── src/                    # Core library
│   ├── IDPG.jl            # Main module entry point
│   └── modules/           # 8 implementation files
├── scripts/
│   ├── scaling_laws/      # Validation simulations (fig1, fig2, fig6)
│   └── heat_maps/         # Paper figure generation (A-E)
├── test/                  # Unit tests
├── docs/                  # Documentation and manuscript
├── output/
│   ├── scaling_laws/      # Simulation outputs
│   └── heat_maps/         # Figure outputs
└── archive/               # Old/superseded material
```

### Design Principles

1. **Library-first**: Generic functionality lives in `src/modules/`
2. **Self-contained scripts**: Each simulation/example runs independently
3. **Data preservation**: Scripts save both figures and raw data (`.jld2`)
4. **Two-mode execution**: Many scripts support `--plot` flag to regenerate figures from saved data

---

## Library Modules

### Location: `src/modules/`

Each module focuses on one aspect of the IDPG framework. Files are included (not submodules) in `src/IDPG.jl`.

---

### 1. LatentSpace.jl (~274 lines)

**Purpose**: Geometry of B^d_+ (the non-negative unit ball in d dimensions).

**Key Types**:
- `LatentPoint{d}` = `SVector{d, Float64}` - a point in B^d_+

**Key Functions**:
- `in_Bd_plus(x)` - check if point is in B^d_+
- `project_to_Bd_plus(x)` - project point onto B^d_+
- `Bd_plus_volume(d)` - compute volume of B^d_+
- `random_point_Bd_plus(d)` - sample uniformly from B^d_+
- `spherical_to_cartesian`, `cartesian_to_spherical` - coordinate transforms

**Review Checklist**:
- [ ] Boundary handling at r=1 and coordinate axes
- [ ] Numerical stability for high dimensions (d > 4)
- [ ] Correct normalization for volume calculations

---

### 2. Intensity.jl (~662 lines)

**Purpose**: Intensity functions ρ(g,r) on G×R space. The heart of IDPG theory.

**Key Types**:
- `AbstractIntensity{d}` - base type for all intensities
- `BdPlusMixture{d}` - mixture of radial-angular distributions
- `ProductIntensity{d}` - factorized ρ(g,r) = ρ_G(g)·ρ_R(r)
- `MixtureOfProductIntensities{d}` - sum of species-specific products
- `ScaledProductEdgeIntensity{d}` - for edge-centric sampling with C_edge = E[N]/2

**Key Functions**:
- `(ρ::AbstractIntensity)(g, r)` - evaluate intensity at point
- `marginal_stats(ρ)` - compute c_G, c_R, μ_G, μ_R, E[N], E[|E|], E[L]
- `total_intensity(ρ)` - integrate over all space
- `marginal_G(ρ)`, `marginal_R(ρ)` - extract marginals

**Review Checklist**:
- [ ] Normalization: product intensities should integrate to Λ = c_G · c_R
- [ ] Centroids μ̃_G, μ̃_R computed correctly (weighted averages)
- [ ] Edge-centric formula: E[L] = (E[N]/2) · p̄ where p̄ = μ̃_G · μ̃_R
- [ ] Mixture weights sum correctly in MixtureOfProductIntensities

---

### 3. PDEEvolution.jl (~443 lines)

**Purpose**: Hand-coded finite difference PDE solver on B^d_+.

**Key Types**:
- `BdPlusGrid{d}` - finite difference grid with boundary info

**Key Functions**:
- `create_Bd_plus_grid(n, d)` - construct grid
- `evolve_diffusion!(ρ, grid, D, dt, n_steps)` - heat equation evolution
- `evolve_advection!(ρ, grid, v, dt, n_steps)` - transport equation
- `evolve_reaction_diffusion!(...)` - combined dynamics
- `apply_reflecting_bc!`, `apply_absorbing_bc!` - boundary conditions

**Review Checklist**:
- [ ] CFL condition: dt ≤ dx²/(2·D·d) for stability
- [ ] Mass conservation under reflecting BC
- [ ] Correct Neumann (reflecting) vs Dirichlet (absorbing) implementation
- [ ] Grid points correctly identified as interior vs boundary

---

### 4. PDESciML.jl (~414 lines)

**Purpose**: SciML-based PDE solver using MethodOfLines.jl.

**Key Functions**:
- `create_Bd_plus_mask_2d()`, `create_Bd_plus_mask_4d()` - domain masks
- `solve_diffusion_mol_2d()`, `solve_diffusion_mol_4d()` - diffusion on B^d_+
- `solve_advection_mol_2d()`, `solve_advection_mol_4d()` - advection

**Review Checklist**:
- [ ] Mask correctly defines B^d_+ region (norm ≤ 1, all coords ≥ 0)
- [ ] Boundary conditions applied at domain edge
- [ ] Initial condition correctly interpolated to PDESystem
- [ ] ODE solver tolerances appropriate for the problem

---

### 5. Sampling.jl (~567 lines)

**Purpose**: Sample sites from Poisson point process with intensity ρ.

**Key Types**:
- `InteractionSite{d}` - a sampled (g, r) pair

**Key Functions**:
- `sample_ppp(ρ, N_expected)` - rejection sampling from general intensity
- `sample_ppp_product(ρ)` - efficient sampling when ρ is product form
- `estimate_max_intensity(ρ)` - find sup for rejection sampling
- `sample_from_grid(ρ_G, ρ_R, grid, n)` - sample from discretized intensity
- `sample_edge_centric(ei)` - sample from ScaledProductEdgeIntensity

**Review Checklist**:
- [ ] Rejection sampling: proposal intensity ≥ true intensity everywhere
- [ ] Number of sites is Poisson-distributed with correct rate
- [ ] Product sampling correctly samples (g,r) independently
- [ ] Grid sampling handles edge cells correctly

---

### 6. GraphGeneration.jl (~400 lines)

**Purpose**: Generate graphs from sampled sites.

**Key Types**:
- `EdgeCentricSample{d}` - collection of edges (connection coords only)
- `FullEdgeCentricSample{d}` - edges with full (g,r) for source and target

**Key Functions**:
- `generate_node_centric(sites)` - sites become nodes, N² pairs tested
- `generate_node_centric_full(sites)` - returns FullEdgeCentricSample
- `generate_edge_centric(sites)` - each site IS an edge (linear in N)

**Review Checklist**:
- [ ] Node-centric: edge probability = g_source · r_target
- [ ] Edge-centric: all "accepted" sites become edges directly
- [ ] Correct graph type: SimpleDiGraph (directed)
- [ ] No self-loops in node-centric generation

---

### 7. Visualization.jl

**Purpose**: Plotting and animation utilities.

**Key Functions**:
- `plot_intensity_Bd_plus(rho, d)` - visualize intensity on B^d_+
- `plot_node_centric_graph(g, sites)` - node-centric graph visualization
- `plot_edge_centric(sample)` - edge-centric sample visualization
- `animate_evolution(...)` - animate PDE dynamics
- `Bd_plus_to_2d(point)` - project high-d point for 2D visualization
- `latex_figure_theme()` - create LaTeX-style Makie theme

**Review Checklist**:
- [ ] Color scales consistent across figures
- [ ] Axis labels use correct notation (g, r)
- [ ] LaTeX rendering works
- [ ] Figure sizes appropriate for paper submission

---

### 8. EcologicalUtils.jl

**Purpose**: Utilities for ecological modeling with guild structure.

**Key Functions**:
- `assign_site_to_guild(site, guild_means)` - classify site by nearest centroid
- `compute_foodweb_matrix(sample, guild_means)` - count guild-guild edges
- `normalize_foodweb_matrix(M; mode)` - normalize edge counts
- `default_trophic_affinity()` - 5x5 trophic affinity matrix
- `construct_food_web_centroids()` - SVD + alternating optimization
- `trophic_layout(n_guilds)` - layout positions for food web plots

**Review Checklist**:
- [ ] Guild assignment uses correct distance metric
- [ ] Centroid optimization converges
- [ ] All centroids remain in B^d_+ after projection

---

## Scaling Laws Scripts

### Location: `scripts/scaling_laws/`

These scripts validate theoretical predictions about quadratic (node-centric) vs linear (edge-centric) scaling with Monte Carlo simulations.

---

### fig1_intensity_scaling.jl

**Purpose**: Verify E[|E|] and E[L] formulas for node- vs edge-centric models.

**Key Results**:
- Node-centric: E[|E|] = E[N]^2 * p_bar (quadratic scaling)
- Edge-centric: E[L] = (E[N]/2) * p_bar (linear scaling)
- Ratio: E[|E|]_node / E[|E|]_edge = 2 * Lambda

**Output**: `output/scaling_laws/fig1_*.png`, `sim1_results.jld2`

**Review Checklist**:
- [ ] 1000+ replications for statistical significance
- [ ] Error bars computed correctly (std/sqrt(n))
- [ ] Theoretical predictions match empirical means within CI

---

### fig2_entity_lifetime.jl

**Purpose**: Demonstrate interpolation between edge-centric (mu -> 0) and node-centric (mu -> infinity) via entity lifetime parameter.

**Key Results**:
- Overlap probability p_overlap = (2/alpha^2)(alpha - 1 + e^{-alpha}) where alpha = W/mu
- Smooth transition between limiting regimes

**Output**: `output/scaling_laws/fig2_*.png`, `sim2_results.jld2`

**Review Checklist**:
- [ ] Lifetime parameter mu correctly controls regime
- [ ] mu -> 0 approaches edge-centric limit
- [ ] mu -> infinity approaches node-centric limit

---

### fig6_temporal_evolution.jl

**Purpose**: Temporal evolution of edge counts when underlying intensity follows PDE dynamics.

**Key Results**:
- Ratio = 2 * Lambda(t) holds at every time point
- Node-centric and edge-centric track intensity evolution differently

**Output**: `output/scaling_laws/fig6_*.png`, `sim6_results.jld2`

**Review Checklist**:
- [ ] Time stepping consistent with PDE solver
- [ ] Intensity evolution tracks expected pattern
- [ ] Edge count evolution follows from intensity

---

## Heat Map Scripts

### Location: `scripts/heat_maps/`

Systematic figure generation for the paper.

---

### construct_centroids.jl

**Purpose**: Find guild centroids mu_i^G and mu_j^R in B^d_+ that approximate a target trophic affinity matrix K*_ij = mu_i^G . mu_j^R. Uses SVD initialization followed by alternating optimization.

**Output**: `output/heat_maps/foodweb_centroids.jld2`

**Review Checklist**:
- [ ] Optimization converges to low reconstruction error
- [ ] All centroids lie in B^d_+ (norm <= 1, coords >= 0)
- [ ] Affinity matrix is non-negative and ecologically reasonable

---

### Figure Scripts

All follow the pattern:
```julia
function compute_figure_X_data()  # Compute and save
function plot_figure_X(data)       # Load and plot
function (@main)(args)             # Entry point with --plot flag
```

| Script | Content | Review Focus |
|--------|---------|--------------|
| figure_A_anatomy.jl | Product intensity basics | Marginals, bound heat, kernel |
| figure_A_prime_asymmetric.jl | Asymmetric variant | Different mu_G, mu_R |
| figure_A_double_prime_nonproduct.jl | Non-product case | 4D raw heat, slices differ |
| figure_B_pde_absorbing.jl | PDE with absorbing BC | Mass decay |
| figure_B_pde_reflecting.jl | PDE with reflecting BC | Mass conservation |
| figure_C_spectral_1d.jl | Spectral decomposition (d=1) | SVD, singular functions |
| figure_C_spectral_4d.jl | Spectral decomposition (d=4) | Higher-d singular values |
| figure_D_foodweb_static.jl | Static food web | Trophic matrix, centroids |
| figure_E_foodweb_dynamic.jl | Dynamic food web | Time evolution, filmstrip |

**Review Checklist for All**:
- [ ] Uses `theme_latexfonts()` for consistent fonts
- [ ] String concatenation (not interpolation) per project convention
- [ ] Colorbar height matches panel height
- [ ] Labels use correct notation (g, r)
- [ ] Output saved via `pkgdir(IDPG)` path resolution

---

## Tests

### Location: `test/`

Run with: `julia --project=. -e 'using Pkg; Pkg.test()'`

| File | Tests |
|------|-------|
| test_latent_space.jl | B^d_+ geometry, projections |
| test_intensity.jl | Intensity evaluation, marginal_stats |
| test_sampling.jl | PPP sampling, estimate_max_intensity |
| test_formulas.jl | E[N], E[|E|], E[L] formulas |
| test_graph_generation.jl | Node/edge-centric generation |
| test_pde.jl | Finite difference PDE solver |
| test_pdesciml_4d.jl | SciML PDE solver |

**Review Checklist**:
- [ ] All tests pass
- [ ] Coverage includes edge cases (d=1, high d, boundary points)
- [ ] Monte Carlo tests have sufficient replications

---

## Review Checklist Summary

### Library Review
- [ ] Types are parameterized by dimension d
- [ ] Callable structs: `(ρ::Intensity)(g, r)` pattern
- [ ] Pure functions where possible (no global state)
- [ ] String concatenation (not interpolation) per project convention

### Script Review
- [ ] Runs successfully: `julia --project=. script.jl`
- [ ] Outputs to correct directory
- [ ] Saves both figure and data
- [ ] Uses library functions (not reimplementing)

### Mathematical Review
- [ ] Node-centric: E[|E|] = E[N]² · p̄ = Λ² · μ̃_G · μ̃_R
- [ ] Edge-centric: E[L] = (E[N]/2) · p̄ = (Λ/2) · μ̃_G · μ̃_R
- [ ] Bound heat: h̄(g,r) = K(g,r) · ρ_G(g) · ρ_R(r)
- [ ] Product intensity: ρ(g,r) = ρ_G(g) · ρ_R(r)
- [ ] Charges: c_G = ∫ρ_G, c_R = ∫ρ_R, Λ = c_G · c_R

---

## Appendix A: Library Function Reference

Key functions extracted to the library for reuse:

### Intensity.jl

| Function | Purpose |
|----------|---------|
| `truncated_gaussian_sigma(x, mu, sigma)` | 1D truncated Gaussian (sigma parameterization) |
| `truncated_gaussian_kappa(x, mu, kappa)` | 1D truncated Gaussian (kappa parameterization) |
| `compute_1d_marginals(grid, mu_G, mu_R, sigma)` | Compute rho_G and rho_R on grid |
| `compute_bound_heat_matrix(grid, rho_G, rho_R)` | Compute bound heat h_bar(g,r) |
| `create_calibrated_product_intensity(lambda)` | Create ProductIntensity with target E[N] |

### EcologicalUtils.jl

| Function | Purpose |
|----------|---------|
| `assign_site_to_guild(site, guild_means)` | Classify site by nearest centroid |
| `assign_point_to_guild(point, guild_means)` | Simpler version for raw points |
| `build_full_guild_means(means_G, means_R)` | Combine G and R centroids |
| `compute_foodweb_matrix(sample, guild_means)` | Count guild-guild edges |
| `normalize_foodweb_matrix(M; mode)` | Normalize edge counts |
| `default_trophic_affinity()` | 5x5 trophic affinity matrix |
| `construct_food_web_centroids()` | SVD + optimization for centroids |
| `trophic_layout(n_guilds)` | Layout positions for food web plots |

### Visualization.jl

| Function | Purpose |
|----------|---------|
| `latex_figure_theme(; fontsize)` | Create LaTeX-style Makie theme |
| `with_latex_theme(f; fontsize)` | Convenience wrapper |
| `draw_Bd_plus_boundary!(ax)` | Draw B^d_+ boundary on axis |
| `Bd_plus_to_2d(point)` | Project high-d point for 2D visualization |

---

## Appendix B: Changelog

### 2025-01-09: Directory Restructuring

**Reorganized:**
- `simulations/` -> `scripts/scaling_laws/`
- `heat_map_sims/phase1/` + `heat_map_sims/phase2/` -> `scripts/heat_maps/`
- `output/simulations/` -> `output/scaling_laws/`

**Renamed scripts to match output figures:**
- `sim1_node_vs_edge.jl` -> `fig1_intensity_scaling.jl`
- `sim2_intermediate.jl` -> `fig2_entity_lifetime.jl`
- `sim6_temporal_comparison.jl` -> `fig6_temporal_evolution.jl`

**Updated path resolution:**
- All scripts now use `pkgdir(IDPG)` instead of relative `@__DIR__` paths

**Archived:**
- `examples/` -> `archive/examples/`

### 2025-01-08: Library Extraction

**Extracted to library:**
- 1D intensity utilities (truncated_gaussian, compute_bound_heat_matrix)
- Ecological utilities (assign_site_to_guild, compute_foodweb_matrix)
- Food web centroid construction (construct_food_web_centroids)
- Calibrated intensity creation (create_calibrated_product_intensity)
- Visualization theme utilities

**Deleted:**
- `simulations/sim_utils.jl` (functionality moved to library)
