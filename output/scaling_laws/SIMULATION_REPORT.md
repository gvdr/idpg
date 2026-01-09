# IDPG Simulation Report

**Date:** January 2026
**Purpose:** Empirical validation of theoretical results from the IDPG manuscript

---

## Executive Summary

Three simulations were conducted to validate the core theoretical predictions of the Intensity Dot Product Graph (IDPG) framework:

1. **Simulation 1:** Verified quadratic vs linear scaling of expected edges under node-centric vs edge-centric rules
2. **Simulation 2:** Demonstrated smooth interpolation between limiting rules via entity lifetime
3. **Simulation 6:** Confirmed that the ratio formula holds under temporal PDE evolution

All theoretical predictions were successfully validated.

---

## Theoretical Background

### The Two Limiting Realization Rules

The IDPG framework defines two limiting rules for generating graphs from a site intensity ρ on Ω = B^d_+ × B^d_+:

**Node-centric (R_∞):** Long-lived entities
- Sample N ~ Poisson(Λ) sites from PPP(ρ)
- All N² ordered pairs are edge opportunities
- Expected edges: E[|E|]_node = Λ² × (μ̃_G · μ̃_R)

**Edge-centric (R_0):** Ephemeral entities
- Sample M ~ Poisson(Λ/2) edge opportunities
- Each opportunity pairs independent source and target
- Expected edges: E[|E|]_edge = (Λ/2) × (μ̃_G · μ̃_R)

**Key ratio:** E[|E|]_node / E[|E|]_edge = 2Λ

### Intermediate Regime (R_μ)

Entities have finite lifetime τ ~ Exp(μ). Two entities interact only if their alive intervals overlap.

**Overlap probability:**
```
p_overlap(μ, W) = (2/α²)(α - 1 + e^{-α})   where α = W/μ
```

**Limiting behavior:**
- μ/W → ∞: p_overlap → 1 (node-centric)
- μ/W → 0: p_overlap → 0 (edge-centric)

---

## Simulation 1: Node-centric vs Edge-centric Comparison

### Configuration

| Parameter | Value |
|-----------|-------|
| Replications per Λ | 1000 |
| Λ values | 10, 25, 50, 100, 200 |
| Dimension d | 2 |
| Mean positions | μ_G = (0.6, 0.4), μ_R = (0.5, 0.5) |
| Concentration κ | 15.0 |

### Results

#### Expected Edges vs Total Intensity

| Λ | E[N]_emp | E[|E|]_node_emp | E[|E|]_edge_emp | Ratio_emp | Ratio_theory |
|---|----------|-----------------|-----------------|-----------|--------------|
| 10 | 9.8 | 50.2 | 2.3 | 22.1 | 20.1 |
| 25 | 24.9 | 301.3 | 5.9 | 51.2 | 50.3 |
| 50 | 49.5 | 1167.1 | 11.6 | 101.0 | 100.6 |
| 100 | 98.6 | 4595.1 | 23.2 | 197.9 | 201.2 |
| 200 | 198.0 | 18480.5 | 46.9 | 394.0 | 402.5 |

#### Log-Log Scaling Analysis

| Model | Expected Slope | Empirical Slope |
|-------|----------------|-----------------|
| Node-centric | 2.0 | **1.97** |
| Edge-centric | 1.0 | **1.01** |

### Key Findings

1. **Quadratic scaling confirmed:** Node-centric edges scale as Λ² (slope ≈ 2)
2. **Linear scaling confirmed:** Edge-centric edges scale as Λ (slope ≈ 1)
3. **Ratio formula verified:** E[|E|]_node / E[|E|]_edge ≈ 2Λ across all tested values

### Generated Figures

- `fig1_1_edges_vs_intensity.png` - Expected edges vs Λ with theoretical curves
- `fig1_2_ratio_verification.png` - Ratio verification (empirical vs 2Λ)
- `fig1_3_loglog_scaling.png` - Log-log plot with fitted slopes
- `fig1_4_sample_graphs.png` - Visual comparison of graph structures
- `fig1_5_degree_distributions.png` - Degree distributions for Λ = 50

---

## Simulation 2: Intermediate Regime (R_η)

### Notation

- η = mean entity lifetime (τ ~ Exp(η))
- α = W/η (dimensionless ratio in overlap formula)
- W = observation window

### Configuration

| Parameter | Value |
|-----------|-------|
| Replications per η | 1000 |
| Observation window W | 1.0 |
| Total intensity Λ | 50.0 |
| η values | 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, ∞ |

### Results

#### Overlap Probability Verification

| η/W | p_overlap (theory) | p_overlap (empirical) | Relative Error |
|-----|-------------------|----------------------|----------------|
| 0.01 | 0.0198 | 0.0199 | 0.5% |
| 0.02 | 0.0392 | 0.0390 | 0.5% |
| 0.05 | 0.0950 | 0.0954 | 0.4% |
| 0.10 | 0.1800 | 0.1812 | 0.7% |
| 0.20 | 0.3205 | 0.3195 | 0.3% |
| 0.50 | 0.5677 | 0.5668 | 0.2% |
| 1.00 | 0.7358 | 0.7348 | 0.1% |
| 2.00 | 0.8522 | 0.8510 | 0.1% |
| 5.00 | 0.9365 | 0.9390 | 0.3% |
| 10.0 | 0.9675 | 0.9671 | 0.04% |
| 50.0 | 0.9934 | 0.9944 | 0.1% |
| ∞ | 1.0000 | 1.0000 | 0.0% |

#### Expected Edges Across Lifetime Spectrum

| η/W | E[opportunities] | E[|E|] | Regime |
|-----|------------------|--------|--------|
| 0.01 | 48.5 | 22.8 | Near edge-centric |
| 0.10 | 438.8 | 205.6 | Transitional |
| 1.00 | 1796.9 | 840.1 | Moderate overlap |
| 10.0 | 2340.9 | 1098.2 | Near node-centric |
| ∞ | 2440.4 | 1140.9 | Node-centric limit |

### Key Findings

1. **Smooth interpolation:** p_overlap increases monotonically from 0 to 1 as η/W increases
2. **Formula verified:** The theoretical formula p_overlap = (2/α²)(α - 1 + e^{-α}) matches empirical data within 1%
3. **Limiting behavior confirmed:**
   - η/W → 0: Approaches edge-centric (few overlaps)
   - η/W → ∞: Approaches node-centric (all pairs overlap)

### Generated Figures

- `fig2_1_opportunities_vs_lifetime.png` - Edge opportunities vs η/W
- `fig2_2_edges_vs_lifetime.png` - Expected edges vs η/W with limit lines
- `fig2_3_overlap_probability.png` - Overlap probability with theoretical curve
- `fig2_4_transition_visualization.png` - Sample graphs at η/W = 0.1, 1.0, 10.0

---

## Simulation 6: Temporal Evolution Comparison

### Configuration

| Parameter | Value |
|-----------|-------|
| Replications per snapshot | 200 |
| Time range | [0, 1] |
| Snapshots | 11 (t = 0, 0.1, ..., 1.0) |
| Grid resolution | 12 (4D) |
| Time step dt | 0.002 |

### PDE Regimes

| Regime | Description | Parameters |
|--------|-------------|------------|
| Static | No evolution | D = 0 |
| Diffusion | Heat equation | D = 0.08 |
| Advection | Transport equation | v = (0.4, 0.3, -0.1, 0) |
| Pursuit-Evasion | Coupled dynamics | pursuit = 0.5, centering = 0.15 |

### Species Configuration (4D Food Web)

| Species | Role | Scale |
|---------|------|-------|
| Producer | Resource-dominated | 300 |
| Herbivore | Intermediate | 250 |
| Carnivore | Intermediate | 200 |
| Apex | Consumer-dominated | 150 |

### Results

#### Total Intensity Evolution

| Regime | Λ(t=0) | Λ(t=1) | Change |
|--------|--------|--------|--------|
| Static | 19.8 | 19.8 | 0% |
| Diffusion | 19.8 | 19.8 | 0% (mass-preserving) |
| Advection | 19.8 | 13.7 | -31% (absorbing boundary) |
| Pursuit-Evasion | 19.8 | 23.4 | +18% |

#### Ratio Verification at Each Time Point

| Regime | MAE (ratio) | Ratio tracks 2Λ(t)? |
|--------|-------------|---------------------|
| Static | 2.23 | ✅ Yes |
| Diffusion | 2.68 | ✅ Yes |
| Advection | 1.86 | ✅ Yes |
| Pursuit-Evasion | 2.10 | ✅ Yes |

*MAE = Mean Absolute Error between empirical ratio and theoretical 2Λ(t)*

#### Detailed Time Series (Advection Regime)

| Time t | Λ(t) | Ratio_theory | Ratio_emp |
|--------|------|--------------|-----------|
| 0.0 | 19.8 | 39.6 | 41.5 |
| 0.2 | 18.7 | 37.4 | 39.7 |
| 0.4 | 18.2 | 36.4 | 39.5 |
| 0.6 | 17.5 | 34.9 | 33.8 |
| 0.8 | 15.8 | 31.7 | 33.9 |
| 1.0 | 13.7 | 27.4 | 30.2 |

### Key Findings

1. **Ratio invariant under dynamics:** The ratio E[|E|]_node / E[|E|]_edge = 2Λ(t) holds at every time point
2. **Tracks changing intensity:** Even as Λ(t) changes due to absorbing boundaries or coupled dynamics, the ratio correctly tracks 2Λ(t)
3. **PDE dynamics validated:** All four regimes (static, diffusion, advection, pursuit-evasion) produce expected behavior

### Generated Figures

- `fig6_1_edges_over_time.png` - 4-panel comparison of E[|E|](t) for each regime
- `fig6_2_ratio_over_time.png` - Ratio evolution across all regimes
- `fig6_3_intensity_over_time.png` - Total intensity Λ(t) evolution
- `fig6_4_logscale_comparison.png` - Log-scale comparison of node vs edge scaling

---

## Conclusions

### Theoretical Predictions Validated

| Prediction | Status | Evidence |
|------------|--------|----------|
| E[|E|]_node = Λ² × (μ̃_G · μ̃_R) | ✅ Verified | Simulation 1: slope = 1.97 |
| E[|E|]_edge = (Λ/2) × (μ̃_G · μ̃_R) | ✅ Verified | Simulation 1: slope = 1.01 |
| Ratio = 2Λ | ✅ Verified | All simulations |
| p_overlap formula | ✅ Verified | Simulation 2: <2% error |
| Ratio invariant under PDE | ✅ Verified | Simulation 6 |

### Implementation Notes

All simulations use:
- Julia 1.10+ with the IDPG package
- Reproducible random seeds for each replication
- Monte Carlo estimation with standard errors
- JLD2 for result persistence

### File Manifest

```
output/simulations/
├── SIMULATION_REPORT.md          (this file)
├── fig1_1_edges_vs_intensity.png
├── fig1_2_ratio_verification.png
├── fig1_3_loglog_scaling.png
├── fig1_4_sample_graphs.png
├── fig1_5_degree_distributions.png
├── fig2_1_opportunities_vs_lifetime.png
├── fig2_2_edges_vs_lifetime.png
├── fig2_3_overlap_probability.png
├── fig2_4_transition_visualization.png
├── fig6_1_edges_over_time.png
├── fig6_2_ratio_over_time.png
├── fig6_3_intensity_over_time.png
├── fig6_4_logscale_comparison.png
├── sim2_results.jld2
└── sim6_results.jld2
```

### Reproducibility

To reproduce these results:

```bash
cd /path/to/idpg
julia --project=. simulations/sim1_node_vs_edge.jl
julia --project=. simulations/sim2_intermediate.jl
julia --project=. simulations/sim6_temporal_comparison.jl
```

---

## Appendix: Key Formulas

### Product Case Statistics

For product intensity ρ(g, r) = ρ_G(g) · ρ_R(r):

- Total intensity: Λ = c_G · c_R
- Normalized means: μ̃_G = μ_G / c_G, μ̃_R = μ_R / c_R
- Average connection probability: p = μ̃_G · μ̃_R

### Expected Edges

| Rule | Formula |
|------|---------|
| Node-centric | E[|E|] = Λ² · p |
| Edge-centric | E[|E|] = (Λ/2) · p |
| Intermediate | E[|E|] = Λ² · p_overlap(η, W) · p |

### Overlap Probability

```
p_overlap(η, W) = (2/α²)(α - 1 + e^{-α})

where α = W/η and η = mean entity lifetime
```

Taylor expansions:
- α → 0 (η → ∞): p_overlap → 1 - α/3 + O(α²) → 1 (node-centric)
- α → ∞ (η → 0): p_overlap → 2/α + O(1/α²) → 0 (edge-centric)
