# Phase 2 Heat Map Simulations: Report

This report summarizes the figures produced in Phase 2 of the IDPG heat map simulations.

---

## Overview

Phase 2 explores the structure of heat maps under various intensity models, PDE dynamics, and ecological applications. The figures demonstrate:

1. **Heat anatomy** under product vs non-product intensity
2. **PDE evolution** of intensities with different boundary conditions
3. **Spectral decomposition** of heat maps
4. **Food web generation** from IDPG (static and dynamic)

---

## Figure A: Heat Map Anatomy (Product Intensity)

**File:** `figure_A_anatomy.png`

**Layout:** 1×3 (single row)
- Panel (a): Marginal intensities ρ_G(g) and ρ_R(r)
- Panel (b): Bound heat h̄(g,r) heatmap
- Panel (c): Kernel K(g,r) = g·r heatmap

**Purpose:** Demonstrate the basic anatomy of heat maps under product intensity (d=1).

### Key Results

| Quantity | Value | Description |
|----------|-------|-------------|
| g₀, r₀ | 0.3, 0.7 | Intensity centers |
| κ | 20.0 | Concentration parameter |
| c_G, c_R | 1.004, 1.004 | Marginal charges |
| Λ = c_G·c_R | 1.008 | Total intensity product |
| μ̃_G, μ̃_R | 0.338, 0.662 | Normalized centroids |
| p̄ = μ̃_G·μ̃_R | 0.224 | Mean connection probability |

### Heat Scaling Relations (Product Intensity)

Under product intensity ρ(g,r) = ρ_G(g)·ρ_R(r):

- **Bound heat:** h̄(g,r) = K(g,r)·ρ_G(g)·ρ_R(r)
- **Raw heat:** H = Λ·H̄ (scales with total intensity)
- **Edge heat:** H⃗ = H̄/2 (half of bound heat)

**Key insight:** All three heats have the **same shape** under product intensity, differing only by scalar factors. The 2D bound heat h̄(g,r) is a complete summary.

---

## Figure A'' (A-double-prime): Non-Product Intensity

**File:** `figure_A_double_prime_nonproduct.png`

**Layout:** 2×3
- Panel (a): Site intensity ρ(g,r) showing correlation structure
- Panel (d): Marginals ρ^(g) and ρ^(r) (identical by symmetry)
- Panels (b,c,e,f): Four slices of 4D raw heat at different (r_s, g_t) values

**Purpose:** Demonstrate why bound heat h̄(g,r) only works under product intensity.

### Setup

Site intensity consists of two off-diagonal Gaussian blobs:
```
ρ(g,r) = N((0.3, 0.7), σ²I) + N((0.7, 0.3), σ²I)
σ = 0.12
```

This is **NOT** a product: ρ(g,r) ≠ ρ^(g)(g)·ρ^(r)(r)

### Key Results

The 4D raw heat h(g_s, r_s, g_t, r_t) = K(g_s, r_t)·ρ(g_s, r_s)·ρ(g_t, r_t) depends on all four coordinates.

| Slice (r_s, g_t) | Total Heat | Shape |
|------------------|------------|-------|
| (0.3, 0.3) | 0.044 | Concentrated at low g_s |
| (0.7, 0.3) | 0.019 | Different pattern |
| (0.3, 0.7) | 0.019 | Different pattern |
| (0.7, 0.7) | 0.008 | Concentrated at high g_s |

**Key insight:** Different slices have **different shapes**. The "inactive" coordinates (r_s, g_t) don't factor out. Under product intensity, all slices would have the same shape (differing only by a scalar), enabling the 2D bound heat summary. With non-product intensity, there is no 2D summary — the full 4D structure matters.

---

## Figure B: PDE Evolution (Reflecting Boundary)

**File:** `figure_B_pde_reflecting.png`

**Purpose:** Show how intensity distributions evolve under diffusion, advection, and pursuit-evasion dynamics with reflecting (Neumann) boundary conditions.

**Layout:** Time series of bound heat snapshots + centroid trajectories

### Dynamics

Three evolution regimes:
1. **Diffusion:** Standard heat equation, intensities spread out
2. **Advection:** Directed transport, intensities shift
3. **Pursuit-Evasion:** Predator-prey dynamics with coupled motion

### Key Features

- Reflecting BC preserves total mass
- Centroid trajectories show μ_G(t) and μ_R(t) over time
- Increased dynamics (T_FINAL=2.0) to show visible evolution

---

## Figure B' (B-prime): PDE Evolution (Absorbing Boundary)

**File:** `figure_B_pde_absorbing.png`

**Purpose:** Show wave-like transients when mass can leave the system through absorbing (Dirichlet) boundary conditions.

**Layout:** 5 time snapshots + mass decay plot

### Parameters
- μ_G = 0.3, μ_R = 0.7, σ = 0.12
- D = 0.05 (diffusion coefficient)
- T_FINAL = 1.0

### Key Results

| Time | c_G | c_R | Notes |
|------|-----|-----|-------|
| 0.0 | 1.0 | 1.0 | Initial |
| 1.0 | ~0.7 | ~0.7 | Mass lost through boundaries |

**Key insight:** Absorbing boundaries cause mass decay over time, visible in both the heat map evolution and the total heat integral plot.

---

## Figure C: Spectral Decomposition

**Files:**
- `figure_C_spectral_1d.png` (d=1)
- `figure_C_spectral_4d.png` (d=4)

**Purpose:** Show spectral decomposition of heat maps using eigenvalue analysis.

### Content

- Guild affinity matrices
- Spectral components
- Reconstruction quality

---

## Figure D: Food Web (Static)

**File:** `figure_D_foodweb_static.png`

**Layout:** 2×2
- Panel (a): Target affinity K*_ij
- Panel (b): Achieved affinity K̂_ij = μ_i^G · μ_j^R
- Panel (c): Expected edges H_ij = Λ²·π_i·π_j·K̂_ij
- Panel (d): Realized food web graph

**Purpose:** Show how IDPG generates realistic trophic networks using guild centroids from Phase 1.

### Parameters
- Λ = 100 (total intensity)
- κ = 30 (concentration)
- 5 guilds: Producer, Small Herbivore, Large Herbivore, Small Predator, Apex

### Key Results

The framework successfully recovers the target affinity structure through optimized guild centroids, producing realistic food web topology with appropriate trophic levels.

---

## Figure E: Food Web (Dynamic)

**File:** `figure_E_foodweb_dynamic.png`

**Purpose:** Show temporal evolution of food web structure as intensities change over time.

**Content:** Time series of food web snapshots showing how guild interactions evolve under PDE dynamics.

---

## Summary of Key Concepts

### Product Intensity (Figure A)

When ρ(g,r) = ρ_G(g)·ρ_R(r):
- Heat maps can be summarized by 2D bound heat h̄(g,r)
- Raw, bound, and edge heats have the same shape
- Scaling relations: H = Λ·H̄, H⃗ = H̄/2

### Non-Product Intensity (Figure A'')

When ρ(g,r) ≠ ρ_G(g)·ρ_R(r):
- No 2D summary exists
- Full 4D heat structure must be considered
- Different slices have genuinely different shapes

### Mixture of Products

For ρ(g,r) = Σ_m ρ_{G,m}(g)·ρ_{R,m}(r):
- Joint intensity is sum of products (not product of sums)
- Each species contributes a product term
- Cross-terms in marginal products don't exist in the mixture

---

## File Inventory

| Figure | File | Size |
|--------|------|------|
| A | figure_A_anatomy.png | 185 KB |
| A'' | figure_A_double_prime_nonproduct.png | 293 KB |
| B | figure_B_pde_reflecting.png | 333 KB |
| B' | figure_B_pde_absorbing.png | 263 KB |
| C (1D) | figure_C_spectral_1d.png | 252 KB |
| C (4D) | figure_C_spectral_4d.png | 266 KB |
| D | figure_D_foodweb_static.png | 365 KB |
| E | figure_E_foodweb_dynamic.png | 396 KB |

**Caption files:**
- figure_A_caption.txt
- figure_A_double_prime_caption.txt

---

## Technical Notes

### LaTeX Notation

All figures use LaTeXStrings.jl for proper mathematical typesetting:
- Bound heat: h̄ (h with overline)
- Marginal charges: c_G, c_R
- Centroids: μ̃_G, μ̃_R (normalized)
- Kernel: K(g,r) = g·r

### Colormap Conventions

- Bound heat: viridis
- Kernel: YlOrRd
- Edge heat: inferno
- Ratio/comparison: RdBu

### Dropped Figures

**Figure A' (Asymmetric Source-Target)** was developed but ultimately dropped. It explored the case where source and target are sampled from different weighted distributions (ρ_S vs ρ_T). The key finding was that for mixture-of-products intensity, the node-centric heat (sum of products) and edge-centric heat (product of weighted sums) have genuinely different spatial distributions.
