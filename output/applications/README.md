# Applications

Ecological network (food web) examples in 2D and 4D latent spaces.

## Plots

### ecological_foodweb.png

**Script:** `examples/ecological_example.jl`

2D food web with bimodal intensity creating trophic levels.

### ecological_4d_comparison.png

**Script:** `examples/ecological_4d_example.jl`

Side-by-side comparison of 2D vs 4D latent spaces:
- Row 1: Food web graph structure
- Row 2: Interaction count heatmaps

**4D Configuration:**
- 5 trophic guilds: Producers, Small Herbivores, Large Herbivores, Small Predators, Apex Predators
- Hyperspherical coordinates with angles phi_small = pi/60, phi_large = 29*pi/60
- Resource radii: 0.95, 0.90, 0.85, 0.75, 0.20 (apex rarely eaten)
- Consumer radii: 0.02, 0.95, 0.90, 0.90, 0.95
- Concentration: kappa = 30, scale = 8000

### ecological_4d_variability.png

**Script:** `examples/ecological_4d_example.jl`

2x2 grid showing variability across 4 independent samples:
- Same intensity function, different realized networks
- Uses nearest-guild assignment (not k-means)

### ecological_4d_detailed.png

**Script:** `examples/ecological_4d_example.jl`

Two-panel analysis:
1. Observed interactions (average over 50 trials)
2. Expected interaction strength (g_i . r_j)

**Parameters:** kappa = 40, scale = 40000, min_edge_count = 20

### ecological_4d_filtered_web.png

**Script:** `examples/ecological_4d_example.jl`

Filtered food web showing only significant edges (avg count >= 20).
Trophic layout: Producers at bottom, Apex at top.
