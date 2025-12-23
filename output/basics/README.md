# Basics

Basic demonstration of IDPG sampling and graph generation.

## Plots

### basic_idpg.png

**Script:** `examples/basic_idpg.jl`

Three-panel visualization showing:
1. Intensity function rho_G (source) on B^2_+
2. Intensity function rho_R (target) on B^2_+
3. Sampled interaction sites (green coordinates)

**Configuration:**
- 2D latent space (d=2)
- 2-component mixture for rho_G: means at (0.8, 0.2) and (0.2, 0.8), concentrations 10.0, weights (0.6, 0.4)
- 2-component mixture for rho_R: means at (0.3, 0.7) and (0.5, 0.5), concentrations (10.0, 5.0), weights (0.5, 0.5)
- Total intensity: c_G = c_R = 50.0

**Expected results:**
- E[N] = 2500 interaction sites
- Average connection probability ~ 0.37
