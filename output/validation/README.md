# Validation

Formula validation plots comparing theoretical predictions to empirical sampling.

## Plots

### product_case_validation.png

**Script:** `examples/product_case.jl`

Three-panel log-log validation:
1. **Node Count:** E[N] theory vs empirical
2. **Node-Centric Edges:** E[|E|] = E[N]^2 * (mu_tilde_G . mu_tilde_R)
3. **Edge-Centric Edges:** E[L] = E[N] * (mu_tilde_G . mu_tilde_R)

Red dashed line = perfect agreement (y = x).

**What to look for:**
- Points should fall on or near the diagonal
- The log-log scale spans ~3 orders of magnitude
- Systematic deviations would indicate formula errors

**Key formulas validated:**
- E[N] = c_G * c_R
- E[|E|] = E[N]^2 * (mu_tilde_G . mu_tilde_R)
- E[L] = E[N] * (mu_tilde_G . mu_tilde_R)
