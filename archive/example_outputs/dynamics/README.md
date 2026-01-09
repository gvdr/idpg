# Dynamics

PDE evolution of intensity functions on B^d_+.

## Plots

### diffusion_evolution.png

**Script:** `examples/diffusion_example.jl`

Multi-panel time evolution showing diffusion on B^2_+:
- Initial: Gaussian-like concentration near (0.8, 0.2) with kappa = 15, c = 50
- Evolution: Heat equation with absorbing boundary

**PDE:**
```
d(rho)/dt = D * nabla^2(rho)
```

**Parameters:**
- Diffusion coefficient: D = 0.1
- Time step: dt = 0.0001
- Final time: t_final = 2.0
- Sample interval: 0.1

**What to look for:**
- Mass spreads from initial concentration
- Total intensity decreases due to absorbing boundary at ||x|| = 1
- Maximum intensity decreases as distribution flattens
