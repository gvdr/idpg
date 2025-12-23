# Temporal Food Web

Temporal evolution of food web structure under different PDE regimes.

## Plots

### regime_comparison_heatmaps.png

**Script:** `examples/temporal_foodweb.jl`

Food web structure at t=0 and t=1 under four regimes:
1. **Static:** No evolution (baseline)
2. **Diffusion:** Niches spread (D = 0.08)
3. **Advection:** Niches shift uniformly (v = [0.4, 0.3, -0.1, 0.0])
4. **Pursuit-Evasion:** Coupled dynamics

### interaction_counts.png

**Script:** `examples/temporal_foodweb.jl`

Total interaction count over time for each regime:
- Static: constant (~72 interactions)
- Diffusion: decreases as niches spread
- Advection: varies with shift direction
- Pursuit-Evasion: dramatic decrease (~25 interactions at t=1)

### final_foodwebs.png

**Script:** `examples/temporal_foodweb.jl`

Food web graphs at t=1 for each regime:
- Trophic layout (Producers at bottom, Apex at top)
- Edge width proportional to interaction frequency

### mean_trajectories.png

**Script:** `examples/temporal_foodweb.jl`

Intensity mean position trajectories (dims 1-2 projection):
- Circle = start, Star = end
- Shows how resource (rho_G) and consumer (rho_R) means move

### pursuit_evasion_filmstrip.png

**Script:** `examples/temporal_foodweb.jl`

Time-lapse of food web under pursuit-evasion (5 snapshots):
- Shows progressive weakening of trophic links
- Resources flee, consumers chase

**Parameters:**
- Grid resolution: 12 (4D)
- dt = 0.002, t_final = 1.0
- Pursuit speed: 0.5
- Centering strength: 0.15
- Center position: (0.25, 0.25, 0.25, 0.25)
