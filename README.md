# IDPG: Intensity Dot Product Graphs

Julia implementation of the Intensity Dot Product Graph (IDPG) framework, accompanying the manuscript:

> **Intensity Dot Product Graphs**
> Giulio Valentino Dalla Riva and Matteo Dalla Riva

IDPG extends Random Dot Product Graphs (RDPG) by treating nodes as samples from a Poisson point process with intensity defined on a continuous latent space. This allows both the number and identity of nodes to emerge from a stochastic process, and enables the study of temporal graph dynamics via partial differential equations.

## Requirements

- Julia 1.10 or later
- Dependencies are specified in `Project.toml`

## Installation

```bash
git clone <repository-url>
cd idpg
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Project Structure

```
idpg/
├── src/                      # Core library
│   ├── IDPG.jl              # Main module
│   └── modules/             # Implementation files
│       ├── LatentSpace.jl   # B^d_+ geometry
│       ├── Intensity.jl     # Intensity functions
│       ├── Sampling.jl      # Poisson point process sampling
│       ├── GraphGeneration.jl
│       ├── PDEEvolution.jl  # Finite difference PDE solver
│       ├── PDESciML.jl      # SciML-based PDE solver
│       ├── Visualization.jl
│       └── EcologicalUtils.jl
├── scripts/
│   ├── scaling_laws/        # Validation simulations (Figures 1, 2, 6)
│   └── heat_maps/           # Paper figure generation (Figures A-E)
├── test/                    # Unit tests
├── docs/                    # Documentation and manuscript
└── output/                  # Generated figures and data
```

## Running Scripts

Scaling law validation (node-centric vs edge-centric comparison):
```bash
julia --project=. scripts/scaling_laws/fig1_intensity_scaling.jl
```

Heat map figures:
```bash
julia --project=. scripts/heat_maps/construct_centroids.jl
julia --project=. scripts/heat_maps/figure_A_anatomy.jl
```

Many figure scripts support a `--plot` flag to regenerate figures from previously saved data:
```bash
julia --project=. scripts/heat_maps/figure_A_anatomy.jl --plot
```

## Key Concepts

- **Heat map**: In classical RDPG, interaction structure is captured by a discrete probability matrix between given nodes. In IDPG, nodes emerge from a continuous intensity landscape, and the natural analog is the heat map: a measure-theoretic object that captures interaction structure between regions of site space rather than between discrete nodes.

- **Node-centric realization**: Sites sampled from intensity become nodes; edges form between all pairs with probability g . r. Expected edges scale as E[N]^2.

- **Edge-centric realization**: Each sampled site represents an edge opportunity. Expected edges scale as E[N]/2.

- **Product intensity**: Factorized form rho(g,r) = rho_G(g) * rho_R(r) where source and target distributions are independent.

- **Mixture of products**: Sum of species-specific products for ecological modeling with guild structure.

- **Temporal dynamics**: The intensity can evolve according to PDEs (diffusion, advection), inducing time-varying graph distributions.

## Documentation

- `CODE_REVIEW_GUIDE.md` - Detailed code review guide
