# Archive

This folder contains material that has been archived during codebase cleanup. The archived content is preserved for reference but is not directly used in the manuscript.

## Directory Structure

```
archive/
├── examples/           # Educational example scripts (10 files)
├── example_outputs/    # Outputs from example scripts
│   ├── basics/
│   ├── dynamics/
│   ├── graphs/
│   ├── mesoscale/
│   ├── applications/
│   ├── temporal_foodweb/
│   └── validation/
├── old_outputs/        # Superseded figures
└── old_logs/           # Old log files
```

## Archived Content

### examples/ (10 scripts, ~3,600 lines)

Educational examples demonstrating IDPG library usage:

| File | Purpose |
|------|---------|
| `basic_idpg.jl` | Minimal working example |
| `product_case.jl` | Formula verification |
| `edge_centric_comparison.jl` | Two edge-centric interpretations |
| `ecological_example.jl` | 2D food web |
| `ecological_4d_node_example.jl` | 4D guild structure (838 lines) |
| `temporal_foodweb.jl` | PDE dynamics + ecology (786 lines) |
| `idpg_vs_erdos_renyi.jl` | IDPG vs ER random graphs |
| `mesoscale_metrics.jl` | Network statistics |
| `graph_visualization.jl` | Graph plotting |
| `diffusion_example.jl` | Simple PDE evolution |

### example_outputs/

Outputs organized by topic:

- **basics/** - Basic IDPG demonstration
- **dynamics/** - Diffusion evolution
- **graphs/** - Graph structure visualizations
- **mesoscale/** - Network metrics (clustering, reciprocity, assortativity)
- **applications/** - Ecological 4D food web examples
- **temporal_foodweb/** - Temporal evolution under different PDE regimes
- **validation/** - Product case and edge-centric formula validation

### old_outputs/

- `ecological_foodweb.png` - Old 2D ecological figure (superseded by 4D examples)
- `ecological_mc_comparison.png` - Old Monte Carlo comparison
- `diffusion_evolution.mp4` - Orphaned video from early experiments

### old_logs/

- `test_output.log` - Captured test output from December 2024

## Cleanup Log

**2025-01-08: Major cleanup**

Archived to streamline the repository for manuscript preparation:
1. Moved all 10 example scripts to `archive/examples/`
2. Moved 7 example output directories to `archive/example_outputs/`
3. Previously archived old figures and logs

The active codebase now contains only:
- `src/` - Core IDPG library (8 modules)
- `simulations/` - Validation simulations for manuscript (sim1, sim2, sim6)
- `heat_map_sims/` - Figure generation scripts (A-E)
- `output/simulations/` and `output/heat_maps/` - Manuscript figures
- `test/` - Unit tests
- `docs/` - Documentation

## Restoring Examples

To restore examples to their original location:

```bash
mv archive/examples/* examples/
mv archive/example_outputs/* output/
```

## See Also

- `output/README.md` - Current manuscript outputs
- `docs/CODE_REVIEW_GUIDE.md` - Full documentation including archived content
