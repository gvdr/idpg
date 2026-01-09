# Graphs

Directed graph structure visualizations.

## Plots

### graph_visualization.png

**Script:** `examples/graph_visualization.jl`

Six-panel visualization of directed graph structure:
1. Node positions in G-space (source propensity)
2. Node positions in R-space (target propensity)
3. Directed graph with spring layout
4. Directed graph embedded in G-space
5. Why edges form: P(i->j) = g_i . r_j
6. Edge-centric sample

### graph_directed.png

**Script:** `examples/graph_visualization.jl`

Single directed graph with nodes colored by ||g|| (source strength):
- Brighter yellow = higher sending propensity
- Spring layout reveals community structure
- Arrows show edge directions

### graph_simple.png

Simple graph visualization for quick inspection.
