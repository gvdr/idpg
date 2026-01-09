# Mesoscale

Network metrics and IDPG vs Erdos-Renyi comparison.

## Plots

### mesoscale_metrics.png

**Script:** `examples/mesoscale_metrics.jl`

Nine-panel analysis of network structure:
- Reciprocity distribution (expected: avg_conn_prob)
- Degree independence test (expected: ~0)
- Clustering distribution
- Giant SCC fraction
- Degree assortativity
- PageRank inequality (Gini)
- In-degree and out-degree distributions
- In vs Out degree scatter

### mesoscale_comparison.png

**Script:** `examples/mesoscale_metrics.jl`

Reciprocity across three connectivity regimes:
1. Aligned (high connectivity): G and R means aligned
2. Orthogonal (medium): G and R means nearly orthogonal
3. Very orthogonal (low): G and R means maximally orthogonal

### idpg_vs_er_comparison.png

**Script:** `examples/idpg_vs_erdos_renyi.jl`

Grid of histograms comparing IDPG (blue) vs matched ER (orange):
- Clustering, Degree Var/Mean, In-degree Std, Assortativity
- Three rows = three intensity configurations

**Key finding:** Degree Var/Mean shows complete separation (IDPG >> ER).

**Configurations tested:**
1. Concentrated (orthogonal): mu_G = (0.65, 0.15), mu_R = (0.15, 0.65), kappa = 40, c = 80
2. Spread (aligned): mu = (0.5, 0.5), kappa = 5, c = 40
3. Bimodal G: mu_G = {(0.8, 0.1), (0.1, 0.8)}, mu_R = (0.5, 0.5), kappa = 30, c = 60

### idpg_vs_er_degrees.png

**Script:** `examples/idpg_vs_erdos_renyi.jl`

Single-realization degree distribution comparison:
- Out-degree distribution
- In-degree distribution
- In vs Out scatter
