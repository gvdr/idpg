# Sampling Analysis: ecological_4d_example vs temporal_foodweb

## Key Findings

### 1. `generate_edge_centric_full` (used in ecological_4d_example)

Location: `src/modules/GraphGeneration.jl:160-178`

```julia
function generate_edge_centric_full(sites::Vector{InteractionSite{d}}; rng) where d
    n = length(sites)
    for i in 1:n
        for j in 1:n
            p = connection_probability(sites[i].g, sites[j].r)
            if rand(rng) < p
                push!(source_sites, sites[i])  # Full (g, r) of source
                push!(target_sites, sites[j])  # Full (g, r) of target
            end
        end
    end
end
```

**Algorithm:**
1. Takes N pre-sampled sites as input (each site has coupled (g, r))
2. Iterates over ALL pairs (i, j) of sites - O(N²) pairs
3. Each pair forms an edge with probability g_i · r_j
4. Preserves FULL (g, r) for BOTH source site i AND target site j

**Key insight:** This is actually the **NODE-CENTRIC model** (sites become nodes, edges form between all pairs), just returning edges with full site info rather than a graph structure. The docstring even says: "This is equivalent to the node-centric model but returns edges with full site info."

### 2. ecological_4d_example Sampling Pipeline

```julia
# Step 1: Sample sites from MixtureOfProductIntensities
labeled_sites = sample_ppp_mixture(ρ; n_samples=n_mc, rng=...)
sites = [site for (_, site) in labeled_sites]

# Step 2: Generate edges using node-centric-style pairing
interactions = generate_edge_centric_full(sites; rng=...)
```

**What `sample_ppp_mixture` does:**
1. Compute γ_m = c_{G,m} · c_{R,m} for each species
2. Sample N ~ Poisson(C) where C = Σ γ_m
3. For each of N sites:
   - Sample species m with probability γ_m / C
   - Sample g from ρ_{G,m}, r from ρ_{R,m} (COUPLED within species)

**Result:** N sites, each with species-coupled (g, r), then O(N²) pairs considered.

### 3. temporal_foodweb Current Implementation (sample_from_species_grids)

```julia
function sample_from_species_grids(...)
    N = rand(rng, Poisson(C))  # C = Σ γ_m

    for _ in 1:N
        # Sample SOURCE species and position
        m_src = sample(rng, 1:n_species, Weights(probs))
        source_site = InteractionSite(g from ρ_{G,m_src}, r from ρ_{R,m_src})

        # Sample TARGET species and position (INDEPENDENTLY!)
        m_tgt = sample(rng, 1:n_species, Weights(probs))
        target_site = InteractionSite(g from ρ_{G,m_tgt}, r from ρ_{R,m_tgt})

        # Accept with probability g_source · r_target
        if rand(rng) < connection_probability(source_site.g, target_site.r)
            push!(source_sites, source_site)
            push!(target_sites, target_site)
        end
    end
end
```

**What this does:**
1. Sample N ~ Poisson(C) edge OPPORTUNITIES directly
2. For EACH opportunity, independently sample source AND target
3. Accept with prob g_source · r_target

**This is TRUE EDGE-CENTRIC sampling** where each opportunity is an independent edge draw.

### 4. Critical Differences

| Aspect | ecological_4d | temporal_foodweb |
|--------|--------------|------------------|
| Model | Node-centric (sites→nodes, all pairs) | True edge-centric (each draw = edge) |
| Complexity | O(N²) pairs from N sites | O(N) edge opportunities |
| Expected edges | N² × E[g·r] | N × E[g·r] |
| Source-target coupling | Same pool of sites | Independent draws |
| Full (g,r) preserved | Yes for both source & target | Yes for both source & target |

### 5. Is ecological_4d_example "Edge-Centric"?

**No, not really.** Despite the function name `generate_edge_centric_full`, the ecological_4d_example uses the NODE-CENTRIC interpretation:
- Sites are sampled first (become "nodes")
- ALL pairs of sites are considered for edges
- This is why we see N² potential interactions from N sites

The "full" in the name refers to preserving full (g, r) coordinates for clustering, NOT to an edge-centric sampling model.

### 6. Verification Needed

To confirm that ecological_4d_example correctly implements the mathematical model:

**Expected behavior for guild i → guild j interactions:**
- Number of species-i sites: n_i
- Number of species-j sites: n_j
- Pairs considered: n_i × n_j
- Connection prob per pair: g_i · r_j (based on guild means)
- Expected edges: n_i × n_j × (g_i · r_j)

This should produce the asymmetric trophic structure where:
- Producers (high g in dim 1, r in dim 4) are eaten by herbivores (r in dim 1)
- Herbivores are eaten by predators (r in dim 2)
- etc.

### 7. What Should temporal_foodweb Do?

To match ecological_4d_example, temporal_foodweb should:
1. Sample N sites from per-species grids (with coupled g,r per species)
2. Call `generate_edge_centric_full(sites)` to create edges

Current implementation samples edge opportunities directly, which is a different (pure edge-centric) model.

### 8. Questions for Clarification

1. Is the node-centric pairing model (O(N²) pairs) the intended behavior?
2. If so, temporal needs restructuring to sample sites first, then pair them
3. The heatmap asymmetry difference might stem from this model mismatch
