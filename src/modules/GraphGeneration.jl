# Graph generation from sampled interaction sites

"""
    generate_node_centric(sites::Vector{InteractionSite{d}}; rng=Random.default_rng())

Generate a graph under the node-centric (long-lived entities) interpretation.

Given sampled interaction sites, create a directed graph where each site
becomes a node, and edge i→j exists with probability g_i · r_j.

# Returns
- `graph::SimpleDiGraph`: The generated graph
- `sites::Vector{InteractionSite{d}}`: The input sites (for reference)
"""
function generate_node_centric(sites::Vector{InteractionSite{d}};
                               rng::AbstractRNG=Random.default_rng()) where d
    n = length(sites)
    graph = SimpleDiGraph(n)

    for i in 1:n
        for j in 1:n
            g_i = sites[i].g
            r_j = sites[j].r
            p = connection_probability(g_i, r_j)

            if rand(rng) < p
                add_edge!(graph, i, j)
            end
        end
    end

    return graph, sites
end

"""
    EdgeCentricSample{d}

Edge-centric interpretation: each sampled point IS an edge.

In the edge-centric view, each sampled interaction site (g, r) represents
a directed edge from a source "at g" to a target "at r". These are ephemeral
entities that exist only for a single interaction.

# Fields
- `sources::Vector{LatentPoint{d}}`: Green coordinates (source positions)
- `targets::Vector{LatentPoint{d}}`: Red coordinates (target positions)

Note: This structure only stores the coordinates used for connection probability.
For full site information (both g and r for source and target), use FullEdgeCentricSample.
"""
struct EdgeCentricSample{d}
    sources::Vector{LatentPoint{d}}
    targets::Vector{LatentPoint{d}}
end

"""
    FullEdgeCentricSample{d}

Edge-centric sample with full site information preserved.

Each edge is formed between a source site (g_s, r_s) and a target site (g_t, r_t).
The connection probability is g_s · r_t, but ALL coordinates are preserved:
- r_s: the source's receiving propensity (available for clustering)
- g_t: the target's sending propensity (available for clustering)

This allows clustering entities using their full (g, r) signature rather than
just the coordinate that contributed to the connection probability.

# Fields
- `source_sites::Vector{InteractionSite{d}}`: Full (g, r) for each source entity
- `target_sites::Vector{InteractionSite{d}}`: Full (g, r) for each target entity
"""
struct FullEdgeCentricSample{d}
    source_sites::Vector{InteractionSite{d}}
    target_sites::Vector{InteractionSite{d}}
end

function Base.length(sample::FullEdgeCentricSample)
    return length(sample.source_sites)
end

function Base.show(io::IO, sample::FullEdgeCentricSample{d}) where d
    print(io, "FullEdgeCentricSample{" * string(d) * "} with " * string(length(sample)) * " edges (full site info)")
end

# Convenience accessors for FullEdgeCentricSample
"""Get source g coordinates (the ones used for connection probability)"""
source_g(sample::FullEdgeCentricSample) = [s.g for s in sample.source_sites]

"""Get source r coordinates (not used for connection, but available for clustering)"""
source_r(sample::FullEdgeCentricSample) = [s.r for s in sample.source_sites]

"""Get target g coordinates (not used for connection, but available for clustering)"""
target_g(sample::FullEdgeCentricSample) = [t.g for t in sample.target_sites]

"""Get target r coordinates (the ones used for connection probability)"""
target_r(sample::FullEdgeCentricSample) = [t.r for t in sample.target_sites]

"""Convert to legacy EdgeCentricSample (loses r_source and g_target)"""
function to_edge_centric(sample::FullEdgeCentricSample{d}) where d
    return EdgeCentricSample{d}(source_g(sample), target_r(sample))
end

function Base.length(sample::EdgeCentricSample)
    return length(sample.sources)
end

function Base.show(io::IO, sample::EdgeCentricSample{d}) where d
    print(io, "EdgeCentricSample{" * string(d) * "} with " * string(length(sample)) * " edges")
end

"""
    generate_edge_centric(sites::Vector{InteractionSite{d}}; rng=Random.default_rng()) -> EdgeCentricSample{d}

Generate edges under the edge-centric (ephemeral entities) interpretation.

Each sampled site (g, r) becomes an edge with probability g · r.
This implements the two-stage interpretation: first sites are sampled from ρ,
then each site converts to an actual edge with probability given by the dot product.

# Returns
EdgeCentricSample containing the accepted edges.
"""
function generate_edge_centric(sites::Vector{InteractionSite{d}};
                               rng::AbstractRNG=Random.default_rng()) where d
    sources = LatentPoint{d}[]
    targets = LatentPoint{d}[]

    for site in sites
        p = connection_probability(site.g, site.r)
        if rand(rng) < p
            push!(sources, site.g)
            push!(targets, site.r)
        end
    end

    return EdgeCentricSample{d}(sources, targets)
end

"""
    generate_edge_centric_full(sites::Vector{InteractionSite{d}}; rng=Random.default_rng()) -> FullEdgeCentricSample{d}

Generate edges with full site information preserved.

For each pair of sites (i, j), creates an edge from i to j with probability g_i · r_j.
Unlike `generate_edge_centric`, this preserves ALL coordinates:
- Source site: both g_i (used for connection) and r_i (available for clustering)
- Target site: both r_j (used for connection) and g_j (available for clustering)

This is equivalent to the node-centric model but returns edges with full site info
rather than a graph structure.

# Arguments
- `sites`: Vector of InteractionSite samples from PPP
- `rng`: Random number generator

# Returns
FullEdgeCentricSample with complete (g, r) for both source and target of each edge.
"""
function generate_edge_centric_full(sites::Vector{InteractionSite{d}};
                                     rng::AbstractRNG=Random.default_rng()) where d
    source_sites = InteractionSite{d}[]
    target_sites = InteractionSite{d}[]

    n = length(sites)
    for i in 1:n
        for j in 1:n
            # Connection probability uses g_i (source's sending) and r_j (target's receiving)
            p = connection_probability(sites[i].g, sites[j].r)
            if rand(rng) < p
                push!(source_sites, sites[i])  # Full (g, r) of source
                push!(target_sites, sites[j])  # Full (g, r) of target
            end
        end
    end

    return FullEdgeCentricSample{d}(source_sites, target_sites)
end

"""
    discretize_edge_centric(sample::EdgeCentricSample{d}, n_clusters::Int;
                            max_iter=100, rng=Random.default_rng()) -> SimpleDiGraph

Discretize an edge-centric sample into a graph by clustering source and target
positions and counting edges between clusters.

This recovers a traditional graph structure from the continuous edge-centric
representation by aggregating edges whose endpoints fall into the same clusters.

# Arguments
- `sample`: Edge-centric sample to discretize
- `n_clusters`: Number of clusters for source and target spaces
- `max_iter`: Maximum iterations for k-means
- `rng`: Random number generator

# Returns
- `graph::SimpleDiGraph`: Weighted directed graph between clusters
- `source_assignments::Vector{Int}`: Cluster assignments for source points
- `target_assignments::Vector{Int}`: Cluster assignments for target points
"""
function discretize_edge_centric(sample::EdgeCentricSample{d}, n_clusters::Int;
                                 max_iter::Int=100,
                                 rng::AbstractRNG=Random.default_rng()) where d
    n_edges = length(sample)

    if n_edges == 0
        return SimpleDiGraph(n_clusters), Int[], Int[]
    end

    # Convert to matrices for clustering (Clustering.jl expects d×n)
    source_matrix = hcat([Vector(s) for s in sample.sources]...)
    target_matrix = hcat([Vector(t) for t in sample.targets]...)

    # Use Clustering.jl's k-means
    source_result = kmeans(source_matrix, n_clusters; maxiter=max_iter)
    target_result = kmeans(target_matrix, n_clusters; maxiter=max_iter)

    source_assignments = source_result.assignments
    target_assignments = target_result.assignments

    # Build graph with edge counts
    graph = SimpleDiGraph(n_clusters)
    edge_counts = zeros(Int, n_clusters, n_clusters)

    for i in 1:n_edges
        s = source_assignments[i]
        t = target_assignments[i]
        edge_counts[s, t] += 1
    end

    # Add edges where count > 0
    for s in 1:n_clusters
        for t in 1:n_clusters
            if edge_counts[s, t] > 0
                add_edge!(graph, s, t)
            end
        end
    end

    return graph, source_assignments, target_assignments
end

# Using Clustering.jl for k-means and DBSCAN

"""
    discretize_edge_centric_joint(sample::EdgeCentricSample{d};
                                   eps=0.15, min_samples=3,
                                   rng=Random.default_rng())

Discretize an edge-centric sample using joint clustering of (g, r) pairs.
Uses DBSCAN density clustering to automatically determine number of clusters.

Each interaction is represented as a 2d-dimensional point (g, r), and interactions
with similar (g, r) pairs are grouped together.

# Arguments
- `sample`: Edge-centric sample to discretize
- `eps`: Maximum distance between points in the same cluster
- `min_samples`: Minimum points to form a dense region

# Returns
- `graph::SimpleDiGraph`: Directed graph between interaction types
- `edge_weights::Matrix{Int}`: Count of interactions between clusters
- `assignments::Vector{Int}`: Cluster assignment for each interaction (0 = noise)
- `n_clusters::Int`: Number of clusters found
"""
function discretize_edge_centric_joint(sample::EdgeCentricSample{d};
                                        eps::Float64=0.15,
                                        min_samples::Int=3,
                                        rng::AbstractRNG=Random.default_rng()) where d
    n_edges = length(sample)

    if n_edges == 0
        return SimpleDiGraph(0), zeros(Int, 0, 0), Int[], 0
    end

    # Concatenate (g, r) into 2d-dimensional points
    joint_points = Matrix{Float64}(undef, 2d, n_edges)
    for i in 1:n_edges
        joint_points[1:d, i] = sample.sources[i]
        joint_points[d+1:2d, i] = sample.targets[i]
    end

    # Use Clustering.jl's DBSCAN
    result = dbscan(joint_points, eps, min_neighbors=min_samples)

    # Convert to assignments vector (0 = noise)
    assignments = zeros(Int, n_edges)
    for (cluster_id, cluster) in enumerate(result.clusters)
        for idx in cluster.core_indices
            assignments[idx] = cluster_id
        end
        for idx in cluster.boundary_indices
            assignments[idx] = cluster_id
        end
    end

    n_clusters = length(result.clusters)

    if n_clusters == 0
        # All points are noise - fall back to treating each as its own cluster
        return SimpleDiGraph(1), ones(Int, 1, 1) * n_edges, ones(Int, n_edges), 1
    end

    # Build weighted graph
    # Simpler: show each cluster as a node with self-weight = count
    graph = SimpleDiGraph(n_clusters)
    edge_weights = zeros(Int, n_clusters, n_clusters)

    # Count interactions per cluster
    for i in 1:n_edges
        c = assignments[i]
        if c > 0
            edge_weights[c, c] += 1
        end
    end

    return graph, edge_weights, assignments, n_clusters
end

"""
    discretize_with_weights(sample::EdgeCentricSample{d}, n_clusters::Int;
                            rng=Random.default_rng())

Discretize with separate source/target clustering but return edge weights.
Uses Clustering.jl's k-means implementation.
"""
function discretize_with_weights(sample::EdgeCentricSample{d}, n_clusters::Int;
                                  rng::AbstractRNG=Random.default_rng()) where d
    n_edges = length(sample)

    if n_edges == 0
        return SimpleDiGraph(n_clusters), zeros(Int, n_clusters, n_clusters), Int[], Int[]
    end

    # Convert to matrices for clustering (Clustering.jl expects d×n)
    source_matrix = hcat([Vector(s) for s in sample.sources]...)
    target_matrix = hcat([Vector(t) for t in sample.targets]...)

    # Use Clustering.jl's k-means
    source_result = kmeans(source_matrix, n_clusters; maxiter=100)
    target_result = kmeans(target_matrix, n_clusters; maxiter=100)

    source_assignments = source_result.assignments
    target_assignments = target_result.assignments

    # Build weighted adjacency matrix
    edge_weights = zeros(Int, n_clusters, n_clusters)

    for i in 1:n_edges
        s = source_assignments[i]
        t = target_assignments[i]
        edge_weights[s, t] += 1
    end

    # Build graph
    graph = SimpleDiGraph(n_clusters)
    for s in 1:n_clusters
        for t in 1:n_clusters
            if edge_weights[s, t] > 0
                add_edge!(graph, s, t)
            end
        end
    end

    return graph, edge_weights, source_assignments, target_assignments
end

"""
    edge_count(graph::SimpleDiGraph) -> Int

Count the number of edges in a directed graph.
"""
edge_count(graph::SimpleDiGraph) = ne(graph)

"""
    node_count(graph::SimpleDiGraph) -> Int

Count the number of nodes in a directed graph.
"""
node_count(graph::SimpleDiGraph) = nv(graph)
