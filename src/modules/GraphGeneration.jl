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
"""
struct EdgeCentricSample{d}
    sources::Vector{LatentPoint{d}}
    targets::Vector{LatentPoint{d}}
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

    # Convert to matrices for clustering
    source_matrix = hcat([Vector(s) for s in sample.sources]...)  # d × n_edges
    target_matrix = hcat([Vector(t) for t in sample.targets]...)  # d × n_edges

    # Cluster sources and targets
    # Using simple k-means assignment based on nearest centroid
    source_centroids = init_centroids(source_matrix, n_clusters, rng)
    target_centroids = init_centroids(target_matrix, n_clusters, rng)

    source_assignments = assign_clusters(source_matrix, source_centroids)
    target_assignments = assign_clusters(target_matrix, target_centroids)

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

"""
Initialize k-means centroids by randomly selecting k points.
"""
function init_centroids(data::Matrix{Float64}, k::Int, rng::AbstractRNG)
    n = size(data, 2)
    indices = sample(rng, 1:n, k, replace=false)
    return data[:, indices]
end

"""
Assign each point to nearest centroid.
"""
function assign_clusters(data::Matrix{Float64}, centroids::Matrix{Float64})
    n = size(data, 2)
    k = size(centroids, 2)
    assignments = zeros(Int, n)

    for i in 1:n
        min_dist = Inf
        for j in 1:k
            dist = sum((data[:, i] .- centroids[:, j]).^2)
            if dist < min_dist
                min_dist = dist
                assignments[i] = j
            end
        end
    end

    return assignments
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
