# Ecological modeling utilities for IDPG
# Functions for guild assignment and food web matrix construction

"""
    assign_site_to_guild(site::InteractionSite, guild_full_means)

Assign a site to the nearest guild based on its full (g, r) signature.

The guild assignment uses Euclidean distance in the concatenated (g, r) space.
This is appropriate when guilds are defined by their position in both the
source (g) and target (r) trait spaces.

# Arguments
- `site`: An `InteractionSite{d}` with fields `g` and `r` (d-dimensional vectors)
- `guild_full_means`: Vector of 2d-dimensional vectors, where each element is
  the concatenation `[means_G[i]; means_R[i]]` representing guild i's centroid

# Returns
Integer index of the nearest guild (1-indexed).

# Example
```julia
# Define 3 guilds in 2D latent space (so 4D full means)
guild_means = [
    [0.8, 0.8, 0.2, 0.2],  # Guild 1: high g, low r
    [0.2, 0.2, 0.8, 0.8],  # Guild 2: low g, high r
    [0.5, 0.5, 0.5, 0.5],  # Guild 3: intermediate
]
guild_idx = assign_site_to_guild(site, guild_means)
```
"""
function assign_site_to_guild(site::InteractionSite, guild_full_means::AbstractVector)
    site_full = vcat(Vector(site.g), Vector(site.r))  # 2d-dimensional
    min_dist = Inf
    best_guild = 1
    for (i, full_mean) in enumerate(guild_full_means)
        dist = norm(site_full .- full_mean)
        if dist < min_dist
            min_dist = dist
            best_guild = i
        end
    end
    return best_guild
end

"""
    assign_point_to_guild(point, guild_means)

Assign a d-dimensional point to the nearest guild centroid.

This is a simpler version of `assign_site_to_guild` that works with raw
d-dimensional points rather than `InteractionSite` objects.

# Arguments
- `point`: A d-dimensional vector (or any iterable)
- `guild_means`: Vector of d-dimensional vectors representing guild centroids

# Returns
Integer index of the nearest guild (1-indexed).
"""
function assign_point_to_guild(point::AbstractVector, guild_means::AbstractVector)
    min_dist = Inf
    best_guild = 1
    for (i, mean) in enumerate(guild_means)
        dist = norm(point .- mean)
        if dist < min_dist
            min_dist = dist
            best_guild = i
        end
    end
    return best_guild
end

"""
    build_full_guild_means(means_G, means_R)

Combine separate G and R guild centroids into full (g, r) guild means.

# Arguments
- `means_G`: Vector of d-dimensional source (g) centroids for each guild
- `means_R`: Vector of d-dimensional target (r) centroids for each guild

# Returns
Vector of 2d-dimensional vectors where each is `[means_G[i]; means_R[i]]`.

# Example
```julia
means_G = [[0.8, 0.8], [0.2, 0.2]]
means_R = [[0.2, 0.2], [0.8, 0.8]]
full_means = build_full_guild_means(means_G, means_R)
# Returns: [[0.8, 0.8, 0.2, 0.2], [0.2, 0.2, 0.8, 0.8]]
```
"""
function build_full_guild_means(means_G::AbstractVector, means_R::AbstractVector)
    @assert length(means_G) == length(means_R) "means_G and means_R must have same length"
    return [vcat(Vector(means_G[i]), Vector(means_R[i])) for i in 1:length(means_G)]
end

"""
    compute_foodweb_matrix(sample::FullEdgeCentricSample, guild_full_means)

Compute food web interaction matrix from edge samples and guild centroids.

Each edge (source → target) is counted into the guild-guild interaction matrix
based on which guilds the source and target sites belong to.

# Arguments
- `sample`: A `FullEdgeCentricSample{d}` containing source and target sites
- `guild_full_means`: Vector of 2d-dimensional guild centroids

# Returns
Matrix where `M[i,j]` counts edges from guild i (as source) to guild j (as target).

# Example
```julia
sample = generate_node_centric_full(sites)
guild_means = build_full_guild_means(means_G, means_R)
fw_matrix = compute_foodweb_matrix(sample, guild_means)
```
"""
function compute_foodweb_matrix(sample::FullEdgeCentricSample, guild_full_means::AbstractVector)
    n_guilds = length(guild_full_means)
    edge_weights = zeros(n_guilds, n_guilds)

    for k in 1:length(sample)
        src_guild = assign_site_to_guild(sample.source_sites[k], guild_full_means)
        tgt_guild = assign_site_to_guild(sample.target_sites[k], guild_full_means)
        edge_weights[src_guild, tgt_guild] += 1
    end

    return edge_weights
end

"""
    compute_foodweb_matrix(sample::FullEdgeCentricSample, means_G, means_R)

Convenience method that builds full guild means from separate G and R centroids.

# Arguments
- `sample`: A `FullEdgeCentricSample{d}` containing source and target sites
- `means_G`: Vector of d-dimensional source (g) centroids for each guild
- `means_R`: Vector of d-dimensional target (r) centroids for each guild

# Returns
Matrix where `M[i,j]` counts edges from guild i to guild j.
"""
function compute_foodweb_matrix(sample::FullEdgeCentricSample, means_G::AbstractVector, means_R::AbstractVector)
    guild_full_means = build_full_guild_means(means_G, means_R)
    return compute_foodweb_matrix(sample, guild_full_means)
end

"""
    compute_foodweb_matrix(sample::EdgeCentricSample, means_G, means_R)

Compute food web matrix from EdgeCentricSample using source g and target r only.

Note: EdgeCentricSample only stores the connection-relevant coordinates (g_source, r_target),
not the full site signatures. Guild assignment uses nearest match in the respective spaces.

# Arguments
- `sample`: An `EdgeCentricSample{d}` with g_source and r_target coordinates
- `means_G`: Vector of d-dimensional source (g) centroids
- `means_R`: Vector of d-dimensional target (r) centroids

# Returns
Matrix where `M[i,j]` counts edges from guild i to guild j.
"""
function compute_foodweb_matrix(sample::EdgeCentricSample, means_G::AbstractVector, means_R::AbstractVector)
    n_guilds = length(means_G)
    edge_weights = zeros(n_guilds, n_guilds)

    for k in 1:sample.n_edges
        src_guild = assign_point_to_guild(sample.g_source[k], means_G)
        tgt_guild = assign_point_to_guild(sample.r_target[k], means_R)
        edge_weights[src_guild, tgt_guild] += 1
    end

    return edge_weights
end

"""
    normalize_foodweb_matrix(M; mode=:row)

Normalize food web matrix rows or columns.

# Arguments
- `M`: Raw edge count matrix
- `mode`: `:row` normalizes rows (outgoing from each guild),
         `:col` normalizes columns (incoming to each guild),
         `:total` divides by total edge count

# Returns
Normalized matrix.
"""
function normalize_foodweb_matrix(M::AbstractMatrix; mode::Symbol=:row)
    if mode == :row
        row_sums = sum(M, dims=2)
        return M ./ max.(row_sums, 1.0)
    elseif mode == :col
        col_sums = sum(M, dims=1)
        return M ./ max.(col_sums, 1.0)
    elseif mode == :total
        return M ./ max(sum(M), 1.0)
    else
        error("Unknown normalization mode: " * string(mode) * ". Use :row, :col, or :total")
    end
end

# ============================================================================
# Guild Sampling
# ============================================================================

"""
    sample_guild_position(μ, κ; rng=Random.default_rng())

Sample a position from a Gaussian centered at μ with concentration κ,
projected onto B^d_+.

# Arguments
- `μ`: Mean position (d-dimensional vector)
- `κ`: Concentration parameter (higher = more peaked); σ = 1/√κ
- `rng`: Random number generator

# Returns
A d-dimensional vector in B^d_+.

# Example
```julia
μ = [0.7, 0.3, 0.0, 0.0]
pos = sample_guild_position(μ, 50.0)
```
"""
function sample_guild_position(μ::AbstractVector, κ::Real; rng::AbstractRNG=Random.default_rng())
    d = length(μ)
    σ = 1.0 / sqrt(κ)
    # Sample from Gaussian
    x = μ .+ σ .* randn(rng, d)
    # Project to B^d_+
    return Vector(project_to_Bd_plus(x))
end

# ============================================================================
# Expected Edge Computation
# ============================================================================

"""
    compute_expected_guild_edges(K_hat, π, Λ)

Compute expected number of edges between guilds under node-centric model.

# Arguments
- `K_hat`: Affinity matrix where K_hat[i,j] = μ_G[i] ⋅ μ_R[j]
- `π`: Guild mixing weights (probabilities summing to 1)
- `Λ`: Total intensity (expected number of sites)

# Returns
Matrix E[i,j] = expected edges from guild i to guild j.

# Formula
E[i,j] = Λ² × π[i] × π[j] × K_hat[i,j]

This follows from the node-centric model where:
- Expected sites from guild i: Λ × π[i]
- Each pair (from i, to j) connects with probability K_hat[i,j]
"""
function compute_expected_guild_edges(K_hat::AbstractMatrix, π::AbstractVector, Λ::Real)
    n = length(π)
    E = zeros(n, n)
    for i in 1:n
        for j in 1:n
            E[i, j] = Λ^2 * π[i] * π[j] * K_hat[i, j]
        end
    end
    return E
end

"""
    compute_guild_affinity(M_G, M_R)

Compute affinity matrix K_hat[i,j] = M_G[i,:] ⋅ M_R[j,:] from guild centroids.

# Arguments
- `M_G`: Matrix of green (source) centroids, each row is a guild (n_guilds × d)
- `M_R`: Matrix of red (target) centroids, each row is a guild (n_guilds × d)

# Returns
Matrix K_hat where K_hat[i,j] is the dot product of guild i's source centroid
with guild j's target centroid.
"""
function compute_guild_affinity(M_G::AbstractMatrix, M_R::AbstractMatrix)
    n = size(M_G, 1)
    K_hat = zeros(n, n)
    for i in 1:n
        for j in 1:n
            K_hat[i, j] = dot(M_G[i, :], M_R[j, :])
        end
    end
    return K_hat
end

# ============================================================================
# Trophic Layout for Visualization
# ============================================================================

"""
    trophic_layout(n_guilds)

Generate node positions for trophic-level visualization.

Returns (x, y) coordinates with producers at bottom, apex predators at top.
Standard 4-guild layout places: Producer (bottom), Herbivore, Carnivore, Apex (top).

# Arguments
- `n_guilds`: Number of guilds/trophic levels

# Returns
Tuple (x, y) of coordinate vectors.
"""
function trophic_layout(n_guilds::Int)
    if n_guilds == 4
        # Standard 4-guild layout: Producer at bottom, Apex at top
        x = [0.0, -0.3, 0.3, 0.0]
        y = [0.0, 0.4, 0.7, 1.0]
    elseif n_guilds == 5
        # 5-guild food web: Prod(TL0), SmallHerb(TL1), LargeHerb(TL1), SmallPred(TL2), Apex(TL3)
        x = [0.5, 0.25, 0.75, 0.5, 0.5]
        y = [0.0, 1.0, 1.0, 2.0, 3.0]
    else
        # Generic layout: spread horizontally, stack vertically
        x = collect(range(-0.5, 0.5, length=n_guilds))
        y = collect(range(0.0, 1.0, length=n_guilds))
    end
    return x, y
end

# ============================================================================
# Food Web Centroid Construction
# ============================================================================

"""
    default_trophic_affinity()

Return a default 5×5 target trophic affinity matrix K*.

Rows = resource (prey), Columns = consumer (predator)
Guilds: 1=Producers, 2=Small Herbivores, 3=Large Herbivores, 4=Small Predators, 5=Apex

Design principles:
- Column 1 ≈ 0: Producers don't consume
- Diagonal small: Limited intra-guild predation
- Clear trophic flow: Prod → Herb → Pred → Apex

# Returns
5×5 Matrix{Float64} representing target guild-guild affinities.
"""
function default_trophic_affinity()
    #        C1(Prod) C2(S.Herb) C3(L.Herb) C4(S.Pred) C5(Apex)
    K = [0.0   0.9       0.9        0.1        0.05;   # R1 (Prod)
         0.0   0.05      0.05       0.8        0.5;    # R2 (S.Herb)
         0.0   0.05      0.05       0.6        0.8;    # R3 (L.Herb)
         0.0   0.0       0.0        0.05       0.9;    # R4 (S.Pred)
         0.0   0.0       0.0        0.05       0.05]   # R5 (Apex)
    return K
end

"""
    svd_initialize_centroids(K_star, d)

Initialize guild centroids via truncated SVD of target affinity matrix.

# Arguments
- `K_star`: Target affinity matrix (n_guilds × n_guilds)
- `d`: Latent dimension for centroids

# Returns
Tuple (M_G, M_R, singular_values) where M_G * M_R' ≈ K_star (rank-d approximation).
"""
function svd_initialize_centroids(K_star::AbstractMatrix, d::Int)
    F = svd(K_star)

    # Truncate to rank d
    r = min(d, length(F.S))
    U_r = F.U[:, 1:r]
    S_r = F.S[1:r]
    V_r = F.V[:, 1:r]

    # Form initial centroids: M_G = U_r * sqrt(Σ), M_R = V_r * sqrt(Σ)
    sqrt_S = sqrt.(S_r)
    M_G = U_r * Diagonal(sqrt_S)
    M_R = V_r * Diagonal(sqrt_S)

    # Pad with zeros if d > rank(K_star)
    if r < d
        n_G, n_R = size(K_star)
        M_G = hcat(M_G, zeros(n_G, d - r))
        M_R = hcat(M_R, zeros(n_R, d - r))
    end

    return M_G, M_R, F.S
end

"""
    count_negative_mass(M_G, M_R)

Count sum of negative entries in centroid matrices (measure of constraint violation).
"""
function count_negative_mass(M_G::AbstractMatrix, M_R::AbstractMatrix)
    neg_G = sum(min.(M_G, 0.0))
    neg_R = sum(min.(M_R, 0.0))
    return neg_G + neg_R
end

"""
    flip_signs_for_positivity!(M_G, M_R)

Flip signs of columns to maximize positive entries.
This preserves dot products since we flip both M_G and M_R.
"""
function flip_signs_for_positivity!(M_G::AbstractMatrix, M_R::AbstractMatrix)
    d = size(M_G, 2)
    for k in 1:d
        pos_G = sum(M_G[:, k] .> 0)
        neg_G = sum(M_G[:, k] .< 0)
        pos_R = sum(M_R[:, k] .> 0)
        neg_R = sum(M_R[:, k] .< 0)

        if (neg_G + neg_R) > (pos_G + pos_R)
            M_G[:, k] .*= -1
            M_R[:, k] .*= -1
        end
    end
end

"""
    apply_givens_rotation!(M_G, M_R, i, j, θ)

Apply Givens rotation in the (i,j) plane by angle θ to both centroid matrices.
"""
function apply_givens_rotation!(M_G::AbstractMatrix, M_R::AbstractMatrix, i::Int, j::Int, θ::Real)
    c, s = cos(θ), sin(θ)
    for M in (M_G, M_R)
        col_i = M[:, i] .* c .- M[:, j] .* s
        col_j = M[:, i] .* s .+ M[:, j] .* c
        M[:, i] .= col_i
        M[:, j] .= col_j
    end
end

"""
    optimize_rotation!(M_G, M_R; n_iter=100, n_angles=20)

Iteratively apply Givens rotations to minimize negative entries in centroids.
"""
function optimize_rotation!(M_G::AbstractMatrix, M_R::AbstractMatrix; n_iter::Int=100, n_angles::Int=20)
    d = size(M_G, 2)
    angles = range(-π/2, π/2, length=n_angles)

    for _ in 1:n_iter
        improved = false
        for i in 1:d
            for j in (i+1):d
                current_neg = count_negative_mass(M_G, M_R)
                best_θ = 0.0
                best_neg = current_neg

                for θ in angles
                    M_G_copy = copy(M_G)
                    M_R_copy = copy(M_R)
                    apply_givens_rotation!(M_G_copy, M_R_copy, i, j, θ)
                    neg = count_negative_mass(M_G_copy, M_R_copy)
                    if neg > best_neg
                        best_neg = neg
                        best_θ = θ
                    end
                end

                if best_neg > current_neg
                    apply_givens_rotation!(M_G, M_R, i, j, best_θ)
                    improved = true
                end
            end
        end
        if !improved
            break
        end
    end
end

"""
    project_rows_to_Bd_plus!(M; scale=0.99)

Project each row of matrix M to B^d_+: clamp negatives to 0, scale to unit ball.

# Arguments
- `M`: Matrix to modify in-place
- `scale`: Maximum norm for rows (default 0.99 to stay interior)
"""
function project_rows_to_Bd_plus!(M::AbstractMatrix; scale::Real=0.99)
    n = size(M, 1)
    for i in 1:n
        M[i, :] .= max.(M[i, :], 0.0)
        norm_i = norm(M[i, :])
        if norm_i > scale
            M[i, :] .*= scale / norm_i
        end
    end
end

"""
    solve_constrained_ls(A, b; max_norm=0.99, max_iter=100, lr=0.1)

Solve min_x ||Ax - b||² s.t. x ≥ 0, ||x|| ≤ max_norm via projected gradient descent.
"""
function solve_constrained_ls(A::AbstractMatrix, b::AbstractVector;
                               max_norm::Real=0.99, max_iter::Int=100, lr::Real=0.1)
    d = size(A, 2)
    x = zeros(d)

    for _ in 1:max_iter
        residual = A * x - b
        grad = 2.0 * A' * residual
        x .-= lr * grad
        x .= max.(x, 0.0)
        nrm = norm(x)
        if nrm > max_norm
            x .*= max_norm / nrm
        end
    end

    return x
end

"""
    alternating_optimization!(M_G, M_R, K_star; max_iter=100, tol=1e-6, verbose=false)

Refine centroids via constrained alternating least squares.
"""
function alternating_optimization!(M_G::AbstractMatrix, M_R::AbstractMatrix,
                                    K_star::AbstractMatrix; max_iter::Int=100, tol::Real=1e-6, verbose::Bool=false)
    n_G, n_R = size(K_star)
    prev_error = Inf

    for iter in 1:max_iter
        # Fix M_R, optimize M_G
        for i in 1:n_G
            M_G[i, :] .= solve_constrained_ls(M_R, K_star[i, :])
        end

        # Fix M_G, optimize M_R
        for j in 1:n_R
            M_R[j, :] .= solve_constrained_ls(M_G, K_star[:, j])
        end

        # Check convergence
        K_approx = M_G * M_R'
        error = norm(K_star - K_approx, 2) / norm(K_star, 2)

        if abs(prev_error - error) < tol
            verbose && println("  Converged at iteration " * string(iter) * " with relative error " * string(round(error, digits=4)))
            break
        end
        prev_error = error

        if verbose && iter % 20 == 0
            println("  Iteration " * string(iter) * ": relative error = " * string(round(error, digits=4)))
        end
    end
end

"""
    verify_centroids(M_G, M_R, K_star)

Report approximation quality and constraint satisfaction for guild centroids.

# Returns
NamedTuple with fields:
- `K_approx`: Achieved approximation M_G * M_R'
- `abs_error`, `rel_error`: Frobenius norm errors
- `all_positive_G`, `all_positive_R`: Whether all entries are non-negative
- `all_in_ball_G`, `all_in_ball_R`: Whether all rows have norm ≤ 1
- `norms_G`, `norms_R`: Row norms
"""
function verify_centroids(M_G::AbstractMatrix, M_R::AbstractMatrix, K_star::AbstractMatrix)
    K_approx = M_G * M_R'

    abs_error = norm(K_star - K_approx, 2)
    rel_error = abs_error / norm(K_star, 2)

    all_positive_G = all(M_G .>= 0)
    all_positive_R = all(M_R .>= 0)
    norms_G = [norm(M_G[i, :]) for i in 1:size(M_G, 1)]
    norms_R = [norm(M_R[j, :]) for j in 1:size(M_R, 1)]
    all_in_ball_G = all(norms_G .<= 1.0 + 1e-10)
    all_in_ball_R = all(norms_R .<= 1.0 + 1e-10)

    return (
        K_approx = K_approx,
        abs_error = abs_error,
        rel_error = rel_error,
        all_positive_G = all_positive_G,
        all_positive_R = all_positive_R,
        all_in_ball_G = all_in_ball_G,
        all_in_ball_R = all_in_ball_R,
        norms_G = norms_G,
        norms_R = norms_R,
    )
end

"""
    construct_food_web_centroids(K_star=default_trophic_affinity(); d=4, verbose=true)

Construct guild centroids in B^d_+ that approximate a target affinity matrix.

Uses a multi-step algorithm:
1. SVD initialization for rank-d approximation
2. Givens rotations to minimize negative entries
3. Projection to B^d_+ (positive orthant of unit ball)
4. Alternating optimization to refine the approximation

# Arguments
- `K_star`: Target affinity matrix (default: 5-guild trophic structure)
- `d`: Latent dimension (default: 4)
- `verbose`: Print progress information

# Returns
NamedTuple with fields:
- `M_G`: n_guilds × d matrix of source (green/proposing) centroids
- `M_R`: n_guilds × d matrix of target (red/accepting) centroids
- `K_star`: The target affinity matrix
- `K_approx`: The achieved approximation
- `singular_values`: Singular values of K_star
- `rel_error`: Relative Frobenius error
- `guild_names`: Names of the guilds

# Example
```julia
result = construct_food_web_centroids()
println("Relative error: ", result.rel_error)
println("M_G centroids: ", result.M_G)
```
"""
function construct_food_web_centroids(K_star::AbstractMatrix=default_trophic_affinity();
                                       d::Int=4, verbose::Bool=true)
    verbose && println("=" ^ 60)
    verbose && println("Food Web Centroid Construction")
    verbose && println("=" ^ 60)

    n_guilds = size(K_star, 1)
    verbose && println("Target matrix: " * string(n_guilds) * "×" * string(n_guilds))
    verbose && println("Latent dimension: d = " * string(d))

    # Step 1: SVD initialization
    verbose && println("\nStep 1: SVD initialization...")
    M_G, M_R, singular_values = svd_initialize_centroids(K_star, d)
    verbose && println("  Singular values of K*: " * string(round.(singular_values, digits=4)))

    # Step 2: Rotation for positivity
    verbose && println("\nStep 2: Rotation for positivity...")
    flip_signs_for_positivity!(M_G, M_R)
    neg_before = count_negative_mass(M_G, M_R)
    optimize_rotation!(M_G, M_R)
    neg_after = count_negative_mass(M_G, M_R)
    verbose && println("  Negative mass: " * string(round(neg_before, digits=4)) * " → " * string(round(neg_after, digits=4)))

    # Step 3: Projection
    verbose && println("\nStep 3: Projection to B^d_+...")
    project_rows_to_Bd_plus!(M_G)
    project_rows_to_Bd_plus!(M_R)

    # Step 4: Alternating optimization
    verbose && println("\nStep 4: Alternating optimization...")
    alternating_optimization!(M_G, M_R, K_star; verbose=verbose)

    # Step 5: Verification
    verbose && println("\nStep 5: Verification...")
    v = verify_centroids(M_G, M_R, K_star)

    verbose && println("  Relative error: " * string(round(v.rel_error, digits=4)))
    verbose && println("  All positive (G): " * string(v.all_positive_G))
    verbose && println("  All positive (R): " * string(v.all_positive_R))
    verbose && println("  All in ball (G): " * string(v.all_in_ball_G))
    verbose && println("  All in ball (R): " * string(v.all_in_ball_R))

    return (
        M_G = M_G,
        M_R = M_R,
        K_star = K_star,
        K_approx = v.K_approx,
        singular_values = singular_values,
        rel_error = v.rel_error,
        guild_names = ["Producers", "Small Herbivores", "Large Herbivores", "Small Predators", "Apex Predators"],
    )
end
