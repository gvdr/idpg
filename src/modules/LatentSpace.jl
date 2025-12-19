# Latent space geometry utilities for the non-negative unit ball B^d_+
# B^d_+ = {x in R^d : x >= 0, ||x|| <= 1}

"""
Type alias for points in the non-negative unit ball B^d_+ in R^d.
"""
const LatentPoint{d} = SVector{d, Float64}

"""
    in_Bd_plus(x::AbstractVector; tol=1e-10) -> Bool

Check if point lies in the non-negative unit ball B^d_+ (within tolerance).
A point is in B^d_+ if all coordinates are non-negative and ||x|| <= 1.
"""
function in_Bd_plus(x::AbstractVector; tol::Float64=1e-10)
    return all(xi -> xi >= -tol, x) && norm(x) <= 1.0 + tol
end

"""
    on_Bd_plus_boundary(x::AbstractVector; tol=1e-10) -> Bool

Check if point lies on the boundary of B^d_+.
Boundary consists of: (1) coordinate faces where x_k = 0, (2) sphere where ||x|| = 1.
"""
function on_Bd_plus_boundary(x::AbstractVector; tol::Float64=1e-10)
    on_coordinate_face = any(xi -> abs(xi) < tol, x)
    on_sphere = abs(norm(x) - 1.0) < tol
    return on_coordinate_face || on_sphere
end

"""
    Bd_plus_outward_normal(x::AbstractVector; tol=1e-10) -> Vector{Float64}

Compute outward normal at a boundary point of B^d_+.
- At coordinate face x_k = 0: normal is -e_k
- At sphere ||x|| = 1: normal is x/||x||
"""
function Bd_plus_outward_normal(x::AbstractVector; tol::Float64=1e-10)
    d = length(x)
    # Check coordinate faces first (they take priority at corners)
    for k in 1:d
        if abs(x[k]) < tol
            n = zeros(d)
            n[k] = -1.0  # outward from B^d_+ is -e_k
            return n
        end
    end
    # Otherwise on sphere: normal is x/||x||
    return x ./ norm(x)
end

"""
    project_to_Bd_plus(x::AbstractVector) -> LatentPoint

Project a point in R^d onto the non-negative unit ball B^d_+.
Two-step projection: (1) clamp to non-negative, (2) scale if norm > 1.
"""
function project_to_Bd_plus(x::AbstractVector{T}) where T
    d = length(x)
    # Step 1: Project onto non-negative orthant
    x_pos = max.(x, 0.0)
    # Step 2: Project onto unit ball
    n = norm(x_pos)
    result = n > 1.0 ? x_pos ./ n : x_pos
    return SVector{d, Float64}(result)
end

"""
    uniform_Bd_plus_sample(d::Int; rng=Random.default_rng()) -> LatentPoint{d}

Sample uniformly from the non-negative unit ball B^d_+ using rejection sampling.
"""
function uniform_Bd_plus_sample(d::Int; rng::AbstractRNG=Random.default_rng())
    while true
        # Sample from [0, 1]^d
        x = rand(rng, d)
        if norm(x) <= 1.0
            return SVector{d, Float64}(x)
        end
    end
end

"""
    Bd_plus_volume(d::Int) -> Float64

Compute the volume of the non-negative unit ball B^d_+.
This is 1/2^d times the volume of the full d-ball.
Volume of d-ball is pi^(d/2) / Gamma(d/2 + 1).
"""
function Bd_plus_volume(d::Int)
    # Volume of unit ball in R^d
    ball_volume = pi^(d/2) / gamma(d/2 + 1)
    # B^d_+ is 1/2^d of the full ball
    return ball_volume / 2^d
end

"""
    connection_probability(g::LatentPoint, r::LatentPoint) -> Float64

Compute the connection probability (dot product) for two latent points.
This is the probability of an edge from a node at `g` to a node at `r`.

Note: For points in B^d_+, the dot product is guaranteed to be in [0, 1]
since g, r >= 0 and ||g||, ||r|| <= 1 implies g·r <= ||g|| ||r|| <= 1.
"""
function connection_probability(g::LatentPoint, r::LatentPoint)
    p = dot(g, r)
    # Clamp to [0, 1] for numerical safety
    return clamp(p, 0.0, 1.0)
end

"""
    radial_coordinate(x::AbstractVector) -> Float64

Return the radial coordinate (norm) of a point.
"""
function radial_coordinate(x::AbstractVector)
    return norm(x)
end

"""
    angular_coordinates(x::AbstractVector) -> Vector{Float64}

Return the angular coordinates (direction) of a point.
Returns x/||x|| if ||x|| > 0, otherwise returns zeros.
"""
function angular_coordinates(x::AbstractVector)
    n = norm(x)
    return n > 0 ? x ./ n : zeros(length(x))
end
