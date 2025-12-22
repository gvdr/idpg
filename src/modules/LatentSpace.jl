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

# =============================================================================
# Hyperspherical Coordinates for B^d_+
# =============================================================================
#
# Hyperspherical coordinates generalize polar (2D) and spherical (3D) to d dimensions.
# A point x ∈ R^d is represented as (r, φ₁, φ₂, ..., φ_{d-1}) where:
#   - r = ||x|| is the radius
#   - φ₁, ..., φ_{d-1} are angles
#
# For B^d_+ (non-negative orthant), all angles are in [0, π/2].
#
# Conversion formulas:
#   x₁ = r cos(φ₁)
#   x₂ = r sin(φ₁) cos(φ₂)
#   x₃ = r sin(φ₁) sin(φ₂) cos(φ₃)
#   ...
#   x_{d-1} = r sin(φ₁) sin(φ₂) ... sin(φ_{d-2}) cos(φ_{d-1})
#   x_d = r sin(φ₁) sin(φ₂) ... sin(φ_{d-1})

"""
    hyperspherical_to_cartesian(r::Real, angles::AbstractVector) -> LatentPoint

Convert hyperspherical coordinates (r, φ₁, ..., φ_{d-1}) to Cartesian coordinates.

# Arguments
- `r`: radius (distance from origin), must be ≥ 0
- `angles`: vector of d-1 angles in radians

# Returns
- `LatentPoint{d}`: Cartesian coordinates as a static vector

# Example
```julia
# 2D: polar coordinates (r=0.8, θ=π/4)
hyperspherical_to_cartesian(0.8, [π/4])  # ≈ [0.566, 0.566]

# 3D: spherical coordinates (r=1, θ=π/4, φ=π/3)
hyperspherical_to_cartesian(1.0, [π/4, π/3])  # ≈ [0.707, 0.354, 0.612]

# 4D point
hyperspherical_to_cartesian(0.9, [π/6, π/4, π/3])
```
"""
function hyperspherical_to_cartesian(r::Real, angles::AbstractVector{<:Real})
    d = length(angles) + 1
    x = zeros(d)

    if r == 0
        return SVector{d, Float64}(x)
    end

    # Precompute cumulative product of sines
    sin_prod = r  # running product: r * sin(φ₁) * sin(φ₂) * ...

    for i in 1:d-1
        x[i] = sin_prod * cos(angles[i])
        sin_prod *= sin(angles[i])
    end
    x[d] = sin_prod  # last coordinate is just the remaining product

    return SVector{d, Float64}(x)
end

"""
    cartesian_to_hyperspherical(x::AbstractVector) -> Tuple{Float64, Vector{Float64}}

Convert Cartesian coordinates to hyperspherical coordinates (r, φ₁, ..., φ_{d-1}).

# Arguments
- `x`: Cartesian coordinates in R^d

# Returns
- `r`: radius (norm of x)
- `angles`: vector of d-1 angles in radians

# Example
```julia
r, angles = cartesian_to_hyperspherical([0.5, 0.5, 0.707])
# r ≈ 1.0, angles ≈ [π/3, π/4]
```
"""
function cartesian_to_hyperspherical(x::AbstractVector{<:Real})
    d = length(x)
    r = norm(x)

    if r == 0 || d == 1
        return r, Float64[]
    end

    angles = zeros(d - 1)
    sin_prod = r  # running product: r * sin(φ₁) * sin(φ₂) * ...

    for i in 1:d-1
        if sin_prod > 0
            # cos(φᵢ) = xᵢ / (r * sin(φ₁) * ... * sin(φᵢ₋₁))
            cos_val = clamp(x[i] / sin_prod, -1.0, 1.0)
            angles[i] = acos(cos_val)
            sin_prod *= sin(angles[i])
        else
            # Degenerate case: set remaining angles to 0
            angles[i] = 0.0
        end
    end

    return r, angles
end

"""
    Bd_plus_from_hyperspherical(r::Real, angles::AbstractVector) -> LatentPoint

Create a point in B^d_+ from hyperspherical coordinates with validation.

Enforces B^d_+ constraints:
- r ∈ [0, 1]
- all angles ∈ [0, π/2] (to ensure non-negative Cartesian coordinates)

# Example
```julia
# Create a 4D point with radius 0.9 and three angles
p = Bd_plus_from_hyperspherical(0.9, [π/6, π/4, π/3])
```
"""
function Bd_plus_from_hyperspherical(r::Real, angles::AbstractVector{<:Real})
    @assert 0 <= r <= 1 "Radius must be in [0, 1] for B^d_+, got r=" * string(r)
    @assert all(0 .<= angles .<= π/2) "All angles must be in [0, π/2] for B^d_+, got angles=" * string(angles)

    return hyperspherical_to_cartesian(r, angles)
end

"""
    Bd_plus_to_hyperspherical(x::LatentPoint) -> Tuple{Float64, Vector{Float64}}

Convert a point in B^d_+ to hyperspherical coordinates.
Alias for `cartesian_to_hyperspherical` with B^d_+ semantics.

# Returns
- `r`: radius in [0, 1]
- `angles`: vector of d-1 angles, all in [0, π/2]
"""
function Bd_plus_to_hyperspherical(x::AbstractVector{<:Real})
    @assert in_Bd_plus(x) "Point must be in B^d_+"
    return cartesian_to_hyperspherical(x)
end
