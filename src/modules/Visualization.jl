# Visualization utilities for IDPG using Makie
# Visualizations for B^d_+ (non-negative unit ball)

"""
    Bd_plus_to_2d(point::LatentPoint{2}) -> Point2f

Convert a 2D point in B^2_+ to 2D coordinates (identity for d=2).
"""
function Bd_plus_to_2d(point::LatentPoint{2})
    return Point2f(point[1], point[2])
end

"""
    Bd_plus_to_2d(point::LatentPoint{3}) -> Point2f

Convert a 3D point in B^3_+ to 2D coordinates by projecting onto the xy plane.
The z-coordinate can be encoded in color/size.
"""
function Bd_plus_to_2d(point::LatentPoint{3})
    return Point2f(point[1], point[2])
end

function Bd_plus_to_2d(point::AbstractVector)
    if length(point) == 2
        return Point2f(point[1], point[2])
    elseif length(point) >= 3
        return Point2f(point[1], point[2])
    else
        error("Point must have at least 2 dimensions")
    end
end

"""
    draw_Bd_plus_boundary!(ax; n_arc_points=50, color=:black, linewidth=2)

Draw the boundary of B^2_+ on an axis: two straight edges along axes and
a quarter-circle arc.
"""
function draw_Bd_plus_boundary!(ax; n_arc_points::Int=50, color=:black, linewidth::Real=2)
    # Draw axes (edges of B^2_+)
    lines!(ax, [Point2f(0, 0), Point2f(1, 0)], color=color, linewidth=linewidth)
    lines!(ax, [Point2f(0, 0), Point2f(0, 1)], color=color, linewidth=linewidth)

    # Draw quarter-circle arc
    θ = range(0, π/2, length=n_arc_points)
    arc_points = [Point2f(cos(t), sin(t)) for t in θ]
    lines!(ax, arc_points, color=color, linewidth=linewidth)

    return ax
end

"""
    plot_intensity_Bd_plus(ρ; resolution=50, colormap=:viridis, kwargs...)

Plot an intensity function on B^2_+ (2D non-negative unit ball).

# Arguments
- `ρ`: Intensity function that takes a 2-vector and returns a scalar
- `resolution`: Number of grid points for sampling
- `colormap`: Colormap to use

# Returns
A Makie Figure.
"""
function plot_intensity_Bd_plus(ρ; resolution::Int=50, colormap=:viridis, kwargs...)
    fig = Figure(size=(600, 550))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Intensity on B²₊")

    plot_intensity_Bd_plus!(ax, ρ; resolution=resolution, colormap=colormap, kwargs...)

    return fig
end

"""
    plot_intensity_Bd_plus!(ax, ρ; resolution=50, colormap=:viridis, kwargs...)

Plot intensity on B^2_+ onto existing axis.
"""
function plot_intensity_Bd_plus!(ax, ρ; resolution::Int=50, colormap=:viridis, kwargs...)
    h = 1.0 / (resolution - 1)
    points_2d = Point2f[]
    values = Float64[]

    for i in 0:(resolution-1)
        for j in 0:(resolution-1)
            x = i * h
            y = j * h

            # Check if in B^2_+
            if x^2 + y^2 <= 1.0
                push!(points_2d, Point2f(x, y))
                push!(values, ρ(SVector{2, Float64}(x, y)))
            end
        end
    end

    # Plot as scatter with color
    scatter!(ax, points_2d, color=values, colormap=colormap, markersize=8; kwargs...)

    # Draw B^2_+ boundary
    draw_Bd_plus_boundary!(ax)

    # Add labels
    text!(ax, Point2f(-0.05, 0), text="x₁", align=(:right, :center), fontsize=12)
    text!(ax, Point2f(0, -0.05), text="x₂", align=(:center, :top), fontsize=12)

    return ax
end

"""
    plot_node_centric_graph(sites::Vector{InteractionSite{d}}, graph::SimpleDiGraph;
                            kwargs...) where d

Plot a node-centric graph with nodes positioned in B^d_+ projected to 2D.

# Arguments
- `sites`: Interaction sites (nodes)
- `graph`: The generated graph

# Returns
A Makie Figure.
"""
function plot_node_centric_graph(sites::Vector{InteractionSite{d}}, graph::SimpleDiGraph;
                                 kwargs...) where d
    fig = Figure(size=(800, 400))
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Green Space (G)")
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Red Space (R)")

    plot_node_centric_graph!(ax1, ax2, sites, graph; kwargs...)

    return fig
end

"""
    plot_node_centric_graph!(ax_g, ax_r, sites, graph; kwargs...)

Plot node-centric graph onto existing axes.
"""
function plot_node_centric_graph!(ax_g, ax_r, sites::Vector{InteractionSite{d}},
                                  graph::SimpleDiGraph; kwargs...) where d
    n = length(sites)

    # Extract positions
    g_positions = [Bd_plus_to_2d(site.g) for site in sites]
    r_positions = [Bd_plus_to_2d(site.r) for site in sites]

    # Draw B^d_+ boundaries
    for ax in [ax_g, ax_r]
        draw_Bd_plus_boundary!(ax)
    end

    # Plot nodes in green space
    scatter!(ax_g, g_positions, color=:green, markersize=10, label="Nodes")

    # Plot nodes in red space
    scatter!(ax_r, r_positions, color=:red, markersize=10, label="Nodes")

    # Draw edges in green space
    for e in edges(graph)
        i, j = src(e), dst(e)
        p1 = g_positions[i]
        p2 = g_positions[j]
        lines!(ax_g, [p1, p2], color=(:blue, 0.3), linewidth=0.5)
    end

    return ax_g, ax_r
end

"""
    plot_edge_centric(sample::EdgeCentricSample{d}; kwargs...) where d

Plot edge-centric sample showing source and target distributions.

# Arguments
- `sample`: Edge-centric sample

# Returns
A Makie Figure with two panels.
"""
function plot_edge_centric(sample::EdgeCentricSample{d}; kwargs...) where d
    fig = Figure(size=(900, 400))
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Sources (G)")
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Targets (R)")

    plot_edge_centric!(ax1, ax2, sample; kwargs...)

    return fig
end

"""
    plot_edge_centric!(ax_sources, ax_targets, sample; kwargs...)

Plot edge-centric sample onto existing axes.
"""
function plot_edge_centric!(ax_sources, ax_targets, sample::EdgeCentricSample{d}; kwargs...) where d
    # Convert to 2D
    source_2d = [Bd_plus_to_2d(s) for s in sample.sources]
    target_2d = [Bd_plus_to_2d(t) for t in sample.targets]

    # Draw B^d_+ boundaries
    for ax in [ax_sources, ax_targets]
        draw_Bd_plus_boundary!(ax)
    end

    # Plot sources and targets
    if !isempty(source_2d)
        scatter!(ax_sources, source_2d, color=:green, markersize=6, alpha=0.6)
    end
    if !isempty(target_2d)
        scatter!(ax_targets, target_2d, color=:red, markersize=6, alpha=0.6)
    end

    # Add count annotation
    text!(ax_sources, Point2f(0.5, -0.1),
          text="n = " * string(length(sample.sources)),
          align=(:center, :top), fontsize=14)

    return ax_sources, ax_targets
end

"""
    animate_evolution(results; filename="evolution.mp4", framerate=10)

Create an animation of intensity evolution on B^2_+.

# Arguments
- `results`: Output from `evolve_and_track`
- `filename`: Output filename
- `framerate`: Frames per second
"""
function animate_evolution(results; filename::String="evolution.mp4", framerate::Int=10)
    times = results.times
    ρ_history = results.ρ_history
    grid = results.grid

    fig = Figure(size=(600, 550))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="t = 0.0")

    # Get value range for consistent colormap
    all_values = vcat(ρ_history...)
    vmin, vmax = extrema(all_values)
    if vmax - vmin < 1e-10
        vmax = vmin + 1.0
    end

    # Initial plot
    points_2d = [Bd_plus_to_2d(p) for p in grid.points]

    # Draw B^2_+ boundary
    draw_Bd_plus_boundary!(ax)

    # Create scatter plot
    sc = scatter!(ax, points_2d, color=ρ_history[1], colormap=:viridis,
                  colorrange=(vmin, vmax), markersize=8)
    Colorbar(fig[1, 2], sc, label="Intensity")

    # Animate
    record(fig, filename, 1:length(times); framerate=framerate) do frame
        ax.title = "t = " * string(round(times[frame], digits=2))
        sc.color = ρ_history[frame]
    end

    return filename
end

"""
    plot_formula_validation(theoretical::Vector{Float64}, empirical::Vector{Float64};
                            quantity::String="edges", kwargs...)

Plot theoretical vs empirical values for formula validation.

# Arguments
- `theoretical`: Theoretical expected values
- `empirical`: Empirical mean values
- `quantity`: Label for what's being compared

# Returns
A Makie Figure.
"""
function plot_formula_validation(theoretical::Vector{Float64}, empirical::Vector{Float64};
                                 quantity::String="edges", kwargs...)
    fig = Figure(size=(500, 500))
    ax = Axis(fig[1, 1],
              xlabel="Theoretical E[" * quantity * "]",
              ylabel="Empirical mean",
              title="Formula Validation")

    # Plot data points
    scatter!(ax, theoretical, empirical, color=:blue, markersize=10, label="Data")

    # Plot y=x line
    lims = (min(minimum(theoretical), minimum(empirical)),
            max(maximum(theoretical), maximum(empirical)))
    lines!(ax, [lims[1], lims[2]], [lims[1], lims[2]],
           color=:red, linestyle=:dash, linewidth=2, label="y = x")

    axislegend(ax, position=:lt)

    return fig
end

"""
    plot_sites_Bd_plus(sites::Vector{InteractionSite{d}}; kwargs...) where d

Plot interaction sites as points in B^d_+, showing both G and R coordinates.
"""
function plot_sites_Bd_plus(sites::Vector{InteractionSite{d}}; kwargs...) where d
    fig = Figure(size=(800, 400))
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Green coordinates (G)")
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Red coordinates (R)")

    # Draw B^d_+ boundaries
    for ax in [ax1, ax2]
        draw_Bd_plus_boundary!(ax)
    end

    # Plot G coordinates
    g_points = [Bd_plus_to_2d(site.g) for site in sites]
    scatter!(ax1, g_points, color=:green, markersize=8; kwargs...)

    # Plot R coordinates
    r_points = [Bd_plus_to_2d(site.r) for site in sites]
    scatter!(ax2, r_points, color=:red, markersize=8; kwargs...)

    return fig
end

