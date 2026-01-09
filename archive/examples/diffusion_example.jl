# Diffusion Example
# Demonstrates intensity evolution via the diffusion equation on B^d_+.
# Shows how an initially concentrated intensity spreads over time.

using IDPG
using Random
using CairoMakie
using Statistics
using StaticArrays

rng = MersenneTwister(42)

println("=" ^ 60)
println("IDPG Diffusion Evolution Example")
println("=" ^ 60)

# Create a B^d_+ grid for PDE evolution
d = 2
resolution = 20
grid = create_Bd_plus_grid(d, resolution)
println("Created grid with ", length(grid.points), " points")

# Initial condition: Gaussian-like concentration near (0.8, 0.2)
# We'll use a BdPlusMixture as initial intensity
initial_mixture = BdPlusMixture([1.0], [[0.8, 0.2]], [15.0], 50.0)

function ρ_initial(p)
    return initial_mixture(p)
end

# Set up diffusion parameters
D = 0.1      # Diffusion coefficient (slower diffusion for better visualization)
dt = 0.0001  # Time step (needs to be small for stability)
t_final = 2.0      # Longer simulation to see full evolution
sample_interval = 0.1  # Snapshots

println("\nDiffusion parameters:")
println("  D = ", D)
println("  dt = ", dt)
println("  t_final = ", t_final)

# Evolve and track
println("\nEvolving intensity...")
results = evolve_and_track(
    ρ_initial, grid;
    pde_type=:diffusion,
    D=D,
    dt=dt,
    t_final=t_final,
    sample_interval=sample_interval
)

println("Evolution complete. ", length(results.times), " snapshots recorded.")

# Print statistics over time
println("\n--- Evolution Statistics ---")
println("Time\tTotal Intensity\tMax Intensity")
for (i, t) in enumerate(results.times)
    total = results.total_intensity[i]
    max_val = maximum(results.ρ_history[i])
    println(round(t, digits=3), "\t", round(total, digits=2), "\t\t", round(max_val, digits=2))
end

# Create visualization
println("\nGenerating plots...")

n_snapshots = length(results.times)
n_cols = min(n_snapshots, 4)
n_rows = ceil(Int, n_snapshots / n_cols)

fig = Figure(size=(300 * n_cols, 300 * n_rows))

for (i, t) in enumerate(results.times)
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1

    # Per-frame color range
    frame_min = minimum(results.ρ_history[i])
    frame_max = maximum(results.ρ_history[i])

    ax = Axis(fig[row, col],
        aspect=DataAspect(),
        title="t = " * string(round(t, digits=2)) * " (max=" * string(round(frame_max, digits=1)) * ")")

    # Draw B^d_+ boundary
    draw_Bd_plus_boundary!(ax)

    # Plot intensity values with per-frame scaling
    points_2d = [Bd_plus_to_2d(p) for p in grid.points]
    scatter!(ax, points_2d, color=results.ρ_history[i],
             colormap=:viridis, colorrange=(frame_min, frame_max), markersize=10)

    hidedecorations!(ax)
end

save("output/dynamics/diffusion_evolution.png", fig)
println("Saved diffusion_evolution.png")

# Create animation
println("\nCreating animation...")
try
    animate_evolution(results; filename="diffusion_evolution.mp4", framerate=5)
    println("Saved diffusion_evolution.mp4")
catch e
    println("Animation failed (may need GLMakie for video): ", e)
end

# Demonstrate effect on expected graph statistics
println("\n" * "=" ^ 60)
println("Effect on Graph Statistics")
println("=" ^ 60)

println("\nNote: As intensity diffuses and spreads to the boundary,")
println("the absorbing boundary condition causes total mass to decrease,")
println("which reduces E[N] over time.")

println("\nTime\tE[Total] (proxy for E[N])")
for (i, t) in enumerate(results.times)
    total = results.total_intensity[i]
    println(round(t, digits=3), "\t", round(total, digits=2))
end

println("\n" * "=" ^ 60)
println("Done!")
