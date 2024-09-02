using FiniteElements

include("cantilever_beam_meshes.jl")
include("plotting.jl")

## FEM Simulation
beam_length = 1000.0 # Length of the beam
beam_height = 100.0 # Height of the beam
nx = 20 # Number of elements in x-direction
ny = 2 # Number of elements in y-direction
applied_force = 1e4  # Total force applied along the right side in -y direction
nodes, elements, forces, fixed_dofs = create_cantilever_beam_quadratic(beam_length,
                                                                       beam_height, nx, ny,
                                                                       applied_force)
material = Material(1e5, 0.3, 6)
u, K, f = solve_fem_quadratic(nodes, elements, material, forces, fixed_dofs)
println("-- results quadratic elements --")
println("number of elements: $(size(elements, 2))")
println("minimum displacement: $(minimum(u[2, :]))")

## Post-Processing
fig = Figure(size=(800, 300))
ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="cantilever beam quadratic")
deformed_nodes = nodes .+ 2 .* u
plot_mesh!(ax1, deformed_nodes, elements; color=u[2, :])
plot_edges_quadratic!(ax1, deformed_nodes, elements; color=:black)
plot_nodes!(ax1, deformed_nodes; color=:black, markersize=7)
save(joinpath("img", "cantilever_beam_quadratic.png"), fig; px_per_unit=2)
fig
