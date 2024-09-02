using FiniteElements

# Problem setup
L = 1.0        # Length of the element (in meters)
H = 1.0        # Height of the element (in meters)
applied_force = 1e8  # Total force applied in the right direction (in Newtons)
E = 200e9      # Young's modulus (e.g., steel) in Pascals (N/m^2)
thickness = 0.1  # Thickness of the structure (in meters)

# Define the nodes
nodes = [0.0  L    L  0.0  L/2  L    L/2  0.0  L/2;
         0.0  0.0  H  H    0.0  H/2  H    H/2  H/2]

# Define elements (single quadrilateral element)
elements = [1; 2; 3; 4; 5; 6; 7; 8; 9;;]

# Material properties
material = Material(E, 0.3, thickness)  # Steel-like material with thickness

# Fixed DOFs (nodes on the left)
# Fix nodes 1 (x and y directions) and 4 (x direction) and 8 (x direction)
fixed_dofs = Dict(1 => :xy, 4 => :x, 8 => :x)

# Distributed force on the right side (nodes 2 and 3)
force_per_node = applied_force / 3
forces = Dict(
    2 => (force_per_node, 0.0),
    3 => (force_per_node, 0.0),
    6 => (force_per_node, 0.0),
)

u, K, f = solve_fem_quadratic(nodes, elements, material, forces, fixed_dofs)

# Analytical displacement for comparison
analytical_displacement(F, L, E, A) = F * L / (E * A)

# Analytical solution
A = thickness * H  # Cross-sectional area (in square meters)
u_analytical = analytical_displacement(applied_force, L, E, A)
println("Numerical displacement at node 2: ", u[1, 2])
println("Numerical displacement at node 3: ", u[1, 3])
println("Analytical displacement: ", u_analytical)

##
include("plotting.jl")
fig = Figure(size=(300, 300))
ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="single element test")
deformed_nodes = nodes .+ 10 .* u
plot_edges_quadratic!(ax1, deformed_nodes, elements; color=:black)
plot_nodes!(ax1, deformed_nodes; color=:black, markersize=7)
save(joinpath("img", "single_element_test_quadratic.png"), fig; px_per_unit=2)
fig
