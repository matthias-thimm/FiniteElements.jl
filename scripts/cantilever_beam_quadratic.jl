using FiniteElements

# Function to create a cantilever beam mesh with quadratic 8-node elements
function create_cantilever_beam_quadratic(beam_length, beam_height, nx, ny, applied_force)
    # Parameters:
    # beam_length: Length of the beam in the x-direction
    # beam_height: Height of the beam in the y-direction
    # nx: Number of elements along the x-direction
    # ny: Number of elements along the y-direction
    # applied_force: Total force applied along the right side of the beam in the -y direction

    # Generate node coordinates (quadratic elements)
    n_nodes_x = 2 * nx + 1  # Number of nodes along the x-direction (quadratic elements have mid-side nodes)
    n_nodes_y = 2 * ny + 1  # Number of nodes along the y-direction
    dx = beam_length / (2 * nx)    # Element size in x direction (account for mid-side nodes)
    dy = beam_height / (2 * ny)    # Element size in y direction

    nodes = zeros(2, n_nodes_x * n_nodes_y)
    node_id = 1
    for j in 0:2*ny
        for i in 0:2*nx
            nodes[:, node_id] = [i * dx, j * dy]
            node_id += 1
        end
    end

    # Generate element connectivity for 8-node quadratic elements
    elements = zeros(Int, 8, nx * ny)
    element_id = 1
    for j in 1:ny
        for i in 1:nx
            n1 = (2*j - 2) * n_nodes_x + 2*i - 1
            n2 = n1 + 2
            n3 = n1 + 2 * n_nodes_x + 2
            n4 = n1 + 2 * n_nodes_x
            n5 = n1 + 1
            n6 = n1 + n_nodes_x + 2
            n7 = n1 + 2 * n_nodes_x + 1
            n8 = n1 + n_nodes_x
            elements[:, element_id] = [n1, n2, n3, n4, n5, n6, n7, n8]
            element_id += 1
        end
    end

    # Generate fixed DOFs (nodes on the left side of the beam)
    left_side_nodes = findall(x -> isapprox(x, 0, atol=1e-6), nodes[1,:])
    fixed_nodes = Dict{Int,Symbol}()
    for node in left_side_nodes
        fixed_nodes[node] = :xy
    end

    # Generate forces (distributed along the right side nodes)
    right_side_nodes = findall(x -> isapprox(x, beam_length, atol=1e-6), nodes[1,:])
    force_per_node = applied_force / length(right_side_nodes)
    forces = Dict()
    for node in right_side_nodes
        forces[node] = (0.0, -force_per_node)
    end

    return nodes, elements, forces, fixed_nodes
end

# Example usage:
beam_length = 1000.0  # Length of the beam
beam_height = 100.0   # Height of the beam
nx = 25          # Number of elements in x-direction
ny = 5           # Number of elements in y-direction
applied_force = 1e4  # Total force applied along the right side in -y direction

nodes, elements, forces, fixed_dofs = create_cantilever_beam_quadratic(beam_length, beam_height, nx, ny, applied_force)

material = Material(1e5, 0.3, 6)

# Solve FEM problem for quadratic elements
u, K, f = solve_fem_quadratic(nodes, elements, material, forces, fixed_dofs)

# Visualize the deformed shape using GLMakie
fig = Figure()
ax1 = Axis(fig[1,1]; aspect=DataAspect())
scatter!(ax1, nodes .+ 2 .* u, color=u[2,:], markersize=5)
fig
