using FiniteElements

function create_cantilever_beam_linear(beam_length, beam_height, nx, ny, applied_force)
    # Parameters:
    # length: Length of the beam in the x-direction
    # height: Height of the beam in the y-direction
    # nx: Number of elements along the x-direction
    # ny: Number of elements along the y-direction
    # applied_force: Total force applied along the right side of the beam in the -y direction

    # Generate node coordinates
    n_nodes_x = nx + 1  # Number of nodes along the x-direction
    n_nodes_y = ny + 1  # Number of nodes along the y-direction
    dx = beam_length / nx    # Element size in x direction
    dy = beam_height / ny    # Element size in y direction

    nodes = zeros(2, n_nodes_x * n_nodes_y)
    node_id = 1
    for j in 0:ny
        for i in 0:nx
            nodes[:, node_id] = [i * dx, j * dy]
            node_id += 1
        end
    end

    # Generate element connectivity
    elements = zeros(Int, 4, nx * ny)
    element_id = 1
    for j in 1:ny
        for i in 1:nx
            n1 = (j - 1) * n_nodes_x + i
            n2 = n1 + 1
            n3 = n1 + n_nodes_x + 1
            n4 = n1 + n_nodes_x
            elements[:, element_id] = [n1, n2, n3, n4]
            element_id += 1
        end
    end

    # Generate fixed DOFs (nodes on the left side of the beam)
    fixed_dofs = Int[]
    for j in 1:n_nodes_y
        node_id = (j - 1) * n_nodes_x + 1  # Left side nodes (x = 0)
        push!(fixed_dofs, 2*node_id - 1)  # Fix x direction
        push!(fixed_dofs, 2*node_id)      # Fix y direction
    end

    # Generate forces (distributed along the right side nodes)
    right_side_nodes = [(n_nodes_x * (j - 1) + n_nodes_x) for j in 1:n_nodes_y]
    # force_per_node = applied_force / length * dy
    force_per_node = applied_force / length(right_side_nodes)
    forces = Dict()
    for node in right_side_nodes
        forces[node] = (0.0, -force_per_node)
    end

    return nodes, elements, forces, fixed_dofs
end

# Example usage:
beam_length = 1000.0  # Length of the beam
beam_height = 100.0   # Height of the beam
nx = 100         # Number of elements in x-direction
ny = 10         # Number of elements in y-direction
applied_force = 1e4  # Total force applied along the right side in -y direction

nodes, elements, forces, fixed_dofs = create_cantilever_beam_linear(beam_length, beam_height, nx, ny, applied_force)

material = Material(1e5, 0.3, 6)

u, K, f = solve_fem_linear(nodes, elements, material, forces, fixed_dofs)

##
using GLMakie

fig = Figure()
ax1 = Axis(fig[1,1]; aspect=DataAspect())
scatter!(ax1, nodes .+ 2 .* u, color=u[2,:], markersize=5)
fig
