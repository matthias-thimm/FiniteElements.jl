
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
    force_per_node = applied_force / length(right_side_nodes)
    forces = Dict()
    for node in right_side_nodes
        forces[node] = (0.0, -force_per_node)
    end

    return nodes, elements, forces, fixed_dofs
end

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
    for y in 1:n_nodes_y
        for x in 1:n_nodes_x
            node_id = (y - 1) * n_nodes_x + x
            nodes[1, node_id] = dx * (x - 1)
            nodes[2, node_id] = dy * (y - 1)
        end
    end

    # Generate element connectivity for 9-node quadratic elements
    elements = zeros(Int, 9, nx * ny)
    for y = 1:ny
        for x = 1:nx
            element_id = (y - 1) * nx + x
            elements[1, element_id] = 2 * x - 1 + (y - 1) * (2 * nx + 1) * 2
            elements[2, element_id] = 2 * x + 1 + (y - 1) * (2 * nx + 1) * 2
            elements[3, element_id] = 2 * x + 1 + y * (2 * nx + 1) * 2
            elements[4, element_id] = 2 * x - 1 + y * (2 * nx + 1) * 2
            elements[5, element_id] = 2 * x + (y - 1) * (2 * nx + 1) * 2
            elements[6, element_id] = 2 * x + 1 + (y - 0.5) * (2 * nx + 1) * 2
            elements[7, element_id] = 2 * x + y * (2 * nx + 1) * 2
            elements[8, element_id] = 2 * x - 1 + (y - 0.5) * (2 * nx + 1) * 2
            elements[9, element_id] = 2 * x + (y - 0.5) * (2 * nx + 1) * 2
        end
    end

    # Generate fixed DOFs (nodes on the left side of the beam)
    left_side_nodes = findall(x -> isapprox(x, 0, atol=1e-6), nodes[1, :])
    fixed_nodes = Dict{Int,Symbol}()
    for node in left_side_nodes
        fixed_nodes[node] = :xy
    end

    # Generate forces (distributed along the right side nodes)
    right_side_nodes = findall(x -> isapprox(x, beam_length, atol=1e-6), nodes[1, :])
    force_per_node = applied_force / length(right_side_nodes)
    forces = Dict()
    for node in right_side_nodes
        forces[node] = (0.0, -force_per_node)
    end

    return nodes, elements, forces, fixed_nodes
end
