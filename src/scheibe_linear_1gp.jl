# using LinearAlgebra, SparseArrays

# Define the material properties (Plane stress)
struct Material
    E::Float64  # Young's modulus
    ν::Float64  # Poisson's ratio
end

function plane_stress_stiffness(material::Material)
    E, ν = material.E, material.ν
    factor = E / (1.0 - ν^2)
    D = factor * [1.0    ν      0.0;
                  ν      1.0    0.0;
                  0.0    0.0    (1.0 - ν) / 2.0]
    return D
end

# Function to calculate the element stiffness matrix for a quadrilateral element
function element_stiffness(nodes, element, D)
    # Nodes for the current element
    x1, y1 = nodes[1, element[1]], nodes[2, element[1]]
    x2, y2 = nodes[1, element[2]], nodes[2, element[2]]
    x3, y3 = nodes[1, element[3]], nodes[2, element[3]]
    x4, y4 = nodes[1, element[4]], nodes[2, element[4]]

    # Shape function derivatives (evaluated at the center for simplicity)
    B = [
        y2 - y3  0       y3 - y1  0       y1 - y2  0       y2 - y4  0;
        0       x3 - x2  0       x1 - x3  0       x2 - x1  0       x4 - x2;
        x3 - x2  y2 - y3 x1 - x3  y3 - y1 x2 - x1  y1 - y2 x4 - x2  y2 - y4
    ]

    # Assuming area is constant for simplicity
    area = 0.5 * abs((x1*y2 + x2*y3 + x3*y4 + x4*y1) - (y1*x2 + y2*x3 + y3*x4 + y4*x1))

    # Element stiffness matrix
    Ke = area * B' * D * B
    return Ke
end

# Function to assemble the global stiffness matrix
function assemble_global_stiffness(nodes, elements, material)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness(nodes, element, D)
        # Global node indices
        dofs = vcat(2*element .- 1, 2*element) |> vec  # Degrees of freedom

        # Add element stiffness to global stiffness matrix
        for i in 1:8, j in 1:8
            K[dofs[i], dofs[j]] += Ke[i, j]
        end
    end

    return K
end

# Function to apply boundary conditions
function apply_boundary_conditions!(K, f, fixed_dofs)
    for dof in fixed_dofs
        K[dof, :] .= 0.0
        K[:, dof] .= 0.0
        K[dof, dof] = 1.0
        f[dof] = 0.0
    end
end

# Main solver function
function solve_fem(nodes, elements, material, forces, fixed_dofs)
    # Assemble global stiffness matrix
    K = assemble_global_stiffness(nodes, elements, material)

    # Force vector
    f = zeros(2 * size(nodes, 2))
    for (node, force) in forces
        f[2 * node - 1] = force[1]
        f[2 * node] = force[2]
    end

    # Apply boundary conditions
    apply_boundary_conditions!(K, f, fixed_dofs)

    # Solve for displacements
    u = K \ f

    return reshape(u, 2, :)
end

# Example usage
nodes = [0.0  1.0  1.0  0.0;
         0.0  0.0  1.0  1.0]
elements = [1; 2; 3; 4;;]
material = Material(200e9, 0.3)  # Steel properties (E, ν)
forces = Dict(3 => (0.0, -1000.0))  # Apply a downward force at node 3
fixed_dofs = [1, 2, 7, 8]  # Fix nodes 1 and 4 (x and y directions)

# Solve the problem
u = solve_fem(nodes, elements, material, forces, fixed_dofs)
println("Displacements: ", u)
