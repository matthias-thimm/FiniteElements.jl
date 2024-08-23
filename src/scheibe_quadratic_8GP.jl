using LinearAlgebra, SparseArrays

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

# 3x3 Gaussian quadrature points and weights for quadratic elements
const GAUSS_POINTS_3x3 = [-sqrt(3/5), 0.0, sqrt(3/5)]
const GAUSS_WEIGHTS_3x3 = [5/9, 8/9, 5/9]

# Shape functions for quadratic quadrilateral elements
function shape_functions_quadratic(ξ, η)
    N = [
        ξ*(ξ-1)*η*(η-1)/4,             # Node 1
        ξ*(ξ+1)*η*(η-1)/4,             # Node 2
        ξ*(ξ+1)*η*(η+1)/4,             # Node 3
        ξ*(ξ-1)*η*(η+1)/4,             # Node 4
        (1-ξ^2)*η*(η-1)/2,             # Node 5 (midpoint of side 1-2)
        ξ*(ξ+1)*(1-η^2)/2,             # Node 6 (midpoint of side 2-3)
        (1-ξ^2)*η*(η+1)/2,             # Node 7 (midpoint of side 3-4)
        ξ*(ξ-1)*(1-η^2)/2              # Node 8 (midpoint of side 4-1)
    ]
    return N
end

# Shape function derivatives for quadratic quadrilateral elements
function shape_function_derivatives_quadratic(ξ, η)
    dN_dξ = [
        (2*ξ-1)*η*(η-1)/4,             # Node 1
        (2*ξ+1)*η*(η-1)/4,             # Node 2
        (2*ξ+1)*η*(η+1)/4,             # Node 3
        (2*ξ-1)*η*(η+1)/4,             # Node 4
        -2*ξ*η*(η-1)/2,                # Node 5
        (2*ξ+1)*(1-η^2)/2,             # Node 6
        -2*ξ*η*(η+1)/2,                # Node 7
        (2*ξ-1)*(1-η^2)/2              # Node 8
    ]

    dN_dη = [
        ξ*(ξ-1)*(2*η-1)/4,             # Node 1
        ξ*(ξ+1)*(2*η-1)/4,             # Node 2
        ξ*(ξ+1)*(2*η+1)/4,             # Node 3
        ξ*(ξ-1)*(2*η+1)/4,             # Node 4
        (1-ξ^2)*(2*η-1)/2,             # Node 5
        ξ*(ξ+1)*(-2*η)/2,              # Node 6
        (1-ξ^2)*(2*η+1)/2,             # Node 7
        ξ*(ξ-1)*(-2*η)/2               # Node 8
    ]

    return dN_dξ, dN_dη
end

# Function to calculate the element stiffness matrix for a quadratic quadrilateral element using Gaussian quadrature
function element_stiffness_quadratic(nodes, element, D)
    Ke = zeros(16, 16)  # Initialize element stiffness matrix (16x16 for a quadratic quadrilateral element)

    for ξ in GAUSS_POINTS_3x3, η in GAUSS_POINTS_3x3
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives_quadratic(ξ, η)

        # Calculate the Jacobian matrix and its determinant
        J = [sum(dN_dξ[i] * nodes[1, element[i]] for i in 1:8)  sum(dN_dξ[i] * nodes[2, element[i]] for i in 1:8);
             sum(dN_dη[i] * nodes[1, element[i]] for i in 1:8)  sum(dN_dη[i] * nodes[2, element[i]] for i in 1:8)]

        detJ = det(J)  # Determinant of the Jacobian
        invJ = inv(J)  # Inverse of the Jacobian

        # Derivatives of shape functions with respect to global coordinates (x, y)
        dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
        dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

        # Construct the B matrix (strain-displacement matrix)
        B = zeros(3, 16)
        for i in 1:8
            B[1, 2*i-1] = dN_dx[i]
            B[2, 2*i]   = dN_dy[i]
            B[3, 2*i-1] = dN_dy[i]
            B[3, 2*i]   = dN_dx[i]
        end

        # Weighting factor (Jacobian determinant and Gauss weights)
        weight = detJ * GAUSS_WEIGHTS_3x3[findfirst(GAUSS_POINTS_3x3 .== ξ)] * GAUSS_WEIGHTS_3x3[findfirst(GAUSS_POINTS_3x3 .== η)]

        # Contribution to the element stiffness matrix
        Ke += weight * B' * D * B
    end

    return Ke
end

# Function to assemble the global stiffness matrix for quadratic elements
function assemble_global_stiffness_quadratic(nodes, elements, material)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness_quadratic(nodes, element, D)
        # Global node indices
        dofs = vcat(2*element .- 1, 2*element) |> vec  # Degrees of freedom

        # Add element stiffness to global stiffness matrix
        for i in 1:16, j in 1:16
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

# Main solver function for quadratic elements
function solve_fem_quadratic(nodes, elements, material, forces, fixed_dofs)
    # Assemble global stiffness matrix
    K = assemble_global_stiffness_quadratic(nodes, elements, material)

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

# Example usage for quadratic elements
nodes = [0.0  1.0  2.0  0.0  1.0  2.0  0.0  1.0;
         0.0  0.0  0.0  1.0  1.0  1.0  2.0  2.0]
elements = [1; 2; 3; 8; 4; 5; 6; 7;;]
material = Material(200e9, 0.3)  # Steel properties (E, ν)
forces = Dict(6 => (0.0, -1000.0))  # Apply a downward force at node 6
fixed_dofs = [1, 2, 13, 14]  # Fix nodes 1 and 7 (x and y directions)

# Solve the problem for quadratic elements
u = solve_fem_quadratic(nodes, elements, material, forces, fixed_dofs)
println("Displacements for quadratic elements: ", u)
