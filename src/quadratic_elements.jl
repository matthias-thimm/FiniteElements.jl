# 3x3 Gaussian quadrature points and weights for quadratic elements
const GAUSS_POINTS_3x3 = ((5 / 9, -sqrt(3 / 5)), (8 / 9, 0.0), (5 / 9, sqrt(3 / 5)))

# Shape functions for quadratic quadrilateral elements in natural coordinates (ξ, η):
#    N1 = 0.25 * ξ * η * (ξ - 1) * (η - 1)
#    N2 = 0.25 * ξ * η * (ξ + 1) * (η - 1)
#    N3 = 0.25 * ξ * η * (ξ + 1) * (η + 1)
#    N4 = 0.25 * ξ * η * (ξ - 1) * (η + 1)
#    N5 = -0.5 * (ξ + 1) * (ξ - 1) * η * (η - 1)
#    N6 = -0.5 * ξ * (ξ + 1) * (η + 1) * (η - 1)
#    N7 = -0.5 * (ξ + 1) * (ξ - 1) * η * (η + 1)
#    N8 = -0.5 * ξ * (ξ - 1) * (η + 1) * (η - 1)
#    N9 = (ξ + 1) * (ξ - 1) * (η + 1) * (η - 1)
#
#    ◍----◍----◍
#   4|   7    3|
#    ◍    ◍    ◍
#   8|   9    6|
#    ◍----◍----◍
#   1    5    2
function shape_function_derivatives_quadratic(ξ, η)
    dN_dξ = SVector{9}(
        0.25 * η * (η - 1) * (2 * ξ - 1),                    # dN1/dξ
        0.25 * η * (η - 1) * (2 * ξ + 1),                    # dN2/dξ
        0.25 * η * (η + 1) * (2 * ξ + 1),                    # dN3/dξ
        0.25 * η * (η + 1) * (2 * ξ - 1),                    # dN4/dξ
        -η * (η - 1) * ξ,                                    # dN5/dξ (Mid-side on edge 1-2)
        -0.5 * (η - 1) * (η + 1) * (2 * ξ + 1),              # dN6/dξ (Mid-side on edge 2-3)
        -η * (η + 1) * ξ,                                    # dN7/dξ (Mid-side on edge 3-4)
        -0.5 * (η - 1) * (η + 1) * (2 * ξ - 1),              # dN8/dξ (Mid-side on edge 4-1)
        2 * ξ * (η + 1) * (η - 1)                            # dN9/dξ (Center node)
    )
    dN_dη = SVector{9}(
        0.25 * ξ * (ξ - 1) * (2 * η - 1),                    # dN1/dη
        0.25 * ξ * (ξ + 1) * (2 * η - 1),                    # dN2/dη
        0.25 * ξ * (ξ + 1) * (2 * η + 1),                    # dN3/dη
        0.25 * ξ * (ξ - 1) * (2 * η + 1),                    # dN4/dη
        -0.5 * (ξ - 1) * (ξ + 1) * (2 * η - 1),              # dN5/dη (Mid-side on edge 1-2)
        -ξ * (ξ + 1) * η,                                    # dN6/dη (Mid-side on edge 2-3)
        -0.5 * (ξ - 1) * (ξ + 1) * (2 * η + 1),              # dN7/dη (Mid-side on edge 3-4)
        -ξ * (ξ - 1) * η,                                    # dN8/dη (Mid-side on edge 4-1)
        2 * η * (ξ + 1) * (ξ - 1)                            # dN9/dη (Center node)
    )
    return dN_dξ, dN_dη
end

# Function to calculate the element stiffness matrix for a quadratic quadrilateral element
# using Gaussian quadrature
function element_stiffness_quadratic(nodes, element, D, thickness)
    # Initialize element stiffness matrix (18x18 for a quadratic quadrilateral element)
    Ke = zeros(MMatrix{18,18})
    node_positions_transposed = @views SMatrix{9,2}(nodes[:, element]')

    for (wξ, ξ) in GAUSS_POINTS_3x3, (wη, η) in GAUSS_POINTS_3x3
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives_quadratic(ξ, η)

        # Calculate the Jacobian matrix and its determinant
        J = SMatrix{2,9}(dN_dξ[1], dN_dη[1],
                         dN_dξ[2], dN_dη[2],
                         dN_dξ[3], dN_dη[3],
                         dN_dξ[4], dN_dη[4],
                         dN_dξ[5], dN_dη[5],
                         dN_dξ[6], dN_dη[6],
                         dN_dξ[7], dN_dη[7],
                         dN_dξ[8], dN_dη[8],
                         dN_dξ[9], dN_dη[9]) * node_positions_transposed
        detJ = det(J)  # Determinant of the Jacobian
        invJ = inv(J)  # Inverse of the Jacobian

        # Derivatives of shape functions with respect to global coordinates (x, y)
        dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
        dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

        # Construct the B matrix (strain-displacement matrix)
        B = SMatrix{3,18}(dN_dx[1], 0.0, dN_dy[1],
                          0.0, dN_dy[1], dN_dx[1],
                          dN_dx[2], 0.0, dN_dy[2],
                          0.0, dN_dy[2], dN_dx[2],
                          dN_dx[3], 0.0, dN_dy[3],
                          0.0, dN_dy[3], dN_dx[3],
                          dN_dx[4], 0.0, dN_dy[4],
                          0.0, dN_dy[4], dN_dx[4],
                          dN_dx[5], 0.0, dN_dy[5],
                          0.0, dN_dy[5], dN_dx[5],
                          dN_dx[6], 0.0, dN_dy[6],
                          0.0, dN_dy[6], dN_dx[6],
                          dN_dx[7], 0.0, dN_dy[7],
                          0.0, dN_dy[7], dN_dx[7],
                          dN_dx[8], 0.0, dN_dy[8],
                          0.0, dN_dy[8], dN_dx[8],
                          dN_dx[9], 0.0, dN_dy[9],
                          0.0, dN_dy[9], dN_dx[9])

        # Weighting factor (Jacobian determinant and Gauss weights)
        weight = detJ * wξ * wη

        # Contribution to the element stiffness matrix
        Ke += thickness * weight * B' * D * B
    end

    return Ke
end

# Function to assemble the global stiffness matrix for quadratic elements
function assemble_global_stiffness_quadratic(nodes, elements, material)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness_quadratic(nodes, element, D, material.thickness)

        # Global node indices
        dofs = get_dofs_quadratic(element)

        # Add element stiffness to global stiffness matrix
        for i in 1:18, j in 1:18
            K[dofs[i], dofs[j]] += Ke[i, j]
        end
    end

    return K
end

function get_dofs_quadratic(element)
    return SVector{18}(get_dof(element[1], :x), get_dof(element[1], :y),
                       get_dof(element[2], :x), get_dof(element[2], :y),
                       get_dof(element[3], :x), get_dof(element[3], :y),
                       get_dof(element[4], :x), get_dof(element[4], :y),
                       get_dof(element[5], :x), get_dof(element[5], :y),
                       get_dof(element[6], :x), get_dof(element[6], :y),
                       get_dof(element[7], :x), get_dof(element[7], :y),
                       get_dof(element[8], :x), get_dof(element[8], :y),
                       get_dof(element[9], :x), get_dof(element[9], :y))
end

# Main solver function for quadratic elements
function solve_fem_quadratic(nodes, elements, material, forces, fixed_dofs)
    # Assemble global stiffness matrix
    K = assemble_global_stiffness_quadratic(nodes, elements, material)

    # Force vector
    f = get_force_vector(nodes, forces)

    # Apply boundary conditions
    apply_boundary_conditions!(K, f, fixed_dofs)

    # Solve for displacements
    u = K \ f

    return reshape(u, 2, :), K, f
end
