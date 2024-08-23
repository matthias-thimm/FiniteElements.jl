# 3x3 Gaussian quadrature points and weights for quadratic elements
const GAUSS_POINTS_3x3 = ((5/9, -sqrt(3/5)), (8/9, 0.0), (5/9, sqrt(3/5)))

# Shape functions for quadratic quadrilateral elements
function shape_functions_quadratic(ξ, η)
    N = SVector{8}(
        ξ*(ξ-1)*η*(η-1)/4,             # Node 1
        ξ*(ξ+1)*η*(η-1)/4,             # Node 2
        ξ*(ξ+1)*η*(η+1)/4,             # Node 3
        ξ*(ξ-1)*η*(η+1)/4,             # Node 4
        (1-ξ^2)*η*(η-1)/2,             # Node 5 (midpoint of side 1-2)
        ξ*(ξ+1)*(1-η^2)/2,             # Node 6 (midpoint of side 2-3)
        (1-ξ^2)*η*(η+1)/2,             # Node 7 (midpoint of side 3-4)
        ξ*(ξ-1)*(1-η^2)/2              # Node 8 (midpoint of side 4-1)
    )
    return N
end

# Shape function derivatives for quadratic quadrilateral elements
function shape_function_derivatives_quadratic(ξ, η)
    dN_dξ = SVector{8}(
        (2*ξ-1)*η*(η-1)/4,             # Node 1
        (2*ξ+1)*η*(η-1)/4,             # Node 2
        (2*ξ+1)*η*(η+1)/4,             # Node 3
        (2*ξ-1)*η*(η+1)/4,             # Node 4
        -2*ξ*η*(η-1)/2,                # Node 5
        (2*ξ+1)*(1-η^2)/2,             # Node 6
        -2*ξ*η*(η+1)/2,                # Node 7
        (2*ξ-1)*(1-η^2)/2              # Node 8
    )
    dN_dη = SVector{8}(
        ξ*(ξ-1)*(2*η-1)/4,             # Node 1
        ξ*(ξ+1)*(2*η-1)/4,             # Node 2
        ξ*(ξ+1)*(2*η+1)/4,             # Node 3
        ξ*(ξ-1)*(2*η+1)/4,             # Node 4
        (1-ξ^2)*(2*η-1)/2,             # Node 5
        ξ*(ξ+1)*(-2*η)/2,              # Node 6
        (1-ξ^2)*(2*η+1)/2,             # Node 7
        ξ*(ξ-1)*(-2*η)/2               # Node 8
    )

    return dN_dξ, dN_dη
end

# Function to calculate the element stiffness matrix for a quadratic quadrilateral element using Gaussian quadrature
function element_stiffness_quadratic(nodes, element, D, thickness)
    Ke = zeros(MMatrix{16,16})  # Initialize element stiffness matrix (16x16 for a quadratic quadrilateral element)
    node_positions_transposed = SMatrix{8,2}(nodes[:, element]')

    for (wξ, ξ) in GAUSS_POINTS_3x3, (wη, η) in GAUSS_POINTS_3x3
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives_quadratic(ξ, η)

        # Calculate the Jacobian matrix and its determinant
        # J = [sum(dN_dξ[i] * nodes[1, element[i]] for i in 1:8)  sum(dN_dξ[i] * nodes[2, element[i]] for i in 1:8);
        #      sum(dN_dη[i] * nodes[1, element[i]] for i in 1:8)  sum(dN_dη[i] * nodes[2, element[i]] for i in 1:8)]
        # J = SMatrix{2,2}(sum(dN_dξ[i] * nodes[1, element[i]] for i in 1:9),
        #                  sum(dN_dξ[i] * nodes[2, element[i]] for i in 1:9),
        #                  sum(dN_dη[i] * nodes[1, element[i]] for i in 1:9),
        #                  sum(dN_dη[i] * nodes[2, element[i]] for i in 1:9))
        J = SMatrix{2,8}(dN_dξ[1], dN_dη[1],
                         dN_dξ[2], dN_dη[2],
                         dN_dξ[3], dN_dη[3],
                         dN_dξ[4], dN_dη[4],
                         dN_dξ[5], dN_dη[5],
                         dN_dξ[6], dN_dη[6],
                         dN_dξ[7], dN_dη[7],
                         dN_dξ[8], dN_dη[8]) * node_positions_transposed

        detJ = det(J)  # Determinant of the Jacobian
        invJ = inv(J)  # Inverse of the Jacobian

        # Derivatives of shape functions with respect to global coordinates (x, y)
        dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
        dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

        # Construct the B matrix (strain-displacement matrix)
        # B = zeros(MMatrix{3,16})
        # for i in 1:8
        #     B[1, 2*i-1] = dN_dx[i]
        #     B[2, 2*i]   = dN_dy[i]
        #     B[3, 2*i-1] = dN_dy[i]
        #     B[3, 2*i]   = dN_dx[i]
        # end
        B = SMatrix{3,16}(dN_dx[1], 0.0, dN_dy[1],
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
                          0.0, dN_dy[8], dN_dx[8])

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
        dofs = zeros(Int, 16)
        dofs[1] = get_dof(element[1], :x)
        dofs[2] = get_dof(element[1], :y)
        dofs[3] = get_dof(element[2], :x)
        dofs[4] = get_dof(element[2], :y)
        dofs[5] = get_dof(element[3], :x)
        dofs[6] = get_dof(element[3], :y)
        dofs[7] = get_dof(element[4], :x)
        dofs[8] = get_dof(element[4], :y)
        dofs[9] = get_dof(element[5], :x)
        dofs[10] = get_dof(element[5], :y)
        dofs[11] = get_dof(element[6], :x)
        dofs[12] = get_dof(element[6], :y)
        dofs[13] = get_dof(element[7], :x)
        dofs[14] = get_dof(element[7], :y)
        dofs[15] = get_dof(element[8], :x)
        dofs[16] = get_dof(element[8], :y)

        # Add element stiffness to global stiffness matrix
        for i in 1:16, j in 1:16
            K[dofs[i], dofs[j]] += Ke[i, j]
        end
    end

    return K
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
