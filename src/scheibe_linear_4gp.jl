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

# 2x2 Gaussian quadrature points and weights
const GAUSS_POINTS = [-1/sqrt(3), 1/sqrt(3)]
const GAUSS_WEIGHTS = [1.0, 1.0]

# Shape functions and derivatives in natural coordinates (ξ, η)
function shape_functions(ξ, η)
    N = [(1 - ξ)*(1 - η)/4,
         (1 + ξ)*(1 - η)/4,
         (1 + ξ)*(1 + η)/4,
         (1 - ξ)*(1 + η)/4]
    return N
end

function shape_function_derivatives(ξ, η)
    dN_dξ = [- (1 - η) / 4,  (1 - η) / 4,  (1 + η) / 4, - (1 + η) / 4]
    dN_dη = [- (1 - ξ) / 4, - (1 + ξ) / 4,  (1 + ξ) / 4,  (1 - ξ) / 4]
    return dN_dξ, dN_dη
end

# Function to calculate the element stiffness matrix for a quadrilateral element using Gaussian quadrature
function element_stiffness(nodes, element, D)
    Ke = zeros(8, 8)  # Initialize element stiffness matrix (8x8 for a quadrilateral element)

    for ξ in GAUSS_POINTS, η in GAUSS_POINTS
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives(ξ, η)

        # Calculate the Jacobian matrix and its determinant
        J = [sum(dN_dξ[i] * nodes[1, element[i]] for i in 1:4)  sum(dN_dξ[i] * nodes[2, element[i]] for i in 1:4);
             sum(dN_dη[i] * nodes[1, element[i]] for i in 1:4)  sum(dN_dη[i] * nodes[2, element[i]] for i in 1:4)]

        detJ = det(J)  # Determinant of the Jacobian
        invJ = inv(J)  # Inverse of the Jacobian

        # Derivatives of shape functions with respect to global coordinates (x, y)
        dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
        dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

        # Construct the B matrix (strain-displacement matrix)
        B = zeros(3, 8)
        for i in 1:4
            B[1, 2*i-1] = dN_dx[i]
            B[2, 2*i]   = dN_dy[i]
            B[3, 2*i-1] = dN_dy[i]
            B[3, 2*i]   = dN_dx[i]
        end

        # Weighting factor (Jacobian determinant and Gauss weights)
        weight = detJ * GAUSS_WEIGHTS[1] * GAUSS_WEIGHTS[2]

        # Contribution to the element stiffness matrix
        Ke += weight * B' * D * B
    end

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

## Example usage
nodes = [0.0  1.0  1.0  0.0;
         0.0  0.0  1.0  1.0]
nodes *= 0.001
elements = [1; 2; 3; 4;;]
material = Material(200e9, 0.3)  # Steel properties (E, ν)
forces = Dict(2 => (1e5, 0.0), 3 => (1e5, 0.0))
fixed_dofs = [1, 2, 7, 8]  # Fix nodes 1 and 4 (x and y directions)

# Solve the problem
u = solve_fem(nodes, elements, material, forces, fixed_dofs)
println("Displacements: ", u)
