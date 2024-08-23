using LinearAlgebra, SparseArrays

# Define the material properties (Plane stress)
struct Material
    E::Float64  # Young's modulus
    ν::Float64  # Poisson's ratio
    thickness::Float64  # Thickness of the structure
end

function plane_stress_stiffness(material::Material)
    (; E, ν) = material
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
function element_stiffness(nodes, element, D, thickness)
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

    Ke *= thickness

    return Ke
end

# Function to assemble the global stiffness matrix
function assemble_global_stiffness(nodes, elements, material)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness(nodes, element, D, material.thickness)
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

# Create nodes, elements, and boundary conditions for a single-element tension test
function create_single_element_test(L, H, applied_force, E, thickness)
    # Define nodes (only one element, so four nodes)
    nodes = [0.0  L  L  0.0;
             0.0  0.0  H  H]

    # Define elements (single quadrilateral element)
    elements = [1; 2; 3; 4]

    # Material properties
    material = Material(E, 0.3, thickness)  # Steel-like material with thickness

    # Fixed DOFs (nodes on the left)
    fixed_dofs = [1, 2, 7, 8]  # Fix nodes 1 and 4 (x and y directions)

    # Distributed force on the right side (nodes 2 and 3)
    force_per_node = applied_force / 2.0
    forces = Dict(2 => (force_per_node, 0.0),
                  3 => (force_per_node, 0.0))

    return nodes, elements, forces, fixed_dofs, material
end

# Analytical displacement for comparison
function analytical_displacement(F, L, E, A)
    return F * L / (E * A)
end

# Main validation function
function validate_single_element_tension()
    # Problem setup
    L = 1.0        # Length of the element (in meters)
    H = 1.0        # Height of the element (in meters)
    applied_force = 100000000.0  # Total force applied in the right direction (in Newtons)
    E = 200e9      # Young's modulus (e.g., steel) in Pascals (N/m^2)
    thickness = 0.1  # Thickness of the structure (in meters)

    # Create simulation data
    nodes, elements, forces, fixed_dofs, material = create_single_element_test(L, H, applied_force, E, thickness)

    # Solve the FEM problem
    u = solve_fem(nodes, elements, material, forces, fixed_dofs)

    # Analytical solution
    A = thickness * H  # Cross-sectional area (in square meters)
    u_analytical = analytical_displacement(applied_force, L, E, A)

    # Check numerical displacement
    u_numeric_2 = u[:, 2]  # Displacement of node 2 (right-bottom)
    u_numeric_3 = u[:, 3]  # Displacement of node 3 (right-top)

    println("Numerical displacement at node 2: ", u_numeric_2)
    println("Numerical displacement at node 3: ", u_numeric_3)
    println("Analytical displacement: ", u_analytical)

    # Validate that both right nodes have the same displacement (tension test)
    @assert isapprox(u_numeric_2[1], u_numeric_3[1], atol=1e-6) "Displacements at nodes 2 and 3 are not equal"

    # Validate that the numerical displacement matches the analytical displacement
    @assert isapprox(u_numeric_2[1], u_analytical, atol=1e-6) "Numerical displacement does not match analytical solution"

    println("Validation successful: The numerical solution matches the analytical solution.")
end

# Run the validation
validate_single_element_tension()
