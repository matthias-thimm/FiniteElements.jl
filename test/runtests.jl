using FiniteElements
using Test

@testset "FiniteElements.jl" begin
    @testset "one element linear quads" begin
        # Problem setup
        L = 1.0        # Length of the element (in meters)
        H = 1.0        # Height of the element (in meters)
        applied_force = 1e8  # Total force applied in the right direction (in Newtons)
        E = 200e9      # Young's modulus (e.g., steel) in Pascals (N/m^2)
        thickness = 0.1  # Thickness of the structure (in meters)

        # Define the nodes
        nodes = [0.0  L  L  0.0;
                 0.0  0.0  H  H]

        # Define elements (single quadrilateral element)
        elements = [1; 2; 3; 4]

        # Material properties
        material = Material(E, 0.3, thickness)  # Steel-like material with thickness

        # Fixed DOFs (nodes on the left)
        # Fix nodes 1 (x and y directions) and 4 (x direction)
        fixed_dofs = Dict(1 => :xy, 4 => :x)

        # Distributed force on the right side (nodes 2 and 3)
        force_per_node = applied_force / 2.0
        forces = Dict(2 => (force_per_node, 0.0),
                      3 => (force_per_node, 0.0))

        u, K, f = solve_fem_linear(nodes, elements, material, forces, fixed_dofs)

        # Analytical displacement for comparison
        analytical_displacement(F, L, E, A) = F * L / (E * A)

        # Analytical solution
        A = thickness * H  # Cross-sectional area (in square meters)
        u_analytical = analytical_displacement(applied_force, L, E, A)

        @test u_analytical ≈ u[1, 2]
        @test u_analytical ≈ u[1, 3]
    end
end
