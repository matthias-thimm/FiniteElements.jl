module FiniteElements

using LinearAlgebra, SparseArrays, StaticArrays, GLMakie

export Material
export solve_fem_linear, solve_fem_quadratic, solve_fem_linear_stress, solve_fem_linear_stress_dynamic

include("autoinfiltrator.jl")
include("../scripts/plotting.jl")
include("material.jl")
include("utils.jl")
include("boundary_conditions.jl")
include("linear_elements.jl")
include("quadratic_elements.jl")
include("stress_calculation.jl")

end