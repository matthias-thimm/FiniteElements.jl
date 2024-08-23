module FiniteElements

using LinearAlgebra, SparseArrays, StaticArrays

export Material
export solve_fem_linear, solve_fem_quadratic

include("material.jl")
include("utils.jl")
include("boundary_conditions.jl")
include("linear_elements.jl")
include("quadratic_elements.jl")

end
