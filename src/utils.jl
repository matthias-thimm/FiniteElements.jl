@inline get_dof(node::Int, dim::Symbol) = get_dof(node, Val(dim))
@inline function get_dof(node::Int, dim::Int)
    if dim == 1
        dim_symb = :x
    elseif dim == 2
        dim_symb = :y
    else
        throw(ArgumentError("unknown dimension number $(dim)!\n"))
    end
    return get_dof(node::Int, dim_symb)
end
@inline get_dof(node::Int, ::Val{:x}) = 2 * node - 1
@inline get_dof(node::Int, ::Val{:y}) = 2 * node
function get_dof(node::Int, ::Val{N}) where {N}
    return throw(ArgumentError("unknown Symbol $(N)!\n"))
end
