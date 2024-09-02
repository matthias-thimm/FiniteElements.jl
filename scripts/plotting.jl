using GLMakie
import GLMakie.GeometryBasics: QuadFace, normal_mesh

function plot_mesh!(ax, nodes, elements; kwargs...)
    node_position = [Point2f(n[1], n[2]) for n in eachcol(nodes)]
    quad_connectivity = [QuadFace(e[1], e[2], e[3], e[4]) for e in eachcol(elements)]
    triangle_mesh = normal_mesh(node_position, quad_connectivity)
    mesh!(ax, triangle_mesh; kwargs...)
    return nothing
end

function plot_edges_linear!(ax, nodes, elements; kwargs...)
    node_position = [Point2f(n[1], n[2]) for n in eachcol(nodes)]
    edges = Vector{Point2f}()
    for i in axes(elements, 2)
        nid1, nid2, nid3, nid4 = elements[1,i], elements[2,i], elements[3,i], elements[4,i]
        n1, n2 = node_position[nid1], node_position[nid2]
        n3, n4 = node_position[nid3], node_position[nid4]
        push!(edges, n1, n2, n2, n3, n3, n4, n4, n1)
    end
    linesegments!(ax, edges; kwargs...)
    return nothing
end

function plot_edges_quadratic!(ax, nodes, elements; kwargs...)
    node_position = [Point2f(n[1], n[2]) for n in eachcol(nodes)]
    edges = Vector{Point2f}()
    for i in axes(elements, 2)
        nid1, nid2, nid3, nid4 = elements[1,i], elements[2,i], elements[3,i], elements[4,i]
        nid5, nid6, nid7, nid8 = elements[5,i], elements[6,i], elements[7,i], elements[8,i]
        n1, n2 = node_position[nid1], node_position[nid2]
        n3, n4 = node_position[nid3], node_position[nid4]
        n5, n6 = node_position[nid5], node_position[nid6]
        n7, n8 = node_position[nid7], node_position[nid8]
        push!(edges, n1, n5, n5, n2, n2, n6, n6, n3, n3, n7, n7, n4, n4, n8, n8, n1)
    end
    linesegments!(ax, edges; kwargs...)
    return nothing
end

function plot_nodes!(ax, nodes; kwargs...)
    scatter!(ax, nodes; kwargs...)
    return nothing
end
