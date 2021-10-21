mutable struct OpinionGraph
    nodeLab
    g
    nodeVal
end

using LightGraphs
using LinearAlgebra
using DelimitedFiles
using Random

function readGraph(gname)
    wjm = "../data/" * gname * ".txt"
    G = OpinionGraph(nothing, nothing, nothing)
    G.nodeLab = Dict()
    n = 0
    E = readdlm(wjm, ' ', Int, '\n')

    function add_node(x)
        if !haskey(G.nodeLab, x)
            n += 1
            G.nodeLab[x] = n
        end
    end

    for i = 1 : size(E, 1)
        u,v = E[i,1],E[i,2]
        add_node(u)
        add_node(v)
    end

    G.g = SimpleGraph(n)

    for i = 1 : size(E, 1)
        u,v = G.nodeLab[E[i,1]],G.nodeLab[E[i,2]]
        add_edge!(G.g, u, v)
    end

    return G
end

function readVal(G, gname)
    wjm = "../data/" * gname * "-val.txt"
    A = readdlm(wjm, ' ', Int, '\n')
    G.nodeVal = zeros(nv(G.g))
    for i = 1 : size(A, 1)
        if !haskey(G.nodeLab, A[i,1]); continue; end
        id, x = G.nodeLab[A[i,1]], A[i,2]
        G.nodeVal[id] = x
    end
end

function randVal(G)
    G.nodeVal = rand([-1.0, 0.0, 1.0], nv(G.g))
end

function dis(G, ans)
    x = deepcopy(G.nodeVal)
    for id in ans
        x[id] = 0
    end
    n = nv(G.g)
    W = inv(laplacian_matrix(G.g) + I + zeros(n, n))
    z = W * x
    return dot(x,z)
end