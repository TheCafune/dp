include("Util.jl")

using LightGraphs
using Laplacians: approxchol_sddm
using LinearAlgebra
using SparseArrays
using Random

function BestPageRank(G, k)
    pr = pagerank(G.g)
    A = [(i, pr[i]) for i = 1 : nv(G.g)]
    sort!(A, by=x->-x[2])
    return [A[i][1] for i = 1 : k]
end

function RandomSelect(G, k)
    Random.seed!(0)
    A = [(i, rand()) for i = 1 : nv(G.g)]
    sort!(A, by=x->-x[2])
    return [A[i][1] for i = 1 : k]
end

function BOMP(G, k)
    n = nv(G.g)
    S = Diagonal(G.nodeVal)
    Q = inv(laplacian_matrix(G.g) + I + zeros(n, n))
    R = Q * S
    tmp = ones(n)
    tar = R * tmp
    ans = []
    for rep = 1 : k
        id,bst = -1,-1e10
        for i = 1 : n
            if tmp[i]<0.5; continue; end
            t = dot(R[:,i], tar)
            if t>bst
                id,bst = i,t
            end
        end
        push!(ans, id)
        tmp[id] = 0
        tar = R * tmp
    end
    return ans
end

function FastGreedy(G, k; eps = 0.5)
    s = deepcopy(G.nodeVal)
    n = nv(G.g)
    f = approxchol_sddm(1.0*laplacian_matrix(G.g)+I)
    ans = []
    diag = zeros(n)
    kkk = round(Int, log2(n)/eps^2)
    for rep = 1 : kkk
        y = randn(n)
        y2 = f(y)
        for i = 1 : n
            diag[i] += y2[i]^2
        end
    end
    diag ./= kkk

    for ii = 1 : k
        s2 = f(f(s))
        id,bst = -1,-1e10
        for i = 1 : n
            if i in ans; continue; end
            t = 2*s[i]*s2[i] - s[i]^2*diag[i]
            if t>bst; id,bst = i,t; end
        end
        push!(ans, id)
        s[id] = 0.0
    end
    return ans
end