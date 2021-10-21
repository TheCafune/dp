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

function ExtGreedy(G, k)
    n,m = nv(G.g),ne(G.g)
    W = inv(laplacian_matrix(G.g) + I + zeros(n, n))
    s = deepcopy(G.nodeVal)
    ret = []
    dd = Dict()
    for rep = 1 : k
        z = W * s
        id,bst = -1,-1e10
        for i = 1 : n
            if haskey(dd, i); continue; end
            t = 2*s[i]*z[i] - s[i]*s[i]*W[i,i]
            if t>bst
                id,bst = i,t
            end
        end
        push!(ret, id)
        s[id] = 0
    end
    return ret
end

function FastGreedy(G, k; eps = 0.5)
    println(" INIT ")
    s = deepcopy(G.nodeVal)
    n,m = nv(G.g),ne(G.g)
    f = approxchol_sddm(1.0*laplacian_matrix(G.g)+I, tol=1e-6)
    f2 = approxchol_sddm(1.0*laplacian_matrix(G.g)+I, tol=1e-6)
    B = incidence_matrix(G.g, Float64, oriented=true)
    ans = []
    diag = zeros(n)
    kkk = round(Int, log2(n)/eps^2)
    println(" STAGE 1 $kkk")
    for rep = 1 : kkk
        println("  $rep --> 1")
        y = randn(n)
        println("  $rep --> 2")
        y2 = f(y)
        println("  $rep --> 3")
        tyy = randn(m)
        println(size(B), " $m")
        yy = B * tyy
        println("  $rep --> 4")
        yy2 = f2(yy)
        println("  $rep --> 5")
        for i = 1 : n
            diag[i] += y2[i]^2 + yy2[i]^2
        end
    end
    diag ./= kkk

    println(" STAGE 2 ")
    for ii = 1 : k
        s2 = f(s)
        id,bst = -1,-1e10
        for i = 1 : n
            if i in ans; continue; end
            t = 2*s[i]*s2[i] - s[i]^2*diag[i]
            if t>bst; id,bst = i,t; end
        end
        println(" $id $bst")
        push!(ans, id)
        s[id] = 0.0
    end
    println(" FINISHED ")
    return ans
end