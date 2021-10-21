include("Util.jl")
include("Alg.jl")

function do_exp(G, maxk, gname, method)
    wjm = "../results/"*gname*method*".txt"
    ori = pol(G, [])
    out = open(wjm, "w")
    println(out, "0 0")
    for k = 1 : maxk
        ans = []
        if method=="Random"
            ans = RandomSelect(G, k)
        elseif method=="PageRank"
            ans = BestPageRank(G, k)
        elseif method=="BOMP"
            ans = BOMP(G, k)
        else
            ans = FastGreedy(G, k)
        end
        y = pol(G, ans)-ori
        println(out, "$k $y")
    end
    close(out)
end

function run_exp_1(gname, maxk)
    G = readGraph(gname)
    readVal(G, gname)
    for method in ["Random", "PageRank", "BOMP", "FastGreedy"]
        do_exp(G, maxk, gname, method)
    end
end

function run_exp_2(gname; p=0.01)
    G = readGraph(gname)
    randVal(G)
    n,m = nv(G.g),ne(G.g)
    k = round(Int, n*p)
    lg = open("../results/log.txt", "a")
    println(lg, "$gname $n $m")
    if n>50000
        t1 = time()
        ans_fast = FastGreedy(G, k)
        t1 = time()-t1
        println("FastGreedy Time : $t1")
    else
        ori = pol(G, [])
        t1 = time()
        ans_fast = FastGreedy(G, k)
        t1 = time()-t1
        score_fast = pol(G, ans_fast)-ori
        t2 = time()
        ans_bomp = BOMP(G, k)
        t2 = time()-t2
        score_bomp = pol(G, ans_bomp)-ori
        err = abs(score_fast-score_bomp)/abs(score_bomp) * 100
        println(lg, "FastGreedy Time : $t1")
        println(lg, "FastGreedy Score : $score_fast")
        println(lg, "BOMP Time : $t2")
        println(lg, "BOMP Score : $score_bomp")
        println(lg, "Error : $err")
    end
    close(lg)
end

#run_exp_1("karate", 20)
#run_exp_1("books", 60)
#run_exp_1("Elections", 300)
#run_exp_1("polblogs", 200)

nameList = ["EmailUniv", "EmailUniv", "Yeast", "Hamster", "GrQc", "Erdos992", "PagesGovernment", "AstroPh", "CondMat", "Gplus", "GemsecRO",
"WikiTalk", "Gowalla", "GooglePlus", "MathSciNet", "Flickr", "IMDB", "YoutubeSnap", "Flixster"]
foreach(x -> run_exp_2(x), nameList)