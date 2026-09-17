#
# Fast tests for the ap/ acyclic-partition machinery (and
# plots/plot_morse_strata.jl's non-drawing return values). Deliberately
# self-contained -- builds its own small complexes inline rather than
# depending on morse_vector_strata.jl's paper-specific example_* helpers,
# so this file can be copied as-is into ConleyDynamics's own test suite.
# Everything here runs in well under a second; the slow, large sweeps
# (disk13's 40232-element AP(X), any full unrestricted sweep on it)
# stay in construct_apx's driver/verify_examples, never here.
#

@testset "Acyclic partition construction" begin
    #
    # The five-cell path v1-e1-v2-e2-v3: f(X) = (3,2), beta(X) = (1,0).
    # |AP(X)| = 16, |AP^d(X)| = 45 (16 of which are connected).
    #
    labels    = ["v1", "v2", "v3"]
    simplices = [[1, 2], [2, 3]]
    lc = create_simplicial_complex(labels, simplices)

    apc = construct_ap_space(lc)                    # connected=true default
    apf = construct_ap_space(lc; connected=false)

    @test length(apc) == 16
    @test length(apf) == 45
    @test all(w -> is_connected_partition(lc, w), apc)
    @test !all(w -> is_connected_partition(lc, w), apf)

    X  = collect(1:lc.ncells)
    M2 = conley_index(lc, X)
    @test M2 == [1, 0]

    stratac = stratum_partition(lc, apc)
    strataf = stratum_partition(lc, apf)
    @test length(strataf[M2]) == 9
    @test length(stratac[M2]) < length(strataf[M2])

    adjf = stratum_adjacency(lc, apf)
    M1 = [3, 2]   # f(X): the finest stratum, every cell its own block
    e  = only(filter(x -> x.M1 == M1 && x.M2 == M2, adjf))
    @test e.weight == 2
    @test e.ntotal - e.ncovered == 5   # 5 of 9 uncovered
end

@testset "is_atomic_refinement and atomic_distances" begin
    labels    = ["v1", "v2", "v3"]
    simplices = [[1, 2], [2, 3]]
    lc = create_simplicial_complex(labels, simplices)
    v1, v2, v3 = lc.indices["v1"], lc.indices["v2"], lc.indices["v3"]
    e1, e2     = lc.indices["v1v2"], lc.indices["v2v3"]

    V2 = [collect(1:lc.ncells)]         # {X}, coarsest
    V1 = [[v1, v2, v3], [e1, e2]]       # one atomic split of X

    @test is_atomic_refinement(lc, V1, V2)
    @test !is_atomic_refinement(lc, V2, V1)

    ap = construct_ap_space(lc; connected=false)
    A  = atomic_distances(lc, ap)
    i2 = findfirst(w -> convert_mvf_partition(lc, w) == convert_mvf_partition(lc, V2), ap)
    i1 = findfirst(w -> convert_mvf_partition(lc, w) == convert_mvf_partition(lc, V1), ap)
    @test i1 in sparse_get_nz_row(A, i2)
end

@testset "Morse vector and rank bookkeeping" begin
    #
    # The 3-ball example: f(X) = (2,3,3,1),
    # rho(X) = (1,2,1), beta(X) = (1,0,0,0).
    #
    defcellbnd = [
        ["v1", 0], ["v2", 0],
        ["e1", 1, "v2", "v1"], ["e2", 1, "v2", "v1"], ["e3", 1, "v2", "v1"],
        ["t1", 2, "e1", "e2"], ["t2", 2, "e2", "e3"], ["t3", 2, "e3", "e1"],
        ["u",  3, "t1", "t2", "t3"],
    ]
    lc = create_lefschetz_gf2(defcellbnd)
    X  = collect(1:lc.ncells)

    @test block_rho(lc, X)[2:4] == [1, 2, 1]
    @test morse_vector(lc, [X]) == [1, 0, 0, 0]

    v1, v2 = lc.indices["v1"], lc.indices["v2"]
    e1     = lc.indices["e1"]
    @test conley_index(lc, [v1, e1]) == [0, 0, 0, 0]   # a regular pair: non-critical
end

@testset "Connectivity" begin
    labels    = ["v1", "v2", "v3"]
    simplices = [[1, 2], [2, 3]]
    lc = create_simplicial_complex(labels, simplices)
    v1, v2, v3 = lc.indices["v1"], lc.indices["v2"], lc.indices["v3"]
    e1, e2     = lc.indices["v1v2"], lc.indices["v2v3"]

    @test is_connected_block(lc, [v1, e1, v2])
    @test !is_connected_block(lc, [v1, v3])   # not adjacent, no shared edge
    @test is_connected_partition(lc, [[v1, e1], [v2, e2]])
    @test !is_connected_partition(lc, [[v1, v3]])
end

@testset "Split spectrum" begin
    labels    = ["v1", "v2", "v3"]
    simplices = [[1, 2], [2, 3]]
    lc = create_simplicial_complex(labels, simplices)
    X  = collect(1:lc.ncells)

    deltas = split_spectrum(lc, X)
    @test [2, 2] in deltas   # the {P,Q} split realizes 2(e0+e1)

    cdeltas = connected_split_deltas(lc, X)
    @test all(r -> is_connected_block(lc, r.A) && is_connected_block(lc, r.B), cdeltas)
    @test length(cdeltas) <= length(block_split_deltas(lc, X))
end

@testset "construct_ap_stratum matches construct_ap_space" begin
    #
    # Mixed-degree six-cell complex: X = Y ⊔ Z, an edge
    # union a disjoint 2/3-cell pair. f(X) = (2,1,1,1), |AP^d(X)| = 47,
    # |AP(X)| well under that. Small enough to fully enumerate both
    # ways and cross-check construct_ap_stratum against the ground truth
    # for two different targets: the global-minimum Morse vector (where
    # pruning provably cannot help, per construct_ap_stratum's docstring)
    # and a non-minimal one (where it should).
    #
    defcellbnd = [["v1", 0], ["v2", 0], ["e", 1, "v2", "v1"],
                  ["c2", 2], ["c3", 3, "c2"]]
    lc = create_lefschetz_gf2(defcellbnd)
    X  = collect(1:lc.ncells)

    for connected in (true, false)
        ap = construct_ap_space(lc; connected=connected)
        strata = stratum_partition(lc, ap)

        for M in keys(strata)
            truth = Set(convert_mvf_partition(lc, ap[i]) for i in strata[M])
            got   = Set(convert_mvf_partition(lc, w)
                        for w in construct_ap_stratum(lc, M; connected=connected))
            @test got == truth
        end
    end
end

@testset "Eventual reachability (stratum_reachability)" begin
    #
    # Same 5-cell path as above; the unrestricted AP^d(X) (45 elements) has
    # more room for multi-step chains than the connected AP(X) (16).
    #
    labels    = ["v1", "v2", "v3"]
    simplices = [[1, 2], [2, 3]]
    lc = create_simplicial_complex(labels, simplices)

    ap = construct_ap_space(lc; connected=false)
    A  = atomic_distances(lc, ap)
    Ms = map(t -> morse_vector(lc, t), ap)

    reach, edges = stratum_reachability(lc, ap; A=A)

    @test length(reach) == length(ap)
    @test all(i -> Ms[i] in reach[i], eachindex(ap))   # reflexive

    # Every direct one-step successor's Morse vector must be included
    @test all(i -> all(w -> Ms[w] in reach[i], sparse_get_nz_row(A, i)), eachindex(ap))

    # Independent brute-force fixed-point closure, cross-checked against
    # the DP-based implementation
    brute = [Set([Ms[i]]) for i in eachindex(ap)]
    changed = true
    while changed
        changed = false
        for i in eachindex(ap), w in sparse_get_nz_row(A, i)
            for M in brute[w]
                if !(M in brute[i])
                    push!(brute[i], M)
                    changed = true
                end
            end
        end
    end
    @test reach == brute

    # Every stratum_adjacency (one-step) edge must also show up as
    # eventually-covered, with ncovered no smaller
    adj = stratum_adjacency(lc, ap; A=A)
    for e in adj
        er = only(filter(x -> x.M1 == e.M1 && x.M2 == e.M2, edges))
        @test er.ncovered >= e.ncovered
    end

    # The transitive reduction used by plot_morse_reachability must drop
    # only genuinely redundant edges: the reachability closure computed
    # from the reduced edge set must match the closure from the full set.
    reduced = ConleyDynamics._stratum_transitive_reduction(edges)
    @test length(reduced) <= length(edges)

    function _closure(es)
        succ = Dict{Vector{Int},Set{Vector{Int}}}()
        for e in es
            push!(get!(() -> Set{Vector{Int}}(), succ, e.M2), e.M1)
        end
        nodes = union(Set(e.M1 for e in es), Set(e.M2 for e in es))
        clo = Dict(M => Set{Vector{Int}}([M]) for M in nodes)
        changed = true
        while changed
            changed = false
            for M in nodes, w in get(succ, M, Set{Vector{Int}}())
                for N in clo[w]
                    if !(N in clo[M])
                        push!(clo[M], N)
                        changed = true
                    end
                end
            end
        end
        return clo
    end

    @test _closure(reduced) == _closure(edges)
end

@testset "_stratum_transitive_reduction" begin
    A, B, C, D = [0, 0], [1, 0], [2, 0], [0, 1]

    # A->B->C plus a redundant direct A->C: the direct edge must be dropped
    edges = [(M1=B, M2=A), (M1=C, M2=B), (M1=C, M2=A)]
    red = ConleyDynamics._stratum_transitive_reduction(edges)
    @test Set((e.M2, e.M1) for e in red) == Set([(A, B), (B, C)])

    # Diamond A->B->D, A->C->D with no direct A->D edge: nothing is
    # redundant here, so every edge must be kept
    edges2 = [(M1=B, M2=A), (M1=D, M2=B), (M1=C, M2=A), (M1=D, M2=C)]
    red2 = ConleyDynamics._stratum_transitive_reduction(edges2)
    @test Set((e.M2, e.M1) for e in red2) == Set((e.M2, e.M1) for e in edges2)
end

@testset "Top layer of a Morse stratum" begin
    #
    # Cross-check is_top_layer, for every element of a small ap, against
    # the global (poset-level) definition: v is in the top layer of its
    # own stratum iff none of its one-step atomic refinements (from
    # atomic_distances) shares its Morse vector.
    #
    labels    = ["v1", "v2", "v3"]
    simplices = [[1, 2], [2, 3]]
    lc = create_simplicial_complex(labels, simplices)

    for connected in (true, false)
        ap = construct_ap_space(lc; connected=connected)
        A  = atomic_distances(lc, ap)
        Ms = map(t -> morse_vector(lc, t), ap)

        for i in eachindex(ap)
            brute = !any(j -> Ms[j] == Ms[i], sparse_get_nz_row(A, i))
            @test is_top_layer(lc, ap[i]; connected=connected) == brute
        end

        strata = stratum_partition(lc, ap)
        stl    = stratum_top_layer(lc, ap; connected=connected)
        for M in keys(strata)
            @test stl[M] == filter(i -> is_top_layer(lc, ap[i]; connected=connected), strata[M])
        end

        # construct_ap_stratum_top matches stratum_top_layer for one target
        M0 = first(keys(strata))
        ctop = Set(convert_mvf_partition(lc, w)
                   for w in construct_ap_stratum_top(lc, M0; connected=connected))
        stop = Set(convert_mvf_partition(lc, ap[i]) for i in stl[M0])
        @test ctop == stop
    end

    # Every Forman vector field is a facet of its stratum, hence in the
    # top layer (Proposition on invariants of the top layer)
    lcf, mvf = example_forman1d()
    @test is_top_layer(lcf, mvf)

    #
    # Twisted 4-cycle over GF(3): u1,u2 (dim 0), w1,w2 (dim 1),
    # dw1=u1+u2, dw2=u1-u2. beta(X)=(0,0), and AP_(0,0)(X) = L_(0,0) =
    # {{X}} -- a single, genuinely rigid partition that is not a Forman
    # vector field (its mvf_length is 1, not the card bound 2).
    #
    labels4 = ["u1", "u2", "w1", "w2"]
    dims4   = [0, 0, 1, 1]
    r4 = [1, 2, 1, 2]
    c4 = [3, 3, 4, 4]
    v4 = Rational{Int}[1, 1, 1, -1]
    bnd4 = sparse_from_lists(4, 4, 0, Rational{Int}(0), Rational{Int}(1), r4, c4, v4)
    lc4 = lefschetz_gfp_conversion(LefschetzComplex(labels4, dims4, bnd4), 3)

    X4 = collect(1:4)
    @test conley_index(lc4, X4) == [0, 0]

    ap00 = construct_ap_stratum(lc4, [0, 0])
    @test length(ap00) == 1
    @test is_top_layer(lc4, ap00[1])
    @test mvf_length(lc4, ap00[1]) < (lc4.ncells + sum([0, 0])) / 2   # not Forman

    #
    # Hex-arrow example: u1,u2,u3,x (dim 0), w1,w2,w3,y (dim 1),
    # dw1=u1-u2, dw2=u2-u3, dw3=u3+u1, dy=x-u3, over GF(3). The twisted
    # hexagon H={u1,u2,u3,w1,w2,w3} is regular and genuinely rigid, and
    # the top layer of AP_(1,1)(X) contains partitions of two different
    # sizes: {H,{x},{y}} (3 blocks) and Forman vector fields (5 blocks).
    #
    labels8 = ["u1", "u2", "u3", "x", "w1", "w2", "w3", "y"]
    dims8   = [0, 0, 0, 0, 1, 1, 1, 1]
    r8 = [1, 2, 2, 3, 3, 1, 4, 3]
    c8 = [5, 5, 6, 6, 7, 7, 8, 8]
    v8 = Rational{Int}[1, -1, 1, -1, 1, 1, 1, -1]
    bnd8 = sparse_from_lists(8, 8, 0, Rational{Int}(0), Rational{Int}(1), r8, c8, v8)
    lc8 = lefschetz_gfp_conversion(LefschetzComplex(labels8, dims8, bnd8), 3)

    H = [1, 2, 3, 5, 6, 7]
    @test conley_index(lc8, H) == [0, 0]
    @test is_connected_block(lc8, H)

    top11 = construct_ap_stratum_top(lc8, [1, 1])
    sizes = sort(unique(mvf_length.(Ref(lc8), top11)))
    @test sizes == [3, 5]
    @test any(w -> convert_mvf_partition(lc8, w) == convert_mvf_partition(lc8, [H]), top11)
end
