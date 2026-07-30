#
# Fast tests for the ap/ acyclic-partition machinery (and
# plots/plot_morse_strata.jl's non-drawing return values). Deliberately
# self-contained -- builds its own small complexes inline rather than
# depending on morse_vector_strata.jl's paper-specific example_* helpers,
# so this file can be copied as-is into ConleyDynamics's own test suite.
# Everything here runs in well under a second; the slow, large sweeps
# (disk13's 40232-element AP^c(X), any full unrestricted sweep on it)
# stay in construct_apx's driver/verify_examples, never here.
#

@testset "Acyclic partition construction" begin
    #
    # The five-cell path v1-e1-v2-e2-v3 (doc/morse-vector-strata.md
    # Section 6.2): f(X) = (3,2), beta(X) = (1,0). |AP^c(X)| = 16,
    # |AP(X)| = 45 (16 of which are connected).
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
    @test e.ntotal - e.ncovered == 5   # 5 of 9 uncovered (Section 6.2)
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
    # The 3-ball example (Section 6.5): f(X) = (2,3,3,1),
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
    @test [2, 2] in deltas   # the {P,Q} split of Section 6.2 realizes 2(e0+e1)

    cdeltas = connected_split_deltas(lc, X)
    @test all(r -> is_connected_block(lc, r.A) && is_connected_block(lc, r.B), cdeltas)
    @test length(cdeltas) <= length(block_split_deltas(lc, X))
end

@testset "construct_ap_stratum matches construct_ap_space" begin
    #
    # Mixed-degree six-cell complex (Section 6.3): X = Y ⊔ Z, an edge
    # union a disjoint 2/3-cell pair. f(X) = (2,1,1,1), |AP(X)| = 47,
    # |AP^c(X)| well under that. Small enough to fully enumerate both
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
