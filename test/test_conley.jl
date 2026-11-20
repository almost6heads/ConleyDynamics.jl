@testset "Connection matrix algorithms agree" begin
    # Use gradient variant 1 of example_three_cm
    lc, mvf  = example_three_cm(1)

    cm_dlms  = connection_matrix(lc, mvf, algorithm="DLMS")
    cm_dhl   = connection_matrix(lc, mvf, algorithm="DHL")
    cm_hms   = connection_matrix(lc, mvf, algorithm="HMS")
    cm_pm    = connection_matrix(lc, mvf, algorithm="pmorse")

    # All algorithms must produce the same Conley indices
    @test cm_dlms.conley == cm_dhl.conley
    @test cm_dlms.conley == cm_hms.conley
    @test cm_dlms.conley == cm_pm.conley

    # All algorithms must produce the same number of Morse sets
    @test length(cm_dlms.morse) == length(cm_dhl.morse)
    @test length(cm_dlms.morse) == length(cm_hms.morse)
    @test length(cm_dlms.morse) == length(cm_pm.morse)
end

@testset "Connection matrix correctness" begin
    # example_three_cm(1): gradient variant with known connection matrix
    lc, mvf = example_three_cm(1)
    cm = connection_matrix(lc, mvf, algorithm="DHL")

    # In gradient variant 1, D-DE and E-CE are regular pairs, so DE is NOT critical;
    # the critical cells are A, C, AC, BD, CD, DF, ABC, EFG
    @test "A"   in cm.labels
    @test "C"   in cm.labels
    @test "CD"  in cm.labels
    @test "ABC" in cm.labels
    @test "EFG" in cm.labels
    @test !("DE" in cm.labels)

    cmfull = full_from_sparse(cm.matrix)
    @test sum(abs.(cmfull)) > 0   # gradient with connections → nonzero matrix

    # Tutorial logo example (already in Tutorial testset, used as cross-check)
    labels = ["A","B","C","D"]
    simplices = [["A","B","C"],["B","C","D"]]
    sclogo = create_simplicial_complex(labels, simplices)
    mvflogo = [["A","AB"],["C","AC"],["B","BC","BD","BCD"]]
    cmlogo = connection_matrix(sclogo, mvflogo)
    @test full_from_sparse(cmlogo.matrix) == [0 0 0; 0 0 1; 0 0 0]
end

@testset "Conley index" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # conley_index agrees with relative_homology via clomo pair
    for cell in ["F","DF","DEF"]
        cl, mo = lefschetz_clomo_pair(sc, [cell])
        @test conley_index(sc, [cell]) == relative_homology(sc, cl, mo)
    end

    # Conley index of larger isolated invariant set
    cl_tri, mo_tri = lefschetz_clomo_pair(sc, ["AB","AC","BC","A","B","C"])
    @test conley_index(sc, ["AB","AC","BC","A","B","C"]) == relative_homology(sc, cl_tri, mo_tri)
end

@testset "Forman cell types" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)
    formanvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]

    # forman_critical_cells returns the same count as mvf_information["N critical"]
    critical = forman_critical_cells(sc, formanvf)
    info = mvf_information(sc, formanvf)
    @test length(critical) == info["N critical"]

    # forman_all_cell_types: every cell appears exactly once across c, s, t
    c, s, t = forman_all_cell_types(sc, formanvf)
    all_classified = sort(vcat(c, s, t))
    @test all_classified == sort(sc.labels)
end

@testset "Conley index algorithms and fields" begin
    # The chaotic Morse set of example_torsion_chaos(3,p) has trivial
    # index unless p = 3, where it has the index of a periodic orbit
    for p in (2, 3, 5)
        sc, vfG, vfC = example_torsion_chaos(3, p)
        ms = morse_sets(sc, vfC)
        @test length.(ms) == [1, 26]
        @test conley_index(sc, ms[2]) == (p == 3 ? [0, 1, 1] : [0, 0, 0])
        @test [conley_index(sc, m) for m in morse_sets(sc, vfG)] ==
              [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    end

    # All algorithms agree with relative homology of the clomo pair,
    # for index and label input, with and without validation
    l1, m1, l2, m2 = example_moebius(0)
    for (lc, mvf) in ((example_three_cm()...,), (l2, m2),
                      (lefschetz_gfp_conversion(l2, 2), m2),
                      (lefschetz_gfp_conversion(l2, 7), m2))
        mvfI = convert_cellsubsets(lc, mvf)
        sets = vcat(mvfI, convert_cellsubsets(lc, morse_sets(lc, mvf)),
                    [collect(1:lc.ncells)])
        for s in sets
            cl, mo = lefschetz_clomo_pair(lc, s)
            ref = relative_homology(lc, cl, mo)
            for alg in ("rank", "matrix", "morse", "pmorse")
                @test conley_index(lc, s, algorithm=alg) == ref
                @test conley_index(lc, lc.labels[s], algorithm=alg) == ref
                @test conley_index(lc, s, algorithm=alg, validate=false) == ref
            end
        end
        # Single multivectors in no Morse set have trivial index
        msI = convert_cellsubsets(lc, morse_sets(lc, mvf))
        for mv in mvfI
            if !any(m -> issubset(mv, m), msI)
                cl, mo = lefschetz_clomo_pair(lc, mv)
                @test sum(relative_homology(lc, cl, mo)) == 0
            end
        end
    end

    # Errors for sets which are not locally closed and unknown algorithms
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)
    @test_throws ErrorException conley_index(sc, ["D","DEF"])
    @test_throws ErrorException conley_index(sc, ["D","DEF"], algorithm="rank")
    @test_throws ErrorException conley_index(sc, ["DEF"], algorithm="unknown")
    @test conley_index(sc, Int[]) == [0, 0, 0]
    @test conley_index(sc, ["F"], algorithm="rank") == [1, 0, 0]
    @test conley_index(sc, ["DF"], algorithm="rank") == [0, 1, 0]
    @test conley_index(sc, ["DEF"], algorithm="rank") == [0, 0, 1]
    @test conley_index(sc, ["AB","AC","BC","A","B","C"], algorithm="rank") == [1, 1, 0]
end

@testset "Automatic Conley index algorithm" begin
    # The threshold is the one described in the documentation
    @test ConleyDynamics.CONLEY_INDEX_RANK_MAX == 100
    # Locally closed sets on both sides of the threshold: closures
    # of growing groups of squares, and these minus their boundaries
    for p in (0, 2, 3)
        lc, coords = create_cubical_rectangle(12, 12, p=p)
        squares = findall(isequal(2), lc.dimensions)
        bnd = lefschetz_closure(lc, manifold_boundary(lc))
        for nsq in (1, 4, 9, 16, 25, 36, 64, 144)
            closed = lefschetz_closure(lc, squares[1:nsq])
            for s in (closed, setdiff(closed, bnd), setdiff(1:lc.ncells, closed))
                isempty(s) && continue
                ci = conley_index(lc, s, algorithm="morse")
                @test conley_index(lc, s) == ci
                @test conley_index(lc, s, algorithm="auto") == ci
                @test conley_index(lc, s, algorithm="rank") == ci
                @test conley_index(lc, lc.labels[s]) == ci
            end
        end
    end
    # Morse sets of a field with one large multivector
    lc, coords = create_simplicial_rectangle(10, 10)
    interior = setdiff(1:lc.ncells, lefschetz_closure(lc, manifold_boundary(lc)))
    @test length(interior) > ConleyDynamics.CONLEY_INDEX_RANK_MAX
    @test conley_index(lc, interior) == conley_index(lc, interior, algorithm="rank")
    # Cells outside the multivector are critical singletons
    ms = morse_sets(lc, [interior])
    @test interior in ms
    @test all(m -> m == interior || length(m) == 1, ms)
end

@testset "Morse reductions with many Morse sets" begin
    # Planar flow whose number of equilibria grows with the grid size
    for p in (0, 2, 3)
        lc, coords = create_cubical_rectangle(24, 24, p=p)
        vf(x) = [sin(3.0*x[1]), sin(3.0*x[2])]
        coordsN = convert_planar_coordinates(coords, [0.1, 0.1], [2pi-0.1, 2pi-0.1])
        mvf = create_planar_mvf(lc, coordsN, vf)
        cm_dlms = connection_matrix(lc, mvf, algorithm="DLMS")
        for alg in ("HMS", "pmorse")
            cm = connection_matrix(lc, mvf, algorithm=alg)
            @test cm.morse == cm_dlms.morse
            @test cm.conley == cm_dlms.conley
        end
        @test length(cm_dlms.morse) > 10
        # Persistent homology with the Morse algorithms
        filt = 1 .+ (create_random_filter(lc) .÷ 5)
        ph = persistent_homology(lc, filt, algorithm="matrix")
        @test persistent_homology(lc, filt, algorithm="morse") == ph
        @test persistent_homology(lc, filt, algorithm="pmorse") == ph
    end
end

@testset "Morse set poset" begin
    # Reference Hasse diagram: compute the reachability relation R between
    # the Morse sets by depth-first search, and keep the covers R \ R∘R
    function hasse_reference(lc, mvfI, msI)
        bnd  = full_from_sparse(lc.boundary)
        succ = [findall(!iszero, bnd[:, c]) for c in 1:lc.ncells]
        for mv in mvfI, m in mv
            append!(succ[m], mv)
        end
        nms = length(msI)
        R = fill(false, nms, nms)
        for k in 1:nms
            seen = fill(false, lc.ncells)
            stack = copy(msI[k])
            seen[stack] .= true
            while !isempty(stack)
                for w in succ[pop!(stack)]
                    if !seen[w]
                        seen[w] = true
                        push!(stack, w)
                    end
                end
            end
            for j in 1:nms
                R[k, j] = (j != k) && seen[msI[j][1]]
            end
        end
        return permutedims(R .& .!((R * R) .> 0))
    end

    function samplevf(x::Vector{Float64})
        x1, x2 = x
        return [x1 * (1.0 - x1*x1 - 3.0*x2*x2), x2 * (1.0 - 3.0*x1*x1 - x2*x2)]
    end

    cases = Any[example_three_cm(), example_three_cm(1), example_clorenz(),
                example_multiflow(), example_torsion_chaos(3, 3)[[1, 2]],
                example_torsion_chaos(3, 3)[[1, 3]], example_dunce_chaos()[[1, 2]],
                example_dunce_chaos()[[1, 3]]]
    # Every cell critical: the poset is the face poset of the complex
    for lc in (create_cubical_rectangle(4, 3)[1], create_simplicial_rectangle(3, 2)[1])
        push!(cases, (lc, [[k] for k in 1:lc.ncells]))
    end
    sc, coords = create_simplicial_rectangle(12, 12)
    coordsN = convert_planar_coordinates(coords, [-2.0, -2.0], [2.0, 2.0])
    push!(cases, (sc, create_planar_mvf(sc, coordsN, samplevf)))

    for (lc, mvf) in cases
        mvfI = mvf isa Vector{Vector{Int}} ? mvf : convert_cellsubsets(lc, mvf)
        msI, hasse = morse_sets(lc, mvfI, poset=true)
        @test msI == morse_sets(lc, mvfI)
        @test hasse == hasse_reference(lc, mvfI, msI)
        msL, hasseL = morse_sets(lc, mvf, poset=true)
        @test hasseL == hasse
    end
end

@testset "Restricted dynamics" begin
    for (lc, mvf) in (example_three_cm(), example_clorenz(),
                      example_torsion_chaos(3, 3)[[1, 3]])
        mvfI = convert_cellsubsets(lc, mvf)
        ms = morse_sets(lc, mvfI)
        mi = morse_interval(lc, mvfI, ms)
        sub = sort(mi)
        # Reference: intersect the multivectors with the subcomplex
        refL = [lc.labels[intersect(mv, sub)] for mv in mvfI]
        refL = refL[findall(x -> length(x) > 1, refL)]
        lcr, mvfr = restrict_dynamics(lc, mvfI, mi)
        @test lcr.labels == lc.labels[sub]
        @test convert_cellsubsets(lcr, mvfr) == refL
        lcrL, mvfrL = restrict_dynamics(lc, mvf, lc.labels[mi])
        @test mvfrL == refL
        lcrN, mvfrN = restrict_dynamics(lc, mvfI, mi, validate=false)
        @test mvfrN == mvfr
        @test full_from_sparse(lcrN.boundary) == full_from_sparse(lcr.boundary)
        # The Morse sets of the restricted field are the same
        @test sort(convert_cellsubsets(lcr, morse_sets(lcr, mvfr))) ==
              sort(convert_cellsubsets(lc, ms))
    end
    # Sets which are not locally closed raise an error
    lc, mvf = example_three_cm()
    nlc = [lc.indices["A"], lc.indices["ABC"]]
    @test !lefschetz_is_locally_closed(lc, nlc)
    @test_throws ErrorException restrict_dynamics(lc, mvf, nlc)
end

@testset "Cells without arcs in the flow digraph" begin
    # The cells a and bc of a triangle form a Lefschetz complex without
    # any facet relations. With the trivial multivector field, no cell is
    # incident to an arc of the flow digraph, but both are Morse sets.
    lc = create_simplicial_complex(["a","b","c"], [["a","b","c"]])
    i = lc.indices
    lcr, mvfr = restrict_dynamics(lc, Vector{Vector{Int}}(), [i["a"], i["bc"]])
    @test lcr.ncells == 2 && isempty(mvfr)
    @test sort(morse_sets(lcr, mvfr)) == [[1], [2]]
    @test sort(morse_sets(lcr, [[1], [2]])) == [[1], [2]]
    ms, hasse = morse_sets(lcr, mvfr; poset=true)
    @test length(ms) == 2 && size(hasse) == (2, 2)
    @test morse_interval(lcr, mvfr, [[2]]) == [2]
    cm = connection_matrix(lcr, mvfr)
    @test sort(cm.morse) == [["a"], ["bc"]]            # Morse sets are returned as labels
    @test sum(cm.conley) == [1, 1]
    # A single open edge
    lce, mvfe = restrict_dynamics(lc, Vector{Vector{Int}}(), [i["ab"]])
    @test morse_sets(lce, mvfe) == [[1]]
    @test connection_matrix(lce, mvfe).conley == [[0, 1]]
    # The last cell (the triangle) is only incident to its facet arcs, and
    # a field which pairs it with an edge must not lose any cell
    mvf = [[i["bc"], i["abc"]]]
    @test sort(reduce(vcat, morse_sets(lc, mvf))) ==
          sort(reduce(vcat, morse_sets(lc, vcat(mvf, [[c] for c in 1:lc.ncells if !(c in mvf[1])]))))
end
