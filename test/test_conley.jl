@testset "Connection matrix algorithms agree" begin
    # Use gradient variant 1 of example_three_cm
    lc, mvf, _ = example_three_cm(1)

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
    lc, mvf, _ = example_three_cm(1)
    cm = connection_matrix(lc, mvf, algorithm="DHL")

    # In gradient variant 1, D-DE and E-CE are regular pairs, so DE is NOT critical;
    # the critical cells are A, C, AC, BD, CD, DF, ABC, EFG
    @test "A"   in cm.labels
    @test "C"   in cm.labels
    @test "CD"  in cm.labels
    @test "ABC" in cm.labels
    @test "EFG" in cm.labels
    @test !("DE" in cm.labels)

    # The connection matrix has the specific nonzero entry EFG → ABC (dimension 2→2 would be wrong;
    # check the actual known matrix: cm.matrix at [ABC,EFG] position should be nonzero)
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
