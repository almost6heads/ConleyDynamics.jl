@testset "MVF acyclicity" begin
    # Single gradient pair on a one-edge complex: a→ab, b is critical → acyclic
    sc_edge = create_simplicial_complex(["a","b"], [["a","b"]])
    @test mvf_is_acyclic(sc_edge, [["a","ab"]])

    # Gradient variant from example_three_cm(1): no periodic orbits → acyclic
    lc_grad, mvf_grad = example_three_cm(1)
    @test mvf_is_acyclic(lc_grad, mvf_grad)

    # Periodic orbit example from example_three_cm(0) is NOT acyclic
    lc_three, mvf_periodic = example_three_cm(0)
    @test !mvf_is_acyclic(lc_three, mvf_periodic)
end

@testset "MVF information" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)
    formanvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]

    info = mvf_information(sc, formanvf)

    # Required keys exist
    @test haskey(info, "N mv")
    @test haskey(info, "N regular")
    @test haskey(info, "N critical")
    @test haskey(info, "Lengths regular")
    @test haskey(info, "Lengths critical")

    # Total multivectors = regular + critical
    @test info["N mv"] == info["N regular"] + info["N critical"]
end

@testset "Morse sets" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)
    formanvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]

    ms = morse_sets(sc, formanvf)

    # The gradient field has critical cells: F, DF, EF, DEF, and the connected component {A,B,C,...}
    # Number of nontrivial Morse sets equals number in cm.morse from Tutorial testset
    cm = connection_matrix(sc, formanvf)
    @test length(ms) == length(cm.morse)

    # Periodic orbit example: one Morse set containing the orbit
    lc_three, mvf_periodic = example_three_cm(0)
    ms_per = morse_sets(lc_three, mvf_periodic)
    @test length(ms_per) >= 1
end
