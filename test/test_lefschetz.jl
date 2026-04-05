@testset "Topological operations" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Closure of a single edge includes its two vertices
    cl_AB = lefschetz_closure(sc, ["AB"])
    @test sort(cl_AB) == sort(["A","B","AB"])

    # Closure of a triangle includes all faces
    cl_DEF = lefschetz_closure(sc, ["DEF"])
    @test sort(cl_DEF) == sort(["D","E","F","DE","DF","EF","DEF"])

    # Closure of a vertex is just the vertex
    @test lefschetz_closure(sc, ["A"]) == ["A"]

    # Skeleton: dimension 0 returns only vertices (as integer indices)
    sk0 = lefschetz_skeleton(sc, 0)
    @test length(sk0) == 6
    @test all(sc.dimensions[v] == 0 for v in sk0)

    # Skeleton: dimension 1 returns only edges, not vertices (as integer indices)
    sk1 = lefschetz_skeleton(sc, 1)
    @test all(sc.dimensions[v] == 1 for v in sk1)
    @test length(sk1) == 7
    @test !(sc.indices["A"] in sk1)

    # Neighbors of vertex F: cells in closure or openhull of {F} but not F itself
    nb_F = lefschetz_neighbors(sc, ["F"])
    @test "DF" in nb_F
    @test "EF" in nb_F
    @test "DEF" in nb_F
    @test !("F" in nb_F)
end

@testset "Open and closed predicates" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Closure is closed
    @test lefschetz_is_closed(sc, lefschetz_closure(sc, ["AB"]))

    # Single edge without vertices is not closed
    @test !lefschetz_is_closed(sc, ["AB"])

    # Vertices alone are closed
    @test lefschetz_is_closed(sc, ["A","B","C"])

    # Full complex is closed
    @test lefschetz_is_closed(sc, sc.labels)

    # Empty set is closed
    @test lefschetz_is_closed(sc, String[])

    # A closed set is always locally closed
    @test lefschetz_is_locally_closed(sc, lefschetz_closure(sc, ["DEF"]))

    # Single edge (without vertices) is locally closed (mouth = {A,B} is closed)
    @test lefschetz_is_locally_closed(sc, ["AB"])
end

@testset "Subcomplex extraction" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Closed subcomplex of triangle DEF: 7 cells
    cl_DEF = lefschetz_closure(sc, ["DEF"])
    sub = lefschetz_subcomplex(sc, cl_DEF)
    @test sub.ncells == 7
    # Filled triangle is contractible: [1, 0, 0]
    @test homology(sub) == [1, 0, 0]

    # Extract the boundary circle of the triangle: 3 vertices + 3 edges
    circle = setdiff(cl_DEF, ["DEF"])
    sub2 = lefschetz_subcomplex(sc, circle)
    @test sub2.ncells == 6
    # Circle has homology [1, 1]
    @test homology(sub2) == [1, 1]
end

@testset "Surfaces" begin
    # Torus: H = [1, 2, 1] over any field
    @test homology(simplicial_torus(0)) == [1, 2, 1]
    @test homology(simplicial_torus(2)) == [1, 2, 1]

    # Klein bottle: H = [1, 1, 0] over rationals, [1, 2, 1] over GF(2)
    @test homology(simplicial_klein_bottle(0)) == [1, 1, 0]
    @test homology(simplicial_klein_bottle(2)) == [1, 2, 1]

    # Projective plane: H = [1, 0, 0] over rationals, [1, 1, 1] over GF(2)
    @test homology(simplicial_projective_plane(0)) == [1, 0, 0]
    @test homology(simplicial_projective_plane(2)) == [1, 1, 1]
end
