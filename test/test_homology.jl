@testset "Homology algorithm variants" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # All three algorithms should agree on the same complex
    @test homology(sc, algorithm="matrix") == [1, 1, 0]
    @test homology(sc, algorithm="morse")  == [1, 1, 0]
    @test homology(sc, algorithm="pmorse") == [1, 1, 0]

    # Torus: all algorithms agree
    torus = simplicial_torus(2)
    @test homology(torus, algorithm="matrix") == [1, 2, 1]
    @test homology(torus, algorithm="morse")  == [1, 2, 1]
    @test homology(torus, algorithm="pmorse") == [1, 2, 1]
end

@testset "Relative homology" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Empty subcomplex gives same as absolute homology
    @test relative_homology(sc, String[]) == homology(sc)

    # Relative homology of the full complex over itself is all zeros
    @test all(relative_homology(sc, sc.labels) .== 0)

    # clomo pair and relative_homology agree with conley_index
    cl1, mo1 = lefschetz_clomo_pair(sc, ["F"])
    @test relative_homology(sc, cl1, mo1) == conley_index(sc, ["F"])

    cl2, mo2 = lefschetz_clomo_pair(sc, ["DEF"])
    @test relative_homology(sc, cl2, mo2) == conley_index(sc, ["DEF"])
end

@testset "Persistent homology" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    filtration = [1,1,1,2,2,2,1,1,1,3,2,2,2,4]
    phsingles, phpairs = persistent_homology(sc, filtration)

    # Each dimension must have a list of singles and pairs
    @test length(phsingles) == sc.dim + 1
    @test length(phpairs)   == sc.dim + 1

    # All birth/death pairs must satisfy birth < death
    for dim_pairs in phpairs
        for (b, d) in dim_pairs
            @test b < d
        end
    end
end

@testset "Homology edge cases" begin
    # Single vertex
    lc_pt = create_simplicial_complex(["v"], Vector{Vector{String}}())
    @test homology(lc_pt) == [1]

    # Two disconnected vertices
    lc_2pt = create_simplicial_complex(["u","v"], Vector{Vector{String}}())
    @test homology(lc_2pt) == [2]
end
