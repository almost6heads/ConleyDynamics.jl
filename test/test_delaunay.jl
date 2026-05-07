@testset "DelaunayTriangulation extension" begin
    #
    # Three points forming a single triangle — the simplest deterministic case.
    # Expected complex: 3 vertices + 3 edges + 1 triangle = 7 cells.
    # Expected homology of a filled triangle: [1, 0, 0].
    #
    pts = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    tri = triangulate(pts)

    sc = delaunay_to_simplicial(tri)
    @test sc isa EuclideanComplex
    @test sc.ncells == 7
    @test homology(sc) == [1, 0, 0]

    sc3 = delaunay_to_simplicial(tri; p=3)
    @test sc3 isa EuclideanComplex
    @test sc3.ncells == 7
    @test homology(sc3) == [1, 0, 0]
end
