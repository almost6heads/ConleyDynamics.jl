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

@testset "delaunay_points_bnd_rectangle" begin
    #
    # The rectangle [-1,1]x[-1,1] should yield exactly 4 corner points and a
    # single closed boundary curve [1,2,3,4,1].
    #
    points, bndcurve = delaunay_points_bnd_rectangle([-1, -1], [1, 1])

    @test length(points) == 4
    @test points[1] == (-1.0, -1.0)
    @test points[2] == ( 1.0, -1.0)
    @test points[3] == ( 1.0,  1.0)
    @test points[4] == (-1.0,  1.0)
    @test bndcurve == [[1, 2, 3, 4, 1]]

    # Integer inputs should be accepted and converted to Float64.
    points2, _ = delaunay_points_bnd_rectangle([0, 0], [3, 2])
    @test eltype(points2) == Tuple{Float64, Float64}
end

@testset "delaunay_points_add_segment — closed curve" begin
    #
    # A closed triangular constraint: [0,0]->[1,0]->[0,1]->[0,0].
    # Three distinct points, 3 edges, no duplicate stored.
    #
    points, _ = delaunay_points_bnd_rectangle([-2, -2], [2, 2])
    @test length(points) == 4

    triangle = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]
    points, segments = delaunay_points_add_segment(points, triangle)

    # 4 boundary corners + 3 new points (duplicate endpoint not stored)
    @test length(points) == 7
    # 3 edges forming the closed triangle
    @test length(segments) == 3
    @test (5, 6) in segments
    @test (6, 7) in segments
    @test (7, 5) in segments
end

@testset "delaunay_points_add_segment — open curve" begin
    #
    # An open line with 3 points: 2 edges, all 3 points stored.
    #
    points, _ = delaunay_points_bnd_rectangle([-2, -2], [2, 2])

    line = [[-1.0, 0.0], [0.0, 1.0], [1.0, 0.0]]
    points, segments = delaunay_points_add_segment(points, line)

    # 4 boundary corners + 3 new points
    @test length(points) == 7
    # 2 edges
    @test length(segments) == 2
    @test (5, 6) in segments
    @test (6, 7) in segments
end

@testset "delaunay_points_add_segment — chained calls" begin
    #
    # Two successive calls accumulate points and edges correctly.
    #
    points, _ = delaunay_points_bnd_rectangle([0, 0], [10, 10])

    seg1 = [[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 1.0]]   # closed triangle
    seg2 = [[5.0, 5.0], [6.0, 5.0], [6.0, 6.0], [5.0, 5.0]]   # another closed triangle

    points, segments = delaunay_points_add_segment(points, seg1)
    points, segments = delaunay_points_add_segment(points, segments, seg2)

    # 4 + 3 + 3 = 10 points; 3 + 3 = 6 edges
    @test length(points) == 10
    @test length(segments) == 6
end

@testset "delaunay_points — full triangulation round-trip" begin
    #
    # Build a constrained triangulation and convert it to a EuclideanComplex.
    # The square [-1,1]^2 with a triangular constraint should triangulate
    # without error and produce a contractible complex (H = [1,0,0]).
    #
    points, bndcurve = delaunay_points_bnd_rectangle([-1, -1], [1, 1])
    triangle = [[0.0, 0.0], [0.5, 0.0], [0.0, 0.5], [0.0, 0.0]]
    points, segments = delaunay_points_add_segment(points, triangle)

    tri = triangulate(points; boundary_nodes = bndcurve, segments = segments)
    sc  = delaunay_to_simplicial(tri)

    @test sc isa EuclideanComplex
    @test homology(sc) == [1, 0, 0]
end
