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

# =============================================================================
# Tests for split_segments
# =============================================================================

@testset "split_segments — simple X crossing" begin
    # Two diagonals of the unit square cross at (1,1).
    segs = [[[0.0, 0.0], [2.0, 2.0]],
            [[0.0, 2.0], [2.0, 0.0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 5     # 4 corners + crossing at (1,1)
    @test length(edg) == 4     # each diagonal split into 2 sub-segments
end

@testset "split_segments — T-intersection" begin
    # Horizontal segment touched at its midpoint by a vertical one.
    segs = [[[0.0, 1.0], [2.0, 1.0]],
            [[1.0, 0.0], [1.0, 1.0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 4
    @test length(edg) == 3
end

@testset "split_segments — collinear partial overlap" begin
    segs = [[[0.0, 0.0], [3.0, 0.0]],
            [[1.0, 0.0], [4.0, 0.0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 4     # x = 0, 1, 3, 4
    @test length(edg) == 3
end

@testset "split_segments — disjoint segments unchanged" begin
    segs = [[[0.0, 0.0], [1.0, 0.0]],
            [[2.0, 0.0], [3.0, 0.0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 4
    @test length(edg) == 2
end

@testset "split_segments — shared endpoint only" begin
    segs = [[[0.0, 0.0], [1.0, 1.0]],
            [[1.0, 1.0], [2.0, 0.0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 3
    @test length(edg) == 2
end

@testset "split_segments — star (four segments through origin)" begin
    segs = [[[-1.0,  0.0], [1.0,  0.0]],
            [[ 0.0, -1.0], [0.0,  1.0]],
            [[-1.0, -1.0], [1.0,  1.0]],
            [[-1.0,  1.0], [1.0, -1.0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 9     # 8 tips + origin
    @test length(edg) == 8
end

@testset "split_segments — integer inputs accepted" begin
    segs = [[[0, 0], [2, 2]],
            [[0, 2], [2, 0]]]
    pts, edg = split_segments(segs)
    @test length(pts) == 5
    @test length(edg) == 4
    @test eltype(pts) == Vector{Float64}
end

# =============================================================================
# Tests for delaunay_points_add_split_segments
# =============================================================================

@testset "delaunay_points_add_split_segments — crossing diagonals" begin
    points, bndcurve = delaunay_points_bnd_rectangle([-2, -2], [2, 2])
    segs = [[[-1.0, -1.0], [1.0,  1.0]],
            [[-1.0,  1.0], [1.0, -1.0]]]
    points2, segments2 = delaunay_points_add_split_segments(points, segs)

    # 4 boundary corners + 5 split-segment vertices (4 tips + origin)
    @test length(points2) == 9
    # 4 sub-segments as constraints
    @test length(segments2) == 4
end

@testset "delaunay_points_add_split_segments — three-argument form" begin
    points, bndcurve = delaunay_points_bnd_rectangle([-2, -2], [2, 2])
    segs1 = [[[-1.0, -1.0], [1.0,  1.0]],
             [[-1.0,  1.0], [1.0, -1.0]]]
    segs2 = [[[0.0, -1.5], [0.0, 1.5]]]   # vertical segment, no new crossings

    points, segs_set = delaunay_points_add_split_segments(points, segs1)
    points, segs_set = delaunay_points_add_split_segments(points, segs_set, segs2)

    # segs2 is a single segment split by nothing → 2 vertices, 1 edge added
    @test length(points) == 11   # 4 + 5 + 2
    @test length(segs_set) == 5  # 4 + 1
end

@testset "delaunay_points_add_split_segments — triangulation round-trip" begin
    points, bndcurve = delaunay_points_bnd_rectangle([-2, -2], [2, 2])
    segs = [[[-1.0, 0.0], [1.0, 0.0]],
            [[ 0.0,-1.0], [0.0, 1.0]]]
    points, segments = delaunay_points_add_split_segments(points, segs)

    tri = triangulate(points; boundary_nodes = bndcurve, segments = segments)
    sc  = delaunay_to_simplicial(tri)

    @test sc isa EuclideanComplex
    @test homology(sc) == [1, 0, 0]
end
