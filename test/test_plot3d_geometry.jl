const CDG = ConleyDynamics

tri_signed_volume(tris) = sum(sum(a .* CDG._cross3(b, c)) / 6 for (a, b, c) in tris)
tri_area(tris) = sum(sqrt(sum(abs2, CDG._cross3(b .- a, c .- a))) / 2 for (a, b, c) in tris)

@testset "3D plot geometry: simplicial" begin
    #
    # Two tetrahedra sharing a face
    #
    points = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0], [1.0,1.0,1.0]]
    ec  = tetmesh_to_euclidean(points, [[1,2,3,4], [2,3,4,5]])
    cof = CDG._coface_lists(ec)
    all3 = trues(ec.ncells)

    # The hull consists of six outward oriented triangles

    tris, segs, pts = CDG._background3d_geometry(ec, all3, cof)
    @test length(tris) == 6
    @test length(segs) == 9
    @test length(pts) == 5
    @test tri_signed_volume(tris) ≈ 1/6 + 1/3
    tris, segs, pts = CDG._background3d_geometry(ec, all3, cof; hull_only=false)
    @test length(tris) == 7
    tris, segs, pts = CDG._background3d_geometry(ec, all3, cof; pdim=[false,false,true,false])
    @test isempty(tris) && isempty(segs) && isempty(pts)

    # Clipping removes the second tetrahedron, and the visible cells
    # form a closed subcomplex

    keep = CDG._clip_mask(ec, ([1,1,1], 1.2), cof)
    @test keep[ec.indices["1234"]]
    @test !keep[ec.indices["2345"]]
    @test !keep[ec.indices["5"]]
    @test keep[ec.indices["234"]]
    @test all(all(keep[ec.boundary.columns[k]]) for k in 1:ec.ncells if keep[k])
    tris, _, _ = CDG._background3d_geometry(ec, keep, cof)
    @test length(tris) == 4
    @test tri_signed_volume(tris) ≈ 1/6

    # Cell subsets: a tetrahedron with some of its faces, a free edge,
    # and an isolated vertex

    cells = [ec.indices[l] for l in ["1234", "123", "12", "1", "25", "5"]]
    tris, segs, pts = CDG._cells3d_geometry(ec, cells, all3)
    @test length(tris) == 4
    @test length(segs) == 1
    @test isempty(pts)
    cells2 = [ec.indices[l] for l in ["1234", "123", "5"]]
    tris, segs, pts = CDG._cells3d_geometry(ec, cells2, all3)
    @test (length(tris), length(segs), length(pts)) == (4, 0, 1)

    # Cells covered by a clipped tetrahedron are not drawn

    cells3 = [ec.indices[l] for l in ["2345", "345", "35"]]
    tris, segs, pts = CDG._cells3d_geometry(ec, cells3, keep)
    @test isempty(tris) && isempty(segs) && isempty(pts)

    # Triangles of a surface in a subset, and planar complexes

    ecs = tetmesh_to_euclidean(trisurface_sphere([0,0,0], 1.0, 1)...)
    cofs = CDG._coface_lists(ecs)
    tris, segs, _ = CDG._background3d_geometry(ecs, trues(ecs.ncells), cofs)
    @test length(tris) == 80
    @test length(segs) == 120
    ec2 = create_simplicial_delaunay(1.0, 1.0, 20; euclidean=true)
    @test_throws ArgumentError CDG._check_coords3d(ec2, "test")
    @test CDG._length_scale(ec) ≈ sqrt(2)
end

@testset "3D plot geometry: cubical" begin
    lc, coords = create_cubical_box(2, 3, 2)
    ec  = lefschetz_to_euclidean(lc, convert_spatial_coordinates(coords, [0,0,0], [2,3,2]))
    cof = CDG._coface_lists(ec)

    # Squares are split into two triangles along a diagonal

    sq = findfirst(==(2), ec.dimensions)
    @test tri_area(CDG._face_triangles(ec, sq)) ≈ 1.0

    # The hull of all cubes is the boundary of the box, oriented outward

    tris, _, _ = CDG._background3d_geometry(ec, trues(ec.ncells), cof)
    @test length(tris) == 2 * 2 * (2*3 + 2*2 + 3*2)
    @test tri_signed_volume(tris) ≈ 12.0
    @test tri_area(tris) ≈ 2 * (6 + 4 + 6)

    # The hull of a set of cubes encloses their volume

    cubes = findall(==(3), ec.dimensions)[1:5]
    tris, segs, pts = CDG._cells3d_geometry(ec, cubes, trues(ec.ncells))
    @test tri_signed_volume(tris) ≈ 5.0
    @test isempty(segs) && isempty(pts)
end

@testset "3D plot geometry: input validation" begin
    ec = tetmesh_to_euclidean([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], [[1,2,3,4]])
    @test isnothing(CDG._check_pdim3d([true,true,true,true], "f"))
    @test_throws ArgumentError CDG._check_pdim3d([true,true,true], "f")
    @test CDG._cells3d_indices(ec, [1, 15], "f") == [1, 15]
    @test CDG._cells3d_indices(ec, ["1", "1234"], "f") == [1, ec.indices["1234"]]
    @test CDG._cells3d_indices(ec, Any["2", 3], "f") == [2, 3]
    @test CDG._cells3d_indices(ec, Int[], "f") == Int[]
    @test_throws ArgumentError CDG._cells3d_indices(ec, ["nolabel"], "f")
    @test_throws ArgumentError CDG._cells3d_indices(ec, [0], "f")
    @test_throws ArgumentError CDG._cells3d_indices(ec, [ec.ncells+1], "f")
    @test_throws ArgumentError CDG._cells3d_indices(ec, [1.0], "f")
end

@testset "3D plot geometry: pdim as BitVector" begin
    ec = tetmesh_to_euclidean([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], [[1,2,3,4]])
    cof  = CDG._coface_lists(ec)
    keep = CDG._clip_mask(ec, nothing, cof)
    @test CDG._background3d_geometry(ec, keep, cof; pdim=trues(4)) ==
          CDG._background3d_geometry(ec, keep, cof; pdim=[true,true,true,true])
    @test CDG._cells3d_geometry(ec, [15], keep; pdim=trues(4)) ==
          CDG._cells3d_geometry(ec, [15], keep; pdim=[true,true,true,true])
end

@testset "Planar and spatial functions check the embedding dimension" begin
    #
    # Planar functions reject complexes in three dimensions, and spatial
    # multivector fields need three-dimensional complexes in space
    #
    ec3 = tetmesh_to_euclidean([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], [[1,2,3,4]])
    ecs = tetmesh_to_euclidean(trisurface_sphere([0,0,0], 1.0, 1)...)
    ec2 = tetmesh_to_euclidean([[0,0],[1,0],[0,1]], [[1,2,3]])
    vf3(x) = [-x[2], x[1], 0.0]
    vf2(x) = [-x[2], x[1]]
    @test_throws ArgumentError create_planar_mvf(ec3, vf2)
    @test_throws ArgumentError create_spatial_mvf(ecs, vf3)
    @test_throws ArgumentError create_spatial_mvf(ec2, vf3)
    @test create_spatial_mvf(ec3, vf3) isa Vector
    @test create_planar_mvf(ec2, vf2) isa Vector
    mktempdir() do dir
        @test_throws ArgumentError plot_planar_simplicial(ec3, joinpath(dir, "a.pdf"))
        @test_throws ArgumentError plot_planar_simplicial_morse(ec3, joinpath(dir, "a.pdf"), [[1]])
    end
end

@testset "3D plot geometry: free faces and invalid input" begin
    #
    # Free triangles of a cell subset are drawn, while triangles which
    # are faces of tetrahedra in the subset are not
    #
    points = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0], [1.0,1.0,1.0]]
    ec   = tetmesh_to_euclidean(points, [[1,2,3,4], [2,3,4,5]])
    keep = trues(ec.ncells)
    tris, segs, pts = CDG._cells3d_geometry(ec, [ec.indices["234"], ec.indices["235"]], keep)
    @test length(tris) == 2 && isempty(segs) && isempty(pts)
    tris, _, _ = CDG._cells3d_geometry(ec, [ec.indices["1234"], ec.indices["234"]], keep)
    @test length(tris) == 4
    tris, _, _ = CDG._cells3d_geometry(ec, [ec.indices["234"]], keep;
                                       pdim=[true,true,false,true])
    @test isempty(tris)
    @test_throws ArgumentError CDG._clip_mask(ec, ([1,0], 0.0), CDG._coface_lists(ec))
end
