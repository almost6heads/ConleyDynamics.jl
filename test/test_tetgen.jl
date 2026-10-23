using LinearAlgebra
using Random

tet_volumes(ec) = [abs(det(hcat((c[2:4] .- Ref(c[1]))...))) / 6
                   for (c, d) in zip(ec.coords, ec.dimensions) if d == 3]

@testset "TetGen extension: boxes" begin
    #
    # A box gives a contractible complex whose tetrahedra fill the box
    #
    res = tetgen_mesh(tetgen_domain_box([0,0,0], [2,1,1]); maxvolume=0.01)
    @test res isa TetGenResult
    ec = tetgen_to_simplicial(res)
    @test ec isa EuclideanComplex
    @test ec.dim == 3
    @test homology(ec) == [1,0,0,0]
    vols = tet_volumes(ec)
    @test all(vols .> 0)
    @test abs(sum(vols) - 2.0) < 1e-10
    # TetGen volume constraints are approximate
    @test maximum(vols) <= 1.5 * 0.01
    @test count(>(0.01), vols) <= 0.05 * length(vols)
    @test count(==(0), ec.dimensions) == length(res.points)

    # The boundary cells form a sphere

    bnd = tetgen_marker_cells(res, ec, 0)
    @test homology(lefschetz_subcomplex(ec, bnd)) == [1,0,1]
    @test lefschetz_is_closed(ec, bnd)
    bndtri = tetgen_marker_cells(res, ec, 0; closure=false)
    @test all(ec.dimensions[bndtri] .== 2)
    @test length(bndtri) == count(==(0), res.trifacemarkers)

    # Other fields, and no quality refinement

    ec3 = tetgen_to_simplicial(res; p=3)
    @test homology(ec3) == [1,0,0,0]
    res0 = tetgen_mesh(tetgen_domain_box([0,0,0], [1,1,1]); minratio=nothing)
    @test length(res0.points) == 8
    resd = tetgen_mesh(tetgen_domain_box([0,0,0], [1,1,1]); maxvolume=0.01)
    resq = tetgen_mesh(tetgen_domain_box([0,0,0], [1,1,1]); maxvolume=0.01,
                       minratio=1.2, mindihedral=15)
    @test length(resq.tetrahedra) > length(resd.tetrahedra)
    @test homology(tetgen_to_simplicial(resq)) == [1,0,0,0]
end

@testset "TetGen extension: holes, regions, and solid tori" begin
    dom = tetgen_domain_box([-2,-2,-2], [2,2,2])
    V, T = trisurface_sphere([0,0,0], 1.0, 2)

    # A spherical hole

    res = tetgen_mesh(tetgen_domain_add_surface(dom, V, T; hole=true); maxvolume=0.1)
    ec = tetgen_to_simplicial(res)
    @test homology(ec) == [1,0,1,0]
    @test sum(tet_volumes(ec)) < 64.0 - 4.0

    # An embedded sphere with a refined region inside

    res = tetgen_mesh(tetgen_domain_add_surface(dom, V, T; region=(7, 0.005)); maxvolume=0.1)
    ec = tetgen_to_simplicial(res)
    @test homology(ec) == [1,0,0,0]
    sph = tetgen_marker_cells(res, ec, 1)
    @test homology(lefschetz_subcomplex(ec, sph)) == [1,0,1]
    # Every triangle of the embedded sphere is a union of mesh triangles
    sphtri = filter(k -> ec.dimensions[k] == 2, sph)
    triarea(c) = norm(cross(c[2] .- c[1], c[3] .- c[1])) / 2
    @test abs(sum(triarea(ec.coords[k]) for k in sphtri) -
              sum(triarea(V[t]) for t in T)) < 1e-10
    inner = tetgen_region_cells(res, ec, 7; closure=false)
    @test all(ec.dimensions[inner] .== 3)
    innervols = tet_volumes(lefschetz_subcomplex(ec, lefschetz_closure(ec, inner)))
    @test maximum(innervols) <= 1.5 * 0.005
    @test count(>(0.005), innervols) <= 0.05 * length(innervols)
    @test homology(lefschetz_subcomplex(ec, tetgen_region_cells(res, ec, 7))) == [1,0,0,0]
    @test sort(unique(res.tetattributes)) == [7, 8]

    # A solid torus

    Vt, Tt = trisurface_torus([0,0,0], 2.0, 0.7, 24, 12)
    res = tetgen_mesh(tetgen_domain_from_surface(Vt, Tt); maxvolume=0.05)
    @test homology(tetgen_to_simplicial(res)) == [1,1,0,0]
end

@testset "TetGen extension: curves, nodes, and shared edges" begin
    box = tetgen_domain_box([0,0,0], [1,1,1])

    # A closed constraint curve becomes a circle of edges

    curve = [[0.2,0.2,0.5], [0.8,0.2,0.5], [0.8,0.8,0.5], [0.2,0.8,0.5], [0.2,0.2,0.5]]
    res = tetgen_mesh(tetgen_domain_add_curve(box, curve); maxvolume=0.001)
    ec = tetgen_to_simplicial(res)
    cc = tetgen_marker_cells(res, ec, 1)
    @test homology(lefschetz_subcomplex(ec, cc)) == [1,1]
    @test length(filter(k -> ec.dimensions[k] == 1, cc)) > 4

    # Interior nodes become vertices, outside nodes are ignored

    res = tetgen_mesh(tetgen_domain_add_nodes(box, [[0.5,0.5,0.5], [0.25,0.25,0.25], [3,3,3]]))
    ec = tetgen_to_simplicial(res)
    vcoords = [ec.coords[k][1] for k in 1:ec.ncells if ec.dimensions[k] == 0]
    @test [0.5,0.5,0.5] in vcoords
    @test [0.25,0.25,0.25] in vcoords
    @test !([3.0,3.0,3.0] in vcoords)

    # Two surfaces meeting along shared edges

    p1 = trisurface_plane_patch([0.2,0.2,0.5], [0.6,0,0], [0,0.6,0], 2, 2)
    p2 = trisurface_plane_patch([0.5,0.2,0.2], [0,0.6,0], [0,0,0.6], 2, 2)
    dom = tetgen_domain_add_surface(tetgen_domain_add_surface(box, p1...), p2...)
    res = tetgen_mesh(dom)
    ec = tetgen_to_simplicial(res)
    m1 = tetgen_marker_cells(res, ec, 1)
    m2 = tetgen_marker_cells(res, ec, 2)
    @test homology(lefschetz_subcomplex(ec, m1)) == [1,0,0]
    @test homology(lefschetz_subcomplex(ec, intersect(m1, m2))) == [1,0]
end

@testset "TetGen extension: errors" begin
    box = tetgen_domain_box([0,0,0], [1,1,1])

    # Crossing constraint surfaces

    p1 = trisurface_plane_patch([0.2,0.2,0.5], [0.6,0,0], [0,0.6,0], 3, 3)
    p2 = trisurface_plane_patch([0.47,0.25,0.21], [0,0.6,0], [0,0,0.6], 3, 3)
    dom = tetgen_domain_add_surface(tetgen_domain_add_surface(box, p1...), p2...)
    @test_throws ArgumentError tetgen_mesh(dom)

    # A hole point which is not enclosed by a surface

    @test_throws ArgumentError tetgen_mesh(tetgen_domain_add_hole(box, [0.5,0.5,0.5]))
    @test_throws ArgumentError tetgen_mesh(box; maxvolume=-1.0)
    @test_throws ArgumentError tetgen_region_cells(tetgen_mesh(box), tetgen_to_simplicial(tetgen_mesh(box)), 1)
end

@testset "TetGen extension: create_simplicial_delaunay3d" begin
    ec1 = create_simplicial_delaunay3d([0,0,0], [1,1,1], 0.15, 20; rng=Xoshiro(1))
    ec2 = create_simplicial_delaunay3d([0,0,0], [1,1,1], 0.15, 20; rng=Xoshiro(1))
    @test ec1 isa EuclideanComplex
    @test ec1.labels == ec2.labels
    @test ec1.coords == ec2.coords
    @test homology(ec1) == [1,0,0,0]
    pts = poisson_points_box([0,0,0], [1,1,1], 0.15; rng=Xoshiro(1))
    @test count(==(0), ec1.dimensions) == length(pts)

    ec3 = create_simplicial_delaunay3d([0,0,0], [1,2,1], 100; box_corners=true, rng=Xoshiro(2))
    @test count(==(0), ec3.dimensions) == 108
    @test abs(sum(tet_volumes(ec3)) - 2.0) < 1e-10
    ec4 = create_simplicial_delaunay3d([0,0,0], [1,2,1], 100; p=0, rng=Xoshiro(2))
    @test count(==(0), ec4.dimensions) == 100
    @test homology(ec4) == [1,0,0,0]
end

@testset "TetGen extension: multivector fields" begin
    #
    # The sample vector field from the create_spatial_mvf docstring has
    # three Morse sets and two connecting orbits
    #
    samplevf(x) = [x[1]*(1-x[1]^2), -x[2], -x[3]]
    res = tetgen_mesh(tetgen_domain_box([-1.5,-1,-1], [1.5,1,1]); maxvolume=0.01)
    ec  = tetgen_to_simplicial(res)
    mvf = create_spatial_mvf(ec, samplevf)
    cm  = connection_matrix(ec, mvf)
    @test sort(cm.conley) == [[0,1,0,0], [1,0,0,0], [1,0,0,0]]
    @test sum(full_from_sparse(cm.matrix)) == 2
end

@testset "TetGen extension: robustness" begin
    box = tetgen_domain_box([-2,-2,-2], [2,2,2])

    # Two spheres given as one surface give two holes

    Va, Ta = trisurface_sphere([ 1.0,0,0], 0.4, 1)
    Vb, Tb = trisurface_sphere([-1.0,0,0], 0.4, 1)
    dom = tetgen_domain_add_surface(box, vcat(Va, Vb),
                                    vcat(Ta, [t .+ length(Va) for t in Tb]); hole=true)
    @test homology(tetgen_to_simplicial(tetgen_mesh(dom; maxvolume=0.1))) == [1,0,2,0]

    # Automatic hole and region points avoid surfaces added before, so
    # that a surface around an existing hole bounds the shell

    Vi, Ti = trisurface_sphere([0,0,0], 0.5, 1)
    Vo, To = trisurface_sphere([0,0,0], 1.0, 1)
    dom = tetgen_domain_add_surface(box, Vi, Ti; hole=true)
    dom = tetgen_domain_add_surface(dom, Vo, To; region=(5, -1.0))
    @test 0.5 < sqrt(sum(abs2, dom.regions[1][1])) < 1.0
    res = tetgen_mesh(dom; maxvolume=0.1)
    @test sort(unique(res.tetattributes)) == [5,6]
    ec = tetgen_to_simplicial(res)
    @test homology(lefschetz_subcomplex(ec, tetgen_region_cells(res, ec, 5))) == [1,0,1,0]
    dom = tetgen_domain_add_surface(tetgen_domain_add_surface(box, Vi, Ti), Vo, To; hole=true)
    @test homology(tetgen_to_simplicial(tetgen_mesh(dom; maxvolume=0.1))) == [2,0,1,0]

    # Closed curves and holes in very large and very small domains

    for s in (1e-9, 1e8)
        V, T = trisurface_sphere([0,0,0], s, 1)
        bx = tetgen_domain_box(-2s .* [1,1,1], 2s .* [1,1,1])
        circ = [[1.5s*cos(t), 1.5s*sin(t), 0.0] for t in range(0, 2*pi, length=30)]
        dom = tetgen_domain_add_curve(tetgen_domain_add_surface(bx, V, T; hole=true), circ)
        res = tetgen_mesh(dom; maxvolume=0.1*s^3)
        ec  = tetgen_to_simplicial(res)
        @test homology(ec) == [1,0,1,0]
        @test homology(lefschetz_subcomplex(ec, tetgen_marker_cells(res, ec, 2))) == [1,1]
    end

    # Triangles need three vertices, and -0.0 is identified with 0.0

    @test_throws ArgumentError tetgen_domain_from_surface(Vi, [[t; t[1]] for t in Ti])
    @test_throws ArgumentError tetgen_domain_add_surface(box, Vi, [t[1:2] for t in Ti])
    dom = tetgen_domain_add_curve(box, [[-0.0,-2.0,-2.0], [0.0,0.0,0.0]])
    dom = tetgen_domain_add_curve(dom, [[0.0,-2.0,-2.0], [1.0,1.0,1.0]])
    @test length(dom.points) == 8 + 3

    # Quality bounds which would make TetGen hang or fail

    @test_throws ArgumentError tetgen_mesh(box; minratio=0.8)
    @test_throws ArgumentError tetgen_mesh(box; mindihedral=0)
    @test_throws ArgumentError tetgen_mesh(box; mindihedral=95)
    err = try tetgen_mesh(box; maxvolume=0.01, mindihedral=50); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("quality bounds", err.msg)

    # Failures leave no TetGen debugging files in the current directory

    V, T = trisurface_sphere([0,0,0], 1.0, 1)
    P = trisurface_plane_patch([-1.5,-1.5,0.2], [3,0,0], [0,3,0], 2, 2)
    bad = tetgen_domain_add_surface(tetgen_domain_add_surface(box, V, T), P...)
    mktempdir() do dir
        cd(dir) do
            touch("tetgen-tmpfile_mine.txt")
            err = try tetgen_mesh(bad); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin("intersect", err.msg)
            @test readdir(dir) == ["tetgen-tmpfile_mine.txt"]
        end
    end

    # Too few sample points for a Delaunay triangulation

    @test_throws ArgumentError create_simplicial_delaunay3d([0,0,0], [1,1,1], 5.0, 20)
    ec = create_simplicial_delaunay3d([0,0,0], [1,1,1], 5.0, 20; box_corners=true)
    @test count(==(0), ec.dimensions) == 9

    # The edge walk along constraint curves agrees with the full search

    circ = [[1.5*cos(s), 1.5*sin(s), 0.1] for s in range(0, 2*pi, length=41)]
    res = tetgen_mesh(tetgen_domain_add_curve(box, circ); maxvolume=0.02)
    ec  = tetgen_to_simplicial(res)
    for (a, b, m) in res.curvesegments
        walk  = ConleyDynamics._tetgen_curve_edges(ec, [(a, b)])
        brute = ConleyDynamics._tetgen_curve_edges_bruteforce(ec, collect(a), collect(b))
        @test sort(walk) == sort(brute)
        @test !isempty(walk)
    end
    @test homology(lefschetz_subcomplex(ec, tetgen_marker_cells(res, ec, 1))) == [1,1]
end

@testset "TetGen extension: options and input validation" begin
    #
    # Switch strings, boundary refinement, and invalid arguments
    #
    box = tetgen_domain_box([-1.5,-1,-1], [1.5,1,1])
    r1 = tetgen_mesh(box; maxvolume=0.01, minratio=1.2)
    r2 = tetgen_mesh(box; switches="pq1.2a0.01Q")
    @test r1.points == r2.points
    @test r1.tetrahedra == r2.tetrahedra
    r3 = tetgen_mesh(box; maxvolume=0.01, steiner_on_boundary=false)
    @test length(r3.points) == 8
    @test length(r3.trifaces) == 12
    empty = TetDomain(NTuple{3,Float64}[], Vector{Int}[], Int[])
    @test_throws ArgumentError tetgen_mesh(empty)
    @test_throws ArgumentError TetDomain(box.points, box.facets, [0])
    @test_throws ArgumentError tetgen_domain_add_curve(box, [[0,0,0]])
    @test_throws ArgumentError create_simplicial_delaunay3d([0,0], [1,1], 10)
    @test_throws ArgumentError create_simplicial_delaunay3d([0,0,0], [1,1,1], 3)
    @test_throws ArgumentError create_simplicial_delaunay3d([0,0,0], [1,0,1], 10)
end

@testset "TetGen extension: consistency checks" begin
    box = tetgen_domain_box([-2,-2,-2], [2,2,2])
    V, T = trisurface_sphere([0,0,0], 1.0, 1)
    res1 = tetgen_mesh(tetgen_domain_add_surface(box, V, T; region=(1, 0.05)); maxvolume=0.5)
    res2 = tetgen_mesh(box; maxvolume=0.1)
    ec1 = tetgen_to_simplicial(res1)
    ec2 = tetgen_to_simplicial(res2)

    # A complex from a different mesh is rejected

    @test_throws ArgumentError tetgen_marker_cells(res1, ec2, 0)
    @test_throws ArgumentError tetgen_marker_cells(res2, ec1, 0)
    @test_throws ArgumentError tetgen_region_cells(res1, ec2, 1)
    @test !isempty(tetgen_marker_cells(res1, ec1, 1))
    @test tetgen_marker_cells(res1, tetgen_to_simplicial(res1; p=3), 1) ==
          tetgen_marker_cells(res1, ec1, 1)

    # Constraint curves outside of the domain are reported

    dom = tetgen_domain_add_curve(box, [[3,3,3], [4,4,4]])
    res = @test_logs (:warn, r"outside") tetgen_mesh(dom)
    @test isempty(tetgen_marker_cells(res, tetgen_to_simplicial(res), 1))
    @test_logs tetgen_mesh(tetgen_domain_add_curve(box, [[0,0,0], [1,1,1]]))
end
