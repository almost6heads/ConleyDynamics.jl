using Random

euler_char(lc) = sum((-1)^d for d in lc.dimensions)

@testset "Vertex labels" begin
    #
    # The label helper has to reproduce the zero-padded format used by
    # create_simplicial_delaunay and delaunay_to_simplicial
    #
    for n in [1, 9, 10, 11, 99, 100, 101, 1000, 1234]
        labwidth = Int(ceil(log(0.5+n)/log(10)))
        labformat = "%0" * string(labwidth) * "d"
        oldlabels = [ConleyDynamics.Printf.format(ConleyDynamics.Printf.Format(labformat), k)
                     for k in 1:n]
        @test ConleyDynamics._vertex_labels(n) == oldlabels
    end

    sc, coords = create_simplicial_delaunay(2.0, 1.0, 50)
    @test sc.labels[1:50] == ConleyDynamics._vertex_labels(50)
    ec = create_simplicial_delaunay(2.0, 1.0, 0.2, 20; euclidean=true)
    nv = count(==(0), ec.dimensions)
    @test ec.labels[1:nv] == ConleyDynamics._vertex_labels(nv)
end

@testset "tetmesh_to_euclidean" begin
    #
    # Two tetrahedra sharing a face, plus an unused vertex
    #
    points = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0],
              [0.0,0.0,1.0], [1.0,1.0,1.0], [9.0,9.0,9.0]]
    tets = [[1,2,3,4], [2,3,4,5]]

    ec = tetmesh_to_euclidean(points, tets)
    @test ec isa EuclideanComplex
    @test ec.ncells == 23
    @test count(==(0), ec.dimensions) == 5
    @test count(==(3), ec.dimensions) == 2
    @test homology(ec) == [1,0,0,0]
    @test "1234" in ec.labels && "2345" in ec.labels
    @test ec.coords[ec.indices["2345"]] == points[2:5]

    ecorph = tetmesh_to_euclidean(points, tets; droporphans=false)
    @test count(==(0), ecorph.dimensions) == 6
    @test homology(ecorph) == [2,0,0,0]

    ec3 = tetmesh_to_euclidean(points, tets; p=3)
    @test homology(ec3) == [1,0,0,0]

    # Matrix form, unsorted cells, and tuple points give the same complex

    ecm = tetmesh_to_euclidean(reduce(hcat, points[1:5]), [4 5; 3 4; 2 3; 1 2])
    @test ecm.labels == ec.labels
    @test ecm.coords == ec.coords
    ect = tetmesh_to_euclidean([Tuple(pt) for pt in points], tets)
    @test ect.labels == ec.labels

    # Surface meshes give two-dimensional complexes in three dimensions

    V, T = trisurface_sphere([0,0,0], 1.0, 1)
    ecs = tetmesh_to_euclidean(V, T)
    @test ecs.dim == 2
    @test euler_char(ecs) == 2
    @test homology(ecs) == [1,0,1]

    # Invalid input

    @test_throws ArgumentError tetmesh_to_euclidean(points, [[1,2,3,7]])
    @test_throws ArgumentError tetmesh_to_euclidean(points, [[1,2,2,3]])
    @test_throws ArgumentError tetmesh_to_euclidean([[0.0,0.0],[1.0,0.0,0.0]], [[1,2]])
    @test_throws ArgumentError tetmesh_to_euclidean(points, Vector{Int}[])
end

@testset "poisson_points_box" begin
    for (bmin, bmax, pdist) in (([0,0], [2,1], 0.05), ([-1,0,0], [1,1,0.5], 0.1))
        pts = poisson_points_box(bmin, bmax, pdist; rng=Xoshiro(123))
        @test length(pts) > 100
        @test all(length(x) == length(bmin) for x in pts)
        @test all(all(bmin .<= x .<= bmax) for x in pts)
        mindist = minimum(sqrt(sum(abs2, pts[i] .- pts[j]))
                          for i in 1:length(pts) for j in i+1:length(pts))
        @test mindist >= pdist
        @test pts == poisson_points_box(bmin, bmax, pdist; rng=Xoshiro(123))
    end

    # The box should be filled: no empty ball of radius 2*pdist

    pts = poisson_points_box([0,0,0], [1,1,1], 0.1; rng=Xoshiro(7))
    rng = Xoshiro(8)
    probes = [0.2 .+ 0.6 .* rand(rng, 3) for _ in 1:200]
    @test all(minimum(sqrt(sum(abs2, q .- x)) for x in pts) < 0.2 for q in probes)

    @test_throws ArgumentError poisson_points_box([0,0], [1,1,1], 0.1)
    @test_throws ArgumentError poisson_points_box([0,0], [1,0], 0.1)
    @test_throws ArgumentError poisson_points_box([0,0], [1,1], 0.0)
end

@testset "trisurface generators" begin
    surfaces = [
        ("sphere",    trisurface_sphere([1,2,3], 2.0, 2),                 true,  2),
        ("ellipsoid", trisurface_ellipsoid([0,0,0], [1,2,3], 1),          true,  2),
        ("cylinder",  trisurface_cylinder([0,0,0], [1,1,1], 0.5, 12, 4),  true,  2),
        ("tube",      trisurface_cylinder([0,0,0], [1,1,1], 0.5, 12, 4;
                                          capped=false),                   false, 0),
        ("torus",     trisurface_torus([0,0,0], 2.0, 0.5, 16, 8;
                                       axis=[1,0,0]),                      true,  0),
        ("plane",     trisurface_plane_patch([0,0,0], [1,0,0], [0,1,0], 3, 4),
                                                                           false, 1),
        ("psphere",   trisurface_parametric((u,v) -> [cos(u)*sin(v), sin(u)*sin(v), cos(v)],
                                            (0,2*pi), (pi,0), 24, 12; uperiodic=true),
                                                                           true,  2),
    ]

    for (name, (V, T), closed, chi) in surfaces
        isclosed, isoriented = ConleyDynamics._trisurface_check(T)
        @test isclosed == closed
        @test isoriented
        @test all(allunique(t) for t in T)
        # No duplicate vertices
        @test length(unique(v -> round.(v, digits=8), V)) == length(V)
        ec = tetmesh_to_euclidean(V, T)
        @test count(==(0), ec.dimensions) == length(V)
        @test euler_char(ec) == chi
        if closed
            @test ConleyDynamics._trisurface_volume(V, T) > 0
        end
    end

    # Sizes and geometry

    V, T = trisurface_sphere([1,2,3], 2.0, 2)
    @test length(V) == 10*4^2 + 2
    @test length(T) == 20*4^2
    @test all(abs(sqrt(sum(abs2, v .- [1,2,3])) - 2.0) < 1e-12 for v in V)

    V, T = trisurface_torus([0,0,0], 2.0, 0.5, 16, 8)
    @test homology(tetmesh_to_euclidean(V, T)) == [1,2,1]
    V, T = trisurface_torus([0,0,0], 2.0, 0.5, 128, 64)
    @test abs(ConleyDynamics._trisurface_volume(V, T) - 2*pi^2*2.0*0.25) < 0.02

    V, T = trisurface_cylinder([0,0,0], [0,0,2], 1.0, 64, 3)
    @test abs(ConleyDynamics._trisurface_volume(V, T) - 2*pi) < 0.02

    # Interior points of closed surfaces

    for (V, T) in (trisurface_sphere([1,2,3], 2.0, 1),
                   trisurface_torus([0,0,0], 2.0, 0.5, 16, 8),
                   trisurface_cylinder([0,0,0], [1,1,1], 0.5, 12, 4))
        q = ConleyDynamics._trisurface_interior_point(V, T)
        @test length(q) == 3
    end
    q = ConleyDynamics._trisurface_interior_point(trisurface_sphere([1,2,3], 2.0, 1)...)
    @test sqrt(sum(abs2, q .- [1,2,3])) < 2.0
    q = ConleyDynamics._trisurface_interior_point(trisurface_torus([0,0,0], 2.0, 0.5, 16, 8)...)
    @test abs(sqrt(q[1]^2 + q[2]^2) - 2.0)^2 + q[3]^2 < 0.25

    # Invalid input

    @test_throws ArgumentError trisurface_sphere([0,0,0], -1.0, 1)
    @test_throws ArgumentError trisurface_torus([0,0,0], 1.0, 2.0, 8, 8)
    @test_throws ArgumentError trisurface_cylinder([0,0,0], [0,0,1], 1.0, 2, 1)
end

@testset "TetDomain builders" begin
    dom = tetgen_domain_box([0,0,0], [2,1,1])
    @test length(dom.points) == 8
    @test length(dom.facets) == 6
    @test dom.markers == zeros(Int, 6)
    @test dom.points[7] == (2.0, 1.0, 1.0)
    # The faces form a closed, outward oriented surface
    tris = [[f[1],f[2],f[3]] for f in dom.facets]
    append!(tris, [[f[1],f[3],f[4]] for f in dom.facets])
    @test ConleyDynamics._trisurface_check(tris) == (true, true)
    @test ConleyDynamics._trisurface_volume(dom.points, tris) ≈ 2.0

    # Nodes, surfaces, curves, holes, and regions

    dom2 = tetgen_domain_add_nodes(dom, [[0.5,0.5,0.5], (1.5,0.5,0.5)])
    @test length(dom2.points) == 10
    @test length(dom.points) == 8

    V, T = trisurface_sphere([1,0.5,0.5], 0.3, 1)
    dom3 = tetgen_domain_add_surface(dom2, V, T; hole=true)
    @test length(dom3.facets) == 6 + length(T)
    @test dom3.markers[end] == 1
    @test length(dom3.holes) == 1
    @test sqrt(sum(abs2, collect(dom3.holes[1]) .- [1,0.5,0.5])) < 0.3

    dom4 = tetgen_domain_add_surface(dom3, trisurface_sphere([0.4,0.5,0.5], 0.2, 1)...;
                                     region=(5, 0.01))
    @test dom4.markers[end] == 2
    @test length(dom4.regions) == 1
    @test dom4.regions[1][2:3] == (5, 0.01)

    curve = [[0.1,0.1,0.1], [0.1,0.9,0.1], [0.1,0.9,0.9], [0.1,0.1,0.1]]
    dom5 = tetgen_domain_add_curve(dom4, curve; marker=9)
    @test dom5.facets[end-2:end] == [[length(dom5.points)-2, length(dom5.points)-1],
                                     [length(dom5.points)-1, length(dom5.points)],
                                     [length(dom5.points), length(dom5.points)-2]]
    @test all(==(9), dom5.markers[end-2:end])
    @test ConleyDynamics._next_marker(dom5) == 10

    # Points shared with existing points are identified

    dom6 = tetgen_domain_add_curve(dom, [[0,0,0], [2,1,1]])
    @test length(dom6.points) == 8
    @test dom6.facets[end] == [1,7]

    dom7 = tetgen_domain_add_hole(tetgen_domain_add_region(dom, [1,1,1], 3), [0,0,0])
    @test dom7.holes == [(0.0,0.0,0.0)]
    @test dom7.regions == [((1.0,1.0,1.0), 3, -1.0)]

    # Domains from closed surfaces

    Vt, Tt = trisurface_torus([0,0,0], 2.0, 0.5, 12, 6)
    domt = tetgen_domain_from_surface(Vt, Tt)
    @test length(domt.facets) == length(Tt)

    # Invalid input

    @test_throws ArgumentError tetgen_domain_box([0,0], [1,1])
    @test_throws ArgumentError tetgen_domain_box([0,0,0], [1,0,1])
    @test_throws ArgumentError tetgen_domain_from_surface(trisurface_plane_patch([0,0,0],[1,0,0],[0,1,0],2,2)...)
    @test_throws ArgumentError tetgen_domain_add_surface(dom, trisurface_plane_patch([0,0,0],[1,0,0],[0,1,0],2,2)...; hole=true)
    @test_throws ArgumentError tetgen_domain_add_nodes(dom, [[0.5,0.5]])
    @test_throws ArgumentError TetDomain([(0.0,0.0,0.0)], [[1,2]], [0])
end

@testset "Meshing input validation" begin
    # Non-finite coordinates and invalid parameters

    @test_throws ArgumentError tetmesh_to_euclidean([[NaN,0,0], [1,0,0], [0,1,0], [0,0,1]],
                                                    [[1,2,3,4]])
    @test_throws ArgumentError tetgen_domain_add_nodes(tetgen_domain_box([0,0,0], [1,1,1]),
                                                       [[Inf,0.5,0.5]])
    @test_throws ArgumentError poisson_points_box([0,0], [1,1], 0.1; attempts=0)

    # The merge tolerance of parametric surfaces has to be positive, and
    # large coordinates must not overflow the merge buckets

    fsph(u,v) = [cos(u)*sin(v), sin(u)*sin(v), cos(v)]
    @test_throws ArgumentError trisurface_parametric(fsph, (0,2*pi), (pi,0), 8, 4;
                                                     uperiodic=true, tol=0)
    @test_throws ArgumentError trisurface_parametric(fsph, (0,2*pi), (pi,0), 8, 4;
                                                     uperiodic=true, tol=-1)
    V1, T1 = trisurface_parametric(fsph, (0,2*pi), (pi,0), 24, 12; uperiodic=true)
    fbig(u,v) = 1e10 .* fsph(u,v)
    @test trisurface_parametric(fbig, (0,2*pi), (pi,0), 24, 12; uperiodic=true) isa Tuple
    V2, T2 = trisurface_parametric(fbig, (0,2*pi), (pi,0), 24, 12; uperiodic=true, tol=1.0)
    @test length(V2) == length(V1) == 266
    @test T2 == T1

    # Holes and regions in surfaces with several components

    Va, Ta = trisurface_sphere([ 1.0,0,0], 0.4, 1)
    Vb, Tb = trisurface_sphere([-1.0,0,0], 0.4, 1)
    VV = vcat(Va, Vb)
    TT = vcat(Ta, [t .+ length(Va) for t in Tb])
    @test length(ConleyDynamics._trisurface_components(TT)) == 2
    box = tetgen_domain_box([-2,-2,-2], [2,2,2])
    dom = tetgen_domain_add_surface(box, VV, TT; hole=true)
    @test length(dom.holes) == 2
    @test sort([round(Int, h[1]) for h in dom.holes]) == [-1, 1]
    @test all(sqrt(sum(abs2, collect(h) .- [round(h[1]),0,0])) < 0.4 for h in dom.holes)
    dom = tetgen_domain_add_surface(box, VV, TT; region=(3, 0.01))
    @test length(dom.regions) == 2
    @test all(r[2] == 3 for r in dom.regions)

    # Holes require consistently oriented surfaces, with either orientation

    V, T = trisurface_sphere([0,0,0], 1.0, 2)
    Tbad = deepcopy(T)
    Tbad[1] = reverse(Tbad[1])
    @test_throws ArgumentError tetgen_domain_add_surface(box, V, Tbad; hole=true)
    dom = tetgen_domain_add_surface(box, V, [reverse(t) for t in T]; hole=true)
    @test sqrt(sum(abs2, collect(dom.holes[1]))) < 1.0

    # Degenerate and repeated curve segments

    dom = tetgen_domain_add_curve(box, [[0,0,0], [1,0,0], [0,0,0]])
    @test length(dom.facets) == 7
    dom = tetgen_domain_add_curve(box, [[0,0,0], [0,0,0], [1,0,0], [1,0,0], [1,1,0]])
    @test length(dom.facets) == 8
    @test all(allunique(f) for f in dom.facets)
    @test_throws ArgumentError tetgen_domain_add_curve(box, [[0,0,0], [0,0,0]])
end

@testset "Meshing inputs: tuples, indices, thin boxes" begin
    box = tetgen_domain_box([-2,-2,-2], [2,2,2])
    V, T = trisurface_sphere([0,0,0], 1.0, 1)
    Vt = [Tuple(v) for v in V]
    Tt = [Tuple(t) for t in T]

    # Surfaces and cells given as tuples

    @test length(tetgen_domain_from_surface(Vt, Tt).facets) == length(T)
    dom = tetgen_domain_add_surface(box, Vt, Tt; hole=true)
    @test length(dom.holes) == 1
    @test sqrt(sum(abs2, collect(dom.holes[1]))) < 1.0
    @test tetmesh_to_euclidean(Vt, Tt).ncells == tetmesh_to_euclidean(V, T).ncells

    # Invalid triangle indices, and holes combined with regions

    @test_throws ArgumentError tetgen_domain_add_surface(box, V, [t .- 1 for t in T])
    @test_throws ArgumentError tetgen_domain_add_surface(box, V[1:10], T)
    @test_throws ArgumentError tetgen_domain_from_surface(V[1:10], T)
    @test_throws ArgumentError tetgen_domain_add_surface(box, V, T; hole=true,
                                                         region=(1, 0.1))

    # Boxes which are thinner than pdist in some directions are filled
    # as lower-dimensional boxes

    mindist(p) = minimum(sqrt(sum(abs2, p[i] .- p[j]))
                         for i in 1:length(p) for j in i+1:length(p))
    n2 = length(poisson_points_box([0,0], [1,1], 0.1; rng=Xoshiro(1)))
    for h in (0.01, 0.05)
        p = poisson_points_box([0,0,0], [1,1,h], 0.1; rng=Xoshiro(1))
        @test length(p) >= 0.8 * n2
        @test mindist(p) >= 0.1
        @test all(all([0,0,0] .<= x .<= [1,1,h]) for x in p)
    end
    p = poisson_points_box([0,0,0], [1,0.01,0.01], 0.1; rng=Xoshiro(1))
    @test length(p) >= 5
    @test mindist(p) >= 0.1
end

@testset "Meshing at extreme scales" begin
    #
    # Interior points and closed curves do not depend on the scale
    #
    for s in (1e-12, 1e-8, 1.0, 1e8, 1e12)
        V, T = trisurface_sphere([s,0,0], s, 1)
        x = ConleyDynamics._trisurface_interior_point(V, T)
        @test sqrt(sum(abs2, x .- [s,0,0])) < 0.9 * s
        dom = tetgen_domain_box(-2s .* [1,1,1], 2s .* [1,1,1])
        dom = tetgen_domain_add_surface(dom, V, T; hole=true)
        @test length(dom.holes) == 1
        circ = [[1.5s*cos(t), 1.5s*sin(t), 0.0] for t in range(0, 2*pi, length=30)]
        dom = tetgen_domain_add_curve(tetgen_domain_box(-2s .* [1,1,1], 2s .* [1,1,1]), circ)
        @test count(f -> length(f) == 2, dom.facets) == 29
        @test length(dom.points) == 8 + 29
    end
    # Endpoints of a tiny open curve are not identified
    dom = tetgen_domain_add_curve(tetgen_domain_box([0,0,0], [1,1,1]),
                                  [[0.5,0.5,0.5], [0.5,0.5+1e-10,0.5], [0.5+1e-10,0.5,0.5]])
    @test count(f -> length(f) == 2, dom.facets) == 2
end
