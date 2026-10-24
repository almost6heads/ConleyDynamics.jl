function gmsh_box_mesh(h)
    gmsh.model.add("box")
    gmsh.model.occ.addBox(0, 0, 0, 2, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(3)
end

@testset "Gmsh extension: conversion" begin
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    # A box

    gmsh_box_mesh(0.3)
    ec = gmsh_to_simplicial()
    @test ec isa EuclideanComplex
    @test ec.dim == 3
    @test homology(ec) == [1,0,0,0]
    tets = [c for (c, d) in zip(ec.coords, ec.dimensions) if d == 3]
    vol(c) = abs(sum((c[2] .- c[1]) .* ConleyDynamics._cross3(c[3] .- c[1], c[4] .- c[1]))) / 6
    @test abs(sum(vol.(tets)) - 2.0) < 1e-10
    @test homology(gmsh_to_simplicial(dim=2)) == [1,0,1]
    @test homology(gmsh_to_simplicial(p=0)) == [1,0,0,0]

    # Boundary cells via an entity and via a physical group

    surfaces = [t for (d, t) in gmsh.model.getEntities(2)]
    pg = gmsh.model.addPhysicalGroup(2, surfaces)
    bnd = gmsh_marker_cells(ec, 2, pg)
    @test homology(lefschetz_subcomplex(ec, bnd)) == [1,0,1]
    face1 = gmsh_marker_cells(ec, 2, surfaces[1]; physical=false)
    @test homology(lefschetz_subcomplex(ec, face1)) == [1,0,0]
    @test issubset(face1, bnd)
    face1tri = gmsh_marker_cells(ec, 2, surfaces[1]; physical=false, closure=false)
    @test all(ec.dimensions[face1tri] .== 2)

    # Physical group of the volume

    pv = gmsh.model.addPhysicalGroup(3, [1])
    @test gmsh_to_simplicial(physical=(3, pv)).ncells == ec.ncells
    @test_throws ArgumentError gmsh_to_simplicial(physical=(2, pg))

    # Second order elements are rejected

    gmsh.model.mesh.setOrder(2)
    @test_throws ArgumentError gmsh_to_simplicial()
    gmsh.model.remove()

    # A box with a spherical hole

    gmsh.model.add("hole")
    box = gmsh.model.occ.addBox(-2, -2, -2, 4, 4, 4)
    sph = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
    gmsh.model.occ.cut([(3, box)], [(3, sph)])
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.6)
    gmsh.model.mesh.generate(3)
    @test homology(gmsh_to_simplicial()) == [1,0,1,0]
    gmsh.model.remove()

    # No mesh

    gmsh.model.add("empty")
    @test_throws ArgumentError gmsh_to_simplicial()
    gmsh.model.remove()

    gmsh.finalize()
    @test_throws ArgumentError gmsh_to_simplicial()
end

@testset "Gmsh extension: intersecting surfaces" begin
    #
    # The fragment example from the manual: a box cut by a plane, with
    # an embedded sphere which intersects the plane in a circle
    #
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("intersecting")

    box    = gmsh.model.occ.addBox(-2, -2, -2, 4, 4, 4)
    sphere = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
    plane  = gmsh.model.occ.addRectangle(-2, -2, 0.3, 4, 4)
    gmsh.model.occ.synchronize()
    sphsurf = gmsh.model.getBoundary([(3, sphere)], false, false)
    gmsh.model.occ.remove([(3, sphere)])
    out, outmap = gmsh.model.occ.fragment([(3, box)], vcat([(2, plane)], sphsurf))
    gmsh.model.occ.synchronize()

    plane_pieces  = [tag for (d, tag) in outmap[2] if d == 2]
    sphere_pieces = [tag for m in outmap[3:end] for (d, tag) in m if d == 2]
    pg_plane  = gmsh.model.addPhysicalGroup(2, plane_pieces)
    pg_sphere = gmsh.model.addPhysicalGroup(2, sphere_pieces)

    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)
    gmsh.model.mesh.generate(3)

    ec = gmsh_to_simplicial()
    @test homology(ec) == [1,0,0,0]
    plane_cells  = gmsh_marker_cells(ec, 2, pg_plane)
    sphere_cells = gmsh_marker_cells(ec, 2, pg_sphere)
    @test homology(lefschetz_subcomplex(ec, plane_cells)) == [1,0,0]
    @test homology(lefschetz_subcomplex(ec, sphere_cells)) == [1,0,1]
    circle = intersect(plane_cells, sphere_cells)
    @test homology(lefschetz_subcomplex(ec, circle)) == [1,1]
    @test all(abs(ec.coords[k][1][3] - 0.3) < 1e-12 for k in circle if ec.dimensions[k] == 0)

    gmsh.finalize()
end

@testset "Gmsh extension: files" begin
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh_box_mesh(0.3)
    ec = gmsh_to_simplicial()
    dir = mktempdir()
    mshfile = joinpath(dir, "box.msh")
    gmsh.write(mshfile)
    gmsh.finalize()

    # Mesh file, with Gmsh not running

    ecf = gmsh_file_to_simplicial(mshfile)
    @test gmsh.isInitialized() == 0
    @test ecf.ncells == ec.ncells
    @test homology(ecf) == [1,0,0,0]

    # Geometry file, which is meshed automatically

    geofile = joinpath(dir, "ball.geo")
    write(geofile, """
        SetFactory("OpenCASCADE");
        Sphere(1) = {0, 0, 0, 1.0};
        Mesh.MeshSizeMax = 0.4;
        """)
    ecb = gmsh_file_to_simplicial(geofile)
    @test homology(ecb) == [1,0,0,0]
    @test homology(gmsh_file_to_simplicial(geofile; dim=2)) == [1,0,1]
    @test_throws ArgumentError gmsh_file_to_simplicial(joinpath(dir, "missing.msh"))
end

@testset "Gmsh extension: non-conforming meshes" begin
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)

    # Two touching boxes without fragment have duplicate nodes

    gmsh.model.add("touching")
    gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.addBox(1, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    ec = @test_logs (:warn, r"identical coordinates") gmsh_to_simplicial()
    @test homology(ec) == [2,0,0,0]
    gmsh.model.remove()

    # With fragment, the mesh is conforming

    gmsh.model.add("fragmented")
    b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    b2 = gmsh.model.occ.addBox(1, 0, 0, 1, 1, 1)
    gmsh.model.occ.fragment([(3, b1)], [(3, b2)])
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    ec = @test_logs gmsh_to_simplicial()
    @test homology(ec) == [1,0,0,0]
    @test homology(gmsh_to_simplicial(dim=1)) == [1,9]

    # A volume without elements, as after a partially failed meshing

    gmsh.model.mesh.clear([(3, b2)])
    ec = @test_logs (:warn, r"contain no mesh elements") gmsh_to_simplicial()
    @test homology(ec) == [1,0,0,0]
    tetvol(c) = abs(sum((c[2] .- c[1]) .* ConleyDynamics._cross3(c[3] .- c[1], c[4] .- c[1]))) / 6
    @test sum(tetvol(c) for (c, d) in zip(ec.coords, ec.dimensions) if d == 3) ≈ 1.0

    gmsh.finalize()
end

@testset "Gmsh extension: common refinements" begin
    #
    # Intersecting constraint surfaces, curves, and volumes are split
    # along their intersections by gmsh_common_refinement
    #
    sub(ec, c) = homology(lefschetz_subcomplex(ec, c))
    verts(ec, c) = [ec.coords[k][1] for k in c if ec.dimensions[k] == 0]
    function spheresurf(x, y, z, r)
        s = gmsh.model.occ.addSphere(x, y, z, r)
        gmsh.model.occ.synchronize()
        bnd = gmsh.model.getBoundary([(3, s)], false, false)
        gmsh.model.occ.remove([(3, s)])
        return bnd
    end
    function meshmodel(h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.model.mesh.generate(3)
        return @test_logs gmsh_to_simplicial()
    end

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    # Two intersecting spheres, and a line through both of them

    gmsh.model.add("spheres")
    box = gmsh.model.occ.addBox(-2, -2, -2, 4, 4, 4)
    sa  = spheresurf(-0.5, 0, 0, 1.0)
    sb  = spheresurf( 0.5, 0, 0, 1.0)
    ln  = gmsh.model.occ.addLine(gmsh.model.occ.addPoint(-1.9, 0, 0),
                                 gmsh.model.occ.addPoint( 1.9, 0, 0))
    pgs = gmsh_common_refinement([(3, box), sa, sb, (1, ln)];
                                 names=["box", "left", "right", "axis"])
    @test length(pgs) == 4
    @test gmsh.model.getPhysicalName(2, pgs[3]) == "right"
    ec = meshmodel(0.35)
    @test homology(ec) == [1,0,0,0]
    @test length(gmsh_marker_cells(ec, 3, pgs[1])) == ec.ncells
    ca = gmsh_marker_cells(ec, 2, pgs[2])
    cb = gmsh_marker_cells(ec, 2, pgs[3])
    cl = gmsh_marker_cells(ec, 1, pgs[4])
    @test sub(ec, ca) == [1,0,1]
    @test sub(ec, cb) == [1,0,1]
    @test sub(ec, union(ca, cb)) == [1,0,3]
    circle = intersect(ca, cb)
    @test sub(ec, circle) == [1,1]
    @test all(abs(v[1]) < 1e-10 && abs(v[2]^2 + v[3]^2 - 0.75) < 1e-8 for v in verts(ec, circle))
    @test sub(ec, cl) == [1,0]
    xs = sort([v[1] for v in verts(ec, intersect(cl, union(ca, cb)))])
    @test length(xs) == 4
    @test maximum(abs.(xs .- [-1.5, -0.5, 0.5, 1.5])) < 1e-10
    gmsh.model.remove()

    # Three coordinate planes meet in a single vertex

    gmsh.model.add("planes")
    box = gmsh.model.occ.addBox(-1, -1, -1, 2, 2, 2)
    p1  = gmsh.model.occ.addRectangle(-1, -1, 0, 2, 2)
    p2  = gmsh.model.occ.addRectangle(-1, -1, 0, 2, 2)
    p3  = gmsh.model.occ.addRectangle(-1, -1, 0, 2, 2)
    gmsh.model.occ.rotate([(2, p2)], 0, 0, 0, 1, 0, 0, pi/2)
    gmsh.model.occ.rotate([(2, p3)], 0, 0, 0, 0, 1, 0, pi/2)
    pgs = gmsh_common_refinement([(3, box), (2, p1), (2, p2), (2, p3)])
    ec = meshmodel(0.5)
    c = [gmsh_marker_cells(ec, 2, t) for t in pgs[2:4]]
    @test all(sub(ec, ci) == [1,0,0] for ci in c)
    @test all(sub(ec, intersect(c[i], c[j])) == [1,0] for (i, j) in ((1,2), (1,3), (2,3)))
    center = intersect(c...)
    @test length(center) == 1
    @test maximum(abs.(verts(ec, center)[1])) < 1e-12
    gmsh.model.remove()

    # A plane which extends beyond the box, and a ball as a volume region

    gmsh.model.add("region")
    box  = gmsh.model.occ.addBox(-2, -2, -2, 4, 4, 4)
    pl   = gmsh.model.occ.addRectangle(-3, -3, 0.2, 6, 6)
    ball = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
    pgbox, pgpl, pgball = gmsh_common_refinement([(3, box), (2, pl), (3, ball)])
    ec = meshmodel(0.4)
    cpl   = gmsh_marker_cells(ec, 2, pgpl)
    cball = gmsh_marker_cells(ec, 3, pgball)
    @test sub(ec, cpl) == [1,0,0]
    @test all(maximum(abs.(v)) <= 2 + 1e-12 for v in verts(ec, cpl))
    @test sub(ec, cball) == [1,0,0,0]
    @test sub(ec, intersect(cpl, cball)) == [1,0,0]
    @test all(sum(abs2, sum(ec.coords[k]) ./ 4) < 1 for k in cball if ec.dimensions[k] == 3)
    gmsh.model.remove()

    # Errors and a single object

    gmsh.model.add("errors")
    b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    r1 = gmsh.model.occ.addRectangle(0, 0, 0.5, 1, 1)
    @test_throws ArgumentError gmsh_common_refinement([])
    @test_throws ArgumentError gmsh_common_refinement([[(3, b1), (2, r1)]])
    @test_throws ArgumentError gmsh_common_refinement([(3, b1)]; names=["a", "b"])
    @test_throws ArgumentError gmsh_common_refinement([Tuple{Int,Int}[]])
    @test length(gmsh_common_refinement([(3, b1)])) == 1

    gmsh.finalize()
    @test_throws ArgumentError gmsh_common_refinement([(3, 1)])
end

@testset "Gmsh extension: triangulated surfaces" begin
    #
    # Triangulated surfaces are added as flat OpenCASCADE surfaces by
    # gmsh_add_trisurface, and refined together with other objects
    #
    sub(ec, c) = homology(lefschetz_subcomplex(ec, c))
    verts(ec, c) = [ec.coords[k][1] for k in c if ec.dimensions[k] == 0]
    function meshmodel(h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.model.mesh.generate(3)
        return @test_logs gmsh_to_simplicial()
    end

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    # A paraboloid patch, a sphere, and a rectangle

    fpar(x, y) = -0.6 + 0.5*(x^2 + y^2)
    V, T = trisurface_parametric((x, y) -> [x, y, fpar(x, y)],
                                 (-1.5, 1.5), (-1.5, 1.5), 12, 12)
    gmsh.model.add("paraboloid")
    box  = gmsh.model.occ.addBox(-2, -2, -2, 4, 4, 4)
    sph  = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
    gmsh.model.occ.synchronize()
    ssph = gmsh.model.getBoundary([(3, sph)], false, false)
    gmsh.model.occ.remove([(3, sph)])
    rect = gmsh.model.occ.addRectangle(-1.6, -1.6, 0, 3.2, 3.2)
    gmsh.model.occ.rotate([(2, rect)], 0, 0, 0, 0, 1, 0, pi/2)
    gmsh.model.occ.translate([(2, rect)], 0.4, 0, 0)
    spar = gmsh_add_trisurface(V, T)
    @test length(spar) == length(T)
    pgs = gmsh_common_refinement([(3, box), ssph, (2, rect), spar])
    ec = meshmodel(0.35)
    @test homology(ec) == [1,0,0,0]
    cs, cr, cp = (gmsh_marker_cells(ec, 2, t) for t in pgs[2:4])
    @test sub(ec, cp) == [1,0,0]
    @test sub(ec, intersect(cs, cp)) == [1,1]
    @test sub(ec, intersect(cr, cp)) == [1,0]
    @test all(0 <= v[3] - fpar(v[1], v[2]) < 0.05 for v in verts(ec, cp))
    triple = verts(ec, intersect(cs, cr, cp))
    @test length(triple) == 2
    @test all(abs(v[1] - 0.4) < 1e-10 && abs(sum(abs2, v) - 1) < 1e-8 for v in triple)
    gmsh.model.remove()

    # A closed surface given as a triangle soup, with separate copies of
    # the vertices of every triangle

    Vs, Ts = trisurface_sphere([0, 0, 0], 1.0, 1)
    Vsoup = [Vs[k] for t in Ts for k in t]
    Tsoup = [[3k-2, 3k-1, 3k] for k in 1:length(Ts)]
    gmsh.model.add("soup")
    box  = gmsh.model.occ.addBox(-2, -2, -2, 4, 4, 4)
    soup = gmsh_add_trisurface(Vsoup, Tsoup)
    pgbox, pgsoup = gmsh_common_refinement([(3, box), soup])
    ec = meshmodel(0.5)
    csoup = gmsh_marker_cells(ec, 2, pgsoup)
    @test sub(ec, csoup) == [1,0,1]
    @test length(verts(ec, csoup)) >= length(Vs)
    gmsh.model.remove()

    # Errors

    gmsh.model.add("errors")
    Vt = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [2, 0, 0], [0, 0, 0]]
    @test_throws ArgumentError gmsh_add_trisurface(Vt, Vector{Int}[])
    @test_throws ArgumentError gmsh_add_trisurface(Vt, [[1, 2]])
    @test_throws ArgumentError gmsh_add_trisurface(Vt, [[1, 2, 6]])
    @test_throws ArgumentError gmsh_add_trisurface(Vt, [[1, 2, 5]])
    @test_throws ArgumentError gmsh_add_trisurface(Vt, [[1, 2, 4]])
    @test_throws ArgumentError gmsh_add_trisurface(Vt, [[1, 2, 3], [3, 1, 2]])
    @test_throws ArgumentError gmsh_add_trisurface([[0, 0], [1, 0], [0, 1]], [[1, 2, 3]])
    @test_throws ArgumentError gmsh_add_trisurface([[0, 0, NaN], [1, 0, 0], [0, 1, 0]],
                                                   [[1, 2, 3]])

    gmsh.finalize()
    @test_throws ArgumentError gmsh_add_trisurface(Vt, [[1, 2, 3]])
end
