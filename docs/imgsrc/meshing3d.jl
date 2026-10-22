#
# Figures for the manual sections on three-dimensional meshes
# (extensions.md) and three-dimensional plotting (plotting.md).
#
# The script needs an environment which contains ConleyDynamics,
# TetGen, Gmsh, and CairoMakie. It writes the PNG files into the
# directory docs/src/man/img. Run it from the docs/imgsrc directory:
#
#    julia --project=<env> -t 4 meshing3d.jl
#
# All random point sets use fixed seeds, and the meshes created by
# TetGen and Gmsh are deterministic, so the figures can be regenerated.
# The whole script takes a few minutes.
#

using Random
using TetGen
using Gmsh: gmsh
using CairoMakie
using ConleyDynamics

imgdir = joinpath(@__DIR__, "..", "src", "man", "img")
savefig3d(name, fig) = save(joinpath(imgdir, name), fig; px_per_unit=1.5)

function panel(fig, pos, title)
    #
    # An Axis3 with equal aspect ratios and without decorations
    #
    ax = Axis3(fig[pos...]; aspect=:data, title=title, titlesize=18)
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

##############################################################################
#
# Extensions: triangulated surfaces
#
##############################################################################

Vs, Ts = trisurface_sphere([0,0,0], 1.0, 2)
Ve, Te = trisurface_ellipsoid([0,0,0], [2.0,1.0,0.5], 2)
Vc, Tc = trisurface_cylinder([0,0,-1], [0,0,1], 0.5, 24, 6; capped=false)
Vt, Tt = trisurface_torus([0,0,0], 2.0, 0.6, 32, 12; axis=[1,0,1])
Vp, Tp = trisurface_plane_patch([-1,-1,0], [2,0,0], [0.5,2,0.5], 8, 8)
fmoebius(u,v) = [(1+v/2*cos(u/2))*cos(u), (1+v/2*cos(u/2))*sin(u), v/2*sin(u/2)]
Vm, Tm = trisurface_parametric(fmoebius, (0,2*pi), (-1,1), 40, 4)

fig = Figure(size=(1200,800))
for (k, (V, T, title)) in enumerate([(Vs, Ts, "trisurface_sphere"),
                                     (Ve, Te, "trisurface_ellipsoid"),
                                     (Vc, Tc, "trisurface_cylinder"),
                                     (Vt, Tt, "trisurface_torus"),
                                     (Vp, Tp, "trisurface_plane_patch"),
                                     (Vm, Tm, "trisurface_parametric")])
    ax = panel(fig, (div(k-1,3)+1, mod1(k,3)), title)
    plot_simplicial3d!(ax, tetmesh_to_euclidean(V, T); alpha=0.6)
end
savefig3d("meshing3d_surfaces.png", fig)

##############################################################################
#
# Extensions: TetGen box
#
##############################################################################

dom = tetgen_domain_box([-1.5,-1,-1], [1.5,1,1])
res = tetgen_mesh(dom; maxvolume=0.01)
ecbox = tetgen_to_simplicial(res)

fig = Figure(size=(1200,500))
ax1 = panel(fig, (1,1), "Boundary surface")
plot_simplicial3d!(ax1, ecbox)
ax2 = panel(fig, (1,2), "Clipped at y = 0")
plot_simplicial3d!(ax2, ecbox; clip=([0,1,0], 0.0), alpha=0.6)
savefig3d("meshing3d_box.png", fig)

##############################################################################
#
# Extensions: constraint surfaces and curves
#
##############################################################################

dom2 = tetgen_domain_box([-2,-2,-2], [2,2,2])
V, T = trisurface_sphere([0,0,0], 1.0, 2)
dom2 = tetgen_domain_add_surface(dom2, V, T)
square = [[-1.8,-1.8,1.5], [1.8,-1.8,1.5], [1.8,1.8,1.5],
          [-1.8,1.8,1.5], [-1.8,-1.8,1.5]]
dom2 = tetgen_domain_add_curve(dom2, square)
spoke = [[0.0,0.0,-1.8], [0.0,0.0,-1.3]]
dom2 = tetgen_domain_add_curve(dom2, spoke; marker=7)
res2 = tetgen_mesh(dom2; maxvolume=0.05)
ec2  = tetgen_to_simplicial(res2)

sphere_cells = tetgen_marker_cells(res2, ec2, 1)
square_cells = tetgen_marker_cells(res2, ec2, 2)
spoke_cells  = tetgen_marker_cells(res2, ec2, 7)

fig = Figure(size=(1200,550))
ax1 = panel(fig, (1,1), "Mesh clipped at x = 0")
plot_simplicial3d!(ax1, ec2; clip=([1,0,0], 0.0), alpha=0.25)
plot_cells3d!(ax1, ec2, sphere_cells; clip=([1,0,0], 0.0), color=:darkorange)
ax2 = panel(fig, (1,2), "Constraint cells")
plot_simplicial3d!(ax2, ec2; alpha=0.05, linealpha=0.1)
plot_cells3d!(ax2, ec2, sphere_cells; color=:darkorange, alpha=0.8)
plot_cells3d!(ax2, ec2, square_cells; color=:red)
plot_cells3d!(ax2, ec2, spoke_cells; color=:green)
savefig3d("meshing3d_constraints.png", fig)

##############################################################################
#
# Extensions: holes
#
##############################################################################

dom3 = tetgen_domain_box([-2,-2,-2], [2,2,2])
dom3 = tetgen_domain_add_surface(dom3, V, T; hole=true)
ec3  = tetgen_to_simplicial(tetgen_mesh(dom3; maxvolume=0.05))

V1, T1 = trisurface_sphere([0,0,0], 1.0, 2)
V2, T2 = trisurface_sphere([0,0,0], 0.5, 1)
dom4 = tetgen_domain_box([-2,-2,-2], [2,2,2])
dom4 = tetgen_domain_add_surface(dom4, V1, T1)
dom4 = tetgen_domain_add_surface(dom4, V2, T2)
dom4 = tetgen_domain_add_hole(dom4, [0.75,0,0])
ec4  = tetgen_to_simplicial(tetgen_mesh(dom4; maxvolume=0.05))

Vh, Th = trisurface_torus([0,0,0], 1.2, 0.4, 32, 12)
dom5 = tetgen_domain_box([-2,-2,-2], [2,2,2])
dom5 = tetgen_domain_add_surface(dom5, Vh, Th; hole=true)
ec5  = tetgen_to_simplicial(tetgen_mesh(dom5; maxvolume=0.05))

fig = Figure(size=(1500,500))
for (k, (ec, title)) in enumerate([(ec3, "Spherical hole"),
                                   (ec4, "Spherical shell removed"),
                                   (ec5, "Toroidal hole")])
    ax = panel(fig, (1,k), title)
    plot_simplicial3d!(ax, ec; clip=([0,1,0], 0.0), alpha=0.6)
end
savefig3d("meshing3d_holes.png", fig)

##############################################################################
#
# Extensions: regions and domains bounded by surfaces
#
##############################################################################

dom6 = tetgen_domain_box([-2,-2,-2], [2,2,2])
dom6 = tetgen_domain_add_surface(dom6, V, T; region=(1, 0.002))
res6 = tetgen_mesh(dom6; maxvolume=0.1)
ec6  = tetgen_to_simplicial(res6)
inner = tetgen_region_cells(res6, ec6, 1)
outer = tetgen_region_cells(res6, ec6, 2)

fig = Figure(size=(1200,550))
ax1 = panel(fig, (1,1), "Refined region, clipped at y = 0")
plot_simplicial3d!(ax1, ec6; clip=([0,1,0], 0.0), alpha=0.6)
ax2 = panel(fig, (1,2), "tetgen_region_cells")
plot_cells3d!(ax2, ec6, inner; clip=([0,1,0], 0.0), color=:darkorange)
plot_cells3d!(ax2, ec6, outer; clip=([0,1,0], 0.0), color=:steelblue, alpha=0.3)
savefig3d("meshing3d_regions.png", fig)

Vt, Tt = trisurface_torus([0,0,0], 2.0, 0.8, 48, 16)
domt = tetgen_domain_from_surface(Vt, Tt)
ect  = tetgen_to_simplicial(tetgen_mesh(domt; maxvolume=0.02))

Ve, Te = trisurface_ellipsoid([0,0,0], [2.0,1.0,0.5], 3)
dome = tetgen_domain_from_surface(Ve, Te)
dome = tetgen_domain_add_surface(dome, trisurface_sphere([ 0.8,0,0], 0.3, 2)...; hole=true)
dome = tetgen_domain_add_surface(dome, trisurface_sphere([-0.8,0,0], 0.3, 2)...; hole=true)
ece  = tetgen_to_simplicial(tetgen_mesh(dome; maxvolume=0.01))

fig = Figure(size=(1200,500))
ax1 = panel(fig, (1,1), "Solid torus, clipped at y = 0")
plot_simplicial3d!(ax1, ect; clip=([0,1,0], 0.0), alpha=0.6)
ax2 = panel(fig, (1,2), "Ellipsoid with two holes, clipped at y = 0")
plot_simplicial3d!(ax2, ece; clip=([0,1,0], 0.0), alpha=0.6)
savefig3d("meshing3d_fromsurface.png", fig)

##############################################################################
#
# Extensions: Delaunay triangulations of random points
#
##############################################################################

ecd1 = create_simplicial_delaunay3d([0,0,0], [2,1,1], 0.15, 20; rng=Xoshiro(1))
ecd2 = create_simplicial_delaunay3d([0,0,0], [2,1,1], 300; box_corners=true,
                                    rng=Xoshiro(1))

fig = Figure(size=(1200,450))
ax1 = panel(fig, (1,1), "Poisson-disk sample, clipped at y = 0.5")
plot_simplicial3d!(ax1, ecd1; clip=([0,1,0], 0.5), pdim=[true,true,true,true])
ax2 = panel(fig, (1,2), "Uniform sample with box corners, clipped at y = 0.5")
plot_simplicial3d!(ax2, ecd2; clip=([0,1,0], 0.5), pdim=[true,true,true,true])
savefig3d("meshing3d_delaunay.png", fig)

##############################################################################
#
# Extensions: a Conley index example on a TetGen mesh
#
##############################################################################

function samplevf(x::Vector{Float64})
    x1, x2, x3 = x
    return [x1 * (1.0 - x1*x1), -x2, -x3]
end

domc = tetgen_domain_box([-1.5,-1,-1], [1.5,1,1])
for (k, c) in enumerate([[-1.0,0,0], [0.0,0,0], [1.0,0,0]])
    Vk, Tk = trisurface_sphere(c, 0.25, 1)
    global domc = tetgen_domain_add_surface(domc, Vk, Tk; region=(k, 0.0005))
end
resc = tetgen_mesh(domc; maxvolume=0.02)
ecc  = tetgen_to_simplicial(resc)
mvfc = create_spatial_mvf(ecc, samplevf)
cmc  = connection_matrix(ecc, mvfc)
println("TetGen Conley example: ", ecc.ncells, " cells, ", cmc.conley, " ", full_from_sparse(cmc.matrix))

fig = Figure(size=(1200,550))
ax1 = panel(fig, (1,1), "Morse sets")
plot_simplicial3d_morse!(ax1, ecc, cmc.morse; ci=true)
ax2 = panel(fig, (1,2), "Refined mesh, clipped at y = 0")
plot_simplicial3d!(ax2, ecc; clip=([0,1,0], 0.0), alpha=0.6)
savefig3d("meshing3d_conley.png", fig)

##############################################################################
#
# Extensions: Gmsh with intersecting surfaces
#
##############################################################################

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
volumes       = [tag for (d, tag) in outmap[1] if d == 3]
upper = [v for v in volumes if gmsh.model.occ.getCenterOfMass(3, v)[3] > 0.3]
pg_plane  = gmsh.model.addPhysicalGroup(2, plane_pieces)
pg_sphere = gmsh.model.addPhysicalGroup(2, sphere_pieces)
pg_upper  = gmsh.model.addPhysicalGroup(3, upper)

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
gmsh.model.mesh.generate(3)

ecg = gmsh_to_simplicial()
plane_cells  = gmsh_marker_cells(ecg, 2, pg_plane)
sphere_cells = gmsh_marker_cells(ecg, 2, pg_sphere)
circle_cells = intersect(plane_cells, sphere_cells)
ecu = gmsh_to_simplicial(physical=(3, pg_upper))
gmsh.finalize()

fig = Figure(size=(1200,550))
ax1 = panel(fig, (1,1), "Plane and sphere meeting in a circle")
plot_simplicial3d!(ax1, ecg; alpha=0.05, linealpha=0.1)
plot_cells3d!(ax1, ecg, plane_cells; color=:steelblue, alpha=0.5)
plot_cells3d!(ax1, ecg, sphere_cells; color=:darkorange, alpha=0.8)
plot_cells3d!(ax1, ecg, [k for k in circle_cells if ecg.dimensions[k] == 1];
              color=:red)
ax2 = panel(fig, (1,2), "Upper volume, clipped at y = 0")
plot_simplicial3d!(ax2, ecu; clip=([0,1,0], 0.0), alpha=0.6)
savefig3d("meshing3d_gmsh.png", fig)

##############################################################################
#
# Extensions: a periodic orbit in a ball meshed by Gmsh
#
##############################################################################

function cyclevf(x::Vector{Float64})
    x1, x2, x3 = x
    r2 = x1*x1 + x2*x2
    return [x1 - 0.5*x2 - x1*r2, 0.5*x1 + x2 - x2*r2, -x3]
end

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.occ.addSphere(0, 0, 0, 2.0)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
gmsh.model.mesh.generate(3)
ecb = gmsh_to_simplicial()
gmsh.finalize()

mvfb = create_spatial_mvf(ecb, cyclevf)
cmb  = connection_matrix(ecb, mvfb)
println("Gmsh ball: ", ecb.ncells, " cells, ", cmb.conley)

fig = Figure(size=(1200,550))
ax1 = panel(fig, (1,1), "Morse sets")
plot_simplicial3d_morse!(ax1, ecb, cmb.morse; ci=true)
ax2 = panel(fig, (1,2), "Morse sets, clipped at z = 0")
plot_simplicial3d_morse!(ax2, ecb, cmb.morse; ci=true, clip=([0,0,1], 0.0))
savefig3d("meshing3d_cycle.png", fig)

##############################################################################
#
# Plotting: basic plots of a complex
#
##############################################################################

ecp = tetgen_to_simplicial(tetgen_mesh(tetgen_domain_box([-1.5,-1,-1], [1.5,1,1]);
                                       maxvolume=0.05))

fig = Figure(size=(1500,500))
ax1 = panel(fig, (1,1), "Default")
plot_simplicial3d!(ax1, ecp)
ax2 = panel(fig, (1,2), "Clipped, with vertices")
plot_simplicial3d!(ax2, ecp; clip=([1,1,1], 0.0), pdim=[true,true,true,true],
                   alpha=0.6)
ax3 = panel(fig, (1,3), "All faces, custom colors")
plot_simplicial3d!(ax3, ecp; hull_only=false, alpha=0.05,
                   facecolor=:gold, edgecolor=:gray40, linewidth=0.5)
savefig3d("makie3d_basic.png", fig)

##############################################################################
#
# Plotting: Morse sets, multivectors, and multivector fields
#
##############################################################################

mvfp = create_spatial_mvf(ecp, samplevf)
msp  = morse_sets(ecp, mvfp)

fig = Figure(size=(1500,500))
ax1 = panel(fig, (1,1), "Default colors")
plot_simplicial3d_morse!(ax1, ecp, msp)
ax2 = panel(fig, (1,2), "Conley index colors")
plot_simplicial3d_morse!(ax2, ecp, msp; ci=true)
ax3 = panel(fig, (1,3), "Custom colors and opacity")
plot_simplicial3d_morse!(ax3, ecp, msp; colors=[:purple, :gold, :seagreen],
                         alpha=0.5, bgalpha=0.02, pdim=[false,false,false,true])
savefig3d("makie3d_morse.png", fig)

mvlargest = mvfp[argmax(length.(mvfp))]

fig = Figure(size=(1500,500))
ax1 = panel(fig, (1,1), "Largest multivector")
plot_simplicial3d_mv!(ax1, ecp, mvlargest)
ax2 = panel(fig, (1,2), "Multivector field, clipped")
plot_simplicial3d_mvf!(ax2, ecp, mvfp; clip=([0,0,1], 0.0))
ax3 = panel(fig, (1,3), "Multivectors by Conley index, clipped")
plot_simplicial3d_mvf!(ax3, ecp, mvfp; ci=true, clip=([0,0,1], 0.0))
savefig3d("makie3d_mvf.png", fig)

##############################################################################
#
# Plotting: cubical complexes
#
##############################################################################

lcq, cq = create_cubical_box(6, 4, 4)
ecq = lefschetz_to_euclidean(lcq, convert_spatial_coordinates(cq, [-1.5,-1,-1],
                                                                  [1.5,1,1]))
mvfq = create_spatial_mvf(ecq, samplevf)

function allencahn3d(x::Vector{Float64})
    lambda = 3.0 * pi
    c      = lambda / pi
    x1, x2, x3 = x
    y1 = (lambda-1)*x1 - 1.5*c * (x1*x1*x1-x1*x1*x3+x2*x2*x3+2*x1*(x2*x2+x3*x3))
    y2 = (lambda-4)*x2 - 1.5*c * x2 * (2*x1*x1+x2*x2+2*x1*x3+2*x3*x3)
    y3 = (lambda-9)*x3 + 0.5*c * (x1*(x1*x1-3*x2*x2)-3*x3*(2*x1*x1+2*x2*x2+x3*x3))
    return [y1, y2, y3]
end

bmax = [1.8, 1.5, 1.0]
lca, ca = create_cubical_box(25, 25, 25)
eca  = lefschetz_to_euclidean(lca, convert_spatial_coordinates(ca, -bmax, bmax))
mvfa = create_spatial_mvf(eca, allencahn3d)
msa  = morse_sets(eca, mvfa)
println("Allen-Cahn: ", [conley_index(eca, m) for m in msa])
println("Sample box: ", [conley_index(ecp, m) for m in msp], " ", length(mvfp), " ", length(mvlargest))

fig = Figure(size=(1500,500))
ax1 = panel(fig, (1,1), "Cubical box")
plot_cubical3d!(ax1, ecq; pdim=[true,true,true,true])
ax2 = panel(fig, (1,2), "Multivector field, clipped")
plot_cubical3d_mvf!(ax2, ecq, mvfq; clip=([0,0,1], 0.0))
ax3 = panel(fig, (1,3), "Allen-Cahn Morse sets")
plot_cubical3d_morse!(ax3, eca, msa; ci=true)
savefig3d("makie3d_cubical.png", fig)

##############################################################################
#
# Plotting: surface complexes and cell subsets
#
##############################################################################

Vt, Tt = trisurface_torus([0,0,0], 2.0, 0.7, 36, 12)
ecs = tetmesh_to_euclidean(Vt, Tt)
meridian = [k for k in 1:ecs.ncells if ecs.dimensions[k] == 1 &&
            all(abs(c[2]) < 1e-12 && c[1] > 0 for c in ecs.coords[k])]
band = [k for k in 1:ecs.ncells if ecs.dimensions[k] == 2 &&
        all(c[1] < -0.5 && abs(c[2]) < 1.2 for c in ecs.coords[k])]

fig = Figure(size=(1200,500))
ax1 = panel(fig, (1,1), "Torus surface")
plot_simplicial3d!(ax1, ecs; alpha=0.5)
ax2 = panel(fig, (1,2), "Added cell subsets")
plot_simplicial3d!(ax2, ecs; alpha=0.15)
plot_cells3d!(ax2, ecs, meridian; color=:red)
plot_cells3d!(ax2, ecs, band; color=:darkorange, alpha=0.9)
plot_cells3d!(ax2, ecs, ecs.labels[1:3]; color=:black)
savefig3d("makie3d_cells.png", fig)
