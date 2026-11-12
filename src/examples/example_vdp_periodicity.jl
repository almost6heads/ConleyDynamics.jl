using Graphs
using LinearAlgebra
using Plots
using ConleyDynamics
using DistMesh
using Random
using DifferentialEquations


Random.seed!(233)

################################################################################
# Helper function
################################################################################
function build_mesh()
        
    function chk(v)
        s = vf_curv(v, vdp_vf_trev, J_vdp_vf_trev)
        q = min(s, 4.5)
        r = (1+q)^-1

        return r/10
    end

    fd(p) = drectangle(p, -7, 7, -7, 7)
    bbox  = ((-7,-7), (7,7))
    pfix = ((-7,-7), (-7,7), (7,-7), (7,7))

    hmin  = 0.05
    fh(p) = chk(p)

    msh = distmesh2d(fd, fh, hmin, bbox, pfix, maxiter=1000)
    pts, tri = msh

    return pts, tri
end


################################################################################
pcolor = ["red", "coral", "firebrick", "orange", "yellow",
          "gold", "green", "lime", "olive", "cyan", "teal",
          "azure", "blue", "indigo", "purple", "pink",
          "magenta", "brown"]
pcolor2 = ["red", "coral", "orange", "green", "blue", "purple"]


################################################################################
# Build mesh using DistMesh
################################################################################
pts, tri = build_mesh()
println("Pts: ", length(pts), " Tri: ", length(tri))

# Translate coords, tris from DistMesh StaticArray
pts_a = [Vector{Float64}() for _ in pts]
for idx in eachindex(pts_a)
    pts_a[idx] = Array(pts[idx])
end

#dt_tri = triangulate(pts_a)
#ec = convert_delaunay_to_simplicial(dt_tri)
tri_a = [Vector{Int}() for _ in tri]
for idx in eachindex(tri_a)
    tri_a[idx] = Array(tri[idx])
end

# Build labels
pad_len = Int(ceil(log10(length(pts_a))))
labels = string.(1:length(pts_a), pad=pad_len)

################################################################################
# Perturbation
################################################################################
p_x = getindex.(pts_a, 1)
p_y = getindex.(pts_a, 2)

A_x_idxs = findall(x -> 4.8 < x < 4.85, p_x)
A_y_idxs = findall(y -> -0.2 < y < -0.1, p_y)
A_idx = intersect(A_x_idxs, A_y_idxs)

B_x_idxs = findall(x -> 4.8 < x < 4.9, p_x)
B_y_idxs = findall(y -> -0.325 < y < -0.25, p_y)
B_idx = intersect(B_x_idxs, B_y_idxs)

pts_a[A_idx][1][2] = -0.2
pts_a[B_idx][1][1] = 4.8


################################################################################
# Create complex
################################################################################
ec = create_simplicial_complex(labels, tri_a, pts_a)

# Find transition space and MVF
sg, _ = create_planar_dts(ec, vdp_vf)
sg = create_mvf(ec, sg, vdp_vf, strong_iso=true)
str_mvf = connected_components(sg)

println("Merge Ratio: ", length(str_mvf)/length(tri_a))


################################################################################
# Create CMVF Digraph and find SCCs
################################################################################
println("Running SCC")
dg = SimpleDiGraph(length(str_mvf))
for idx in eachindex(str_mvf)
    this_cl_mv = lefschetz_closure(ec, str_mvf[idx])
    this_mo = setdiff(this_cl_mv, str_mvf[idx])

    for ele in this_mo
        pi_v = findall(v -> ele in v, str_mvf)

        for mv in pi_v
            add_edge!(dg, idx, mv)
        end
    end
end

# Map SCC (indexed to MVF) back to elements (indexed to EuclideanComplex)
println("Map SCCs to MVs")
scc = strongly_connected_components(dg)
sort!(scc, by=length, rev=true)

scc_ele = [Vector{Int}() for _ in eachindex(scc)]
for idx in eachindex(scc)
    for mv_idx in scc[idx]
        append!(scc_ele[idx], str_mvf[mv_idx])
    end
    unique!(scc_ele[idx])
end
sort!(scc_ele, by=length, rev=true)

red_scc_cl, red_scc_mo = lefschetz_clomo_pair(ec, scc_ele[1])
red_scc_homology = relative_homology(ec, red_scc_cl, red_scc_mo)
println(red_scc_homology)

ylw_scc_cl, ylw_scc_mo = lefschetz_clomo_pair(ec, scc_ele[2])
ylw_scc_homology = relative_homology(ec, ylw_scc_cl, ylw_scc_mo)
println(ylw_scc_homology)


################################################################################
# Get largest SCC
################################################################################
sub_dg, sub_dg_vmap = induced_subgraph(dg, scc[1])
sub_bfs = bfs_tree(sub_dg, 100)
sub_lp_idxs = dag_longest_path(sub_bfs)

# Find coarsened blocks by propagating in the BFS
seed_mvf = sub_lp_idxs[[5, 25, 45]]
pset = digraph_propagate(sub_dg, seed_mvf, 150)

# Get MVF indices from BFS subgraph indices
ppt = [Vector{Int}() for _ in eachindex(pset)]
for idx in eachindex(pset)
    this_ele = str_mvf[sub_dg_vmap[pset[idx]]]
    ppt[idx] = reduce(vcat, this_ele)
end

ppt_cl = [Vector{Int}() for _ in eachindex(ppt)]
for idx in eachindex(pset)
    ppt_cl[idx] = lefschetz_closure(ec, ppt[idx])
end


################################################################################
# Check for 'backward' flow between coarsened blocks
################################################################################
vio13 = intersect(ppt_cl[1], ppt[3])
vio21 = intersect(ppt_cl[2], ppt[1])
vio32 = intersect(ppt_cl[3], ppt[2])

println("1->3 ", vio13)
println("2->1 ", vio21)
println("3->2 ", vio32)


################################################################################
# Get verts in regions adjacent to coarsened blocks
################################################################################
qvert = Vector{Int}()
append!(qvert, lefschetz_skeleton(ec, intersect(ppt_cl[1], ppt_cl[2]), 0))
append!(qvert, lefschetz_skeleton(ec, intersect(ppt_cl[1], ppt_cl[3]), 0))
append!(qvert, lefschetz_skeleton(ec, intersect(ppt_cl[2], ppt_cl[3]), 0))


bnd_xy = ec.coords[qvert]
bnd_uv = [Vector{Float64}() for _ in eachindex(bnd_xy)]
for idx in eachindex(bnd_xy)
    this_vf = vdp_vf(bnd_xy[idx][1])
    bnd_uv[idx] = this_vf/norm(this_vf)
end

bnd_x = [x[1][1] for x in bnd_xy]
bnd_y = [x[1][2] for x in bnd_xy]
bnd_u = [x[1] for x in bnd_uv]
bnd_v = [x[2] for x in bnd_uv]

# Get vec field flow for visualization
vxr = -7.0:0.5:7.0
vyr = -7.0:0.5:7.0
vf_xy = [[xi, yi] for xi in vxr, yi in vyr]

vf_uv = [Vector{Float64}() for _ in eachindex(vf_xy)]
for idx in eachindex(vf_xy)
    this_vf = vdp_vf(vf_xy[idx])
    vf_uv[idx] = this_vf/norm(this_vf)
end

vf_x = [x[1] for x in vf_xy]
vf_y = [x[2] for x in vf_xy]
vf_u = [x[1] for x in vf_uv]
vf_v = [x[2] for x in vf_uv]


# solve for periodic orbit
u0 = [0.0, 2.0]
tspan = (0.0, 100)

ode_prob = ODEProblem(vdp_vf_ode!, u0, tspan)
ode_sol = solve(ode_prob, Tsit5())
u1 = ode_sol(100)
period_ode = ODEProblem(vdp_vf_ode!, u1, tspan)
new_sol = solve(period_ode)

################################################################################
# Plot stuff
################################################################################
pythonplot()

# MVF plot
scale=0.3
p1 = plot_simplicial_mvf(ec, str_mvf, pdim=[false,false,true], mvfcolor=pcolor2)
quiver!(p1,
        vf_x[1:2:end], vf_y[1:2:end],
        quiver=(scale*vf_u[1:2:end], scale*vf_v[1:2:end]),
        lw=0.4, arrow=(Plots.arrow(:closed,:head,0.4,0.2)), color=:magenta)
plot!(p1, new_sol, vars=(1,2), linecolor=:lime)
savefig(p1, "distmesh_vdp_curv_mvf_7x7_233.pdf")

# SCC plot
scale=0.3
p2 = plot_simplicial_mvf(ec, scc_ele, pdim=[false,false,true], mvfcolor=pcolor2)
quiver!(p2,
        vf_x[1:2:end], vf_y[1:2:end],
        quiver=(scale*vf_u[1:2:end], scale*vf_v[1:2:end]),
        lw=0.4, arrow=(Plots.arrow(:closed,:head,0.4,0.2)), color=:magenta)
plot!(p2, new_sol, vars=(1,2), linecolor=:lime)
savefig(p2, "distmesh_vdp_curv_scc_7x7_233.pdf")

# Coarsen Plot
scale = 0.3
p3 = plot_simplicial_mvf(ec, ppt, pdim=[false,false,true], mvfcolor=pcolor2)
quiver!(p3,
        bnd_x, bnd_y,
        quiver=(scale*bnd_u, scale*bnd_v),
        lw=0.4, arrow=(Plots.arrow(:closed,:head,0.4,0.2)), color=:magenta)
plot!(p3, new_sol, vars=(1,2), linecolor=:lime)
savefig(p3, "distmesh_vdp_curv_periodic_7x7_233.pdf")