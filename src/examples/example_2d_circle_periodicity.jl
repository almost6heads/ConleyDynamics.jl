using Graphs
using LinearAlgebra
using Plots
using ConleyDynamics
using DifferentialEquations


################################################################################
pcolor = ["red", "coral", "firebrick", "orange", "yellow",
          "gold", "green", "lime", "olive", "cyan", "teal",
          "azure", "blue", "indigo", "purple", "pink",
          "magenta", "brown"]
pcolor2 = ["red", "coral", "orange", "green", "blue", "purple", "magenta"]
################################################################################

# Create complex
ec = create_simplicial_rectangle(40, 40, euclidean=true)

# Find transition space and MVF
sg, _ = create_planar_dts(ec, circ_vf)
sg = create_mvf(ec, sg, circ_vf, strong_iso=true)
str_mvf = connected_components(sg)

# Create CMVF Digraph and find SCCs
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


################################################################################
# Get largest SCC
################################################################################
max_scc = argmax(length, scc)
sub_dg, sub_dg_vmap = induced_subgraph(dg, max_scc)

sub_bfs = bfs_tree(sub_dg, 1)

sub_lp_idxs = dag_longest_path(sub_bfs)
sub_lp_del = Int(floor(length(sub_lp_idxs)/3))

seed_mvf = sub_lp_idxs[[1, 1+sub_lp_del, 1+2*sub_lp_del]]

# Coarsen SCC
pset = digraph_propagate(sub_dg, seed_mvf, 50)

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
# Get verts in regions adjacent to coarsened blocks
################################################################################
qvert = Vector{Int}()
append!(qvert, lefschetz_skeleton(ec, intersect(ppt_cl[1], ppt_cl[2]), 0))
append!(qvert, lefschetz_skeleton(ec, intersect(ppt_cl[1], ppt_cl[3]), 0))
append!(qvert, lefschetz_skeleton(ec, intersect(ppt_cl[2], ppt_cl[3]), 0))


bnd_xy = ec.coords[qvert]
bnd_uv = [Vector{Float64}() for _ in eachindex(bnd_xy)]
for idx in eachindex(bnd_xy)
    this_vf = circ_vf(bnd_xy[idx][1])
    bnd_uv[idx] = this_vf/norm(this_vf)
end

bnd_x = [x[1][1] for x in bnd_xy]
bnd_y = [x[1][2] for x in bnd_xy]
bnd_u = [x[1] for x in bnd_uv]
bnd_v = [x[2] for x in bnd_uv]

# Create vector field for visualization
#vxr = 0.0:4.0:40.0
#vyr = 0.0:4.0:40.0
#vf_xy = [[xi, yi] for xi in vxr, yi in vyr]

vf_xy = ec.coords[1:2:41*41]
vf_uv = [Vector{Float64}() for _ in eachindex(vf_xy)]
for idx in eachindex(vf_xy)
    this_vf = circ_vf(vf_xy[idx][1])
    vf_uv[idx] = this_vf/norm(this_vf)
end

vf_x = [x[1][1] for x in vf_xy]
vf_y = [x[1][2] for x in vf_xy]
vf_u = [x[1] for x in vf_uv]
vf_v = [x[2] for x in vf_uv]


# solve for periodic orbit
u0 = [20.0, 30.0]
tspan = (0.0, 100)

ode_prob = ODEProblem(circ_vf_ode!, u0, tspan)
ode_sol = solve(ode_prob)



################################################################################
# Plot stuff
################################################################################
# # Initial plot
# plot_planar_simplicial(ec, "simrect_40x40_circ_init.pdf",
#                        pdim=[false,true,false])

# MVF plot
qscale = 0.75
p1 = plot_simplicial_mvf(ec, str_mvf, pdim=[false,false,true], mvfcolor=pcolor2)
quiver!(p1, vf_x, vf_y, quiver=(qscale*vf_u, qscale*vf_v),
        lw=0.3, arrow=(Plots.arrow(:closed,:head,0.3,0.15)), color=:magenta)
savefig(p1, "0simrect_40x40_circ_mvf.pdf")

# SCC plot
qscale = 0.75
p2 = plot_simplicial_mvf(ec, scc_ele, pdim=[false,false,true], mvfcolor=pcolor2)
quiver!(p2, vf_x, vf_y, quiver=(qscale*vf_u, qscale*vf_v),
        lw=0.3, arrow=(Plots.arrow(:closed,:head,0.3,0.15)), color=:magenta)
plot!(p2, ode_sol, vars=(1,2), linecolor=:lime)
savefig(p2, "0simrect_40x40_circ_scc.pdf")

qscale = 0.75
p3 = plot_simplicial_mvf(ec, ppt,
                        pdim=[false,false,true], mvfcolor=pcolor2)
quiver!(p3, bnd_x, bnd_y, quiver=(qscale*bnd_u, qscale*bnd_v),
        lw=0.3, arrow=(Plots.arrow(:closed,:head,0.3,0.15)), color=:magenta)
plot!(p3, ode_sol, vars=(1,2), linecolor=:lime)
savefig(p3, "0simrect_40x40_circ_periodic.pdf")