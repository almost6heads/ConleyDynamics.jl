using Graphs
using LinearAlgebra
using Plots
using ConleyDynamics


################################################################################
function edge_coords(ec, ei)
    e_x = Vector{Float64}(undef, 3*length(ei))
    e_y = Vector{Float64}(undef, 3*length(ei))

    for idx in eachindex(ei)
        this_coords = ec.coords[ei[idx]]
        e_x[3*(idx-1)+1] = this_coords[1][1]
        e_x[3*(idx-1)+2] = this_coords[2][1]
        e_x[3*(idx-1)+3] = NaN
        e_y[3*(idx-1)+1] = this_coords[1][2]
        e_y[3*(idx-1)+2] = this_coords[2][2]
        e_y[3*(idx-1)+3] = NaN
    end

    return e_x, e_y
end


################################################################################
pcolor = ["red", "coral", "firebrick", "orange", "yellow",
          "gold", "green", "lime", "olive", "cyan", "teal",
          "azure", "blue", "indigo", "purple", "pink",
          "magenta", "brown"]
pcolor2 = ["red", "coral", "orange", "green", "blue", "purple"]


################################################################################
pythonplot()
scale = 0.25

# Create simplicial complex
ec = create_simplicial_rectangle(4, 5, euclidean=true)

# Find non-transverse edges for plotting
nt_e = planar_nontransverse_edges(ec, ouroboros_vf)
nte_x, nte_y = edge_coords(ec, nt_e)

# Find vector field for plotting
vf_xy = ec.coords[1:50]
vf_uv = [Vector{Float64}() for _ in eachindex(vf_xy)]
for idx in eachindex(vf_xy)
    this_vf = ouroboros_vf(vf_xy[idx][1])
    vf_uv[idx] = this_vf/norm(this_vf)
end

vf_x = [x[1][1] for x in vf_xy]
vf_y = [x[1][2] for x in vf_xy]
vf_u = [x[1] for x in vf_uv]
vf_v = [x[2] for x in vf_uv]


# Create DTS and merge once, plot before and after
sg, _ = create_planar_dts(ec, ouroboros_vf)

sg_cc1 = connected_components(sg)
sg2 = graph_merge(ec, sg)
sg_cc2 = connected_components(sg2)

p1 = plot_simplicial_mvf(ec, sg_cc1, pdim=[false,true,true], mvfcolor=pcolor2)
plot!(p1, nte_x, nte_y, linewidth=3, linestyle=:dot, linecolor=:lime)
quiver!(p1, vf_x, vf_y, quiver=(scale*vf_u, scale*vf_v),
        lw=0.5, arrow=(Plots.arrow(:closed,:head,0.5,0.25)), color=:magenta)
savefig(p1, "simrect_4x5_ouro_init_dts.pdf")

p2 = plot_simplicial_mvf(ec, sg_cc2, pdim=[false,true,true], mvfcolor=pcolor2)
plot!(p2, nte_x, nte_y, linewidth=3, linestyle=:dot, linecolor=:lime)
quiver!(p2, vf_x, vf_y, quiver=(scale*vf_u, scale*vf_v),
        lw=0.5, arrow=(Plots.arrow(:closed,:head,0.5,0.25)), color=:magenta)
savefig(p2, "simrect_4x5_ouro_fwd_merge.pdf")


# Create dual space and plot
sg_dual = SimpleGraph(ec.ncells)
sg_cc3 = connected_components(sg2)
for cc in sg_cc3
    this_dual = create_dual_set(ec, cc, ouroboros_vf)
    for m = 1:length(this_dual)-1
        add_edge!(sg_dual, this_dual[m], this_dual[m+1])
    end
end

sg_dual_cc1 = connected_components(sg_dual)

p3 = plot_simplicial_mvf(ec, sg_dual_cc1,
                         pdim=[false,true,true], mvfcolor=pcolor2)
plot!(p3, nte_x, nte_y, linewidth=3, linestyle=:dot, linecolor=:lime)
quiver!(p3, vf_x, vf_y, quiver=(-scale*vf_u, -scale*vf_v),
        lw=0.5, arrow=(Plots.arrow(:closed,:head,0.5,0.25)), color=:magenta)
savefig(p3, "simrect_4x5_ouro_init_dual.pdf")


# Merge in dual space and plot
sg3, _ = dual_graph_merge(ec, sg2, ouroboros_vf)

sg_dual2 = SimpleGraph(ec.ncells)
sg_cc4 = connected_components(sg3)
for cc in sg_cc4
    this_dual = create_dual_set(ec, cc, ouroboros_vf)
    for m = 1:length(this_dual)-1
        add_edge!(sg_dual2, this_dual[m], this_dual[m+1])
    end
end

sg_dual_cc2 = connected_components(sg_dual2)

p4 = plot_simplicial_mvf(ec, sg_dual_cc2,
                         pdim=[false,true,true], mvfcolor=pcolor2)
plot!(p4, nte_x, nte_y, linewidth=3, linestyle=:dot, linecolor=:lime)
quiver!(p4, vf_x, vf_y, quiver=(-scale*vf_u, -scale*vf_v),
        lw=0.5, arrow=(Plots.arrow(:closed,:head,0.5,0.25)), color=:magenta)
savefig(p4, "simrect_4x5_ouro_dual_merge.pdf")


# Plot final MVF
sg_cc3 = connected_components(sg3)
p5 = plot_simplicial_mvf(ec, sg_cc3, pdim=[false,true,true], mvfcolor=pcolor2)
plot!(p5, nte_x, nte_y, linewidth=3, linestyle=:dot, linecolor=:lime)
quiver!(p5, vf_x, vf_y, quiver=(scale*vf_u, scale*vf_v),
        lw=0.5, arrow=(Plots.arrow(:closed,:head,0.5,0.25)), color=:magenta)
savefig(p5, "simrect_4x5_ouro_final.pdf")