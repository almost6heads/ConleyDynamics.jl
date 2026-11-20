#
# Figure docs/src/man/img/ouroboros.png for the Ouroboros example in
# docs/src/man/strongmvf.md. Left: the complete accumulating multivector
# field, with the Ouroboros multivector in gold and the multivector in
# front of its head in red. Right: the finest strongly accumulating
# multivector field, in which the two are merged. The arrows show the
# normalized vector field at the edge midpoints.
#
using Plots
using LinearAlgebra
using ConleyDynamics

ec, vf = example_ouroboros()
lc  = euclidean_to_lefschetz(ec)
fwd = planar_transitions(ec, vf)
bwd = planar_transitions(ec, x -> -vf(x))
mvf = create_mvf_hull(lc, fwd)
duals = mvf_dual(lc, mvf, bwd)
smvf, sdual = create_strong_mvf(lc, fwd, bwd)

iouro = 8                                           # Ouroboros multivector
jouro = findfirst(v -> !isempty(intersect(v, mvf[iouro])), smvf)
ihead = [k for k in eachindex(mvf) if k != iouro &&
         !isempty(intersect(mvf[k], smvf[jouro]))][1]  # merged neighbor

# Colors: highlighted multivectors stand out, all others are muted
muted = ["lightsteelblue", "lightblue", "powderblue", "thistle",
         "lightgray", "paleturquoise", "lavender", "lightcyan"]
function mvcolors(n, nall, special)
    cols = [muted[mod1(k, length(muted))] for k in 1:nall]
    for (k, c) in special
        cols[k] = c
    end
    return cols
end
nall(m) = length(m) + count(k -> all(v -> !(k in v), m), 1:ec.ncells)

col1 = mvcolors(length(mvf), nall(mvf),
                [iouro => "gold", ihead => "orangered"])
col3 = mvcolors(length(smvf), nall(smvf), [jouro => "gold"])

# Normalized vector field at the edge midpoints, drawn as small arrows
function arrows!(p; len=0.2, hl=0.07, hw=0.035, color=:black)
    lx, ly, tx, ty = Float64[], Float64[], Float64[], Float64[]
    for k in 1:ec.ncells
        ec.dimensions[k] == 1 || continue
        m = (ec.coords[k][1] + ec.coords[k][2]) / 2
        d = vf(m); d = d / norm(d); n = [-d[2], d[1]]
        a = m - (len/2) * d; b = m + (len/2) * d; c = b - hl * d
        push!(lx, a[1], c[1], NaN); push!(ly, a[2], c[2], NaN)
        q1 = c + hw * n; q2 = c - hw * n
        push!(tx, b[1], q1[1], q2[1], b[1], NaN)
        push!(ty, b[2], q1[2], q2[2], b[2], NaN)
    end
    plot!(p, lx, ly, seriestype=:path, color=color, lw=1.2, label="")
    plot!(p, tx, ty, seriestype=:shape, fillcolor=color, linewidth=0,
          label="")
end

opts = (pdim=[false,true,false], mvfalpha=0.6, tubefac=0.07)
p1 = plot_simplicial_mvf(ec, mvf; mvfcolor=col1, opts...)
arrows!(p1)
p3 = plot_simplicial_mvf(ec, smvf; mvfcolor=col3, opts...)
arrows!(p3)

pall = plot(p1, p3, layout=(1,2), size=(1100, 620), dpi=250)
savefig(pall, "ouroboros.png")
