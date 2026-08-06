export plot_morse_strata

"""
    plot_morse_strata(lc, ap, fname; A=nothing, title="", figw=nothing, figh=nothing, pv=false)

Hasse-style diagram of the Morse-vector strata of AP^d(X). Nodes are
distinct Morse vectors, drawn as circles labeled `(M)` and `n=...` (the
number of AP^d(X) partitions with that Morse vector), placed by *level*
sum(M) (monotone under refinement) -- note that "adjacent"
levels need not differ by 1 in sum(M): every weight-1 jump already
raises it by 2, since e_k + e_{k-1} sums to 2.

Edges are drawn between every pair of strata found adjacent by
`stratum_adjacency`, labeled `|c|=<weight>, <covered>/<total>` (how
many of the coarser stratum's partitions admit an atomic refinement
into the finer one, out of how many total), and styled by the
classification spelled out in the figure's own legend:
  * solid green:  weight |c| = 1 (always uniform)
  * solid gray:   weight |c| >= 2, but every partition happens to be
                  covered for this particular complex anyway
  * dashed red:   weight |c| >= 2 and *not* every partition is covered
                  (a realized counterexample to the naive question)
An edge that skips an intermediate level is bent aside so it doesn't
draw on top of the shorter edges passing through that level.

`figw` and `figh` default to `nothing`, which auto-sizes the canvas: the
full layout (nodes and every bend point) is computed first, then the
canvas is sized and the whole picture shifted to fit it exactly (with
a floor wide enough for the legend text even when there's very little
else to draw), so nothing is ever cut off regardless of how far a bend
ends up needing to reach. Pass explicit values to override -- explicit
`figw` and `figh` are used as-is, not auto-corrected. `title`, if
given, is drawn at the top of the figure.

See also `plot_morse_reachability` for the analogous diagram of *eventual*
(multi-step) reachability between strata, rather than single-step
adjacency.

# Examples

Consider the triangle `ABC` with a pendant edge `CD` attached at `C`:

```julia
labels    = ["A", "B", "C", "D"]
simplices = [[1, 2, 3], [3, 4]]
lc = create_simplicial_complex(labels, simplices)

ap = construct_ap_space(lc)
plot_morse_strata(lc, ap, "strata.pdf", title="Triangle ABC with pendant edge CD")
```

If both `plot_morse_strata` and `plot_morse_reachability` are wanted for the
same complex, precompute `atomic_distances` once and pass it to both, since
it is the expensive step and is otherwise recomputed by each call:

```julia
A = atomic_distances(lc, ap)
plot_morse_strata(lc, ap, "adjacency.pdf"; A=A)
plot_morse_reachability(lc, ap, "reachability.pdf"; A=A)
```
"""
function plot_morse_strata(lc::AbstractComplex,
                           ap::Vector{Vector{Vector{Int}}},
                           fname::String;
                           A=nothing,
                           title::String="",
                           figw::Union{Int,Nothing}=nothing,
                           figh::Union{Int,Nothing}=nothing,
                           pv::Bool=false)
    strata = stratum_partition(lc, ap)
    edges  = stratum_adjacency(lc, ap; A=A)

    edge_style(e) = e.uniform ? (e.weight <= 1 ? ("seagreen", "solid", 2.6) :
                                                  ("gray50",   "solid", 1.8)) :
                                 ("firebrick", "shortdashed", 2.0)

    legend_header = [
        "Node: (Morse vector M)  --  n = number of AP^d(X) partitions with that M",
        "Edge label: |c| = jump weight of M1-M2  --  covered/total = how many of the",
        "coarser stratum's partitions admit an atomic refinement into the finer one",
    ]

    legend_swatches = [
        (color="seagreen", dash="solid", width=2.6,
         caption=["|c|=1 and uniform (guaranteed for the full, unrestricted AP^d(X))"]),
        (color="gray50", dash="solid", width=1.8,
         caption=["|c|>=2, but every partition happens to be covered here"]),
        (color="firebrick", dash="shortdashed", width=2.0,
         caption=["NOT every partition is covered (a realized counterexample -- can happen even",
                   "at |c|=1 if this is a restricted sub-poset like AP(X); the guarantee only covers AP^d(X))"]),
    ]

    pos, figwI, fighI = _plot_strata_diagram(strata, edges, fname;
                                             title=title, figw=figw, figh=figh, pv=pv,
                                             legendh=155.0,
                                             edge_style=edge_style,
                                             legend_header=legend_header,
                                             legend_swatches=legend_swatches)

    return (pos=pos, edges=edges, strata=strata, figw=figwI, figh=fighI)
end

"""
    _plot_strata_diagram(strata, edges, fname; title, figw, figh, pv, legendh,
                         edge_style, legend_header, legend_swatches)

Shared Luxor rendering core behind `plot_morse_strata` and
`plot_morse_reachability`: lays out one node per key of `strata` (a
`Dict{Vector{Int},Vector{Int}}` as returned by `stratum_partition`) by
level sum(M), and draws an edge for every element of `edges`
(`NamedTuple`s shaped like the output of `stratum_adjacency` and
`stratum_reachability`, using `.M1`, `.M2`, `.weight`, `.ncovered`,
`.ntotal`).

The legend is built from two pieces:
  * `legend_header`: plain text lines shown at the top of the legend.
  * `legend_swatches`: one colored line-and-caption block per entry,
    given as a vector of named tuples with fields `color`, `dash`,
    `width`, and `caption` (a `Vector{String}`).

`edge_style(e)` maps an edge to its `(color, dash, width)` -- this is
the one piece of classification logic that differs between callers.
`legendh` is the fixed height reserved for the legend block; callers
size it to their own header/swatch content.

Returns `(pos, figwI, fighI)` for the caller to fold into its own return
value alongside whatever edge/strata data it wants to expose.
"""
function _plot_strata_diagram(strata::Dict{Vector{Int},Vector{Int}},
                              edges::Vector{<:NamedTuple},
                              fname::String;
                              title::String,
                              figw::Union{Int,Nothing},
                              figh::Union{Int,Nothing},
                              pv::Bool,
                              legendh::Float64,
                              edge_style::Function,
                              legend_header::Vector{String},
                              legend_swatches::Vector{<:NamedTuple})
    if !(lowercase(fname[end-3:end]) in [".png", ".pdf", ".eps", ".svg"])
        error("The filename must have one of the following extensions: .png, .pdf, .eps, .svg")
    end

    Ms     = collect(keys(strata))
    ranks  = Dict(M => sum(M) for M in Ms)
    levels = sort(unique(values(ranks)))
    nlevels = length(levels)

    # "Adjacent" strata are judged by *level index* (position in the
    # sorted list of distinct sum(M) values that actually occur), not by
    # raw rank difference: every weight-1 jump e_k+e_{k-1} already raises
    # sum(M) by 2, so consecutive levels routinely differ by 2 or more
    # even when no intermediate stratum is skipped.
    lvidx = Dict(r => i for (i, r) in enumerate(levels))

    bylevel = Dict(r => sort([M for M in Ms if ranks[M] == r]) for r in levels)
    maxrow  = maximum(length(v) for v in values(bylevel))

    # Sizing/layout strategy: rather than *guessing* how much extra room
    # a bent skip-edge will need and hoping it's enough (a fixed guess
    # here previously caused edges to be cut off at the canvas boundary
    # -- the guess just happened to be enough for whichever example was
    # tested most recently), lay everything out first in an
    # unconstrained coordinate system, actually compute every bend
    # point's real position (searching for genuine clearance against
    # every other node), and only *then* measure how wide the result
    # truly is and size/shift the canvas to fit it exactly.

    colspacing = 230.0
    margin     = 110.0
    titleh     = title == "" ? 20.0 : 46.0

    provw = max(950.0, 2margin + max(0, maxrow - 1) * colspacing)
    fh = figh === nothing ? max(560.0, titleh + 2margin + max(0, nlevels - 1) * 155.0 + legendh) :
                             Float64(figh)
    fighI = Int(round(fh))
    plotbottom = fighI - legendh

    yof(r) = nlevels == 1 ? titleh + margin + (plotbottom - titleh - 2margin) / 2 :
             titleh + margin + (plotbottom - titleh - 2margin) *
             (nlevels - lvidx[r]) / (nlevels - 1)

    pos = Dict{Vector{Int},Point}()
    for r in levels
        Mr = bylevel[r]
        n  = length(Mr)
        for (i, M) in enumerate(Mr)
            x = n == 1 ? provw / 2 : margin + (provw - 2margin) * (i - 1) / (n - 1)
            pos[M] = Point(x, yof(r))
        end
    end

    # Edges that skip a level (there is no stratum at an intermediate
    # sum(M) actually reached in between) would be drawn exactly on top
    # of the shorter edges between the levels they pass through if drawn
    # straight, so they're bent into a V-shaped path offset to the side
    # instead. A *fixed* offset can just as easily bend the path
    # directly toward some other node that happens to sit near the
    # straight-line route (this did happen: with a single node per
    # intermediate level, that node is often dead-centered exactly where
    # a naive bend lands), so the offset is chosen by actually checking
    # clearance against every other node, starting from the direction
    # the straight line already slants in and growing until both new
    # segments clear every node they pass near.

    function _point_segment_dist(p::Point, a::Point, b::Point)
        abx, aby = b.x - a.x, b.y - a.y
        ablen2 = abx^2 + aby^2
        t = ablen2 == 0 ? 0.0 : clamp(((p.x - a.x) * abx + (p.y - a.y) * aby) / ablen2, 0.0, 1.0)
        qx, qy = a.x + t * abx, a.y + t * aby
        return hypot(p.x - qx, p.y - qy)
    end

    clearance = 46.0

    function _bend_point(p1::Point, p2::Point, obstacles::Vector{Point})
        mx, my = (p1.x + p2.x) / 2, (p1.y + p2.y) / 2
        preferred = p2.x >= p1.x ? 1 : -1
        for mult in 0:8, side in (preferred, -preferred)
            off = side * (70.0 + 45.0 * mult)
            cand = Point(mx + off, my)
            if all(o -> _point_segment_dist(o, p1, cand) > clearance &&
                       _point_segment_dist(o, cand, p2) > clearance, obstacles)
                return cand
            end
        end
        return Point(mx + preferred * 430.0, my)
    end

    bends = Dict{Int,Point}()
    for (i, e) in enumerate(edges)
        gap = lvidx[ranks[e.M1]] - lvidx[ranks[e.M2]]
        if gap > 1
            obstacles = [q for (N, q) in pos if N != e.M1 && N != e.M2]
            bends[i] = _bend_point(pos[e.M2], pos[e.M1], obstacles)
        end
    end

    # Measure the true horizontal extent of every node and every bend
    # point (plus room for the text labels drawn at each), then size
    # the canvas to fit that exactly and shift everything to sit inside
    # it -- unless the caller passed an explicit figw, which is used
    # as-is (an intentional override, not auto-corrected).

    noderad  = 34.0
    labelpad = 95.0

    xs = Float64[]
    for p in values(pos)
        push!(xs, p.x - noderad - 4.0, p.x + noderad + 4.0)
    end
    for p in values(bends)
        push!(xs, p.x - labelpad, p.x + labelpad)
    end
    xmin, xmax = minimum(xs), maximum(xs)

    if figw === nothing
        # Also enforce a minimum canvas width, wide enough for the
        # legend text regardless of how little else there is to draw
        # (a sparse diagram -- few strata, no bends -- would otherwise
        # size the canvas to its few small nodes and clip the legend).
        # When the content is narrower than that minimum it is centered
        # in the extra space rather than left-aligned against it.
        outerpad = 20.0
        minfigw  = 950.0
        contentw = xmax - xmin
        figwF  = max(contentw + 2 * outerpad, minfigw)
        figwI  = Int(round(figwF))
        shiftx = (figwF - contentw) / 2 - xmin
    else
        figwI  = figw
        shiftx = 0.0
    end

    if shiftx != 0.0
        for M in keys(pos)
            pos[M] = Point(pos[M].x + shiftx, pos[M].y)
        end
        for i in keys(bends)
            bends[i] = Point(bends[i].x + shiftx, bends[i].y)
        end
    end

    Drawing(figwI, fighI, fname)
    background("white")

    if title != ""
        sethue("black")
        fontsize(17)
        Luxor.text(title, Point(figwI / 2, 28), halign=:center)
    end

    # Edges first, so nodes are drawn on top of them.

    for (i, e) in enumerate(edges)
        p1, p2 = pos[e.M2], pos[e.M1]
        gap = lvidx[ranks[e.M1]] - lvidx[ranks[e.M2]]

        color, dash, width = edge_style(e)
        sethue(color)
        setdash(dash)
        setline(width)

        if gap <= 1
            line(p1, p2, action = :stroke)
            mid = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)
        else
            mid = bends[i]
            poly([p1, mid, p2], action = :stroke)
        end

        sethue("black")
        fontsize(11)
        Luxor.text("|c|=$(e.weight), $(e.ncovered)/$(e.ntotal)", mid, halign=:center)
    end

    # Nodes

    for (M, p) in pos
        sethue("white")
        circle(p, 34, action = :fill)
        sethue("royalblue4")
        setline(1.6)
        circle(p, 34, action = :stroke)
        sethue("black")
        fontsize(14)
        Luxor.text(string(Tuple(M)), p + Point(0, -5), halign=:center)
        fontsize(11)
        Luxor.text("n=$(length(strata[M]))", p + Point(0, 14), halign=:center)
    end

    # Legend, explaining what the node and edge labels mean, followed by
    # one colored line-and-caption block per entry of legend_swatches.

    lx, ly = 22.0, plotbottom + 14.0
    sethue("black")
    fontsize(12)
    headerlinesp = 17.0
    for (i, ln) in enumerate(legend_header)
        Luxor.text(ln, Point(lx, ly + headerlinesp * (i - 1)), halign=:left)
    end

    swy = ly + headerlinesp * (length(legend_header) - 1) + 22.0
    capgap    = 4.0
    caplinesp = 16.0
    swatchgap = 20.0
    for sw in legend_swatches
        sethue(sw.color); setline(sw.width); setdash(sw.dash)
        line(Point(lx, swy), Point(lx + 40, swy), action = :stroke)
        sethue("black"); fontsize(12)
        for (j, ln) in enumerate(sw.caption)
            Luxor.text(ln, Point(lx + 50, swy + capgap + caplinesp * (j - 1)), halign=:left)
        end
        swy += swatchgap + caplinesp * max(0, length(sw.caption) - 1)
    end

    finish()
    pv && preview()

    return (pos, figwI, fighI)
end
