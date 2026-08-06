export plot_morse_reachability

"""
    plot_morse_reachability(lc, ap, fname; A=nothing, title="", figw=nothing, figh=nothing, pv=false)

Hasse-style diagram of *eventual* reachability between the Morse-vector
strata of AP(X), i.e. the multi-step generalization of `plot_morse_strata`:
an edge `M2 -> M1` here means some element of `AP_{M2}` reaches `AP_{M1}`
via a chain of atomic refinements of any length (`stratum_reachability`),
not just a single one (`stratum_adjacency`). Nodes are laid out exactly as
in `plot_morse_strata` -- distinct Morse vectors, drawn as circles labeled
`(M)` and `n=...`, placed by level sum(M).

Because eventual reachability is transitive, `stratum_reachability` reports
an edge for *every* comparable pair of strata reachable by some chain --
far denser than the single-step covering relation `plot_morse_strata`
draws. To keep the picture legible, only the *transitive reduction* of
that relation is drawn (`_stratum_transitive_reduction`): an edge `M2 ->
M1` is omitted whenever it is already implied by some other drawn path
`M2 -> ... -> M1`. This means the *absence* of a drawn edge between two
strata does not mean one is unreachable from the other -- only that a
shown path already establishes it. The full (unreduced) edge set and
underlying element-level reachability are both available in the return
value for anything needing the raw data.

Edges are labeled `|c|=<weight>, <covered>/<total>` (how many of the
coarser stratum's partitions *eventually* reach the finer one, out of how
many total), and colored:
  * solid green:  uniform -- every partition of the coarser stratum
                  eventually reaches the finer one
  * dashed red:   not uniform -- at least one partition of the coarser
                  stratum never reaches the finer one
Unlike `plot_morse_strata`, there is no weight-based three-way split:
Theorem 5.4's uniformity guarantee is specifically about single atomic
refinements and does not extend to multi-step chains, so no analogous
"guaranteed" tier exists here to draw a distinction around.

`figw`, `figh`, `title`, and `pv` behave exactly as in `plot_morse_strata`.

# Examples

Consider the triangle `ABC` with a pendant edge `CD` attached at `C`:

```julia
labels    = ["A", "B", "C", "D"]
simplices = [[1, 2, 3], [3, 4]]
lc = create_simplicial_complex(labels, simplices)

ap = construct_ap_space(lc)
plot_morse_reachability(lc, ap, "reachability.pdf",
                        title="Triangle ABC with pendant edge CD")
```

If both `plot_morse_reachability` and `plot_morse_strata` are wanted for the
same complex, precompute `atomic_distances` once and pass it to both, since
it is the expensive step and is otherwise recomputed by each call:

```julia
A = atomic_distances(lc, ap)
plot_morse_reachability(lc, ap, "reachability.pdf"; A=A)
plot_morse_strata(lc, ap, "adjacency.pdf"; A=A)
```
"""
function plot_morse_reachability(lc::AbstractComplex,
                                 ap::Vector{Vector{Vector{Int}}},
                                 fname::String;
                                 A=nothing,
                                 title::String="",
                                 figw::Union{Int,Nothing}=nothing,
                                 figh::Union{Int,Nothing}=nothing,
                                 pv::Bool=false)
    strata = stratum_partition(lc, ap)
    reach, edges = stratum_reachability(lc, ap; A=A)
    edges_drawn  = _stratum_transitive_reduction(edges)

    edge_style(e) = e.uniform ? ("seagreen", "solid", 2.2) : ("firebrick", "shortdashed", 2.0)

    legend_header = [
        "Node: (Morse vector M)  --  n = number of AP(X) partitions with that M",
        "Edge label: |c| = jump weight of M1-M2  --  covered/total = how many of the coarser",
        "stratum's partitions eventually reach the finer one via some chain of atomic refinements",
    ]

    legend_swatches = [
        (color="seagreen", dash="solid", width=2.2,
         caption=["uniform -- every partition of the coarser stratum eventually reaches the finer one"]),
        (color="firebrick", dash="shortdashed", width=2.0,
         caption=["not uniform -- at least one partition of the coarser stratum never reaches it"]),
    ]

    pos, figwI, fighI = _plot_strata_diagram(strata, edges_drawn, fname;
                                             title=title, figw=figw, figh=figh, pv=pv,
                                             legendh=115.0,
                                             edge_style=edge_style,
                                             legend_header=legend_header,
                                             legend_swatches=legend_swatches)

    return (pos=pos, edges=edges, edges_drawn=edges_drawn, reach=reach,
            strata=strata, figw=figwI, figh=fighI)
end
