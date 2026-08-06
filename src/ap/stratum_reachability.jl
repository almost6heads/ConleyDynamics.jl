export stratum_reachability

"""
    stratum_reachability(lc::AbstractComplex, ap::Vector{Vector{Vector{Int}}}; A=nothing)

For every element of `ap`, compute the full set of Morse vectors reachable
by *any* chain of atomic refinements of any length -- not just the single
step `stratum_adjacency` checks. `atomic_distances` already contains every
such one-step edge, including the ones between two elements of the *same*
stratum (`stratum_adjacency` deliberately skips `M1 == M2` pairs when
reporting, but the edges are there); this walks the full poset via those
edges to answer the "can you always eventually get there" question,
regardless of how many intermediate steps -- possibly staying within the
source stratum for a while -- it takes.

Since every atomic refinement increases the total block count
(`mvf_length`) by exactly 1, `AP(X)` (or `AP^d(X)`) is graded by block
count and hence acyclic; reachability is computed with a single dynamic-
programming pass in decreasing block-count order (finest elements, which
have no successors, first), rather than a full pairwise transitive
closure: `reach[v] = {M(v)} \\cup \\bigcup_{v \\to w} reach[w]`.

`A`, if supplied, must be `atomic_distances(lc, ap)`; otherwise computed
internally.

Returns `(reach, edges)`:
- `reach::Vector{Set{Vector{Int}}}`: `reach[i]` is the set of Morse vectors
  reachable from `ap[i]`, including `morse_vector(lc, ap[i])` itself.
- `edges::Vector{NamedTuple}`: one entry per pair of distinct strata
  `M2`, `M1` (`M1 >= M2` componentwise) for which at least one element of
  `AP_{M2}` eventually reaches `AP_{M1}`, in the same field layout as
  `stratum_adjacency`'s output (`M1`, `M2`, `delta`, `weight`, `ncovered`,
  `ntotal`, `uniform`, `uncovered`), except `ncovered`/`uniform` now count
  *eventual* rather than one-step reachability.
"""
function stratum_reachability(lc::AbstractComplex,
                              ap::Vector{Vector{Vector{Int}}};
                              A=nothing)
    #
    # Compute, for every element, the set of Morse vectors reachable via
    # any-length chain of atomic refinements
    #

    nap   = length(ap)
    Amat  = A === nothing ? atomic_distances(lc, ap) : A
    Ms    = map(t -> morse_vector(lc, t), ap)
    aplen = map(t -> mvf_length(lc, t), ap)

    # Process finest elements (largest block count, no successors) first,
    # so that every successor of v has already been resolved by the time
    # v itself is processed.

    order = sortperm(aplen; rev=true)

    reach = Vector{Set{Vector{Int}}}(undef, nap)
    for v in order
        s = Set{Vector{Int}}()
        push!(s, Ms[v])
        for w in sparse_get_nz_row(Amat, v)
            union!(s, reach[w])
        end
        reach[v] = s
    end

    # Collapse to a stratum-level table, mirroring stratum_adjacency

    strata = stratum_partition(lc, ap)
    Mlist  = collect(keys(strata))
    edges  = NamedTuple[]

    for M2 in Mlist, M1 in Mlist
        M1 == M2 && continue
        delta = M1 .- M2
        any(<(0), delta) && continue

        idx2 = strata[M2]
        uncovered = filter(v -> !(M1 in reach[v]), idx2)
        ncovered  = length(idx2) - length(uncovered)

        if ncovered > 0
            push!(edges, (M1=M1, M2=M2, delta=delta,
                          weight=stratum_jump_weight(delta),
                          ncovered=ncovered, ntotal=length(idx2),
                          uniform=isempty(uncovered), uncovered=uncovered))
        end
    end

    return reach, edges
end

"""
    _stratum_transitive_reduction(edges::Vector{<:NamedTuple})

Transitive reduction of the stratum graph described by `edges` (as
returned by `stratum_reachability`, or any `NamedTuple` vector with `.M1`
and `.M2` fields defining a directed arc `M2 -> M1`): drop an arc
`M2 -> M1` whenever some other direct successor `w` of `M2` (`w != M1`)
already reaches `M1` via the remaining arcs, since that makes `M2 -> M1`
redundant for reachability -- it is implied by the two-hop-or-longer path
`M2 -> w -> ... -> M1`. The underlying graph is a DAG (`stratum_reachability`
only ever reports `M1 >= M2` componentwise), so this reduction is unique.

Used by `plot_morse_reachability` to draw only the covering relation of
the reachability order (a proper Hasse diagram) instead of every
comparable pair, which `stratum_reachability`'s raw `edges` -- being
closed under transitivity -- reports for all of them.

Returns the filtered subset of `edges` (same `NamedTuple`s, just fewer).
"""
function _stratum_transitive_reduction(edges::Vector{<:NamedTuple})
    succ = Dict{Vector{Int},Vector{Vector{Int}}}()
    for e in edges
        push!(get!(() -> Vector{Int}[], succ, e.M2), e.M1)
    end

    memo = Dict{Vector{Int},Set{Vector{Int}}}()
    function descendants(u::Vector{Int})
        haskey(memo, u) && return memo[u]
        s = Set{Vector{Int}}()
        memo[u] = s
        for w in get(succ, u, Vector{Int}[])
            push!(s, w)
            union!(s, descendants(w))
        end
        return s
    end

    return filter(edges) do e
        !any(w -> w != e.M1 && e.M1 in descendants(w), get(succ, e.M2, Vector{Int}[]))
    end
end
