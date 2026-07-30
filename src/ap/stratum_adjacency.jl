export stratum_adjacency

"""
    stratum_adjacency(lc::AbstractComplex, ap::Vector{Vector{Vector{Int}}}; A=nothing)

For every pair of distinct strata `M1 > M2` for which at least one element
of `AP_{M2}` admits an atomic refinement into `AP_{M1}`, report the jump
`delta = M1-M2`, its weight `|c|`, how many of `AP_{M2}`'s elements are
covered from `M1`, and whether coverage is uniform.

`A`, if supplied, must be `atomic_distances(lc, ap)`; otherwise it is
computed internally (this is the expensive step, so pass it in if already
available).
"""
function stratum_adjacency(lc::AbstractComplex,
                           ap::Vector{Vector{Vector{Int}}};
                           A=nothing)
    strata = stratum_partition(lc, ap)
    Amat   = A === nothing ? atomic_distances(lc, ap) : A

    edges = NamedTuple[]
    Mlist = collect(keys(strata))

    for M2 in Mlist, M1 in Mlist
        M1 == M2 && continue
        delta = M1 .- M2
        any(<(0), delta) && continue

        idx2, idx1 = strata[M2], strata[M1]
        idx1set = Set(idx1)

        uncovered = Int[]
        for w in idx2
            nz = sparse_get_nz_row(Amat, w)
            if isempty(intersect(nz, idx1set))
                push!(uncovered, w)
            end
        end
        ncovered = length(idx2) - length(uncovered)

        if ncovered > 0
            push!(edges, (M1=M1, M2=M2, delta=delta,
                          weight=stratum_jump_weight(delta),
                          ncovered=ncovered, ntotal=length(idx2),
                          uniform=isempty(uncovered), uncovered=uncovered))
        end
    end

    return edges
end
