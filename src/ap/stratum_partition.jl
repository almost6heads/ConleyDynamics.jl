export stratum_partition, stratum_jump_weight

"""
    stratum_partition(lc::AbstractComplex, ap::Vector{Vector{Vector{Int}}})

Group the enumerated acyclic partitions `ap = construct_ap_space(lc)` by
their Morse vector. Returns `Dict{Vector{Int},Vector{Int}}` mapping `M` to
the indices into `ap` with `morse_vector(lc, ap[i]) == M`.
"""
function stratum_partition(lc::AbstractComplex, ap::Vector{Vector{Vector{Int}}})
    strata = Dict{Vector{Int}, Vector{Int}}()
    for i in eachindex(ap)
        M = morse_vector(lc, ap[i])
        push!(get!(() -> Int[], strata, M), i)
    end
    return strata
end

"""
    stratum_jump_weight(delta::Vector{Int})

The weight `|c|` of the canonical decomposition
`delta = sum_k c_k (e_k+e_{k-1})`.
"""
function stratum_jump_weight(delta::Vector{Int})
    return sum(_rank_decomposition(delta))
end
