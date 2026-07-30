export convert_partition_mvf, convert_mvf_partition

"""
    convert_partition_mvf(sp::Set{Set{Int}})

Convert a set partition into a multivector field.

`sp` is a full partition of a complex's cells into blocks, including
singleton blocks (the representation used internally by `construct_ap_space`
and `construct_ap_stratum`). Returns the corresponding `CellSubsets`-style
multivector field: singleton blocks are dropped (they are implicit critical
cells) and every remaining block is returned as a sorted `Vector{Int}`.
"""
function convert_partition_mvf(sp::Set{Set{Int}})
    #
    # Convert a set partition into a multivector field
    #
    spvec = collect.(collect(sp))
    mvf = Vector{Vector{Int}}()
    for mv in spvec
        if length(mv)>1
            push!(mvf, sort(mv))
        end
    end

    return mvf
end

"""
    convert_mvf_partition(lc::AbstractComplex, mvf::CellSubsets)

Convert a multivector field into a set partition.

The inverse of `convert_partition_mvf`: adds the implicit singleton blocks
(every cell of `lc` not already covered by `mvf`) and returns the result as
a `Set{Set{Int}}`, the full-partition representation used internally by
`construct_ap_space`/`construct_ap_stratum`.
"""
function convert_mvf_partition(lc::AbstractComplex,
                               mvf::CellSubsets)
    #
    # Convert a multivector field into a partition
    #

    if typeof(mvf) == Vector{Vector{String}}
        mvfI = convert_cellsubsets(lc, mvf)
    else
        mvfI = deepcopy(mvf)
    end

    # Add the singletons

    singletons = setdiff(collect(1:lc.ncells), vcat(mvfI...))
    for cc in singletons
        push!(mvfI, [cc,])
    end

    # Create and return the partition

    return Set(Set.(mvfI))
end
