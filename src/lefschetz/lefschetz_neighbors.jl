export lefschetz_neighbors

"""
    lefschetz_neighbors(lc::AbstractComplex, subcomp::Vector{Int})

Find all cells neighboring a Lefschetz complex subset.

The function returns all cells which are not contained in `subcomp`,
but which are either in the closure or the open hull of `subcomp`.
In other words, if one adds a neighbor to the set `subcomp`, then
the number of connected components of the combined set does not
increase. It could, however, decrease.
"""
function lefschetz_neighbors(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Find the neighbors of a Lefschetz complex subset
    #

    # Handle the trivial case first

    if length(subcomp) == 0
        return Vector{Int}()
    end

    # Collect all neighbors and return the result

    subcl = lefschetz_closure(lc, subcomp)
    suboh = lefschetz_openhull(lc, subcomp)

    return sort(setdiff(union(subcl, suboh), subcomp))
end

"""
    lefschetz_neighbors(lc::AbstractComplex, subcomp::Vector{String})

Find all cells neighboring a Lefschetz complex subset.

The function returns all cells which are not contained in `subcomp`,
but which are either in the closure or the open hull of `subcomp`.
In other words, if one adds a neighbor to the set `subcomp`, then
the number of connected components of the combined set does not
increase. It could, however, decrease.
"""
function lefschetz_neighbors(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Find the neighbors of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcnI = lefschetz_neighbors(lc, subcompI)
    lcn  = lc.labels[lcnI]
    return lcn
end

