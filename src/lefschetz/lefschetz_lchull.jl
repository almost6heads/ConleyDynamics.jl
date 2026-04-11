export lefschetz_lchull

"""
    lefschetz_lchull(lc::AbstractComplex, subcomp::Vector{Int})

Compute the locally closed hull of a Lefschetz complex subset.

The locally closed hull is the smallest locally closed set which
contains the given cells. It is the intersection of the closure
and the open hull.
"""
function lefschetz_lchull(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Compute the locally closed hull of a Lefschetz complex subset
    #

    # Determine the open hull and the closure

    openhull = lefschetz_openhull(lc, subcomp)
    closure  = lefschetz_closure(lc, subcomp)

    # Find the intersection and return

    lchull = sort(intersect(openhull, closure))
    return lchull
end

"""
    lefschetz_lchull(lc::AbstractComplex, subcomp::Vector{String})

Compute the locally closed hull of a Lefschetz complex subset.

The locally closed hull is the smallest locally closed set which
contains the given cells. It is the intersection of the closure
and the open hull.
"""
function lefschetz_lchull(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Compute the locally closed hull of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lchullI = lefschetz_lchull(lc, subcompI)
    lchull  = lc.labels[lchullI]
    return lchull
end

