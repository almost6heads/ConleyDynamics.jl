export lefschetz_topboundary

"""
    lefschetz_topboundary(lc::AbstractComplex, subcomp::Vector{Int})

Compute the topological boundary of a Lefschetz complex subset.

In contrast to the algebraic boundary defined via the boundary operator,
this function computes the boundary of the Lefschetz complex subset
specified in `subcomp` if the Lefschetz complex is interpreted as a
finite topological space. In other words, the topological boundary
is the set difference of the closure and the interior of the subset.
"""
function lefschetz_topboundary(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Compute the topological boundary of a Lefschetz complex subset
    #

    # Compute the closure and the interior

    lcclosure  = lefschetz_closure(lc, subcomp)
    lcinterior = lefschetz_interior(lc, subcomp)

    # Compute and return the topological boundary

    return setdiff(lcclosure, lcinterior)
end

"""
    lefschetz_topboundary(lc::AbstractComplex, subcomp::Vector{String})

Compute the topological boundary of a Lefschetz complex subset.

In contrast to the algebraic boundary defined via the boundary operator,
this function computes the boundary of the Lefschetz complex subset
specified in `subcomp` if the Lefschetz complex is interpreted as a
finite topological space. In other words, the topological boundary
is the set difference of the closure and the interior of the subset.
"""
function lefschetz_topboundary(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Compute the topological boundary of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lctopbndI = lefschetz_topboundary(lc, subcompI)
    lctopbnd  = lc.labels[lctopbndI]
    return lctopbnd
end

