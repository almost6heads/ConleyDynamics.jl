export lefschetz_is_locally_closed

"""
    lefschetz_is_locally_closed(lc::AbstractComplex, subcomp::Vector{Int})

Determine whether a Lefschetz complex subset is locally closed.
"""
function lefschetz_is_locally_closed(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Determine whether a Lefschetz complex subset is locally closed
    #

    # Deal with the trivial case first

    if length(subcomp) == 0
        return true
    end

    # Determine the closure-mouth-pair

    clo, mo = lefschetz_clomo_pair(lc, subcomp)

    # Do we have a locally closed set?

    return lefschetz_is_closed(lc, mo)
end

"""
    lefschetz_is_locally_closed(lc::AbstractComplex, subcomp::Vector{String})

Determine whether a Lefschetz complex subset is locally closed.
"""
function lefschetz_is_locally_closed(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Determine whether a Lefschetz complex subset is locally closed
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    return lefschetz_is_locally_closed(lc, subcompI)
end

