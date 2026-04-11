export lefschetz_interior

"""
    lefschetz_interior(lc::AbstractComplex, subcomp::Vector{Int})

Compute the interior of a Lefschetz complex subset.
"""
function lefschetz_interior(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Compute the interior of a Lefschetz complex subset
    #

    # Initialize the interior

    lcinterior = Vector{Int}()

    # For every cell, check whether its open hull is contained in subcomp

    for kc in subcomp
        kcopenhull = lefschetz_openhull(lc, [kc])
        if issubset(kcopenhull, subcomp)
            push!(lcinterior, kc)
        end
    end

    # Return the interior

    return lcinterior
end

"""
    lefschetz_interior(lc::AbstractComplex, subcomp::Vector{String})

Compute the interior of a Lefschetz complex subset.
"""
function lefschetz_interior(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Compute the interior of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcinteriorI = lefschetz_interior(lc, subcompI)
    lcinterior  = lc.labels[lcinteriorI]
    return lcinterior
end

