export lefschetz_is_closed

"""
    lefschetz_is_closed(lc::LefschetzComplex, subcomp::Vector{Int})

Determine whether a Lefschetz complex subset is closed.
"""
function lefschetz_is_closed(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Determine whether a Lefschetz complex subset is closed
    #

    # Deal with the trivial case first

    if length(subcomp) == 0
        return true
    end

    # Compute the closure

    subclosure = lefschetz_closure(lc, subcomp)

    # Determine whether subcomp is closed

    k = 1
    isclosed = true
    while isclosed && (k <= length(subclosure))
        if subclosure[k] in subcomp
            k += 1
        else
            isclosed = false
        end
    end

    # Return the result

    return isclosed
end

"""
    lefschetz_is_closed(lc::LefschetzComplex, subcomp::Vector{String})

Determine whether a Lefschetz complex subset is closed.
"""
function lefschetz_is_closed(lc::LefschetzComplex, subcomp::Vector{String})
    #
    # Determine whether a Lefschetz complex subset is closed
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    return lefschetz_is_closed(lc, subcompI)
end

