export lefschetz_openhull

"""
    lefschetz_openhull(lc::AbstractComplex, subcomp::Vector{Int})

Compute the open hull of a Lefschetz complex subset.
"""
function lefschetz_openhull(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Compute the open hull of a Lefschetz complex subset
    #

    # Extract the boundary matrix
        
    bnd = lc.boundary

    # Create list of cells to consider

    lcsubset = deepcopy(subcomp)
    lcopenhull = Vector{Int}()

    while length(lcsubset)>0
        cell = popfirst!(lcsubset)
        push!(lcopenhull, cell)
        if length(bnd.rows[cell])>0
            append!(lcsubset, bnd.rows[cell])
        end
    end

    # Remove duplicates, sort, and return

    return sort(unique(lcopenhull))
end

"""
    lefschetz_openhull(lc::AbstractComplex, subcomp::Vector{String})

Compute the open hull of a Lefschetz complex subset.
"""
function lefschetz_openhull(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Compute the open hull of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcopenhullI = lefschetz_openhull(lc, subcompI)
    lcopenhull  = lc.labels[lcopenhullI]
    return lcopenhull
end

