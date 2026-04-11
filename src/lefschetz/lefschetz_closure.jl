export lefschetz_closure

"""
    lefschetz_closure(lc::AbstractComplex, subcomp::Vector{Int})

Compute the closure of a Lefschetz complex subset.
"""
function lefschetz_closure(lc::AbstractComplex, subcomp::Vector{Int})
    #
    # Compute the closure of a Lefschetz complex subset
    #

    # Extract the boundary matrix
        
    bnd = lc.boundary

    # Create list of cells to consider

    lcsubset = deepcopy(subcomp)
    lcclosure = Vector{Int}()

    while length(lcsubset)>0
        cell = popfirst!(lcsubset)
        push!(lcclosure, cell)
        if length(bnd.columns[cell])>0
            append!(lcsubset, bnd.columns[cell])
        end
    end

    # Remove duplicates, sort, and return

    return sort(unique(lcclosure))
end

"""
    lefschetz_closure(lc::AbstractComplex, subcomp::Vector{String})

Compute the closure of a Lefschetz complex subset.
"""
function lefschetz_closure(lc::AbstractComplex, subcomp::Vector{String})
    #
    # Compute the closure of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcclosureI = lefschetz_closure(lc, subcompI)
    lcclosure  = lc.labels[lcclosureI]
    return lcclosure
end

