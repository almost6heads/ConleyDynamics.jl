export manifold_boundary

"""
    manifold_boundary(lc::AbstractComplex)

Extract the manifold boundary from a Lefschetz complex.

The function expects a Lefschetz complex which represents a
compact d-dimensional manifold with boundary. It returns a 
list of all cells which lie on the topological boundary of
the manifold, in the form of a `Vector{Int}`.
"""
function manifold_boundary(lc::AbstractComplex)
    #
    # Extract the manifold boundary from a Lefschetz complex
    #

    # Loop through all the cells and find the ones on
    # the boundary of the Lefschetz complex

    boundarylist = Vector{Int}()
    lowdimcells = findall(x -> x<lc.dim, lc.dimensions)

    for cell in lowdimcells
        cellcobnd = lefschetz_coboundary(lc, cell)
        celldim   = lc.dimensions[cell]
        if (celldim < lc.dim-1) & (length(cellcobnd) == 0)
            push!(boundarylist, cell)
        elseif (celldim == lc.dim-1) & (length(cellcobnd) <= 1)
            push!(boundarylist, cell)
        end
    end

    # Determine the closure of the boundary facets

    bndclosure = lefschetz_closure(lc, boundarylist)

    # Return the manifold boundary

    return bndclosure
end

