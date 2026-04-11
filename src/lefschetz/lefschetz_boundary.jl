export lefschetz_boundary, lefschetz_coboundary

"""
    lefschetz_boundary(lc::AbstractComplex, cellI::Int)

Compute the support of the boundary of a Lefschetz complex cell.

This method returns the boundary support as a `Vector{Int}`.
"""
function lefschetz_boundary(lc::AbstractComplex, cellI::Int)
    #
    # Compute the support of the boundary of a Lefschetz complex cell
    #

    # Extract the boundary matrix
        
    bnd = lc.boundary

    # Create the vector of boundary indices

    if length(bnd.columns[cellI]) == 0
        bndvecI = Vector{Int}([])
    else
        bndvecI = deepcopy(bnd.columns[cellI])
    end

    # Return the result

    return bndvecI
end

"""
    lefschetz_boundary(lc::AbstractComplex, cellS::String)

Compute the support of the boundary of a Lefschetz complex cell.

This method returns the boundary support as a `Vector{String}`.
"""
function lefschetz_boundary(lc::AbstractComplex, cellS::String)
    #
    # Compute the support of the boundary of a Lefschetz complex cell
    #

    # Extract the boundary matrix
        
    bnd = lc.boundary

    # Convert to integer format, get the boundary, convert back

    cellI = lc.indices[cellS]
    bndvecI = lefschetz_boundary(lc, cellI)
    bndvecS = lc.labels[bndvecI]

    # Return the result

    return bndvecS
end

"""
    lefschetz_coboundary(lc::AbstractComplex, cellI::Int)

Compute the support of the coboundary of a Lefschetz complex cell.

This method returns the boundary support as a `Vector{Int}`.
"""
function lefschetz_coboundary(lc::AbstractComplex, cellI::Int)
    #
    # Compute the support of the coboundary of a Lefschetz complex cell
    #

    # Extract the boundary matrix
        
    bnd = lc.boundary

    # Create the vector of coboundary indices

    if length(bnd.rows[cellI]) == 0
        cobndvecI = Vector{Int}([])
    else
        cobndvecI = deepcopy(bnd.rows[cellI])
    end

    # Return the result

    return cobndvecI
end

"""
    lefschetz_coboundary(lc::AbstractComplex, cellS::String)

Compute the support of the coboundary of a Lefschetz complex cell.

This method returns the boundary support as a `Vector{String}`.
"""
function lefschetz_coboundary(lc::AbstractComplex, cellS::String)
    #
    # Compute the support of the coboundary of a Lefschetz complex cell
    #

    # Extract the boundary matrix
        
    bnd = lc.boundary

    # Convert to integer format, get the coboundary, convert back

    cellI = lc.indices[cellS]
    cobndvecI = lefschetz_coboundary(lc, cellI)
    cobndvecS = lc.labels[cobndvecI]

    # Return the result

    return cobndvecS
end

