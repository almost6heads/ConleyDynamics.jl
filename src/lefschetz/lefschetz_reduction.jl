export lefschetz_reduction

"""
    lefschetz_reduction(lc::LefschetzComplex, redpairs::Vector{Vector{Int}})

Apply a sequence of elementary reductions to a Lefschetz complex.

The reduction pairs have to be specified in the argument `redpairs`. Each entry
has to be a vector of length two which contains an elementary reduction pair
in index form. In particular, the dimensions of the two cells in the pair have
to differ by one, and once the pair is reached in the reduction sequence, one cell
has to be a face of the other. The function returns a new Lefschetz complex, where
all cells in `redpairs` have been removed.
"""
function lefschetz_reduction(lc::LefschetzComplex, redpairs::Vector{Vector{Int}})
    #
    # Apply a sequence of elementary reductions to a Lefschetz complex
    #
    
    # Make sure the reduction pairs are pairs

    for redpair in redpairs
        if !(length(redpair) == 2)
            error("The reduction pairs have to come in pairs!")
        end
    end

    # Make sure the dimensions of the reduction pairs differ by 1

    for redpair in redpairs
        if !(abs(lc.dimensions[redpair[1]] - lc.dimensions[redpair[2]]) == 1)
            error("The reduction pair dimensions have to differ by 1!")
        end
    end

    # Make sure the reduction pair cells are distinct

    if !(length(unique(reduce(vcat,redpairs))) == 2*length(redpairs))
        error("The reduction pair cells have to be distinct!")
    end

    # Extract the necessary Lefschetz complex information
    
    labels = deepcopy(lc.labels)
    dims   = deepcopy(lc.dimensions)
    bnd    = deepcopy(lc.boundary)

    # Loop through the elementary reductions

    removed_cells = Vector{Int}()

    for redpair in redpairs

        # Extract the two cells

        if dims[redpair[1]] < dims[redpair[2]]
            a, b = redpair
        else
            b, a = redpair
        end

        push!(removed_cells, a)
        push!(removed_cells, b)

        # Extract the incidence coefficient

        lambda = bnd[a,b]

        if lambda == 0
            error("Reduction pairs have to consist of a face and coface!")
        end

        # Modify the columns of the coboundary of a, except for the
        # removed cells. Since we have modified the boundary matrix,
        # we need to get the row of bnd directly.

        for c in setdiff(sparse_get_nz_row(bnd,a), removed_cells)
            sparse_add_column!(bnd, c, b, -bnd[a,c], lambda)
        end
    end

    # Determine the remaining cells

    redcells = setdiff(Vector{Int}(1:lc.ncells), removed_cells)

    redlabels = labels[redcells]
    reddims   = dims[redcells]
    redbnd    = sparse_minor(bnd, redcells, redcells)

    # Create and return the reduced Lefschetz complex

    lcred = LefschetzComplex(redlabels, reddims, redbnd; validate=false)
    return lcred
end

"""
    lefschetz_reduction(lc::LefschetzComplex, redpairs::Vector{Vector{String}})

Apply a sequence of elementary reductions to a Lefschetz complex.

The reduction pairs have to be specified in the argument `redpairs`. Each entry
has to be a vector of length two which contains an elementary reduction pair
in label form. In particular, the dimensions of the two cells in the pair have
to differ by one, and once the pair is reached in the reduction sequence, one cell
has to be a face of the other. The function returns a new Lefschetz complex, where
all cells in `redpairs` have been removed.
"""
function lefschetz_reduction(lc::LefschetzComplex, redpairs::Vector{Vector{String}})
    #
    # Apply a sequence of elementary reductions to a Lefschetz complex
    #

    redpairsI = Vector{Vector{Int}}()
    for rp in redpairs
        push!(redpairsI, [lc.indices[rp[1]], lc.indices[rp[2]]])
    end

    lcred = lefschetz_reduction(lc, redpairsI)
    return lcred
end

"""
    lefschetz_reduction(lc::LefschetzComplex, r1::Int, r2::Int)

Apply a single elementary reduction to a Lefschetz complex.

This method expects that the two cells `r1` and `r2` which form the reduction
pair are given in index form. The function returns the reduced Lefschetz complex.
"""
lefschetz_reduction(lc::LefschetzComplex, r1::Int, r2::Int) = lefschetz_reduction(lc,[[r1,r2]])

"""
    lefschetz_reduction(lc::LefschetzComplex, r1::String, r2::String)

Apply a single elementary reduction to a Lefschetz complex.

This method expects that the two cells `r1` and `r2` which form the reduction
pair are given in label form. The function returns the reduced Lefschetz complex.
"""
lefschetz_reduction(lc::LefschetzComplex, r1::String, r2::String) = lefschetz_reduction(lc,[[r1,r2]])

