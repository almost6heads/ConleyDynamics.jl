export restrict_dynamics

"""
    restrict_dynamics(lc::AbstractComplex, mvf::CellSubsets, lcsub::Cells)

Restrict a multivector field to a Lefschetz subcomplex.

For a given multivector field `mvf` on a Lefschetz complex `lc`, and a subcomplex
which is given by the locally closed set represented by `lcsub`, create the
associated Lefschetz subcomplex `lcreduced` and the induced multivector
field `mvfreduced` on the subcomplex. The multivectors of the new multivector
field are the intersections of the original multivectors and the subcomplex.
"""
function restrict_dynamics(lc::AbstractComplex, mvf::CellSubsets, lcsub::Cells)
    #
    # Restrict a multivector field to a Lefschetz subcomplex
    #

    # Convert the multivector field and cell subset to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
        return_as_labels = false
    else
        mvfI = convert_cellsubsets(lc, mvf)
        return_as_labels = true
    end

    if lcsub isa Vector{Int}
        lcsubI = lcsub
    else
        lcsubI = convert_cells(lc, lcsub)
    end

    # Create the restricted multivector field

    mvfredI = deepcopy(mvfI)
    Threads.@threads for k in eachindex(mvfredI)
        intersect!(mvfredI[k], lcsubI)
    end
    mvfredI = mvfredI[findall(x -> length(x)>1, mvfredI)]
    mvfredL = convert_cellsubsets(lc, mvfredI)

    # Create the reduced Lefschetz complex

    lcred = lefschetz_subcomplex(lc, lcsubI)
 
    # Return the results

    if return_as_labels
        mvfred = mvfredL
    else
        mvfred = convert_cellsubsets(lcred, mvfredL)
    end

    return lcred, mvfred
end

