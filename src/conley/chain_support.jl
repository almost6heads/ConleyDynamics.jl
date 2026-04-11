export chain_support

"""
    chain_support(lc::AbstractComplex, cvec::SparseMatrix; coeff::Bool=false)

Determine the support of a chain given as a sparse vector.

This function returns the support of a chain, given in the form of a 
sparse vector. In the default usage, the function returns a `Vector{String}`
which contains the labels of all the cells which have nonzero coefficients
in the chain. If one passes the optional parameter `coeff=true`, then the
function returns two arguments: In addition to the vector of labels as
above, it also returns the vector of associated coefficients.
"""
function chain_support(lc::AbstractComplex, cvec::SparseMatrix; coeff::Bool=false)
    #
    # Determine the support of a chain given as a sparse vector
    #

    # Create empty return vectors

    clabels = Vector{String}([])

    if coeff
        ccoeffs = Vector{typeof(lc.boundary.one)}([])
    end

    # Loop through the chain and add the cells and coefficients

    nzind = sparse_get_nz_column(cvec,1)
    for k in nzind
        push!(clabels, lc.labels[k])
        if coeff
            push!(ccoeffs, cvec[k,1])
        end
    end

    # Return the results

    if !coeff
        return clabels
    else
        return clabels, ccoeffs
    end
end

