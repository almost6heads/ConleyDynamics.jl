export sparse_inverse

"""
    sparse_inverse(matrix::SparseMatrix)

Compute the inverse of a sparse matrix.

The function automatically performs the computations over the
underlying field of the sparse matrix.
"""
function sparse_inverse(matrix::SparseMatrix)
    #
    # Compute the inverse of a sparse matrix
    #

    # Make sure the matrix is square, copy it, and initialize
    # the identity matrix

    if !(matrix.nrow == matrix.ncol)
        error("The matrix is not square!")
    end

    p     = matrix.char
    tzero = matrix.zero
    tone  = matrix.one
    n     = matrix.nrow
    mwork = deepcopy(matrix)
    minv  = sparse_identity(n, p=p)

    # Reduce mwork to the identity matrix using only
    # column operations

    for m=1:n

        # Make sure there is a nonzero entry at (m,m)

        if mwork[m,m] == tzero
            ki = findfirst(x -> x>m, matrix.rows[m])
            if ki == nothing
                error("The matrix is not invertible!")
            end
            k = matrix.rows[m][ki]
            permr = Vector{Int}(1:n)
            permc = Vector{Int}(1:n)
            permc[k] = m
            permc[m] = k
            mwork = sparse_permute(mwork, permr, permc)
            minv  = sparse_permute(minv,  permr, permc)
        end

        # Scale the m-th column so that mwork[m,m]= 1

        if !(mwork[m,m] == tone)
            lambdainv = scalar_inverse(mwork[m,m], p)
            for k in mwork.columns[m]
                mwork[k,m] = mwork[k,m] * lambdainv
            end
            for k in minv.columns[m]
                minv[k,m]  = minv[k,m]  * lambdainv
            end
        end

        # Zero out off-diagonal entries in m-th row

        nzentries = setdiff(mwork.rows[m], [m])

        for k in nzentries
            lambda = mwork[m,k]
            sparse_add_column!(mwork, k, m, -lambda, tone)
            sparse_add_column!(minv,  k, m, -lambda, tone)
        end
    end

    # Return the answer

    return minv
end

