export sparse_hcat

"""
    sparse_hcat(A::SparseMatrix, B::SparseMatrix)

Concatenate two sparse matrices horizontally.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of rows.
"""
function sparse_hcat(A::SparseMatrix, B::SparseMatrix)
    #
    # Concatenate two sparse matrices horizontally
    #

    # Check whether the concatenation makes sense

    @assert A.nrow == B.nrow "The matrices need to have the same number of rows!"
    @assert A.char == B.char "The matrices have to be over the same field!"

    # Perform the matrix concatenation

    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(A)
    for k=1:B.ncol
        if length(B.columns[k]) > 0
            for m in 1:length(B.columns[k])
                push!(c, k+nc)
                push!(r, B.columns[k][m])
                push!(v, B.entries[k][m])
            end
        end
    end
    nc = nc + B.ncol

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, v)
    return sm
end

"""
    sparse_hcat(v)

Concatenate a vector of sparse matrices horizontally.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of rows.
"""
sparse_hcat(v::Vector) = reduce(sparse_hcat, v)

