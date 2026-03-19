export sparse_hcat
export sparse_vcat

################################
#                              #
#   Horizontal concatenation   #
#                              #
################################

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
    sparse_hcat(S1, S2, S3, ...)

Concatenate several sparse matrices horizontally.

This method of `sparse_hcat` allows for any number of 
arguments. If you want to concatenate a vector `vec`
of sparse matrices, use the call `sparse_hcat(vec...)`.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of rows.
"""
sparse_hcat(vargs::SparseMatrix...) = reduce(sparse_hcat, vargs)

"""
    Base.hcat(smvec::SparseMatrix...)

Concatenate several sparse matrices horizontally.

This method allows one to use the shortform `[S1 S2 S3...]`
for the horizontal concatenation, just like in the matrix case.
"""
Base.hcat(smvec::SparseMatrix...) = sparse_hcat(smvec...)

##############################
#                            #
#   Vertical concatenation   #
#                            #
##############################

"""
    sparse_vcat(A::SparseMatrix, B::SparseMatrix)

Concatenate two sparse matrices vertically.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of columns.
"""
function sparse_vcat(A::SparseMatrix, B::SparseMatrix)
    #
    # Concatenate two sparse matrices vertically
    #

    # Check whether the concatenation makes sense

    @assert A.ncol == B.ncol "The matrices need to have the same number of columns!"
    @assert A.char == B.char "The matrices have to be over the same field!"

    # Perform the matrix concatenation

    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(A)
    for k=1:B.ncol
        if length(B.columns[k]) > 0
            for m in 1:length(B.columns[k])
                push!(c, k)
                push!(r, B.columns[k][m] + nr)
                push!(v, B.entries[k][m])
            end
        end
    end
    nr = nr + B.nrow

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, v)
    return sm
end

"""
    sparse_vcat(S1, S2, S3, ...)

Concatenate several sparse matrices vertically.

This method of `sparse_vcat` allows for any number of 
arguments. If you want to concatenate a vector `vec`
of sparse matrices, use the call `sparse_vcat(vec...)`.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of columns.
"""
sparse_vcat(vargs::SparseMatrix...) = reduce(sparse_vcat, vargs)

"""
    Base.vcat(smvec::SparseMatrix...)

Concatenate several sparse matrices vertically.

This method allows one to use the shortform `[S1; S2; S3; ...]`
for the vertical concatenation, just like in the matrix case.
"""
Base.vcat(smvec::SparseMatrix...) = sparse_vcat(smvec...)

