export sparse_subtract

"""
    sparse_subtract(A::SparseMatrix, B::SparseMatrix)

Subtract two sparse matrices.

The function returns `A - B`. Exceptions are raised if the
matrix difference is not defined or the entry types do not match.
"""
function sparse_subtract(A::SparseMatrix, B::SparseMatrix)
    #
    # Subtract two sparse matrices.
    #

    # Check whether the difference exists

    if !((A.ncol == B.ncol) && (A.nrow == B.nrow))
        error("Matrix difference undefined!")
    end

    if !(typeof(A.zero) == typeof(B.zero))
        error("The sparse matrices need to have the same type!")
    end

    if !(A.char == B.char)
        error("The sparse matrices have to be over the same field!")
    end

    # Perform the matrix difference

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(A.zero)}()
    tzero = A.zero

    for m=1:A.ncol
        diffindex = union(A.columns[m],B.columns[m])
        if length(diffindex) > 0
            for k in diffindex
                entrydiff = A[k,m] - B[k,m]
                if !(entrydiff == tzero)
                    push!(r,k)
                    push!(c,m)
                    push!(v,entrydiff)
                end
            end
        end
    end

    # Construct and return the sparse matrix difference

    sm = sparse_from_lists(A.nrow, A.ncol, A.char, A.zero, A.one, r, c, v)
    return sm
end

"""
    Base.:-(A::SparseMatrix, B::SparseMatrix)

Subtract two sparse matrices.

Exceptions are raised if the matrix difference is not defined 
or the entry types do not match.
"""
Base.:-(A::SparseMatrix,B::SparseMatrix) = sparse_subtract(A,B)

