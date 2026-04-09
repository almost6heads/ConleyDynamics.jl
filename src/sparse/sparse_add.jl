export sparse_add

"""
    sparse_add(A::SparseMatrix, B::SparseMatrix)

Add two sparse matrices.

Exceptions are raised if the matrix sum is not defined 
or the entry types do not match.
"""
function sparse_add(A::SparseMatrix, B::SparseMatrix)
    #
    # Add two sparse matrices.
    #

    # Check whether the sum exists

    if !((A.ncol == B.ncol) && (A.nrow == B.nrow))
        error("Matrix sum undefined!")
    end

    if !(typeof(A.zero) == typeof(B.zero))
        error("The sparse matrices need to have the same type!")
    end

    if !(A.char == B.char)
        error("The sparse matrices have to be over the same field!")
    end

    # Perform the matrix sum

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(A.zero)}()
    tzero = A.zero

    for m=1:A.ncol
        sumindex = union(A.columns[m],B.columns[m])
        if length(sumindex) > 0
            for k in sumindex
                entrysum = A[k,m] + B[k,m]
                if !(entrysum == tzero)
                    push!(r,k)
                    push!(c,m)
                    push!(v,entrysum)
                end
            end
        end
    end

    # Construct and return the sparse matrix sum

    sm = sparse_from_lists(A.nrow, A.ncol, A.char, A.zero, A.one, r, c, v)
    return sm
end

"""
    Base.:+(A::SparseMatrix, B::SparseMatrix)

Add two sparse matrices.

Exceptions are raised if the matrix sum is not defined 
or the entry types do not match.
"""
Base.:+(A::SparseMatrix,B::SparseMatrix) = sparse_add(A,B)

