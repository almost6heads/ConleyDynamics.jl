export sparse_scale

"""
    sparse_scale(sfac, A::SparseMatrix)

Scale a sparse matrix by a scalar.

An exception is raised if the entry types do not match.
"""
function sparse_scale(sfac, A::SparseMatrix)
    #
    # Scale a sparse matrix.
    #

    # Check whether the scalar multiplication is possible

    if !(typeof(A.zero) == typeof(sfac))
        error("The sparse matrix and the scalar have to have the same type!")
    end

    # Perform the scalar multiplication

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(A.zero)}()
    tzero = A.zero

    for m=1:A.ncol
        for k in A.columns[m]
            newentry = sfac * A[k,m]
            if !(newentry == tzero)
                push!(r,k)
                push!(c,m)
                push!(v,newentry)
            end
        end
    end

    # Construct and return the scalar multiple

    sm = sparse_from_lists(A.nrow, A.ncol, A.char, A.zero, A.one, r, c, v)
    return sm
end

"""
    Base.:*(sfac, A::SparseMatrix)

Compute the scalar product of `sfac` and the sparse matrix `A`.

An exception is raised if the scalar `sfac` does not have the
same type as the matrix entries.
"""
Base.:*(sfac,A::SparseMatrix) = sparse_scale(sfac,A)

