export sparse_transpose

"""
    sparse_transpose(A::SparseMatrix)

Create the transpose of a sparse matrix.

This does seem pretty self-explanatory, doesn't it?
"""
function sparse_transpose(A::SparseMatrix)
    #
    # Transpose of a sparse matrix
    #

    # Construct the transpose

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(A.zero)}()
    tzero = A.zero

    for k=1:A.ncol
        if length(A.columns[k]) > 0
            for m in 1:length(A.columns[k])
                push!(r, k)
                push!(c, A.columns[k][m])
                push!(v, A.entries[k][m])
            end
        end
    end

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(A.ncol, A.nrow, A.char, A.zero, A.one, r, c, v)
    return sm
end

"""
    Base.adjoint(A::SparseMatrix)

Create the transpose of a sparse matrix.

This method allows one to use `A'` for the computation of the 
transpose matrix.
"""
Base.adjoint(A::SparseMatrix) = sparse_transpose(A)

