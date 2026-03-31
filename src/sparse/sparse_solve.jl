export sparse_solve

"""
    sparse_solve(A::SparseMatrix, b::SparseMatrix)

Find a solution of a sparse linear system.

This function computes a solution of the sparse linear
system `Ax = b`. It is expected that `A` and `b` are sparse
matrices over the same field, and that both have the same
number of rows. Usually, this function is used with a column
vector `b`. However, if `b` is a matrix, then the respective
matrix equation is solved. The function returns a particular
solution of the system, if one exists. Otherwise it
returns `nothing`.

In order to obtain all solutions of the system, one has to
add an arbitrary linear combination of the basis vectors
determined via `sparse_basis_kernel(A)`.

# Example
```jldoctest
julia> Afull = [1 0 1 2 4 3 2; 3 4 6 2 3 4 2; 2 0 2 4 1 1 1];

julia> A = sparse_from_full(Afull, p=7);

julia> bfull = [1 2 5 3; 3 6 6 4; 2 4 2 1];

julia> b = sparse_from_full(bfull, p=7);

julia> sparse_show(A)
 1 . 1 2 4 3 2
 3 4 6 2 3 4 2
 2 . 2 4 1 1 1

julia> sparse_show(b)
 1 2 5 3
 3 6 6 4
 2 4 2 1

julia> x = sparse_solve(A, b);

julia> sparse_show(x)
 1 2 3 .
 . . 5 .
 . . . .
 . . . .
 . . . .
 . . 3 1
 . . . .

julia> A*x == b
true
```
"""
function sparse_solve(A::SparseMatrix, b::SparseMatrix)
    #
    # Solve a sparse linear system
    #

    @assert A.nrow==b.nrow "A and b have to have the same number of rows!"
    @assert A.char==b.char "A and b have to have the same base field!"

    # Create the extended matrix and find it's rref

    Aext = [A b]
    sparse_rref!(Aext)

    # Deal with the trivial cases first

    if sparse_is_zero(Aext)
        x = sparse_zero(A.ncol, b.ncol, p=A.char)
        for k = 1:b.ncol
            x[1,k] = A.one
        end
        return x
    end

    # Find the number of pivots and the pivot columns

    npivot = length(findall(t -> t > 0, length.(Aext.rows)))
    pivotcols = [Aext.rows[k][1] for k = 1:npivot]

    # Deal with the case of no solution

    if maximum(pivotcols) > A.ncol
        return nothing
    end

    # Construct the answer

    x = sparse_zero(A.ncol, b.ncol, p=A.char)
    for k = 1:npivot
        pc = pivotcols[k]
        for j = 1:b.ncol
            x[pc,j] = Aext[k, A.ncol + j]
        end
    end

    return x
end

