export sparse_rref!
export sparse_rref
export sparse_is_rref

"""
    sparse_rref!(A::SparseMatrix)

Compute the reduced row Echelon form of the given matrix.

This function computes the reduced row Echelon form of a sparse 
matrix in place. The input matrix will be altered!

# Example
```jldoctest
julia> Afull = [1 0 1 2 4 3 2; 3 4 6 2 3 4 2; 2 0 2 4 1 1 1];

julia> A = sparse_from_full(Afull, p=7);

julia> sparse_rref!(A)

julia> sparse_show(A)
 1 . 1 2 4 . 3
 . 1 6 6 3 . 5
 . . . . . 1 2
```
"""
function sparse_rref!(A::SparseMatrix)
    #
    # Compute reduced row Echelon form
    #
    rows  = A.nrow
    cols  = A.ncol
    tone  = A.one
    tzero = A.zero
    p     = A.char
    pivot_row  = 1

    for col in 1:cols
        # Find pivot row: first non-zero entry in this column at or below pivot_row

        nonzero_col = A.columns[col]
        pivot = findfirst(t -> t >= pivot_row, nonzero_col)

        # No pivot in this column — move to next column
        if pivot === nothing
            continue
        end

        # Adjust index to be relative to the full matrix
        pivot = nonzero_col[pivot]

        # If pivot > pivot_row, add pivot_row to pivot
        if pivot > pivot_row
            sparse_add_row!(A, pivot_row, pivot, tone, tone)
        end

        # Scale the pivot row so the pivot element becomes 1
        fac = scalar_inverse(A[pivot_row,col], p)
        for j in A.rows[pivot_row]
            A[pivot_row,j] = scalar_multiply(A[pivot_row,j], fac, p)
        end

        # Eliminate all other entries in this column (above and below)
        nonzero_col = copy(A.columns[col])
        for row in nonzero_col
            if row != pivot_row
                sparse_add_row!(A, row, pivot_row, -A[row,col], tone)
            end
        end

        pivot_row += 1
        if pivot_row > rows
            break
        end
    end
end

"""
    sparse_rref(A::SparseMatrix)

Compute the reduced row Echelon form of the given matrix.

This function computes and returns the reduced row Echelon
form of a sparse matrix. The argument `A` is not altered.

# Example
```jldoctest
julia> Afull = [1 0 1 2 4 3 2; 3 4 6 2 3 4 2; 2 0 2 4 1 1 1];

julia> A = sparse_from_full(Afull, p=2);

julia> sparse_show(A)
 1 . 1 . . 1 .
 1 . . . 1 . .
 . . . . 1 1 1

julia> sparse_show(sparse_rref(A))
 1 . . . . 1 1
 . . 1 . . . 1
 . . . . 1 1 1
```
"""
function sparse_rref(A::SparseMatrix)
    #
    # Compute reduced row Echelon form
    #
    Acopy = deepcopy(A)
    sparse_rref!(Acopy)
    return Acopy
end

"""
    sparse_is_rref(A::SparseMatrix)

Check whether a matrix is in reduced row Echelon format.

This function determines whether a given matrix is in reduced
row Echelon format. It returns the appropriate boolean value.
"""
function sparse_is_rref(A::SparseMatrix)
    #
    # Check whether a matrix has reduced row Echelon form
    #

    # Let's get the trivial case out of the way
    if sparse_is_zero(A)
        return true
    end

    # Check that the nonempty rows are the first ones
    rowlengths = length.(A.rows)
    nzrows = findall(t -> t > 0, rowlengths)

    if !(nzrows == collect(1:length(nzrows)))
        return false
    end

    # Make sure the pivots are strictly increasing
    for k = 1:length(nzrows)-1
        if minimum(A.rows[k]) >= minimum(A.rows[k+1])
            return false
        end
    end

    # Make sure the pivots are 1
    for k = 1:length(nzrows)-1
        if !(A[k,minimum(A.rows[k])] == A.one)
            return false
        end
    end

    # Make sure the pivot columns only contain the pivot
    for k = 1:length(nzrows)-1
        if !(length(A.columns[minimum(A.rows[k])]) == 1)
            return false
        end
    end

    # We have passed all tests, so it is in rref
    return true
end

