export sparse_get_entry, sparse_set_entry!
export sparse_get_column
export sparse_get_nz_column, sparse_get_nz_row
export getindex, setindex!
       
"""
    sparse_get_entry(matrix::SparseMatrix, ri::Int, ci::Int)

Get the sparse matrix entry at location `(ri,ci)`.
"""
function sparse_get_entry(matrix::SparseMatrix, ri::Int, ci::Int)
    #
    # Get matrix[ri,ci]
    #

    # Find the location of (ri,ci) in column ci

    index_in_col = findfirst(x -> x==ri, matrix.columns[ci])

    # Determine the value and return it

    if index_in_col == nothing
        return matrix.zero
    else
        return matrix.entries[ci][index_in_col]
    end
end

"""
    Base.getindex(matrix::SparseMatrix, ri::Int, ci::Int)

Get the sparse matrix entry at location `(ri,ci)`.
"""
function Base.getindex(matrix::SparseMatrix, ri::Int, ci::Int)
    return sparse_get_entry(matrix, ri, ci)
end

"""
    sparse_set_entry!(matrix::SparseMatrix, ri::Int, ci::Int, val)

Set the sparse matrix entry at location `(ri,ci)` to `val`.
"""
function sparse_set_entry!(matrix::SparseMatrix, ri::Int, ci::Int, val)
    #
    # Set matrix[ri,ci] = val
    #

    # Determine the location of ri in the column ci, if it exists

    index_in_col = findfirst(x -> x==ri, matrix.columns[ci])

    # Incorporate the value into the matrix

    if val == matrix.zero
        # If the entry was present, it has to be removed
        if !(index_in_col == nothing)
            sparse_remove!(matrix, ri, ci)
        end
    elseif !(index_in_col == nothing)
        # Entry present, just needs to be updated
        matrix.entries[ci][index_in_col] = val
    else
        # Insert ci into row ri
        insert_index_row = findfirst(x -> x>ci, matrix.rows[ri])
        if insert_index_row == nothing
            push!(matrix.rows[ri], ci)
        else
            insert!(matrix.rows[ri], insert_index_row, ci)
        end

        # Insert ri into column ci and add the new entry
        insert_index_col = findfirst(x -> x>ri, matrix.columns[ci])
        if insert_index_col == nothing
            push!(matrix.columns[ci], ri)
            push!(matrix.entries[ci], val)
        else
            insert!(matrix.columns[ci], insert_index_col, ri)
            insert!(matrix.entries[ci], insert_index_col, val)
        end
    end
    return
end

"""
    sparse_set_entry!(matrix::SparseMatrix{Int}, ri::Int, ci::Int, val)

Set the sparse matrix entry at location `(ri,ci)` to `val`.

This method assumes that the sparse matrix is over the integers, which
means that it is interpreted over a finite field. Thus, the new entry
is automatically reduced via modular arithmetic.
"""
function sparse_set_entry!(matrix::SparseMatrix{Int}, ri::Int, ci::Int, val)
    #
    # Set matrix[ri,ci] = val
    #

    # Determine the location of ri in the column ci, if it exists

    index_in_col = findfirst(x -> x==ri, matrix.columns[ci])

    # Incorporate the value into the matrix

    if !(typeof(val) == typeof(matrix.zero))
        error("The value has to be an integer!")
    end

    val = mod(val, matrix.char)    # Reduce the value modulo p

    if val == matrix.zero
        # If the entry was present, it has to be removed
        if !(index_in_col == nothing)
            sparse_remove!(matrix, ri, ci)
        end
    elseif !(index_in_col == nothing)
        # Entry present, just needs to be updated
        matrix.entries[ci][index_in_col] = val
    else
        # Insert ci into row ri
        insert_index_row = findfirst(x -> x>ci, matrix.rows[ri])
        if insert_index_row == nothing
            push!(matrix.rows[ri], ci)
        else
            insert!(matrix.rows[ri], insert_index_row, ci)
        end

        # Insert ri into column ci and add the new entry
        insert_index_col = findfirst(x -> x>ri, matrix.columns[ci])
        if insert_index_col == nothing
            push!(matrix.columns[ci], ri)
            push!(matrix.entries[ci], val)
        else
            insert!(matrix.columns[ci], insert_index_col, ri)
            insert!(matrix.entries[ci], insert_index_col, val)
        end
    end
    return
end

"""
    Base.setindex!(matrix::SparseMatrix, val, ri::Int, ci::Int)

Set the sparse matrix entry at location `(ri,ci)` to `val`.
"""
function Base.setindex!(matrix::SparseMatrix, val, ri::Int, ci::Int)
    sparse_set_entry!(matrix, ri, ci, val)
    return
end

"""
    sparse_get_column(matrix::SparseMatrix, ci::Int)

Get the ci-th column of the sparse matrix.
"""
function sparse_get_column(matrix::SparseMatrix, ci::Int)
    #
    # Get matrix[:,ci]
    #

    # Initialize the empty column vector

    column = fill(matrix.zero, matrix.nrow)

    # Add the nonzero elements

    for k=1:length(matrix.columns[ci])
        column[matrix.columns[ci][k]] = matrix.entries[ci][k]
    end

    # Return the column

    return column
end

"""
    sparse_get_nz_column(matrix::SparseMatrix, ci::Int)

Get the row indices for the nonzero entries in the ci-th column
of the sparse matrix.
"""
function sparse_get_nz_column(matrix::SparseMatrix, ci::Int)
    #
    # Get all k for which matrix[k,ci]<>0
    #

    nzcol = deepcopy(matrix.columns[ci])
    return nzcol
end

"""
    sparse_get_nz_row(matrix::SparseMatrix, ri::Int)

Get the column indices for the nonzero entries in the ri-th row
of the sparse matrix.
"""
function sparse_get_nz_row(matrix::SparseMatrix, ri::Int)
    #
    # Get all k for which matrix[ri,k]<>0
    #

    nzrow = deepcopy(matrix.rows[ri])
    return nzrow
end

