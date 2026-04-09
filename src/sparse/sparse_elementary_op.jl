export sparse_add_column!, sparse_add_row!

"""
    sparse_add_column!(matrix::SparseMatrix, ci1::Int, ci2::Int, cn, cd)

Replace `column[ci1]` by `column[ci1] + (cn/cd) * column[ci2]`.
"""
function sparse_add_column!(matrix::SparseMatrix, ci1::Int, ci2::Int, cn, cd)
    #
    # Replace column[ci1] by column[ci1] + (cn/cd) * column[ci2]
    #

    # Find all nonzero entries in column[ci2]

    nonzero_rows    = matrix.columns[ci2]
    nonzero_entries = matrix.entries[ci2]

    # Loop though the entries and perform the column operation

    for k in 1:length(nonzero_rows)
        row_index = nonzero_rows[k]
        new_value = cn * nonzero_entries[k] / cd
        new_value += matrix[row_index,ci1]
        matrix[row_index,ci1] = new_value
    end
end

"""
    sparse_add_column!(matrix::SparseMatrix{Int}, ci1::Int, ci2::Int,
                       cn::Int, cd::Int)

Replace `column[ci1]` by `column[ci1] + (cn/cd) * column[ci2]`.

The computation is performed mod p, where the characteristic is
taken from `matrix.char`. An error is thrown if `matrix.char==0`.
"""
function sparse_add_column!(matrix::SparseMatrix{Int}, ci1::Int, ci2::Int,
                            cn::Int, cd::Int)
    #
    # Replace column[ci1] by column[ci1] + (cn/cd) * column[ci2]
    #

    p = matrix.char

    # Find all nonzero entries in column[ci2]

    nonzero_rows    = matrix.columns[ci2]
    nonzero_entries = matrix.entries[ci2]

    # Loop though the entries and perform the column operation

    for k in 1:length(nonzero_rows)
        row_index = nonzero_rows[k]
        new_value = cn * nonzero_entries[k] * invmod(cd, p)
        new_value = mod(new_value + matrix[row_index,ci1], p)
        matrix[row_index,ci1] = new_value
    end
end

"""
    sparse_add_row!(matrix::SparseMatrix, ri1::Int, ri2::Int, cn, cd)

Replace `row[ri1]` by `row[ri1] + (cn/cd) * row[ri2]`.
"""
function sparse_add_row!(matrix::SparseMatrix, ri1::Int, ri2::Int, cn, cd)
    #
    # Replace row[ri1] by row[ri1] + (cn/cd) * row[ri2]
    #

    # Find all nonzero entries in row[ri2]

    nonzero_columns = matrix.rows[ri2]

    # Loop though the entries and perform the column operation

    for kc in nonzero_columns
        new_value = cn * matrix[ri2,kc] / cd
        new_value += matrix[ri1,kc]
        matrix[ri1,kc] = new_value
    end
end

"""
    sparse_add_row!(matrix::SparseMatrix{Int}, ri1::Int, ri2::Int,
                    cn::Int, cd::Int)

Replace `row[ri1]` by `row[ri1] + (cn/cd) * row[ri2]`.

The computation is performed mod p, where the characteristic is
taken from `matrix.char`. An error is thrown if `matrix.char==0`.
"""
function sparse_add_row!(matrix::SparseMatrix{Int}, ri1::Int, ri2::Int,
                         cn::Int, cd::Int)
    #
    # Replace row[ri1] by row[ri1] + (cn/cd) * row[ri2]
    #

    p = matrix.char

    # Find all nonzero entries in row[ri2]

    nonzero_columns = matrix.rows[ri2]

    # Loop though the entries and perform the column operation

    for kc in nonzero_columns
        new_value = cn * matrix[ri2,kc] * invmod(cd, p)
        new_value = mod(new_value + matrix[ri1,kc], p)
        matrix[ri1,kc] = new_value
    end
end

