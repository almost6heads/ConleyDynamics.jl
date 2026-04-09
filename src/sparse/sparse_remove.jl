export sparse_remove!

"""
    sparse_remove!(matrix::SparseMatrix, ri::Int, ci::Int)

Remove the sparse matrix entry at location `(ri,ci)`.
"""
function sparse_remove!(matrix::SparseMatrix, ri::Int, ci::Int)
    #
    # Remove the sparse matrix entry at location (ri,ci).
    #

    # Find the location of (ri,ci) in column ci

    index_in_col = searchsortedfirst(matrix.columns[ci], ri)

    if index_in_col <= length(matrix.columns[ci]) && matrix.columns[ci][index_in_col] == ri
        index_in_row = searchsortedfirst(matrix.rows[ri], ci)
        deleteat!(matrix.rows[ri],index_in_row)
        deleteat!(matrix.columns[ci],index_in_col)
        deleteat!(matrix.entries[ci],index_in_col)
    end
end

