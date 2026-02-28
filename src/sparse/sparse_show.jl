export sparse_show

"""
    sparse_show(sm::SparseMatrix)

Display a sparse matrix in readable format.

This function will print the complete matrix, so be careful
with sparse matrices of large size!
"""
function sparse_show(sm::SparseMatrix)
    #
    # Display a sparse matrix in readable format
    #
    ncols = sm.ncol
    nrows = sm.nrow

    # Determine the maximal field size

    maxl = 1
    for k=1:ncols
        if length(sm.entries[k]) > 0
            maxentry = maximum(length.(string.(sm.entries[k])))
            maxl     = maximum([maxl, maxentry])
        end
    end

    # Print the matrix

    for m=1:nrows
        pstr = ""
        for k=1:ncols
            estr = string(sm[m,k])
            nd   = 1 + maxl - length(estr)
            pstr = pstr * repeat(" ",nd) * estr
        end
        println(pstr)
    end
end

"""
    sparse_show(sm::SparseMatrix, rlabels::Vector{String}, clabels::Vector{String})

Display a sparse matrix in readable format, with labels.

The input parameter `clabels` contains the labels for the columns, 
while `rlabels` give the row labels.

# Examples
```jldoctest
julia> lc, mvf = example_forman1d();

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show(cm.matrix, cm.labels, cm.labels)
  ┆  A CD  F BF DE
--┆---------------
 A┆  0  0  0  0  1
CD┆  0  0  0  0  0
 F┆  0  0  0  0  1
BF┆  0  0  0  0  0
DE┆  0  0  0  0  0
```
"""
function sparse_show(sm::SparseMatrix, rlabels::Vector{String}, clabels::Vector{String})
    #
    # Display a sparse matrix in readable format, with labels
    #

    # Make sure the row and column label vectors have the correct lengths

    ncols = sm.ncol
    nrows = sm.nrow
    @assert ncols == length(clabels) "Number of labels does not equal the number of columns!"
    @assert nrows == length(rlabels) "Number of labels does not equal the number of rows!"

    # Determine the maximal field size

    maxl = maximum([1, maximum(length.(clabels)), maximum(length.(rlabels))])

    for k=1:ncols
        if length(sm.entries[k]) > 0
            maxentry = maximum(length.(string.(sm.entries[k])))
            maxl     = maximum([maxl, maxentry])
        end
    end

    # Print the header

    pstr = repeat(" ",maxl) * "┆"
    for k=1:ncols
        nd   = 1 + maxl - length(clabels[k])
        pstr = pstr * repeat(" ",nd) * string(clabels[k])
    end
    println(pstr)

    pstr = repeat("-",maxl) * "┆"
    pstr = pstr * repeat("-", (maxl+1)*ncols)
    println(pstr)

    # Print the matrix

    for m=1:nrows
        nd   = maxl - length(rlabels[m])
        pstr = repeat(" ",nd) * rlabels[m] * "┆"
        for k=1:ncols
            estr = string(sm[m,k])
            nd   = 1 + maxl - length(estr)
            pstr = pstr * repeat(" ",nd) * estr
        end
        println(pstr)
    end
end

"""
    sparse_show(sm::SparseMatrix, labels::Vector{String})

Display a sparse matrix in readable format, with labels.

This function assumes that the matrix is square, and that the labels
are the same for rows and columns. They are provided in `labels`.
"""
function sparse_show(sm::SparseMatrix, labels::Vector{String})
    #
    # Display a sparse matrix in readable format, with labels
    #
    sparse_show(sm, labels, labels)
end

"""
    sparse_show(cm::ConleyMorseCM)

Display the connection matrix with labels.

This function displays the (sparse) connection matrix including
its labels. It uses `sparse_show(cm.matrix, cm.labels, cm.labels)`.

# Examples
```jldoctest
julia> lc, mvf = example_forman1d();

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show(cm)
  ┆  A CD  F BF DE
--┆---------------
 A┆  0  0  0  0  1
CD┆  0  0  0  0  0
 F┆  0  0  0  0  1
BF┆  0  0  0  0  0
DE┆  0  0  0  0  0
```
"""
function sparse_show(cm::ConleyMorseCM)
    #
    # Display the connection matrix with labels
    #
    sparse_show(cm.matrix, cm.labels, cm.labels)
end

"""
    Base.show(io::IO, ::MIME"text/plain", sm::SparseMatrix)

Display the sparse matrix `sm` when hitting return in REPL.
"""
function Base.show(io::IO, ::MIME"text/plain", sm::SparseMatrix)
    #
    # Display information for a sparse matrix
    #
    Nshow = 50

    # Display the size and type info

    pstr = string(sm.nrow) * "x" * string(sm.ncol) * "-dimensional "
    pstr = pstr * string(typeof(sm)) * ", "
    pstr = pstr * "sparsity=" * string(sparse_sparsity(sm)) * ":"
    println(io, pstr)

    # Only basic information if the matrix is too large

    if (sm.nrow > Nshow) || (sm.ncol > Nshow)
        print(io, "Use `sparse_show` for printing this large matrix!")
        return
    end

    # Otherwise, print the matrix

    ncols = sm.ncol
    nrows = sm.nrow
    maxl = 0
    for k=1:ncols
        if length(sm.entries[k]) > 0
            maxentry = maximum(length.(string.(sm.entries[k])))
            maxl     = maximum([maxl, maxentry])
        end
    end
    for m=1:nrows
        pstr = ""
        for k=1:ncols
            estr = string(sm[m,k])
            nd   = 1 + maxl - length(estr)
            pstr = pstr * repeat(" ",nd) * estr
        end
        if (m < nrows)
            println(io, pstr)
        else
            print(io, pstr)
        end
    end
end

