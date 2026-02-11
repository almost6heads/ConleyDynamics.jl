export sparse_show_nice

"""
    sparse_show_nice(sm::SparseMatrix, rlabels::Vector{String}, clabels::Vector{String})

Display a sparse matrix in readable format, with labels.

The input parameter `clabels` contains the labels for the columns, 
while `rlabels` give the row labels.

# Examples
```jldoctest
julia> lc, mvf = example_forman1d();

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show_nice(cm.matrix, cm.labels, cm.labels)
    A AD  F BF DE
 A  0  0  0  0  1
AD  0  0  0  0  0
 F  0  0  0  0  1
BF  0  0  0  0  0
DE  0  0  0  0  0
```
"""
function sparse_show_nice(sm::SparseMatrix, rlabels::Vector{String}, clabels::Vector{String})
    #
    # Display a sparse matrix in readable format, with labels
    #

    # Make sure the row and column label vectors have the correct lengths

    ncols = sm.ncol
    nrows = sm.nrow
    @assert ncols == length(clabels) "Number of labels does not equal the number of columns!"
    @assert nrows == length(rlabels) "Number of labels does not equal the number of rows!"

    # Determine the maximal field size

    maxl = maximum([maximum(length.(clabels)), maximum(length.(rlabels))])

    for k=1:ncols
        if length(sm.entries[k]) > 0
            maxentry = maximum(length.(string.(sm.entries[k])))
            maxl     = maximum([maxl, maxentry])
        end
    end

    # Print the header

    pstr = repeat(" ",maxl)   # +1 (if using separator below)
    for k=1:ncols
        nd   = 1 + maxl - length(clabels[k])
        pstr = pstr * repeat(" ",nd) * string(clabels[k])
    end
    println(pstr)

    # Print the matrix

    for m=1:nrows
        nd   = maxl - length(rlabels[m])
        pstr = repeat(" ",nd) * rlabels[m]   # * "│"
        for k=1:ncols
            estr = string(sm[m,k])
            nd   = 1 + maxl - length(estr)
            pstr = pstr * repeat(" ",nd) * estr
        end
        println(pstr)
    end
end


