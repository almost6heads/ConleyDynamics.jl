export lefschetz_cell_count

"""
    lefschetz_cell_count(lc::AbstractComplex; bounds::Bool=false)

Returns the number of cells in each dimension.

The function returns the number of cells in each dimension.
The return variable is of type `Vector{Int}`, has length `lc.dim + 1`,
and its k-th entry contains the number of cells in dimension k-1.
If the optional parameter `bounds=true` is passed, then the 
function also returns two integer vectors `lo` and `hi`. These
contain the beginning and end indices of the cells in each 
dimension.
"""
function lefschetz_cell_count(lc::AbstractComplex; bounds::Bool=false)
    #
    # Return the number of cells in each dimension
    #

    cellcount = Vector{Int}()

    # Determine the cell counts

    for k = 0:lc.dim
        dimcount = length(findall(x -> x==k, lc.dimensions))
        push!(cellcount, dimcount)
    end

    # Create the index bounds in each dimension

    if bounds == true
        lo = zeros(Int, 1 + lc.dim)
        hi = zeros(Int, 1 + lc.dim)
        offset = 0

        for k=0:lc.dim
            if !(cellcount[k+1] == 0)
                lo[k+1] = offset + 1
                hi[k+1] = offset + cellcount[k+1]
                offset  = offset + cellcount[k+1]
            end
        end
    end

    # Return the results
    
    if (bounds == false)
        return cellcount
    else
        return cellcount, lo, hi
    end
end

