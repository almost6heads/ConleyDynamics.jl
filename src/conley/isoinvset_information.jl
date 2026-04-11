export isoinvset_information

"""
    isoinvset_information(lc::AbstractComplex, mvf::CellSubsets, iis::Cells)

Compute basic information about an isolated invariant set.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The isolated invariant set is specified in the 
argument `iis`. The function returns the information in the form of a
`Dict{String,Any}`. The command `keys` can be used to see the keyset of
the return dictionary. These describe the following information:
- `"Conley index"` contains the Conley index of the 
  isolated invariant set.
- `"N multivectors"` contains the number of multivectors
  in the isolated invariant set.
"""
function isoinvset_information(lc::AbstractComplex, mvf::CellSubsets, iis::Cells)
    #
    # Compute basic information about an isolated invariant set
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = deepcopy(mvf)
    else
        mvfI = convert_cellsubsets(lc, mvf)
    end

    # Create the cell-to-mv map and check whether it is
    # a valid multivector field

    cell2mv = fill(Int(0), lc.ncells)

    for k = 1:length(mvfI)
        for m in mvfI[k]
            if cell2mv[m] == 0
                cell2mv[m] = k
            else
                error(" The multivector field is not a partition!")
            end
        end
    end

    # Add the critical cells explicitly to mvfI

    ccimplied = findall(x -> x==0, cell2mv)

    for cc in ccimplied
        push!(mvfI, [cc])
    end

    # Update the cell-to-mv map

    cell2mv = fill(Int(0), lc.ncells)

    for k = 1:length(mvfI)
        for m in mvfI[k]
            cell2mv[m] = k
        end
    end

    # Convert the isolated invariant set to integer form

    if iis isa Vector{Int}
        iisI = deepcopy(iis)
    else
        iisI = convert_cells(lc, iis)
    end

    # Create the empty return dictionary

    rdict = Dict{String,Any}()

    # Compute the Conley index

    rdict["Conley index"] = conley_index(lc, iis)

    # Determine the number of multivectors in the set

    rdict["N multivectors"] = length(unique(cell2mv[iisI]))

    # Return the results

    return rdict
end

