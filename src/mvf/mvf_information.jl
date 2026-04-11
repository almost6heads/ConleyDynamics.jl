export mvf_information

"""
    mvf_information(lc::AbstractComplex, mvf::CellSubsets)

Extract basic information about a multivector field.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The function returns the information in the form of
a `Dict{String,Any}`. You can use the command `keys` to see the keyset of
the return dictionary:
- `"N mv"`: Number of multivectors
- `"N critical"`: Number of critcal multivectors
- `"N regular"`: Number of regular multivectors
- `"Lengths critical"`: Length distribution of critical multivectors
- `"Lengths regular"`: Length distribution of regular multivectors
In the last two cases, the dictionary entries are vectors of pairs
`(length,frequency)`, where each pair indicates that there are
`frequency` multivectors of length `length`.
"""
function mvf_information(lc::AbstractComplex, mvf::CellSubsets)
    #
    # Extract basic information about a multivector field
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

    # Create array of multivector lengths

    mvfIlen = length.(mvfI)

    # Determine the Betti number sums of the multivectors

    bnumsum = fill(Int(0), length(mvfI))

    Threads.@threads for k = 1:length(mvfI)
        bnumsum[k] = sum(conley_index(lc, mvfI[k]))
    end

    # Find the indices of the critical and regular multivectors

    indexcritical = findall(x -> x>0,  bnumsum)
    indexregular  = findall(x -> x==0, bnumsum)

    # Create the empty return dictionary

    rdict = Dict{String,Any}()

    # Determine the number of critical and regular multivectors

    rdict["N mv"]       = length(mvfI)
    rdict["N critical"] = length(indexcritical)
    rdict["N regular"]  = length(indexregular)

    # Determine the critical multivector length distributions

    lencritical = mvfIlen[indexcritical]
    distcritical = Vector{Vector{Int}}()

    for k in sort(unique(lencritical))
        kn = length(findall(x -> x==k, lencritical))
        push!(distcritical, [k,kn])
    end

    rdict["Lengths critical"] = distcritical

    # Determine the regular multivector length distributions

    lenregular = mvfIlen[indexregular]
    distregular = Vector{Vector{Int}}()

    for k in sort(unique(lenregular))
        kn = length(findall(x -> x==k, lenregular))
        push!(distregular, [k,kn])
    end

    rdict["Lengths regular"] = distregular

    # Return the results

    return rdict
end

