export extract_multivectors

"""
    extract_multivectors(lc::AbstractComplex, mvf::Vector{Vector{Int}},
                         scells::Vector{Int})

Extract all multivectors containing a provided selection of cells.

The function returns all multivectors which contain at least one of the cells
in the input vector `scells`. The return argument has type `Vector{Vector{Int}}`.
"""
function extract_multivectors(lc::AbstractComplex, mvf::Vector{Vector{Int}},
                              scells::Vector{Int})
    #
    # Extract all multivectors containing the provided cells
    #

    # Create the cell to multivector mapping

    cell2mv = fill(Int(0), lc.ncells)
    for k in eachindex(mvf)
        for mvcell in mvf[k]
            cell2mv[mvcell] = k
        end
    end

    # Extract the needed critical cells and multivectors

    indexc = findall(isequal(0),cell2mv[scells])
    indexm = findall(!isequal(0),cell2mv[scells])

    # Initialize the result with the multivectors

    mvextract = mvf[sort(unique(cell2mv[scells[indexm]]))]
    for k in indexc
        push!(mvextract,[scells[k]])
    end

    # Return the result

    return mvextract
end

"""
    extract_multivectors(lc::AbstractComplex, mvf::Vector{Vector{String}},
                         scells::Vector{String})

Extract all multivectors containing a provided selection of cells.

The function returns all multivectors which contain at least one of the cells
in the input vector `scells`. The return argument has type `Vector{Vector{String}}`.
"""
function extract_multivectors(lc::AbstractComplex, mvf::Vector{Vector{String}},
                              scells::Vector{String})
    #
    # Extract all multivectors containing the provided cells
    #

    mvfI = convert_cellsubsets(lc, mvf)

    scellsI = Vector{Int}()
    for k in scells
        push!(scellsI,lc.indices[k])
    end

    mvextractI = extract_multivectors(lc, mvfI, scellsI)

    mvextract = Vector{Vector{String}}()
    for m in mvextractI
        push!(mvextract, lc.labels[m])
    end

    return mvextract
end

