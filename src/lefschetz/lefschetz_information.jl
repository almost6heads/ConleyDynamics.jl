export lefschetz_information

"""
    lefschetz_information(lc::AbstractComplex)

Extract basic information about a Lefschetz complex.

The input argument `lc` contains the Lefschetz complex. The function returns
the information in the form of a `Dict{String,Any}`. You can use the command
`keys` to see the keyset of the return dictionary:
- `"Dimension"`: Dimension of the Lefschetz complex
- `"Coefficient field"`: Underlying coefficient field
- `"Euler characteristic"`: Euler characteristic of the complex
- `"Homology"`: Betti numbers of the Lefschetz complex
- `"Boundary sparsity"`: Sparsity percentage of the boundary matrix
- `"Number of cells"`: Total number of cells in the complex
- `"Cell counts by dim"`: Cell counts by dimension
In the last case, the dictionary entry is a vector of pairs
`(dimension, cell count)`.
"""
function lefschetz_information(lc::AbstractComplex)
    #
    # Extract basic information about a Lefschetz complex
    #

    # Extract the dimension and the underlying field

    lfield = lefschetz_field(lc)

    # Extract cell information

    dmin = minimum(lc.dimensions)
    dmax = maximum(lc.dimensions)
    cellcount = Vector{Vector{Int}}()
    eulerchar = 0

    for dd in dmin:dmax
        ddcount = length(findall(isequal(dd),lc.dimensions))
        if ddcount > 0
            push!(cellcount, [dd,ddcount])
            eulerchar = eulerchar + (-1)^dd * ddcount
        end
    end

    # Extract sparsity information of the boundary matrix

    bndsparsity = round(sparse_sparsity(lc.boundary) * 10000.0) / 100.0

    # Compute the homology

    lchom = homology(lc)

    # Create the empty return dictionary

    rdict = Dict{String,Any}()

    # Determine the number of critical and regular multivectors

    rdict["Dimension"]            = lc.dim
    rdict["Number of cells"]      = lc.ncells
    rdict["Boundary sparsity"]    = bndsparsity
    rdict["Coefficient field"]    = lfield
    rdict["Euler characteristic"] = eulerchar
    rdict["Homology"]             = lchom
    rdict["Cell counts by dim"]   = cellcount

    # Return the results

    return rdict
end

