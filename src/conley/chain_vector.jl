export chain_vector

"""
    chain_vector(lc::AbstractComplex, chcells::Vector{Int})

Create a sparse vector representing a chain.

This function returns a sparse matrix in the form of a column vector,
which has a `1` in every cell location indicated by the entries in the
input argument `chcells`.
"""
function chain_vector(lc::AbstractComplex, chcells::Vector{Int})
    #
    # Create a sparse vector representing a chain
    #

    # Create an empty sparse vector

    svec = sparse_zero(lc.ncells, 1, p=lc.boundary.char)

    # Enter 1's in the location of the specified cells

    tone = lc.boundary.one
    for k in chcells
        svec[k,1] = tone
    end

    # Return the sparse vector

    return svec
end

"""
    chain_vector(lc::AbstractComplex, chcells::Vector{String})

Create a sparse vector representing a chain.

This function returns a sparse matrix in the form of a column vector,
which has a `1` for every cell indicated by the labels listed in the
input argument `chcells`.
"""
function chain_vector(lc::AbstractComplex, chcells::Vector{String})
    #
    # Create a sparse vector representing a chain
    #

    chcellsI = convert_cells(lc, chcells)
    svec = chain_vector(lc, chcellsI)
    return svec
end

"""
    chain_vector(lc::AbstractComplex, chcell::Int)

Create a sparse vector representing a chain.

This function returns a sparse matrix in the form of a column vector,
which has a `1` at the cell index specified in the second argument.
"""
chain_vector(lc::AbstractComplex, chcell::Int) = chain_vector(lc,[chcell])

"""
    chain_vector(lc::AbstractComplex, chcell::String)

Create a sparse vector representing a chain.

This function returns a sparse matrix in the form of a column vector,
which has a `1` at the cell specified by the second argument.
"""
chain_vector(lc::AbstractComplex, chcell::String) = chain_vector(lc,[chcell])

"""
    chain_vector(lc::AbstractComplex, chcells::Vector{Int}, chcoeff)

Create a sparse vector representing a chain.

This function returns a sparse matrix in the form of a column vector, which
has the entry `chcoeff[k]` in the `chcells[k]`-th row. In other words, it 
constructs the vector representation of the chain consisting of the cells
specified in `chcells`, with coefficients as specified in `chcoeff`.
"""
function chain_vector(lc::AbstractComplex, chcells::Vector{Int}, chcoeff)
    #
    # Create a sparse vector representing a chain
    #

    # Create an empty sparse vector

    svec = sparse_zero(lc.ncells, 1, p=lc.boundary.char)

    # Enter the coefficients in the location of the specified cells

    for k in eachindex(chcells)
        svec[chcells[k],1] = chcoeff[k]
    end

    # Return the sparse vector

    return svec
end

"""
    chain_vector(lc::AbstractComplex, chcells::Vector{String}, chcoeff)

Create a sparse vector representing a chain.

This function returns a sparse matrix in the form of a column vector, which
has the entry `chcoeff[k]` in the row corresponding to the cell with label
`chcells[k]`. In other words, it constructs the vector representation of the
chain consisting of the cells specified in `chcells`, with coefficients as
specified in `chcoeff`.
"""
function chain_vector(lc::AbstractComplex, chcells::Vector{String}, chcoeff)
    #
    # Create a sparse vector representing a chain
    #

    chcellsI = convert_cells(lc, chcells)
    svec = chain_vector(lc, chcellsI, chcoeff)
    return svec
end

"""
    chain_vector(lc::AbstractComplex, chcell::Int, chcoeff)

Create a sparse vector representing a chain.

This short-cut method specifies a chain consiting of one cell and its
coefficient. The cell is given as its index.
"""
chain_vector(lc::AbstractComplex, chcell::Int, chcoeff) =
   chain_vector(lc,[chcell],[chcoeff])

"""
    chain_vector(lc::AbstractComplex, chcell::String, chcoeff)

Create a sparse vector representing a chain.

This short-cut method specifies a chain consiting of one cell and its
coefficient. The cell is given via its label.
"""
chain_vector(lc::AbstractComplex, chcell::String, chcoeff) =
   chain_vector(lc,[chcell],[chcoeff])

