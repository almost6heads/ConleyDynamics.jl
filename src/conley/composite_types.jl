export LefschetzComplex, ConleyMorseCM, Cells, CellSubsets

"""
    LefschetzComplex

Collect the Lefschetz complex information in a struct.

The struct is created via the following fields:
* `labels::Vector{String}`: Vector of labels associated with cell indices
* `dimensions::Vector{Int}`: Vector cell dimensions
* `boundary::SparseMatrix`: Boundary matrix, columns give the cell boundaries
It is expected that the dimensions are given in increasing order, and that
the square of the boundary matrix is zero. Otherwise, exceptions are raised.
In addition, the following fields are created during initialization:
* `ncells::Int`: Number of cells
* `dim::Int`: Dimension of the complex
* `indices::Dict{String,Int}`: Dictionary for finding cell index from label
The coefficient field is specified by the boundary matrix.

!!! warning
    Note that the constructor does not check whether the boundary matrix
    squares to zero. It is the responsibility of the user to ensure that!
"""
struct LefschetzComplex
    #
    # Fields that have to be declared
    #
    labels::Vector{String}
    dimensions::Vector{Int}
    boundary::SparseMatrix
    #
    # Fields that will be created
    #
    ncells::Int
    dim::Int
    indices::Dict{String,Int}
    #
    # Inner constructor
    #
    function LefschetzComplex(labels::Vector{String},
                              dimensions::Vector{Int},
                              boundary::SparseMatrix)
        #
        # Create a Lefschetz complex instance
        #
        
        # Perform basic length checks

        ncells = length(labels)
        if !(ncells == length(dimensions))
            error("Input vectors need to have the same length!")
        end
        if !(sparse_size(boundary,1) == sparse_size(boundary,2))
            error("The boundary matrix has to be square!")
        end
        if !(sparse_size(boundary,1) == ncells)
            error("The boundary matrix size has to be the number of cells!")
        end

        # Make sure the cell dimensions increase

        for k in 1:ncells-1
            if dimensions[k] > dimensions[k+1]
                error("The cells dimensions cannot decrease!")
            end
        end
        dim = dimensions[ncells]

        # Create the label to indices dictionary

        indices = Dict{String,Int}([(labels[k],k) for k in 1:ncells])

        # Create the composite type
        
        new(labels, dimensions, boundary, ncells, dim, indices)
    end
end

"""
    ConleyMorseCM{T}

Collect the connection matrix information in a struct.

The struct has the following fields:
* `matrix::SparseMatrix{T}`: Connection matrix
* `columns::Vector{Int}`: Corresponding columns in the boundary matrix
* `poset::Vector{Int}`: Poset indices for the connection matrix columns
* `labels::Vector{String}`: Labels for the connection matrix columns
* `morse::Vector{Vector{String}}`: Vector of Morse sets in original complex
* `conley::Vector{Vector{Int}}`: Vector of Conley indices for the Morse sets
* `complex::LefschetzComplex`: The Conley complex as a Lefschetz complex
"""
struct ConleyMorseCM{T}
    matrix::SparseMatrix{T}
    columns::Vector{Int}
    poset::Vector{Int}
    labels::Vector{String}
    morse::Vector{Vector{String}}
    conley::Vector{Vector{Int}}
    complex::LefschetzComplex
end

"""
    Cell = Union{Int,String}

A cell of a Lefschetz complex.

This data type is used to represent a cell of a Lefschetz
complex. The cell can be specified either via its index,
or its label.
"""
Cell = Union{Int,String}

"""
    Cells = Union{Vector{Int},Vector{String}}

A list of cells of a Lefschetz complex.

This data type is used to represent subsets of a Lefschetz
complex. It is used for individual isolated invariant sets,
locally closed subsets, and multivectors.
"""
Cells = Union{Vector{Int},Vector{String}}

"""
    CellSubsets = Union{Vector{Vector{Int}},Vector{Vector{String}}}

A collection of cell lists.

This data type is used to represent a collection of subsets of
a Lefschetz complex. It is used for Morse decompositions and
for multivector fields.
"""
CellSubsets = Union{Vector{Vector{Int}},Vector{Vector{String}}}

######################################################################

"""
    Base.show(io::IO, ::MIME"text/plain", lc::LefschetzComplex)

Display Lefschetz complex information when hitting return in REPL.
"""
function Base.show(io::IO, ::MIME"text/plain", lc::LefschetzComplex)
    #
    # Display information for a Lefschetz complex
    #

    # Display the size and type info

    pstr = string(typeof(lc)) * ": struct with the following fields:"
    println(io, pstr)
    println(io, "  labels, indices, dimensions")
    println(io, "  dim:      " * string(lc.dim))
    println(io, "  ncells:   " * string(lc.ncells))
    print(io, "  boundary: field " * lefschetz_field(lc) *
              ", sparsity " * string(sparse_sparsity(lc.boundary)))
end

"""
    Base.show(io::IO, ::MIME"text/plain", cm::ConleyMorseCM)

Display connection matrix information when hitting return in REPL.
"""
function Base.show(io::IO, ::MIME"text/plain", cm::ConleyMorseCM)
    #
    # Display information for a connection matrix
    #

    # Display the size and type info
    
    pstr = string(typeof(cm)) * ": struct with the following fields:"
    println(io, pstr)
    println(io, "  columns, poset, labels, morse, conley, complex")
    pstr = "  matrix: "
    pstr = pstr * string(cm.matrix.nrow) * "x" * string(cm.matrix.ncol)
    pstr = pstr * "-dimensional matrix with sparsity "
    pstr = pstr * string(sparse_sparsity(cm.matrix))
    println(io, pstr)
    print(io,"  inspect the connection matrix with sparse_show(cm)")
end

