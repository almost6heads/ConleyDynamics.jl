export LefschetzComplex, EuclideanComplex, ConleyMorseCM
export AbstractComplex
export Cell, Cells, CellSubsets

"""
    AbstractComplex

Abstract base type for Lefschetz complexes.

Both `LefschetzComplex` and `EuclideanComplex` are subtypes of `AbstractComplex`.
All topology, homology, and dynamics functions accept any `AbstractComplex`.
"""
abstract type AbstractComplex end

"""
    LefschetzComplex
    LefschetzComplex(labels, dimensions, boundary; validate=true)

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

The optional keyword argument `validate` controls whether the constructor
verifies that the boundary matrix squares to zero. By default this check is
enabled. Internal library functions that are guaranteed to produce a valid
complex by construction (simplicial, cubical, restriction, permutation,
basis-change, etc.) pass `validate=false` explicitly to suppress the check,
since for large complexes the matrix multiplication is expensive. Users
constructing a `LefschetzComplex` manually can also call
[`validate_lefschetz_complex`](@ref) to check an existing complex at any point.
"""
struct LefschetzComplex <: AbstractComplex
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
                              boundary::SparseMatrix;
                              validate::Bool=true)
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

        # Optionally verify that the boundary matrix squares to zero

        if validate
            if !sparse_is_zero(sparse_multiply(boundary, boundary))
                error("Boundary matrix does not square to zero: invalid Lefschetz complex.")
            end
        end

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
    EuclideanComplex
    EuclideanComplex(labels, dimensions, boundary, coords; validate=true)

A Lefschetz complex with embedded Euclidean coordinates for every cell.

The struct shares the fields of `LefschetzComplex`:
* `labels::Vector{String}`: Vector of labels associated with cell indices
* `dimensions::Vector{Int}`: Vector of cell dimensions
* `boundary::SparseMatrix`: Boundary matrix, columns give the cell boundaries
* `ncells::Int`: Number of cells
* `dim::Int`: Dimension of the complex (must satisfy `dim ≤ 3`)
* `indices::Dict{String,Int}`: Dictionary for finding cell index from label

In addition it carries:
* `coords::Vector{Vector{Vector{Float64}}}`: Per-cell vertex coordinates.
  `coords[k]` is the list of vertex coordinate vectors needed to draw cell `k`:
  - Vertex: `[[x,y]]` (length 1)
  - Edge: `[[x1,y1],[x2,y2]]` (length 2)
  - Triangle: `[[x1,y1],[x2,y2],[x3,y3]]` (length 3)
  - Cubical 2-cell: `[[x1,y1],[x2,y2],[x3,y3],[x4,y4]]` (length 4)

An `EuclideanComplex` can be created from an existing `LefschetzComplex` using
[`lefschetz_to_euclidean`](@ref), and converted back using
[`euclidean_to_lefschetz`](@ref).
"""
struct EuclideanComplex <: AbstractComplex
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
    # Embedded coordinates
    #
    coords::Vector{Vector{Vector{Float64}}}
    #
    # Inner constructor
    #
    function EuclideanComplex(labels::Vector{String},
                              dimensions::Vector{Int},
                              boundary::SparseMatrix,
                              coords::Vector{Vector{Vector{Float64}}};
                              validate::Bool=true)
        #
        # Create a EuclideanComplex instance
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
        if !(length(coords) == ncells)
            error("The coords vector length has to be the number of cells!")
        end

        # Make sure the cell dimensions increase

        for k in 1:ncells-1
            if dimensions[k] > dimensions[k+1]
                error("The cells dimensions cannot decrease!")
            end
        end
        dim = dimensions[ncells]

        # Check that the dimension does not exceed 3

        if dim > 3
            error("EuclideanComplex only supports complexes of dimension at most 3!")
        end

        # Optionally verify that the boundary matrix squares to zero

        if validate
            if !sparse_is_zero(sparse_multiply(boundary, boundary))
                error("Boundary matrix does not square to zero: invalid complex.")
            end
        end

        # Create the label to indices dictionary

        indices = Dict{String,Int}([(labels[k],k) for k in 1:ncells])

        # Create the composite type

        new(labels, dimensions, boundary, ncells, dim, indices, coords)
    end
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
    Base.show(io::IO, ::MIME"text/plain", ec::EuclideanComplex)

Display EuclideanComplex information when hitting return in REPL.
"""
function Base.show(io::IO, ::MIME"text/plain", ec::EuclideanComplex)
    #
    # Display information for a Euclidean complex
    #

    pstr = string(typeof(ec)) * ": struct with the following fields:"
    println(io, pstr)
    println(io, "  labels, indices, dimensions, coords")
    println(io, "  dim:      " * string(ec.dim))
    println(io, "  ncells:   " * string(ec.ncells))
    print(io, "  boundary: field " * lefschetz_field(ec) *
              ", sparsity " * string(sparse_sparsity(ec.boundary)))
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

