export sparse_size, sparse_low, sparse_identity, sparse_zero
export sparse_fullness, sparse_sparsity, sparse_nz_count
export sparse_is_zero, sparse_is_equal, sparse_is_identity
export scalar_inverse

"""
    sparse_size(matrix::SparseMatrix, dim::Int)

Number of rows (`dim=1`) or columns (`dim=2`) of a sparse matrix.
"""
function sparse_size(matrix::SparseMatrix, dim::Int)
    #
    # Return the number of columns or rows
    #

    if dim == 1
        return matrix.nrow
    elseif dim == 2
        return matrix.ncol
    else
        error("dim has to be 1 or 2!")
    end
end

"""
    sparse_low(matrix::SparseMatrix, col::Int)

Row index of the lowest nonzero matrix entry in column `col`.
"""
function sparse_low(matrix::SparseMatrix, col::Int)
    #
    # Return the location of the lowest nonzero entry in 
    # column col of the matrix
    #

    if length(matrix.columns[col]) == 0
        return 0
    else
        return last(matrix.columns[col])
    end
end

"""
    sparse_identity(n::Int; p::Int=0)

Create a sparse identity matrix with n rows and columns.

The optional argument `p` specifies the field characteristic.
If `p=0` then the sparse matrix is over the rationals, while
if `p>0` is a prime, then the matrix is an integer matrix 
whose entries are interpreted in `GF(p)`.
"""
function sparse_identity(n::Int; p::Int=0)
    #
    # Create a sparse identity matrix with n rows and columns,
    # and with diagonal entry tone.
    #

    # Initialize the variables and lists

    if p == 0
        tzero = Rational{Int}(0)
        tone  = Rational{Int}(1)
    elseif p > 0
        tzero = Int(0)
        tone  = Int(1)
    else
        error("The characteristic cannot be negative!")
    end

    r = Vector{Int}(1:n)
    vals = fill(tone,n)

    # Create the sparse identity and return it

    sm = sparse_from_lists(n, n, p, tzero, tone, r, r, vals)
    return sm
end

"""
    sparse_zero(nr::Int, nc::Int; p::Int=0)

Create a sparse zero matrix with `nr` rows and `nc` columns.

The optional argument `p` specifies the field characteristic.
If `p=0` then the sparse matrix is over the rationals, while
if `p>0` is a prime, then the matrix is an integer matrix 
whose entries are interpreted in `GF(p)`.
"""
function sparse_zero(nr::Int, nc::Int; p::Int=0)
    #
    # Create a sparse zero matrix with nr rows and nc columns
    #

    # Initialize the variables and lists

    if p == 0
        tzero = Rational{Int}(0)
        tone  = Rational{Int}(1)
    elseif p > 0
        tzero = Int(0)
        tone  = Int(1)
    else
        error("The characteristic cannot be negative!")
    end

    r = Vector{Int}()

    # Create the sparse zero matrix and return it

    sm = sparse_from_lists(nr, nc, p, tzero, tone, r, r, r)
    return sm
end

"""
    sparse_fullness(sm::SparseMatrix)

Return the fullness of the sparse matrix `sm`,
which equals the percentage of nonzero elements.
"""
function sparse_fullness(sm::SparseMatrix)
    #
    # Compute the fullness of a sparse matrix
    #

    nz = 0
    for k=1:sm.ncol
        nz += length(sm.columns[k])
    end
    fullness = Float64(nz) / (Float64(sm.ncol) * Float64(sm.nrow))
    
    return fullness
end

"""
    sparse_sparsity(sm::SparseMatrix)

Return the sparsity of the sparse matrix `sm`,
which equals the percentage of zero entries.
"""
function sparse_sparsity(sm::SparseMatrix)
    #
    # Compute the sparsity of a sparse matrix
    #

    return 1.0 - sparse_fullness(sm)
end

"""
    sparse_nz_count(sm::SparseMatrix)

Return the number of nonzero entries of the sparse matrix `sm`.
"""
function sparse_nz_count(sm::SparseMatrix)
    #
    # Compute the number of nonzero entries
    #

    nz = 0
    for k=1:sm.ncol
        nz += length(sm.columns[k])
    end
    
    return nz
end

"""
    sparse_is_zero(sm::SparseMatrix)

Test whether the sparse matrix `sm` is the zero matrix.
"""
function sparse_is_zero(sm::SparseMatrix)
    #
    # Test whether a sparse matrix is zero
    #

    nz = 0
    for k=1:sm.ncol
        nz += length(sm.columns[k])
    end

    if nz == 0
        return true
    else
        return false
    end
end

"""
    sparse_is_equal(A::SparseMatrix, B::SparseMatrix)

Test whether the sparse matrices `A` and `B` are equal.

For equality, the two sparse matrices not only have to have the
same size and the same entries, they also have to be of the same
type and be defined over the same field.
"""
function sparse_is_equal(A::SparseMatrix, B::SparseMatrix)
    #
    # Equality test for sparse matrices
    #
    
    # First check the type and size

    if !(A.nrow == B.nrow); return false end
    if !(A.ncol == B.ncol); return false end
    if !(A.char == B.char); return false end
    if !(typeof(A.one) == typeof(B.one)); return false end

    # Now check the locations of the nonzero entries

    if !(A.columns == B.columns); return false end
    if !(A.rows == B.rows); return false end   # not necessary

    # Finally check the entries

    if !(A.entries == B.entries); return false end

    # They are the same!

    return true
end

"""
    Base.:(==)(A::SparseMatrix,B::SparseMatrix)

Test whether the sparse matrices `A` and `B` are equal.

For equality, the two sparse matrices not only have to have the
same size and the same entries, they also have to be of the same
type and be defined over the same field.
"""
Base.:(==)(A::SparseMatrix,B::SparseMatrix) = sparse_is_equal(A,B)

"""
    sparse_is_identity(A::SparseMatrix)

Test whether the sparse matrix `A` is the identity matrix.
"""
function sparse_is_identity(A::SparseMatrix)
    #
    # Test whether a sparse matrix is the identity
    #
    
    # Create the identity matrix with the same number of rows 
    # as A, and over the same field

    Imatrix = sparse_identity(A.nrow, p=A.char)

    # Return true if this matrix equals A

    return sparse_is_equal(A,Imatrix)
end

"""
    scalar_inverse(s, p::Int)

Compute the inverse of a scalar.
"""
function scalar_inverse(s, p::Int)
    #
    # Compute the inverse of a scalar
    #
    return 1 / s
end

"""
    scalar_inverse(s::Int, p::Int)

Compute the inverse of a scalar.

This function computes the inverse in modular arithmetic
with base `p`.
"""
function scalar_inverse(s::Int, p::Int)
    #
    # Compute the inverse of a scalar
    #
    return invmod(s, p)
end

