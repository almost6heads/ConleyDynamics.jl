export sparse_multiply

"""
    sparse_multiply(A::SparseMatrix, B::SparseMatrix)

Multiply two sparse matrices.

Exceptions are raised if the matrix product is not defined 
or the entry types do not match.
"""
function sparse_multiply(A::SparseMatrix, B::SparseMatrix)
    #
    # Multiply two sparse matrices.
    #

    # Check whether the product exists

    if !(A.ncol == B.nrow)
        error("Matrix product undefined!")
    end

    if !(typeof(A.zero) == typeof(B.zero))
        error("The sparse matrices need to have the same type!")
    end

    if !(A.char == B.char)
        error("The sparse matrices have to be over the same field!")
    end

    # Create row and column index sets

    rset = fill(Set{Int}(), A.nrow)
    cset = fill(Set{Int}(), B.ncol)
    for k = 1:A.nrow
        rset[k] = Set(A.rows[k])
    end
    for k = 1:B.ncol
        cset[k] = Set(B.columns[k])
    end

    # Perform the matrix product

    echannel = Channel{Tuple{Int,Int,typeof(A.zero)}}(Inf)
    tzero    = A.zero

    Threads.@threads for k=1:A.nrow
        for m=1:B.ncol
            dpindex = intersect(rset[k],cset[m])
            if length(dpindex) > 0
                dotprod = tzero
                for j in dpindex
                    dotprod += A[k,j] * B[j,m]
                end
                if !(dotprod == tzero)
                    put!(echannel, (k,m,dotprod))
                end
            end
        end
    end

    # Collect the entries

    close(echannel)
    entries = collect(echannel)
    if length(entries) > 0
        r = [entries[k][1] for k in 1:length(entries)]
        c = [entries[k][2] for k in 1:length(entries)]
        v = [entries[k][3] for k in 1:length(entries)]
    else
        r = Vector{Int}()
        c = Vector{Int}()
        v = Vector{typeof(A.zero)}()
    end

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(A.nrow, B.ncol, A.char, A.zero, A.one, r, c, v)
    return sm
end

"""
    sparse_multiply(alpha::Int, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for finite fields. An exception is raised if
the sparse matrix is not defined over a finite field.
"""
function sparse_multiply(alpha::Int, A::SparseMatrix)
    #
    # Multiply a sparse matrix by a scalar
    #
    p = A.char
    @assert p > 0 "Wrong scalar type!"

    if alpha == 0
        B = sparse_zero(A.nrow, A.ncol, p=p)
    else
        B = deepcopy(A)
    end
    for k in 1:A.ncol
        for m in 1:length(A.columns[k])
            val = A.entries[k][m]
            B.entries[k][m] = scalar_multiply(val, alpha, p)
        end
    end
    return B
end

"""
    sparse_multiply(alpha::Rational{Int}, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for the rational numbers. An exception is raised
if the sparse matrix is not defined over field of rationals.
"""
function sparse_multiply(alpha::Rational{Int}, A::SparseMatrix)
    #
    # Multiply a sparse matrix by a scalar
    #
    p = A.char
    @assert p == 0 "Wrong scalar type!"

    if alpha == 0
        B = sparse_zero(A.nrow, A.ncol, p=p)
    else
        B = deepcopy(A)
    end
    for k in 1:A.ncol
        for m in 1:length(A.columns[k])
            B.entries[k][m] = alpha * A.entries[k][m]
        end
    end
    return B
end

"""
    Base.:*(A::SparseMatrix, B::SparseMatrix)

Multiply two sparse matrices.

Exceptions are raised if the matrix product is not defined 
or the entry types do not match.
"""
Base.:*(A::SparseMatrix,B::SparseMatrix) = sparse_multiply(A,B)

"""
    Base.:*(alpha::Int, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for finite fields. An exception is raised if
the sparse matrix is not defined over a finite field.
"""
Base.:*(alpha::Int, A::SparseMatrix) = sparse_multiply(alpha,A)

"""
    Base.:*(alpha::Rational{Int}, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for the rational numbers. An exception is raised
if the sparse matrix is not defined over field of rationals.
"""
Base.:*(alpha::Rational{Int}, A::SparseMatrix) = sparse_multiply(alpha,A)

