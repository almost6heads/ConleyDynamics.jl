export sparse_from_full, full_from_sparse

"""
    sparse_from_full(matrix::Matrix{Int}; [p::Int=2])

Create sparse matrix from full integer matrix. If the optional
argument `p` is specified and positive, then the returned matrix
is an integer matrix which is interpreted over `GF(p)`. On the
other hand, if `p` is equal to zero, then the return matrix has
rational type.
"""
function sparse_from_full(matrix::Matrix{Int}; p::Int=2)
    #
    # Create sparse matrix from full matrix

    # Initialize the row, column, and entry vectors

    nr = size(matrix,1)
    nc = size(matrix,2)

    if p == 0
        tzero = Rational{Int}(0)
        tone  = Rational{Int}(1)
    elseif p > 0
        tzero = Int(0)
        tone  = Int(1)
    else
        error("The characteristic p has to be 0 or positive!")
    end
    
    # Create the lists for the nonzero entries

    r = Vector{Int}([])
    c = Vector{Int}([])
    vals = Vector{typeof(tone)}([])

    for k=1:nr
        for m=1:nc
            if p == 0
                mvalue = Rational{Int}(matrix[k,m])
            else
                mvalue = mod(matrix[k,m], p)
            end
            if !(mvalue == tzero)
                push!(r,k)
                push!(c,m)
                push!(vals,mvalue)
            end
        end
    end

    # Create the sparse matrix and return it

    sm = sparse_from_lists(nr, nc, p, tzero, tone, r, c, vals)
    return sm
end

"""
    sparse_from_full(vec::Vector{Int}; [p::Int=2])

Create sparse column matrix from a full integer vector. If the
optional argument `p` is specified and positive, then the returned
matrix is an integer matrix which is interpreted over `GF(p)`. On
the other hand, if `p` is equal to zero, then the return matrix has
rational type.
"""
function sparse_from_full(vec::Vector{Int}; p::Int=2)
    #
    # Create sparse matrix from full vector
    #
    return sparse_from_full(hcat(vec), p=p)
end

"""
    full_from_sparse(sm::SparseMatrix)

Create full matrix from sparse matrix.
"""
function full_from_sparse(sm::SparseMatrix)
    #
    # Create full matrix from sparse matrix
    #

    # Initialize the full matrix with zeros

    fm = fill(sm.zero, (sm.nrow, sm.ncol))

    # Extract the nonzero elements

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(sm)

    # Set the nonzero elements

    for k=1:length(r)
        fm[r[k],c[k]] = vals[k]
    end

    # Return the full matrix
    
    return fm
end

