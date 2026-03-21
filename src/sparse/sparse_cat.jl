export sparse_hcat
export sparse_vcat
export sparse_hvcat
export sparse_hvncat
export sparse_cat

################################
#                              #
#   Horizontal concatenation   #
#                              #
################################

"""
    sparse_hcat(A::SparseMatrix, B::SparseMatrix)

Concatenate two sparse matrices horizontally.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of rows.
"""
function sparse_hcat(A::SparseMatrix, B::SparseMatrix)
    #
    # Concatenate two sparse matrices horizontally
    #

    # Check whether the concatenation makes sense

    @assert A.nrow == B.nrow "The matrices need to have the same number of rows!"
    @assert A.char == B.char "The matrices have to be over the same field!"

    # Perform the matrix concatenation

    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(A)
    for k=1:B.ncol
        if length(B.columns[k]) > 0
            for m in 1:length(B.columns[k])
                push!(c, k+nc)
                push!(r, B.columns[k][m])
                push!(v, B.entries[k][m])
            end
        end
    end
    nc = nc + B.ncol

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, v)
    return sm
end

"""
    sparse_hcat(S1, S2, S3, ...)

Concatenate several sparse matrices horizontally.

This method of `sparse_hcat` allows for any number of 
arguments. If you want to concatenate a vector `vec`
of sparse matrices, use the call `sparse_hcat(vec...)`.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of rows.
"""
sparse_hcat(vargs::SparseMatrix...) = reduce(sparse_hcat, vargs)

"""
    Base.hcat(::SparseMatrix...)

Concatenate several sparse matrices horizontally.

This method allows one to use the shortform `[S1 S2 S3...]`
for the horizontal concatenation, just like in the matrix case.
"""
Base.hcat(smvec::SparseMatrix...) = sparse_hcat(smvec...)

##############################
#                            #
#   Vertical concatenation   #
#                            #
##############################

"""
    sparse_vcat(A::SparseMatrix, B::SparseMatrix)

Concatenate two sparse matrices vertically.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of columns.
"""
function sparse_vcat(A::SparseMatrix, B::SparseMatrix)
    #
    # Concatenate two sparse matrices vertically
    #

    # Check whether the concatenation makes sense

    @assert A.ncol == B.ncol "The matrices need to have the same number of columns!"
    @assert A.char == B.char "The matrices have to be over the same field!"

    # Perform the matrix concatenation

    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(A)
    for k=1:B.ncol
        if length(B.columns[k]) > 0
            for m in 1:length(B.columns[k])
                push!(c, k)
                push!(r, B.columns[k][m] + nr)
                push!(v, B.entries[k][m])
            end
        end
    end
    nr = nr + B.nrow

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, v)
    return sm
end

"""
    sparse_vcat(S1, S2, S3, ...)

Concatenate several sparse matrices vertically.

This method of `sparse_vcat` allows for any number of 
arguments. If you want to concatenate a vector `vec`
of sparse matrices, use the call `sparse_vcat(vec...)`.

Exceptions are raised if the matrices are not defined over
the same field, or if they have different numbers of columns.
"""
sparse_vcat(vargs::SparseMatrix...) = reduce(sparse_vcat, vargs)

"""
    Base.vcat(::SparseMatrix...)

Concatenate several sparse matrices vertically.

This method allows one to use the shortform `[S1; S2; S3; ...]`
for the vertical concatenation, just like in the matrix case.
"""
Base.vcat(smvec::SparseMatrix...) = sparse_vcat(smvec...)

#############################################
#                                           #
#   Horizontal and vertical concatenation   #
#                                           #
#############################################

"""
    sparse_hvcat(rblocks::Tuple{Vararg{Int}}, A::SparseMatrix...)

Arrange sparse matrices in rectangular form.

The fist argument has to be a tuple which specifies how many blocks
have to be combined in each row. This means that the total number
of sparse matrices has to be equal to the sum of the entries in
`rblocks`. Exceptions are raised if the dimensions of the matrices
are not compatible, or if the matrices are not defined over the
same field.
"""
function sparse_hvcat(b::Tuple{Vararg{Int}}, Avec::SparseMatrix...)
    #
    # Concatenate sparse matrices in rectangular form
    #

    # Perform some sanity checks
    
    @assert minimum(b) >= 1 "We need at least one block in each row!"
    @assert reduce(+,b) == length(Avec) "Wrong number of matrices!"

    # Perform the horizontal concatenations

    rows = Vector{SparseMatrix}()
    offset = 0
    for k = 1:length(b)
        push!(rows, sparse_hcat(Avec[offset+1:offset+b[k]]...))
        offset = offset + b[k]
    end

    # Perform the vertical concatention and return the result

    return sparse_vcat(rows...)
end

"""
    Base.hvcat(::Tuple{Vararg{Int}}, ::SparseMatrix...)

Arrange sparse matrices in rectangular form.

This method allows one to use the shortform `[S1 S2; S3 S4 S5]`
etc for the concatenation, just like in the matrix case. It
is based on the function `sparse_hvcat`. For more information,
refer to its documentation. The number of blocks in each row
can vary, as long as the dimensions are overall compatible.
"""
Base.hvcat(b::Tuple{Vararg{Int}}, s::SparseMatrix...) = sparse_hvcat(b, s...)

"""
    sparse_hvcat(rblocks::Int, A::SparseMatrix...)

Arrange sparse matrices in rectangular form.

The fist argument specifies how many blocks are contained
in each row. Thus, the number of sparse matrices in `Avec`
has to be divisible by `rblocks`. Exceptions are raised if
the dimensions of the matrices are not compatible, or if
the matrices are not defined over the same field.
"""
function sparse_hvcat(bpr::Int, Avec::SparseMatrix...)
    #
    # Concatenate sparse matrices in rectangular form
    #

    # Reduce to the tuple method

    nrow = div(length(Avec), bpr)
    b = Tuple(fill(bpr, nrow))
    
    return sparse_hvcat(b, Avec...)
end

"""
    Base.hvcat(::Int, ::SparseMatrix...)

Arrange sparse matrices in rectangular form.

This method allows one to use the shortform `[S1 S2; S3 S4]`
etc for the concatenation, just like in the matrix case. It
is based on the function `sparse_hvcat`. For more information,
refer to its documentation.
"""
Base.hvcat(b::Int, s::SparseMatrix...) = sparse_hvcat(b, s...)

"""
    sparse_hvncat(dims::Tuple{Int,Int}, rowfirst::Bool, A::SparseMatrix...)

Arrange sparse matrices in rectangular form.

The function expects dimensions `(d1,d2)`, as well as `d1*d2`
sparse matrices. The boolean argument `rowfirst` indicates whether
the sparse matrices are arranged row-first or column-first. The
matrices are then arranged in an `d1 x d2`-times rectangular array
to create a new sparse matrix. Exceptions are raised if the
dimensions of the matrices are not compatible, or if the
matrices are not defined over the same field.
"""
function sparse_hvncat(dims::Tuple{Int,Int}, rowfirst::Bool, Avec::SparseMatrix...)
    #
    # Arrange sparse matrices in rectangular form
    #

    # Perform some sanity checks

    d1, d2 = dims
    @assert d1 >= 1 "We need at least one row!"
    @assert d2 >= 1 "We need at least one column!"
    @assert d1*d2 == length(Avec) "Wrong number of matrices!"

    # Consider the row-first case

    if rowfirst
        rows = Vector{SparseMatrix}()
        for k = 0:d1-1   # Loop through the rows
            push!(rows, sparse_hcat(Avec[k*d2+1:k*d2+d2]...))
        end
        return sparse_vcat(rows...)
    else
        columns = Vector{SparseMatrix}()
        for m = 0:d2-1   # Loop through the columns
            push!(columns, sparse_vcat(Avec[m*d1+1:m*d1+d1]...))
        end
        return sparse_hcat(columns...)
    end
end

"""
    Base.hvncat(::Tuple{Int,Int}, ::Bool, ::SparseMatrix...)

Arrange sparse matrices in rectangular form.

This method allows one to use the shortform `[S1; S2;; S3; S4]`
etc for the concatenation, just like in the matrix case. It
is based on the function `sparse_hvncat`. For more information,
refer to its documentation. This form does require that the same
number of blocks are used in each row and each column. For varying
numbers please use the functions `hvcat` or `sparse_hvcat`.
"""
Base.hvncat(d::Tuple{Int,Int}, r::Bool, s::SparseMatrix...) = sparse_hvncat(d, r, s...)

"""
    sparse_cat(A::SparseMatrix...; dims=1)

Concatenate sparse matrices along a dimension.

The integer `dims` specifies along which dimension the sparse
matrices are combined. In other words, `dims=1` reduces to `vcat`,
while `dims=2` reduces to `hcat`. In addition, this function
accepts the argument `dims=(1,2)`. In this case, the sparse
matrices are placed as blocks along the diagonal of a large
sparse block matrix, with all remaining blocks being zero.
"""
function sparse_cat(Avec::SparseMatrix...; dims=1)
    #
    # Concatenate sparse matrices along a dimension
    #

    # Deal with the integer arguments

    if dims == 1
        return sparse_vcat(Avec...)
    end

    if dims == 2
        return sparse_hcat(Avec...)
    end

    # Now consider the case dims=(1,2)

    if dims == (1,2)
        nr = [Avec[k].nrow for k in 1:length(Avec)]
        nc = [Avec[k].ncol for k in 1:length(Avec)]
        np = [Avec[k].char for k in 1:length(Avec)]
        rows = Vector{SparseMatrix}()

        for k = 1:length(Avec)
            currentrow = Vector{SparseMatrix}()
            for m = 1:length(Avec)
                if m == k
                    push!(currentrow, Avec[k])
                else
                    push!(currentrow, sparse_zero(nr[k], nc[m], p=np[k]))
                end
            end
            push!(rows, sparse_hcat(currentrow...))
        end
        return sparse_vcat(rows...)
    end

    error("dims needs to be 1, 2, or (1,2)!")
end

