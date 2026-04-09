export sparse_basis_kernel
export sparse_basis_range

"""
    sparse_basis_kernel(A::SparseMatrix)

Compute a kernel basis for a sparse matrix.

This function determines a basis for the kernel of the
given matrix `A`. The result is returned as a vector
of sparse matrices. The latter are column vectors which
span the kernel.

# Example
```jldoctest
julia> Afull = [1 0 1 2 4 3 2; 3 4 6 2 3 4 2; 2 0 2 4 1 1 1];

julia> A = sparse_from_full(Afull, p=7);

julia> sparse_rref!(A)

julia> sparse_show(A)
 1 . 1 2 4 . 3
 . 1 6 6 3 . 5
 . . . . . 1 2

julia> kbasis = sparse_basis_kernel(A);

julia> length(kbasis)
4

julia> sparse_show(sparse_transpose(kbasis[1]))
 6 1 1 . . . .

julia> sparse_show(sparse_transpose(kbasis[2]))
 5 1 . 1 . . .

julia> sparse_show(sparse_transpose(kbasis[3]))
 3 4 . . 1 . .

julia> sparse_show(sparse_transpose(kbasis[4]))
 4 2 . . . 5 1
```
"""
function sparse_basis_kernel(A::SparseMatrix)
    #
    # Compute a kernel basis for a sparse matrix
    #

    basisvecs = Vector{SparseMatrix}()

    # Deal with the trivial case first

    if sparse_is_zero(A)
        n = A.ncol
        for k = 1:n
            bvec = sparse_zero(n, 1, p=A.char)
            bvec[k,1] = A.one
            push!(basisvecs, bvec)
        end
        return basisvecs
    end
    
    # Check whether or not the matrix is already in rref

    if sparse_is_rref(A)
        AA = A
    else
        AA = sparse_rref(A)
    end

    # Find the number of pivots and the pivot columns

    n = AA.ncol
    npivot = length(findall(t -> t > 0, length.(AA.rows)))
    pivotcols = [AA.rows[k][1] for k = 1:npivot]
    
    # Find the columns corresponding to kernel basis vectors

    kernelcols = setdiff(collect(1:n), pivotcols)

    # Create the basis vectors

    for kc in kernelcols
        bvec = sparse_zero(n, 1, p=AA.char)
        bvec[kc,1] = AA.one
        for j = 1:npivot
            bvec[pivotcols[j],1] = -AA[j,kc]
        end
        push!(basisvecs, bvec)
    end

    return basisvecs
end

"""
    sparse_basis_range(A::SparseMatrix)

Compute a basis for the range of a sparse matrix.

This function determines a basis for the range of the
given matrix `A`. The result is returned as a vector
of sparse matrices. The latter are column vectors which
span the range, also known as column space.

# Example
```jldoctest
julia> Afull = [1 0 1 2 4 3 2; 3 4 6 2 3 4 2; 2 0 2 4 1 1 1];

julia> A = sparse_from_full(Afull, p=7);

julia> sparse_show(A)
 1 . 1 2 4 3 2
 3 4 6 2 3 4 2
 2 . 2 4 1 1 1

julia> rbasis = sparse_basis_range(A);

julia> length(rbasis)
3

julia> sparse_show(sparse_transpose(rbasis[1]))
 1 3 2

julia> sparse_show(sparse_transpose(rbasis[2]))
 . 4 .

julia> sparse_show(sparse_transpose(rbasis[3]))
 3 4 1

julia> sparse_show(sparse_rref(A))
 1 . 1 2 4 . 3
 . 1 6 6 3 . 5
 . . . . . 1 2
```
"""
function sparse_basis_range(A::SparseMatrix)
    #
    # Compute a basis for the range of a sparse matrix
    #

    basisvecs = Vector{SparseMatrix}()

    # Deal with the trivial case first

    if sparse_is_zero(A)
        return basisvecs
    end
    
    # Check whether or not the matrix is already in rref

    if sparse_is_rref(A)
        AA = A
    else
        AA = sparse_rref(A)
    end

    # Find the number of pivots and the pivot columns

    n = AA.nrow
    npivot = length(findall(t -> t > 0, length.(AA.rows)))
    pivotcols = [AA.rows[k][1] for k = 1:npivot]

    # Create the basis vectors

    for kc in pivotcols
        bvec = sparse_zero(n, 1, p=AA.char)
        for j = 1:n
            bvec[j,1] = A[j,kc]
        end
        push!(basisvecs, bvec)
    end

    return basisvecs
end

