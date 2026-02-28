export cm_reduce_dlms24!
export cm_reduce_dhl26!

"""
    cm_reduce_dlms24!(matrix::SparseMatrix, psetvec::Vector{Int};
                      [returnbasis::Bool],[returntm::Bool])

Compute the connection matrix.

This function uses the algorithm from the paper Dey, Lipinski,
Mrozek, and Slechta (SIAM Journal on Applied Dynamical Systems, 2024).
It assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
* `basisvecs` (optional): If the argument `returnbasis=true` is given,
  this returns information about the computed basis. The k-th entry
  of `basisvecs` is a vector containing the columns making up the
  k-th basis vector, which corresponds to column `cmatrix_cols[k]`.
* `tmatrix` (optional): If the argument `returntm=true` is given
  in addition to `returnbasis=true`, then instead of `basisvecs`
  the function returns the complete transformation matrix. In this
  case, `basicvecs` is not returned.
"""
function cm_reduce_dlms24!(matrix::SparseMatrix, psetvec::Vector{Int};
                           returnbasis::Bool=false, returntm::Bool=false)
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Create the identity matrix for basis computation

    numcolumns = sparse_size(matrix, 2)
    if returnbasis
        basis = sparse_identity(numcolumns, p=tchar)
    end

    # Initialize the main computation

    for j = 1:numcolumns
        jlow = sparse_low(matrix,j)
        if jlow > 0
            for i = jlow:-1:1
                if !(matrix[i,j] == tzero)
                    s = 1
                    found_s = false
                    while (!found_s) & (s <= numcolumns)
                        if ((!(s == j)) & (sparse_low(matrix,s) == i) & is_homogeneous(matrix,psetvec,s))
                            found_s = true
                        else
                            s += 1
                        end
                    end

                    if found_s
                        gamma1 = matrix[i,j]
                        gamma2 = matrix[i,s]
                        sparse_add_column!(matrix,j,s,-gamma1,gamma2)
                        if returnbasis
                            sparse_add_column!(basis,j,s,-gamma1,gamma2)
                        end
                        sparse_add_row!(matrix,s,j,gamma1,gamma2)
                    end
                end
            end
        end
    end

    cmatrix_cols = cm_columns(matrix, psetvec)
    cmatrix      = sparse_minor(matrix, cmatrix_cols, cmatrix_cols)
    
    if returnbasis
        if returntm
            return cmatrix, cmatrix_cols, basis
        else
            # Create and return a vector of basis vectors
            basisvecs = Vector{Vector{Int}}()
            for k=1:length(cmatrix_cols)
                bvec = sort(sparse_get_nz_column(basis,cmatrix_cols[k]),rev=true)
                push!(basisvecs,bvec)
            end
            return cmatrix, cmatrix_cols, basisvecs
        end
    else
        return cmatrix, cmatrix_cols
    end
end

"""
    cm_reduce_dhl26!(matrix::SparseMatrix, psetvec::Vector{Int})

Compute the connection matrix.

This function uses the algorithm from the paper Dey, Haas, and
Lipinski (SIAM Journal on Applied Dynamical Systems, 2026).
Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
"""
function cm_reduce_dhl26!(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Initialize the main computation
    
    numcolumns = sparse_size(matrix, 2)

    for j = 1:numcolumns
        jlow = sparse_low(matrix,j)
        if jlow > 0
            for i = jlow:-1:1
                if !(matrix[i,j] == tzero)
                    found_s = false
                    if length(matrix.rows[i]) > 0
                        s = 0  # Needed to make it accessible outside for loop
                        for sj in matrix.rows[i]
                            if ((sj < j) && (sparse_low(matrix,sj) == i) && is_homogeneous(matrix,psetvec,sj))
                                found_s = true
                                s = sj
                                break
                            end
                        end
                    end

                    if found_s
                        gamma1 = matrix[i,j]
                        gamma2 = matrix[i,s]
                        sparse_add_column!(matrix,j,s,-gamma1,gamma2)
                    end
                end
            end
        end
    end

#    #
#    # The following uses standard sparse matrix interfaces, which results in
#    # more memory allocation, but the speed seems to be the same. For now I
#    # keep the above version, as it is using less memory.
#    #
#    for j = 1:numcolumns
#        jlow = sparse_low(matrix,j)
#        if jlow > 0
#            for i = jlow:-1:1
#                if !(matrix[i,j] == tzero)
#                    nzrow   = sparse_get_nz_row(matrix, i)
#                    nzrowj  = nzrow[nzrow.<j]
#                    found_s = false
#
#                    if length(nzrowj) > 0
#                        s = 0  # Needed to make it accessible outside for loop
#                        for sj in nzrowj
#                            if ((sparse_low(matrix,sj) == i) && is_homogeneous(matrix,psetvec,sj))
#                                found_s = true
#                                s = sj
#                                break
#                            end
#                        end
#                    end
#
#                    if found_s
#                        gamma1 = matrix[i,j]
#                        gamma2 = matrix[i,s]
#                        sparse_add_column!(matrix,j,s,-gamma1,gamma2)
#                    end
#                end
#            end
#        end
#    end

    cmatrix_cols = cm_columns(matrix, psetvec)
    cmatrix      = sparse_minor(matrix, cmatrix_cols, cmatrix_cols)
    return cmatrix, cmatrix_cols
end

###############################################################################

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

###############################################################################

function cm_columns(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Create a vector of column indices for the connection matrix
    #
    numcolumns = sparse_size(matrix, 2)
    hcols      = homogeneous_columns(matrix, psetvec)
    tcols      = target_columns(matrix, psetvec)
    ccols      = Vector{Int}()

    for j = 1:numcolumns
        if (hcols[j] == false) & (tcols[j] == false)
            push!(ccols,j)
        end
    end

    return ccols
end

###############################################################################

function homogeneous_columns(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Determine which columns are homogenous columns
    #
    numcolumns = sparse_size(matrix, 2)
    hcols      = Vector(zeros(Bool,numcolumns))

    for j = 1:numcolumns
        hcols[j] = is_homogeneous(matrix, psetvec, j)
    end

    return hcols
end

###############################################################################

function is_homogeneous(matrix::SparseMatrix, psetvec::Vector{Int},
                        cindex::Int)
    #
    # Decide whether a column is homogeneous.
    #

    lowc = sparse_low(matrix,cindex)

    if lowc == 0
        return false
    else
        targetindex = lowc
        if psetvec[cindex] == psetvec[targetindex]
            return true
        else
            return false
        end
    end
end

###############################################################################

function target_columns(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Determine which columns are target columns
    #
    numcolumns = sparse_size(matrix, 2)
    tcols      = Vector(zeros(Bool,numcolumns))

    for j = 1:numcolumns
        if is_homogeneous(matrix, psetvec, j)
            tcols[sparse_low(matrix,j)] = true
        end
    end

    return tcols
end

###############################################################################

