export cm_create! 

"""
    cmatrix, cmatrix_cols = cm_create!(matrix, psetvec)
    cmatrix, cmatrix_cols, basisvecs = cm_create!(matrix, psetvec;
                                                  returnbasis=true)

Compute the connection matrix.

Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`. If the optional
argument `returnbasis=true` is given, then the function also
returns information about the computed basis. The k-th entry of
`basisvecs` is a vector containing the columns making up the
k-th basis vector, which corresponds to column `cmatrix_cols[k]`.
"""
function cm_create!(matrix::Matrix, psetvec::Vector{Int};
                    returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Create the identity matrix for basis computation

    numcolumns = size(matrix, 2)
    if returnbasis
        basis = deepcopy(matrix)
        basis = 0 * basis
        for i=1:numcolumns
            basis[i,i] = basis[i,i] + 1
        end
    end

    # Initialize the main computation

    lowvec = zeros(Int,numcolumns)
    update_low!(matrix, lowvec, startindex=1)

    for j = 1:numcolumns
        jlow = lowvec[j]
        if jlow > 0
            for i = jlow:-1:1
                if !(matrix[i,j] == 0)
                    s = 1
                    found_s = false
                    while (!found_s) & (s <= numcolumns)
                        if ((!(s == j)) & (lowvec[s] == i) & is_homogeneous(matrix, lowvec, psetvec, s))
                            found_s = true
                        else
                            s += 1
                        end
                    end

                    if found_s
                        gamma = matrix[i,j] / matrix[i,s]
                        @views matrix[:,j] = matrix[:,j] .- gamma .* matrix[:,s]
                        if returnbasis
                            @views basis[:,j] = basis[:,j] .- gamma .* basis[:,s]
                        end
                        @views matrix[s,:] = matrix[s,:] .+ gamma .* matrix[j,:]
                        update_low!(matrix, lowvec, startindex=j)
                    end
                end
            end
        end
    end

    cmatrix_cols = cm_columns(matrix, lowvec, psetvec)
    cmatrix      = matrix[cmatrix_cols, cmatrix_cols]
    
    if returnbasis

        # Create a vector of basis vectors

        basisvecs = Vector{Vector{Int}}()

        for k=1:length(cmatrix_cols)
            bvec = Vector{Int}()
            for m = size(basis,2):-1:1
                if !(basis[m,cmatrix_cols[k]]==0)
                    push!(bvec,m)
                end
            end
            push!(basisvecs,bvec)
        end

        # Return the results

        return cmatrix, cmatrix_cols, basisvecs
    else
        return cmatrix, cmatrix_cols
    end
end

"""
    cmatrix, cmatrix_cols = cm_create!(matrix, psetvec)
    cmatrix, cmatrix_cols, basisvecs = cm_create!(matrix, psetvec;
                                                  returnbasis=true)

Compute the connection matrix.

Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`. If the optional
argument `returnbasis=true` is given, then the function also
returns information about the computed basis. The k-th entry of
`basisvecs` is a vector containing the columns making up the
k-th basis vector, which corresponds to column `cmatrix_cols[k]`.
"""
function cm_create!(matrix::SparseMatrix, psetvec::Vector{Int};
                    returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tzero = matrix.zero
    tone  = matrix.one

    # Create the identity matrix for basis computation

    numcolumns = sparse_size(matrix, 2)
    if returnbasis
        basis = sparse_identity(numcolumns, tone)
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
                        gamma = matrix[i,j] / matrix[i,s]
                        sparse_add_column!(matrix,j,s,-gamma)
                        if returnbasis
                            sparse_add_column!(basis,j,s,-gamma)
                        end
                        sparse_add_row!(matrix,s,j,gamma)
                    end
                end
            end
        end
    end

    cmatrix_cols = cm_columns(matrix, psetvec)
    cmatrix      = sparse_minor(matrix, cmatrix_cols, cmatrix_cols)
    
    if returnbasis

        # Create a vector of basis vectors

        basisvecs = Vector{Vector{Int}}()

        for k=1:length(cmatrix_cols)
            bvec = sort(sparse_get_nz_column(basis,cmatrix_cols[k]),rev=true)
            push!(basisvecs,bvec)
        end

        # Return the results

        return cmatrix, cmatrix_cols, basisvecs
    else
        return cmatrix, cmatrix_cols
    end
end
