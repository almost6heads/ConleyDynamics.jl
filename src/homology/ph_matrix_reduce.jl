export ph_matrix_reduce

"""
    ph_matrix_reduce(lc::LefschetzComplex, filtration::Vector{Int})

Compute the persistent homology of a Lefschetz complex filtration over
the field associated with the Lefschetz complex boundary matrix.

The function returns the two values
* `phsingles::Vector{Vector{Int}}`
* `phpairs::Vector{Vector{Tuple{Int,Int}}}`
It assumes that the order given by the filtration values is admissible,
i.e., the permuted boundary matrix is strictly upper triangular. The
function returns the starting filtration values for infinite length
persistence intervals in `phsingles`, and the birth- and death-filtration
values for finite length persistence intervals in `phpairs`.
"""
function ph_matrix_reduce(lc::LefschetzComplex, filtration::Vector{Int})
    #
    # Compute the persistent homology of a Lefschetz complex filtration
    #

    # Find an admissible permutation of the cells

    fvals = sort(unique(filtration))
    adperm = Vector{Int}([])
    for kv in fvals
        append!(adperm, findall(t -> t==kv, filtration))
    end

    # Create the permuted boundary matrix and make sure it
    # is strictly upper triangular

    bndperm = sparse_permute(lc.boundary, adperm, adperm)

    if !sparse_is_sut(bndperm)
        error("Filtration error!")
        return
    end

    # Perform the persistence algorithm

    permsingles, permpairs = ph_reduce!(bndperm)

    # Extract the correct persistence intervals based on
    # the original order

    phsingles = [Vector{Int}() for _ in 0:lc.dim]
    phpairs = [Vector{Tuple{Int,Int}}() for _ in 0:lc.dim]

    for k=1:length(permsingles)
        singleindex  = adperm[permsingles[k]]
        singlefilter = filtration[singleindex]
        singledim    = lc.dimensions[singleindex]
        push!(phsingles[1+singledim],singlefilter)
    end

    for k=1:length(permpairs)
        pairindex1  = adperm[permpairs[k][1]]
        pairindex2  = adperm[permpairs[k][2]]
        pairfilter1 = filtration[pairindex1]
        pairfilter2 = filtration[pairindex2]
        if !(pairfilter1 == pairfilter2)
            pairdim = lc.dimensions[pairindex1]
            push!(phpairs[1+pairdim],(pairfilter1,pairfilter2))
        end
    end

    # Return the results

    return phsingles, phpairs
end

## """
##     ph_reduce!(matrix::SparseMatrix; [returnbasis=true])
## 
## Apply the persistence reduction algorithm to the matrix.
## 
## The function returns the values
## * `phsingles::Vector{Vector{Int}}`
## * `phpairs::Vector{Vector{Tuple{Int,Int}}}`
## * `basis::SparseMatrix` (if `returnbasis=true`)
## It assumes that `matrix` is strictly upper triangular. The
## function returns the starting columns for infinite length
## persistence intervals in `phsingles`, and the birth- and
## death-columns for finite length persistence intervals in
## `phpairs`. If the optional argument `returnbasis=true` is
## given, then the function also returns the computed basis
## matrix B with `reduced = matrix * B`.
## """
function ph_reduce!(matrix::SparseMatrix; returnbasis::Bool=false)
    #
    # Apply the persistence reduction algorithm to the matrix
    #

    if !sparse_is_sut(matrix)
        error("The matrix has to be strictly upper triangular!")
        return
    end

    # Extract the zero and one elements

    tzero = matrix.zero
    tone  = matrix.one

    # Create the identity matrix for basis computation

    numcolumns = sparse_size(matrix, 2)
    if returnbasis
        basis = sparse_identity(numcolumns, p=matrix.char)
    end

    # Initialize the main computation

    lowtocolumn = fill(Int(0),numcolumns)
    partofinterval = fill(false,numcolumns)

    for j=1:numcolumns
        keepgoing = true
        while (length(matrix.columns[j]) > 0) & keepgoing
            columnlow = sparse_low(matrix,j)
            if lowtocolumn[columnlow] == 0
                keepgoing = false
                lowtocolumn[columnlow] = j
                partofinterval[columnlow] = true
                partofinterval[j] = true
            else
                s = lowtocolumn[columnlow]
                gamma1 = matrix[columnlow,j]
                gamma2 = matrix[columnlow,s]
                sparse_add_column!(matrix,j,s,-gamma1,gamma2)
                if returnbasis
                    sparse_add_column!(basis,j,s,-gamma1,gamma2)
                end
            end
        end
    end

    # Prepare the return arrays

    phpairs = Vector{Tuple{Int,Int}}()
    phsingles = Vector{Int}()

    for k=1:numcolumns
        if lowtocolumn[k] > 0
            push!(phpairs,(k,lowtocolumn[k]))
        end
        if !partofinterval[k]
            push!(phsingles,k)
        end
    end

    # Return the results

    if returnbasis
        return phsingles, phpairs, basis
    else
        return phsingles, phpairs
    end
end

