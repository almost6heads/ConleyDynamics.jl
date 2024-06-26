export persistent_homology

"""
    persistent_homology(lc::LefschetzComplex, filtration::Vector{Int})

Complete the persistent homology of a Lefschetz complex filtration over
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
function persistent_homology(lc::LefschetzComplex, filtration::Vector{Int})
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

