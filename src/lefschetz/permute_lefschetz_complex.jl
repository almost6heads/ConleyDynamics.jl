export permute_lefschetz_complex

"""
    permute_lefschetz_complex(lc::AbstractComplex,
                              permutation::Vector{Int})

Permute the indices of a Lefschetz complex.

The vector `permutation` contains a permutation of the indices for the
given Lefschetz complex `lc`. If no permutation is specified, or if the
length of the vector is not correct, then a randomly generated one will
be used. Note that the permutation has to respect the ordering of the
cells by dimension, otherwise an error is raised. In other words, 
the permutation has to decompose into permutations within each dimension.
This is automatically done if no permutation is explicitly specified.
"""
function permute_lefschetz_complex(lc::AbstractComplex,
                                   perm::Vector{Int}=Vector{Int}([]))
    #
    # Create a Lefschetz complex which permutes the cell indices
    #
    
    # Create a random permutation if none is set

    if !(length(perm) == lc.ncells)
        perm = Vector{Int}([])
        for k = 0:lc.dim
            crange = findall(x -> x==k, lc.dimensions)
            if length(crange) > 0
                crmin = minimum(crange)
                crmax = maximum(crange)
                crperm = randperm(1 + crmax - crmin) .+ (crmin - 1)
                append!(perm, crperm)
            end
        end
    end

    # Copy the data from the original complex

    ncells1     = lc.ncells
    lcdim1      = lc.dim
    boundary1   = lc.boundary
    labels1     = lc.labels
    dimensions1 = lc.dimensions

    # Create the permuted data

    labels2     = labels1[perm]
    dimensions2 = dimensions1[perm]
    boundary2   = sparse_permute(boundary1, perm, perm)

    # Create the permuted complex (preserving EuclideanComplex if applicable)

    if lc isa EuclideanComplex
        coords2 = lc.coords[perm]
        lc2 = EuclideanComplex(labels2, dimensions2, boundary2, coords2; validate=false)
    else
        lc2 = LefschetzComplex(labels2, dimensions2, boundary2; validate=false)
    end

    # Return the permuted complex

    return lc2
end

