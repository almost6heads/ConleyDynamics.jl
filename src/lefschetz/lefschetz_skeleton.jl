export lefschetz_skeleton

"""
    lefschetz_skeleton(lc::AbstractComplex, subcomp::Vector{Int}, skdim::Int)

Compute the `skdim`-dimensional skeleton of a Lefschetz complex subset.

The computed skeleton is for the closure of the subcomplex given by `subcomp`.
"""
function lefschetz_skeleton(lc::AbstractComplex, subcomp::Vector{Int}, skdim::Int)
    #
    # Compute the skdim-dimensional skeleton of a Lefschetz complex subset
    #

    # Compute the closure of the subcomplex
    
    lcclosure = lefschetz_closure(lc, subcomp)

    # Extract the cells of dimension skdim

    lcskeleton = Vector{Int}()

    for k in lcclosure
        if lc.dimensions[k] == skdim
            push!(lcskeleton,k)
        end
    end

    # Return the skeleton indices

    return lcskeleton
end

"""
    lefschetz_skeleton(lc::AbstractComplex, subcomp::Vector{String}, skdim::Int)

Compute the `skdim`-dimensional skeleton of a Lefschetz complex subset.

The computed skeleton is for the closure of the subcomplex given by `subcomp`.
"""
function lefschetz_skeleton(lc::AbstractComplex, subcomp::Vector{String}, skdim::Int)
    #
    # Compute the skdim-dimensional skeleton  of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcskeletonI = lefschetz_skeleton(lc, subcompI, skdim)
    lcskeleton  = lc.labels[lcskeletonI]
    return lcskeleton
end

"""
    lefschetz_skeleton(lc::AbstractComplex, skdim::Int)

Compute the `skdim`-dimensional skeleton of a Lefschetz complex.

The computed skeleton is for the full Lefschetz complex.
"""
function lefschetz_skeleton(lc::AbstractComplex, skdim::Int)
    #
    # Compute the skdim-dimensional skeleton of a Lefschetz complex
    #

    # Extract the cells of dimension skdim

    lcskeleton = Vector{Int}()

    for k in 1:lc.ncells
        if lc.dimensions[k] == skdim
            push!(lcskeleton,k)
        end
    end

    # Return the skeleton indices

    return lcskeleton
end

