#
# Internal helper functions for multivector fields
#

function cell_to_sets(ncells::Int, sets::Vector{Vector{Int}})
    #
    # For every cell, list the indices of the sets which
    # contain it, in increasing order and without repetitions
    #
    c2s = [Vector{Int}() for _ in 1:ncells]
    for k in 1:length(sets)
        for c in sets[k]
            if isempty(c2s[c]) || (last(c2s[c]) != k)
                push!(c2s[c], k)
            end
        end
    end
    return c2s
end

function sets_containing(c2s::Vector{Vector{Int}}, cells::Vector{Int}, nsets::Int)
    #
    # Return the indices of all sets which contain at least one of
    # the given cells, in increasing order. Here `c2s` is the output
    # of `cell_to_sets`, and cells outside the complex are ignored.
    #
    isset = fill(false, nsets)
    for c in cells
        if 1 <= c <= length(c2s)
            for k in c2s[c]
                isset[k] = true
            end
        end
    end
    return findall(isset)
end

function euclidean_vertex_coords(ec::EuclideanComplex)
    #
    # Return a coordinate vector indexed by cell, in which the entries
    # of the vertices contain their coordinates, and all other entries
    # are empty, as expected by create_planar_mvf and create_spatial_mvf
    #
    coords = [Vector{Float64}() for _ in 1:ec.ncells]
    for k in 1:ec.ncells
        if ec.dimensions[k] == 0
            coords[k] = ec.coords[k][1]
        end
    end
    return coords
end
