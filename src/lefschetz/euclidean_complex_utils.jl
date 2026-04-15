export euclidean_to_lefschetz, lefschetz_to_euclidean, rescale_coords

"""
    euclidean_to_lefschetz(ec::EuclideanComplex) -> LefschetzComplex

Convert a `EuclideanComplex` to a `LefschetzComplex` by stripping the
embedded coordinates. All other fields are shared unchanged.
"""
function euclidean_to_lefschetz(ec::EuclideanComplex)
    #
    # Strip the coords field and return a plain LefschetzComplex
    #
    return LefschetzComplex(ec.labels, ec.dimensions, ec.boundary; validate=false)
end

"""
    lefschetz_to_euclidean(lc::LefschetzComplex,
                           vertex_coords::Vector{<:Vector{<:Real}}) -> EuclideanComplex

Embed a `LefschetzComplex` into Euclidean space by attaching coordinates
to every cell.

The vector `vertex_coords` contains one coordinate vector per vertex of `lc`,
indexed by the cell index of the vertex (i.e., `vertex_coords[k]` gives the
coordinates of the vertex whose cell index is `k`). Since vertices are always
stored first in the cell list (dimensions are non-decreasing), the indices
1 through `nvertices` correspond to the vertices.

The function raises an `ArgumentError` if the complex is not simplicial or
cubical, i.e., if any 1-cell has a vertex count other than 2, or any 2-cell
has a vertex count other than 3 or 4.
"""
function lefschetz_to_euclidean(lc::LefschetzComplex,
                                vertex_coords::Vector{<:Vector{<:Real}})
    #
    # Build per-cell coords by finding each cell's vertices via lefschetz_skeleton
    #

    coords = Vector{Vector{Vector{Float64}}}(undef, lc.ncells)

    for k in 1:lc.ncells
        verts = lefschetz_skeleton(lc, [k], 0)
        n = length(verts)
        cdim = lc.dimensions[k]

        # Validate vertex counts for simplicial/cubical complexes
        if cdim == 1 && n != 2
            throw(ArgumentError(
                "lefschetz_to_euclidean: complex is not simplicial or cubical — " *
                "cell \"$(lc.labels[k])\" has $n vertices (expected 2 for edges)"))
        elseif cdim == 2 && n != 3 && n != 4
            throw(ArgumentError(
                "lefschetz_to_euclidean: complex is not simplicial or cubical — " *
                "cell \"$(lc.labels[k])\" has $n vertices (expected 3 for triangles, " *
                "4 for quads)"))
        end

        coords[k] = [Float64.(vertex_coords[v]) for v in verts]
    end

    return EuclideanComplex(lc.labels, lc.dimensions, lc.boundary, coords; validate=false)
end

"""
    lefschetz_to_euclidean(lc::LefschetzComplex,
                           coords::Vector{<:Vector{<:Vector{<:Real}}})
                           -> EuclideanComplex

Embed a `LefschetzComplex` into Euclidean space by attaching coordinates
to every cell.

The vector `coords` contains coordinates for each of the cells in `lc`,
indexed by the cell index. This function only performs basic check on the
suitability of this coordinate vector. More precisely, the function raises
an `ArgumentError` if any vertex has coordinate count other than 1, any 1-cell
has a coordinate count other than 2, or any 2-cell has a coordinate count
other than 3 or 4.
"""
function lefschetz_to_euclidean(lc::LefschetzComplex,
                                coords::Vector{<:Vector{<:Vector{<:Real}}})
    #
    # Convert a Lefschetz to a Euclidean complex by adding coords
    #

    @assert length(coords)==lc.ncells "I do need coordinates for every cell!"

    # Convert coordinates to Float64 and perform basic checks

    coordF = Vector{Vector{Vector{Float64}}}(undef, lc.ncells)
    embdim = length(coords[1][1])  # Dimension of the ambient space

    for k in 1:lc.ncells
        vertcoords = Vector{Vector{Float64}}()
        cdim = lc.dimensions[k]    # Dimension of the cell
        nc = length(coords[k])     # Number of provided coordinates

        # Validate vertex coordinate counts
        if cdim == 0 && nc != 1
            throw(ArgumentError(
                "lefschetz_to_euclidean: wrong number of vertex coordinates — " *
                "cell \"$(lc.labels[k])\" has $nc vertices (expected 1 for vertices)"))
        elseif cdim == 1 && nc != 2
            throw(ArgumentError(
                "lefschetz_to_euclidean: wrong number of edge coordinates — " *
                "cell \"$(lc.labels[k])\" has $nc vertices (expected 2 for edges)"))
        elseif cdim == 2 && nc != 3 && nc != 4
            throw(ArgumentError(
                "lefschetz_to_euclidean: wrong number of 2-cell coordinates — " *
                "cell \"$(lc.labels[k])\" has $nc vertices (expected 3 for triangles, " *
                "4 for quads)"))
        end

        # Make sure all povided vertex coordinates have the same length

        for m in 1:length(coords[k])
            if !(length(coords[k][m])==embdim)
                throw(ArgumentError(
                "lefschetz_to_euclidean: all vertices need to have the same number of coordinates"))
            end
            push!(vertcoords, Float64.(coords[k][m]))
        end

        coordF[k] = vertcoords
    end
            
    return EuclideanComplex(lc.labels, lc.dimensions, lc.boundary, coordF; validate=false)
end

"""
    rescale_coords(ec::EuclideanComplex,
                   a::Vector{<:Real}, b::Vector{<:Real}) -> EuclideanComplex

Rescale the coordinates of a `EuclideanComplex` so that component `j`
maps from its current range `[min_j, max_j]` to `[a[j], b[j]]`.

`a` and `b` must be vectors of length equal to the spatial dimension
(number of coordinate components). The rescaling is linear and applied
uniformly to all vertex coordinates across all cells.

# Example
```julia
# Normalize x to [0,1] and y to [-1,1]:
ec2 = rescale_coords(ec, [0.0, -1.0], [1.0, 1.0])
```
"""
function rescale_coords(ec::EuclideanComplex,
                        a::Vector{<:Real}, b::Vector{<:Real})
    #
    # Rescale coordinates per component from [min_j, max_j] to [a[j], b[j]]
    #

    # Determine spatial dimension from the first vertex cell
    d = length(a)
    if length(b) != d
        error("rescale_coords: vectors a and b must have the same length!")
    end

    # Find global min and max per component, scanning all vertex coords in all cells

    cmin = fill(Inf, d)
    cmax = fill(-Inf, d)

    for k in 1:ec.ncells
        for cvec in ec.coords[k]
            for j in 1:d
                if cvec[j] < cmin[j]
                    cmin[j] = cvec[j]
                end
                if cvec[j] > cmax[j]
                    cmax[j] = cvec[j]
                end
            end
        end
    end

    # Build the rescaled coords

    new_coords = Vector{Vector{Vector{Float64}}}(undef, ec.ncells)

    for k in 1:ec.ncells
        new_coords[k] = Vector{Vector{Float64}}(undef, length(ec.coords[k]))
        for (m, cvec) in enumerate(ec.coords[k])
            new_vec = Vector{Float64}(undef, d)
            for j in 1:d
                span = cmax[j] - cmin[j]
                if span == 0.0
                    new_vec[j] = (Float64(a[j]) + Float64(b[j])) / 2.0
                else
                    new_vec[j] = Float64(a[j]) +
                                 (cvec[j] - cmin[j]) / span *
                                 (Float64(b[j]) - Float64(a[j]))
                end
            end
            new_coords[k][m] = new_vec
        end
    end

    return EuclideanComplex(ec.labels, ec.dimensions, ec.boundary,
                            new_coords; validate=false)
end

"""
    rescale_coords(ec::EuclideanComplex, a::Real, b::Real) -> EuclideanComplex

Rescale all coordinate components uniformly to the interval `[a, b]`.

This is a convenience overload of [`rescale_coords`](@ref) that applies the
same target interval to every spatial dimension.
"""
function rescale_coords(ec::EuclideanComplex, a::Real, b::Real)
    #
    # Convenience: uniform rescaling to [a, b] for all components
    #

    # Determine spatial dimension from the first vertex entry
    d = 0
    for k in 1:ec.ncells
        if ec.dimensions[k] == 0
            d = length(ec.coords[k][1])
            break
        end
    end
    if d == 0
        error("rescale_coords: could not determine spatial dimension (no vertices found)")
    end

    return rescale_coords(ec, fill(Float64(a), d), fill(Float64(b), d))
end

