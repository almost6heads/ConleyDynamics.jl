export spatial_transitions

"""
    spatial_transitions(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}}, vf;
                        npts::Int=3)

Compute the sampled dynamical transitions of a vector field in three
dimensions.

The function expects a three-dimensional Lefschetz complex `lc`, a
coordinate vector `coords` of coordinates for all the 0-dimensional cells
in the complex, and a vector field `vf(z::Vector{Float64})::Vector{Float64}`
in three dimensions, as in `create_spatial_mvf`. For every vertex, edge,
and face `c` of `lc`, it determines the regions, i.e., the
three-dimensional cells, which the flow can enter from `c`. The result is a
vector `transitions` of integer vectors, where every entry has the form
`[c; targets]`. Cells for which no region is found, such as boundary cells
through which the flow leaves the domain, are omitted, and the regions
themselves are not listed.

The transitions are approximations of the immediate future of the cells,
i.e., of the sets `{c} ∪ Ift(c)`. They are computed from finitely many
samples of the vector field, and are therefore not rigorous. The keyword
`npts` gives the number of sample points in the interior of every edge,
and every face is sampled on a grid with spacing `1/(npts+1)`. Larger
values detect more tangencies, at a higher cost.

The smallest multivector field which contains all transitions is
`create_mvf_hull(lc, transitions)`, which is the result of
`create_spatial_mvf`. The transitions of the reversed flow, i.e., the
approximations of the immediate past, are obtained via
`spatial_transitions(lc, coords, x -> -vf(x))`, which uses the same
sample points.
"""
function spatial_transitions(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}}, vf;
                             npts::Int=3)
    #
    # Compute the sampled dynamical transitions of a spatial vector field
    #

    # Find the indices for the vertices, the edges, and the faces

    vi = findall(isequal(0), lc.dimensions)
    ei = findall(isequal(1), lc.dimensions)
    fi = findall(isequal(2), lc.dimensions)

    # Find the vertex to region map

    v2r = [Vector{Int}() for _ in 1:length(vi)]

    Threads.@threads for k = 1:length(vi)
        v2r[k] = spatial_mvf_vertex2region(vi[k], lc, coords, vf)
    end

    # Find the edge to region map

    e2r = [Vector{Int}() for _ in 1:length(ei)]

    Threads.@threads for k = 1:length(ei)
        e2r[k] = spatial_mvf_edge2region(ei[k], lc, coords, vf, npts=npts)
    end

    # Find the face to region map

    f2r = [Vector{Int}() for _ in 1:length(fi)]

    Threads.@threads for k = 1:length(fi)
        f2r[k] = spatial_mvf_face2region(fi[k], lc, coords, vf, npts=npts)
    end

    # Collect the transitions

    transitions = Vector{Vector{Int}}()

    for (cells, targets) in ((vi, v2r), (ei, e2r), (fi, f2r))
        for k in 1:length(cells)
            if length(targets[k]) > 0
                push!(transitions, vcat([cells[k]], targets[k]))
            end
        end
    end

    # Return the transitions

    return transitions
end

"""
    spatial_transitions(ec::EuclideanComplex, vf; npts::Int=3)

Compute the sampled dynamical transitions of a vector field in three
dimensions on a `EuclideanComplex`.

This method extracts the vertex coordinates from the embedded coordinates
of `ec` and delegates to the standard `spatial_transitions` implementation.
The cell indices in the result refer to `ec`. The complex has to be
three-dimensional and embedded in three dimensions, otherwise an
`ArgumentError` is raised.
"""
function spatial_transitions(ec::EuclideanComplex, vf; npts::Int=3)
    _check_spatial_complex(ec, "spatial_transitions")
    lc = euclidean_to_lefschetz(ec)
    return spatial_transitions(lc, euclidean_vertex_coords(ec), vf, npts=npts)
end
