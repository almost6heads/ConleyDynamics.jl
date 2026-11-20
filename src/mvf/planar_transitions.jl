export planar_transitions

"""
    planar_transitions(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}}, vf;
                       npts::Int=10, nontransverse::Bool=false)

Compute the sampled dynamical transitions of a planar vector field.

The function expects a planar Lefschetz complex `lc`, a coordinate vector
`coords` of coordinates for all the 0-dimensional cells in the complex,
and a vector field `vf(z::Vector{Float64})::Vector{Float64}` on the plane,
as in `create_planar_mvf`. For every vertex and every edge `c` of `lc`, it
determines the regions, i.e., the two-dimensional cells, which the flow can
enter from `c`. The result is a vector `transitions` of integer vectors,
where every entry has the form `[c; targets]`. Cells for which no region
is found, such as boundary cells through which the flow leaves the domain,
are omitted, and the regions themselves are not listed.

The transitions are approximations of the immediate future of the cells,
i.e., of the sets `{c} ∪ Ift(c)`. They are computed from finitely many
samples of the vector field, and are therefore not rigorous. For a vertex,
the direction of the vector field at the vertex is compared with the
adjacent edges. For an edge, the sign of the normal component of the
vector field is evaluated at `npts` equally spaced points in the interior
of the edge, and both adjacent regions are selected if the sign is not
constant. If `nontransverse=true`, then every interior edge which is
flagged by `planar_nontransverse_edges` is in addition joined to both of
its regions. This test is stricter, since it also samples the vertices of
the edge and treats zeros as tangencies.

The smallest multivector field which contains all transitions is
`create_mvf_hull(lc, transitions)`, which is the result of
`create_planar_mvf`. The transitions of the reversed flow, i.e., the
approximations of the immediate past, are obtained via
`planar_transitions(lc, coords, x -> -vf(x))`. Both computations use
the same sample points, so that forward and backward transitions are
consistent with each other. This also holds for `nontransverse=true`,
since an edge is flagged for the reversed flow if and only if it is
flagged for the original one.

This function is based on the function `create_planar_dts` by Justin
Thorpe from pull request #10, which built the corresponding transition
graph, see also his dissertation *Periodic Orbits in the Setting of
Combinatorial Dynamics* (George Mason University, 2026).
"""
function planar_transitions(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}}, vf;
                            npts::Int=10, nontransverse::Bool=false)
    #
    # Compute the sampled dynamical transitions of a planar vector field
    #

    # Find the indices for the vertices and the edges

    vi = findall(isequal(0), lc.dimensions)
    ei = findall(isequal(1), lc.dimensions)

    # Find the vertex to region map

    v2r = [Vector{Int}() for _ in 1:length(vi)]

    Threads.@threads for k = 1:length(vi)
        v2r[k] = planar_mvf_vertex2region(vi[k], lc, coords, vf)
    end

    # Find the edge to region map

    e2r = [Vector{Int}() for _ in 1:length(ei)]

    Threads.@threads for k = 1:length(ei)
        e2r[k] = planar_mvf_edge2region(ei[k], lc, coords, vf, npts=npts)
    end

    # If requested, join every interior edge which is not flow
    # transverse to both of its regions. Boundary edges are left
    # unchanged, since their single region is already selected
    # whenever the flow can enter it.

    if nontransverse
        eindex = zeros(Int, lc.ncells)
        eindex[ei] = 1:length(ei)
        for e in planar_nontransverse_edges(lc, coords, vf)
            regions = lefschetz_coboundary(lc, e)
            if length(regions) == 2
                e2r[eindex[e]] = sort(regions)
            end
        end
    end

    # Collect the transitions

    transitions = Vector{Vector{Int}}()

    for k in 1:length(vi)
        if length(v2r[k]) > 0
            push!(transitions, vcat([vi[k]], v2r[k]))
        end
    end

    for k in 1:length(ei)
        if length(e2r[k]) > 0
            push!(transitions, vcat([ei[k]], e2r[k]))
        end
    end

    # Return the transitions

    return transitions
end

"""
    planar_transitions(ec::EuclideanComplex, vf; npts::Int=10,
                       nontransverse::Bool=false)

Compute the sampled dynamical transitions of a planar vector field on a
`EuclideanComplex`.

This method extracts the vertex coordinates from the embedded coordinates
of `ec` and delegates to the standard `planar_transitions` implementation.
The cell indices in the result refer to `ec`. The complex has to be
embedded in two dimensions, otherwise an `ArgumentError` is raised.
"""
function planar_transitions(ec::EuclideanComplex, vf; npts::Int=10,
                            nontransverse::Bool=false)
    _check_coords2d(ec, "planar_transitions",
                    "For complexes in three dimensions, use spatial_transitions.")
    lc = euclidean_to_lefschetz(ec)
    return planar_transitions(lc, euclidean_vertex_coords(ec), vf, npts=npts,
                              nontransverse=nontransverse)
end
