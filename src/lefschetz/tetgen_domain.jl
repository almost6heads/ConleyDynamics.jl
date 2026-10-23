export TetDomain, tetgen_domain_box, tetgen_domain_from_surface
export tetgen_domain_add_nodes, tetgen_domain_add_surface, tetgen_domain_add_curve
export tetgen_domain_add_hole, tetgen_domain_add_region

"""
    TetDomain

Description of a three-dimensional domain for tetrahedral mesh generation
with `TetGen.jl`, as a piecewise linear complex.

The struct has the following fields:
* `points::Vector{NTuple{3,Float64}}`: The vertices of the domain description.
* `facets::Vector{Vector{Int}}`: The facets, each given as a vector of
  indices into `points`. A facet with at least three vertices is a planar
  polygon, a facet with two vertices is a segment, and a facet with one
  vertex is a single point. All of them are contained in the final mesh.
* `markers::Vector{Int}`: One marker for every facet. By convention, the
  outer boundary has marker `0`, and the `k`-th added constraint surface
  or curve has marker `k`.
* `holes::Vector{NTuple{3,Float64}}`: For every hole, one point inside it.
  The region around this point which is enclosed by facets is not meshed.
* `regions::Vector{Tuple{NTuple{3,Float64},Int,Float64}}`: For every
  region, a point inside it, an integer attribute, and the maximum
  tetrahedron volume in that region. A nonpositive maximum volume means
  that there is no regional volume constraint.

A `TetDomain` is usually created with `tetgen_domain_box` or
`tetgen_domain_from_surface`, and then extended with the functions
`tetgen_domain_add_nodes`, `tetgen_domain_add_surface`,
`tetgen_domain_add_curve`, `tetgen_domain_add_hole`, and
`tetgen_domain_add_region`. None of these functions require `TetGen.jl`.
The domain is then meshed with `tetgen_mesh`.

All facets have to be non-intersecting, except along shared vertices
and edges. In particular, constraint surfaces cannot cross each other
or the outer boundary.
"""
struct TetDomain
    points::Vector{NTuple{3,Float64}}
    facets::Vector{Vector{Int}}
    markers::Vector{Int}
    holes::Vector{NTuple{3,Float64}}
    regions::Vector{Tuple{NTuple{3,Float64},Int,Float64}}
    #
    # Inner constructor
    #
    function TetDomain(points, facets, markers, holes=NTuple{3,Float64}[],
                       regions=Tuple{NTuple{3,Float64},Int,Float64}[])
        #
        # Create a TetDomain instance after basic consistency checks
        #
        if length(facets) != length(markers)
            throw(ArgumentError("TetDomain: need one marker per facet"))
        end
        npts = length(points)
        for f in facets
            if isempty(f) || any(v -> (v < 1) || (v > npts), f)
                throw(ArgumentError("TetDomain: facet $f is empty or contains " *
                                    "an index outside of 1:$npts"))
            end
        end
        return new(points, facets, markers, holes, regions)
    end
end

function Base.show(io::IO, dom::TetDomain)
    print(io, "TetDomain with ", length(dom.points), " points, ",
          length(dom.facets), " facets, markers ", sort(unique(dom.markers)),
          ", ", length(dom.holes), " holes, ", length(dom.regions), " regions")
end

"""
    tetgen_domain_box(bmin::AbstractVector{<:Real}, bmax::AbstractVector{<:Real};
                      markers::AbstractVector{<:Integer}=zeros(Int,6)) -> TetDomain

Create a `TetDomain` for the box `[bmin[1],bmax[1]] x [bmin[2],bmax[2]] x
[bmin[3],bmax[3]]`.

The domain consists of the eight corners of the box and its six faces,
which are quadrilateral facets with outward orientation. By default, all
faces have marker `0`. The optional argument `markers` assigns individual
markers to the faces, in the order `x = bmin[1]`, `x = bmax[1]`,
`y = bmin[2]`, `y = bmax[2]`, `z = bmin[3]`, and `z = bmax[3]`. This
makes it possible to find the cells on some of the faces with
`tetgen_marker_cells`, for example to describe an exit set. This is the
three-dimensional analogue of `delaunay_points_bnd_rectangle`.
"""
function tetgen_domain_box(bmin::AbstractVector{<:Real}, bmax::AbstractVector{<:Real};
                           markers::AbstractVector{<:Integer}=zeros(Int,6))
    #
    # Create the domain description of a box
    #

    if (length(bmin) != 3) || (length(bmax) != 3)
        throw(ArgumentError("tetgen_domain_box: bmin and bmax need to have length 3"))
    end
    if any(bmax .<= bmin)
        throw(ArgumentError("tetgen_domain_box: bmax needs to be larger than " *
                            "bmin in every coordinate"))
    end

    if length(markers) != 6
        throw(ArgumentError("tetgen_domain_box: markers needs one entry for " *
                            "each of the six faces"))
    end

    x0, y0, z0 = Float64.(bmin)
    x1, y1, z1 = Float64.(bmax)
    points = [(x0,y0,z0), (x1,y0,z0), (x1,y1,z0), (x0,y1,z0),
              (x0,y0,z1), (x1,y0,z1), (x1,y1,z1), (x0,y1,z1)]

    # The facets are the faces z = z0, z = z1, y = y0, x = x1, y = y1,
    # and x = x0, so that markers in the documented order are permuted

    facets = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [2,3,7,6], [3,4,8,7], [4,1,5,8]]
    fmarkers = Int.(collect(markers))[[5, 6, 3, 2, 4, 1]]

    return TetDomain(points, facets, fmarkers)
end

"""
    tetgen_domain_from_surface(V::AbstractVector, T::AbstractVector) -> TetDomain

Create a `TetDomain` whose outer boundary is the closed triangulated
surface given by the vertices `V` and the triangles `T`.

The surface is usually created by one of the `trisurface_*` functions,
and it has to be closed, i.e., every edge has to belong to exactly two
triangles. Otherwise an `ArgumentError` is raised. Vertices with the
same coordinates are identified before this check, so the surface can
also be given as a list of separate triangles, each with its own copies
of the vertices, as in STL files. For example, the
output of `trisurface_sphere` gives a ball, and the output of
`trisurface_torus` gives a solid torus. All facets get marker `0`.
"""
function tetgen_domain_from_surface(V::AbstractVector, T::AbstractVector)
    #
    # Create a domain bounded by a closed triangulated surface
    #

    nv = length(V)
    if any(any(v -> !(1 <= v <= nv), t) for t in T)
        throw(ArgumentError("tetgen_domain_from_surface: the triangles contain " *
                            "indices outside of 1:$nv"))
    end
    if any(length(t) != 3 for t in T)
        throw(ArgumentError("tetgen_domain_from_surface: every triangle needs " *
                            "three vertex indices"))
    end

    # Vertices with the same coordinates are merged before the surface
    # is checked, so that triangle soups as read from STL files work

    points, newindex = _merge_points(NTuple{3,Float64}[], [_point3(v) for v in V])
    facets = [newindex[collect(t)] for t in T]
    _check_triangles(facets, "tetgen_domain_from_surface")
    closed, _ = _trisurface_check(facets)
    if !closed
        throw(ArgumentError("tetgen_domain_from_surface: the surface is not closed"))
    end

    return TetDomain(points, facets, zeros(Int, length(facets)))
end

"""
    tetgen_domain_add_nodes(dom::TetDomain, newpoints::AbstractVector) -> TetDomain

Add nodes to a `TetDomain`.

The vector `newpoints` contains the coordinates of the new nodes, as
vectors or tuples of length three. Each of these nodes inside the domain
becomes a vertex of the final mesh. Nodes outside of the domain or in
holes are ignored by the mesh generator, and nodes which coincide with
points of the domain are merged with them. Nodes must not lie on a
facet or a constraint curve without being one of its vertices, since
`TetGen` reports such nodes as self-intersections and fails. The
function returns a new `TetDomain`, and the original domain is not
changed.
"""
function tetgen_domain_add_nodes(dom::TetDomain, newpoints::AbstractVector)
    #
    # Add interior nodes to a domain
    #
    points = vcat(dom.points, [_point3(v) for v in newpoints])
    return TetDomain(points, deepcopy(dom.facets), copy(dom.markers),
                     copy(dom.holes), copy(dom.regions))
end

"""
    tetgen_domain_add_surface(dom::TetDomain, V::AbstractVector, T::AbstractVector;
                              marker::Union{Int,Nothing}=nothing,
                              hole::Bool=false,
                              region::Union{Tuple{Int,Real},Nothing}=nothing)
                              -> TetDomain

Add a triangulated constraint surface to a `TetDomain`.

The surface is given by its vertices `V` and its triangles `T`, for
example as created by one of the `trisurface_*` functions. Its triangles
become facets of the domain, so the final mesh will contain the surface
as a union of triangles. The facets get the marker `marker`, and if no
marker is given, the next free marker `maximum(dom.markers) + 1` is used,
or `1` if this maximum is zero or less. The cells of the final complex
on the surface can be found with `tetgen_marker_cells`.

For closed surfaces, two additional options are available, which
cannot be combined:
* `hole=true` removes the inside of the surface from the domain.
* `region=(attr, maxvol)` marks the inside of the surface as a region
  with integer attribute `attr`, and with maximum tetrahedron volume
  `maxvol` (use a nonpositive value for no volume constraint). The
  tetrahedra in this region can be found with `tetgen_region_cells`.

In both cases, the surface has to be closed and consistently oriented,
and a point inside the surface is determined automatically. This point
avoids the surfaces which are already part of the domain, so that a
region bounded by a sphere around a previously added hole is the shell
between the two surfaces. If the surface consists of several connected
components, such as two disjoint spheres, one point is determined for
every component, so that all of them become holes or regions. Note that for nested components, such as
two concentric spheres, both the inner ball and the shell between the
spheres are removed or marked. If this is not intended, add the
surfaces separately, starting with the inner one, or use
`tetgen_domain_add_hole` or `tetgen_domain_add_region` with explicit
points.

Vertices of the surface which coincide exactly with each other or with
existing points of the domain are identified, before the surface is
checked. Triangles which then have repeated vertices are rejected.
Apart from such shared vertices
and edges, the surface must not intersect other surfaces or the
boundary of the domain. For intersecting surfaces, use `Gmsh.jl`
instead, see the manual. The function returns a new `TetDomain`, and the
original domain is not changed.
"""
function tetgen_domain_add_surface(dom::TetDomain, V::AbstractVector, T::AbstractVector;
                                   marker::Union{Int,Nothing}=nothing,
                                   hole::Bool=false,
                                   region::Union{Tuple{Int,Real},Nothing}=nothing)
    #
    # Add a triangulated constraint surface to a domain
    #

    nv = length(V)
    if any(any(v -> !(1 <= v <= nv), t) for t in T)
        throw(ArgumentError("tetgen_domain_add_surface: the triangles contain " *
                            "indices outside of 1:$nv"))
    end
    if any(length(t) != 3 for t in T)
        throw(ArgumentError("tetgen_domain_add_surface: every triangle needs " *
                            "three vertex indices"))
    end
    if hole && !isnothing(region)
        throw(ArgumentError("tetgen_domain_add_surface: a surface cannot bound " *
                            "both a hole and a region"))
    end

    # Vertices with the same coordinates are merged before the surface
    # is checked, so that triangle soups as read from STL files work

    cmarker = isnothing(marker) ? _next_marker(dom) : marker
    points, newindex = _merge_points(dom.points, [_point3(v) for v in V])
    newfacets = [newindex[collect(t)] for t in T]
    _check_triangles(newfacets, "tetgen_domain_add_surface")

    if hole || !isnothing(region)
        closed, oriented = _trisurface_check(newfacets)
        if !closed
            throw(ArgumentError("tetgen_domain_add_surface: holes and regions " *
                                "require a closed surface"))
        end
        if !oriented
            throw(ArgumentError("tetgen_domain_add_surface: holes and regions " *
                                "require a consistently oriented surface"))
        end
    end

    facets  = vcat(deepcopy(dom.facets), newfacets)
    markers = vcat(dom.markers, fill(cmarker, length(T)))
    holes   = copy(dom.holes)
    regions = copy(dom.regions)

    # One hole or region point inside every connected component. The
    # ray which determines it stops at the new surface and at all
    # polygonal facets of the domain, so that the point does not end
    # up inside a surface which was added before, such as a smaller
    # concentric sphere.

    if hole || !isnothing(region)
        oldtris = [[f[1], f[k], f[k+1]] for f in dom.facets for k in 2:length(f)-1]
        alltris = vcat(newfacets, oldtris)
        for comp in _trisurface_components(newfacets)
            q = _point3(_trisurface_interior_point(points, newfacets[comp], alltris))
            hole && push!(holes, q)
            if !isnothing(region)
                push!(regions, (q, region[1], Float64(region[2])))
            end
        end
    end

    return TetDomain(points, facets, markers, holes, regions)
end

"""
    tetgen_domain_add_curve(dom::TetDomain, curve::AbstractVector;
                            marker::Union{Int,Nothing}=nothing) -> TetDomain

Add a polygonal constraint curve to a `TetDomain`.

The vector `curve` contains the consecutive points of the curve, as
vectors or tuples of length three. If the first and the last point
agree up to `1e-9` times the size of the curve, the curve is closed.
Every segment of the curve becomes a facet of the domain, so the final mesh will
contain the curve as a union of edges. Repeated consecutive points and
segments which occur more than once are ignored. The facets get the marker
`marker`, and if no marker is given, the next free marker is used as in
`tetgen_domain_add_surface`.

Vertices of the curve which coincide exactly with existing points of
the domain are identified with them. Apart from such shared vertices,
the curve must not intersect any facets of the domain, and in
particular it cannot lie inside a face of the outer boundary. Curves
outside of the domain or inside holes are ignored by `TetGen`, and
`tetgen_mesh` prints a warning in this case. The function
returns a new `TetDomain`, and the original domain is not changed.
"""
function tetgen_domain_add_curve(dom::TetDomain, curve::AbstractVector;
                                 marker::Union{Int,Nothing}=nothing)
    #
    # Add a polygonal constraint curve to a domain
    #

    if length(curve) < 2
        throw(ArgumentError("tetgen_domain_add_curve: the curve needs at least two points"))
    end

    # The curve is closed if its endpoints agree up to a tolerance
    # relative to its size, so that rounding errors in the coordinates
    # of large curves do not matter

    cpoints = [_point3(v) for v in curve]
    csize   = maximum(maximum(p[k] for p in cpoints) - minimum(p[k] for p in cpoints)
                      for k in 1:3)
    closed  = (length(cpoints) > 2) &&
              (sqrt(sum(abs2, cpoints[1] .- cpoints[end])) <= 1e-9 * csize)
    if closed
        pop!(cpoints)
    end

    cmarker = isnothing(marker) ? _next_marker(dom) : marker
    points, newindex = _merge_points(dom.points, cpoints)
    nc = length(cpoints)
    allsegs = [(newindex[k], newindex[k+1]) for k in 1:nc-1]
    if closed
        push!(allsegs, (newindex[nc], newindex[1]))
    end

    # Remove degenerate segments (repeated consecutive points) and
    # segments which occur more than once

    segments = Vector{Vector{Int}}()
    seen = Set{Tuple{Int,Int}}()
    for (i, j) in allsegs
        key = (min(i,j), max(i,j))
        if (i != j) && !(key in seen)
            push!(seen, key)
            push!(segments, [i, j])
        end
    end
    if isempty(segments)
        throw(ArgumentError("tetgen_domain_add_curve: the curve needs at least " *
                            "two distinct points"))
    end

    facets  = vcat(deepcopy(dom.facets), segments)
    markers = vcat(dom.markers, fill(cmarker, length(segments)))

    return TetDomain(points, facets, markers, copy(dom.holes), copy(dom.regions))
end

"""
    tetgen_domain_add_hole(dom::TetDomain, point::AbstractVector{<:Real}) -> TetDomain

Add a hole to a `TetDomain`.

The point `point` has to lie inside a region of the domain which is
completely enclosed by facets. This region is then removed from the
mesh. A hole point which is not enclosed by facets removes the whole
domain. The function returns a new `TetDomain`.
"""
function tetgen_domain_add_hole(dom::TetDomain, point::AbstractVector{<:Real})
    #
    # Add a hole point to a domain
    #
    return TetDomain(copy(dom.points), deepcopy(dom.facets), copy(dom.markers),
                     vcat(dom.holes, [_point3(point)]), copy(dom.regions))
end

"""
    tetgen_domain_add_region(dom::TetDomain, point::AbstractVector{<:Real},
                             attr::Int, maxvol::Real=-1.0) -> TetDomain

Add a region to a `TetDomain`.

The region is the part of the domain which contains `point` and which
is bounded by facets. All tetrahedra in this region get the integer
attribute `attr`, and they can be found with `tetgen_region_cells`. If
`maxvol` is positive, the tetrahedra in the region will have volume at
most `maxvol`, up to the tolerance of `TetGen` described in
`tetgen_mesh`. The function returns a new `TetDomain`.
"""
function tetgen_domain_add_region(dom::TetDomain, point::AbstractVector{<:Real},
                                  attr::Int, maxvol::Real=-1.0)
    #
    # Add a region seed to a domain
    #
    return TetDomain(copy(dom.points), deepcopy(dom.facets), copy(dom.markers),
                     copy(dom.holes),
                     vcat(dom.regions, [(_point3(point), attr, Float64(maxvol))]))
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function _point3(v)
    #
    # Convert a point to a 3-tuple of Float64
    #
    if length(v) != 3
        throw(ArgumentError("Points need to have three coordinates, got $v"))
    end
    if !all(isfinite, v)
        throw(ArgumentError("Point coordinates need to be finite, got $v"))
    end
    # Adding 0.0 turns -0.0 into 0.0, so that equal points are merged
    return (Float64(v[1]) + 0.0, Float64(v[2]) + 0.0, Float64(v[3]) + 0.0)
end

function _check_triangles(T::Vector{Vector{Int}}, fname::String)
    #
    # Reject triangles with repeated vertices, which can also arise
    # if distinct vertices with the same coordinates are merged
    #
    for t in T
        if !allunique(t)
            throw(ArgumentError("$fname: the triangle $t has repeated vertices, " *
                                "possibly after merging vertices with the same " *
                                "coordinates"))
        end
    end
end

function _next_marker(dom::TetDomain)
    #
    # Next free facet marker, at least 1
    #
    return max(1, maximum(dom.markers; init=0) + 1)
end

function _merge_points(points::Vector{NTuple{3,Float64}},
                       newpoints::Vector{NTuple{3,Float64}})
    #
    # Append new points, identifying points which coincide exactly
    # with existing ones. Returns the new point list and the indices
    # of the new points in it.
    #
    lookup = Dict{NTuple{3,Float64},Int}(points[k] => k for k in 1:length(points))
    allpoints = copy(points)
    newindex  = Vector{Int}(undef, length(newpoints))
    for (k, pt) in enumerate(newpoints)
        if haskey(lookup, pt)
            newindex[k] = lookup[pt]
        else
            push!(allpoints, pt)
            lookup[pt] = length(allpoints)
            newindex[k] = length(allpoints)
        end
    end
    return allpoints, newindex
end
