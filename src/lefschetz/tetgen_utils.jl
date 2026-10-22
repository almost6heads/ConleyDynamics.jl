export TetGenResult, tetgen_mesh, tetgen_to_simplicial
export tetgen_marker_cells, tetgen_region_cells, create_simplicial_delaunay3d

"""
    TetGenResult

Tetrahedral mesh created by `tetgen_mesh`.

The struct has the following fields:
* `points::Vector{Vector{Float64}}`: The vertex coordinates of the mesh,
  including all points added by the mesh generator.
* `tetrahedra::Vector{Vector{Int}}`: The tetrahedra, as vectors of four
  indices into `points`.
* `trifaces::Vector{Vector{Int}}`: The boundary and constraint triangles,
  as vectors of three indices into `points`.
* `trifacemarkers::Vector{Int}`: The marker of the facet which contains
  each triangle in `trifaces`.
* `tetattributes::Vector{Int}`: The region attribute of every tetrahedron,
  or an empty vector if the domain has no regions.
* `curvesegments::Vector{Tuple{NTuple{3,Float64},NTuple{3,Float64},Int}}`:
  The segments of all constraint curves of the domain, given by their
  endpoints and their marker.

The mesh can be converted into an `EuclideanComplex` with
`tetgen_to_simplicial`, and the cells on constraint surfaces and curves,
or in regions, can be found with `tetgen_marker_cells` and
`tetgen_region_cells`.
"""
struct TetGenResult
    points::Vector{Vector{Float64}}
    tetrahedra::Vector{Vector{Int}}
    trifaces::Vector{Vector{Int}}
    trifacemarkers::Vector{Int}
    tetattributes::Vector{Int}
    curvesegments::Vector{Tuple{NTuple{3,Float64},NTuple{3,Float64},Int}}
end

function Base.show(io::IO, res::TetGenResult)
    print(io, "TetGenResult with ", length(res.points), " points, ",
          length(res.tetrahedra), " tetrahedra, ", length(res.trifaces),
          " boundary and constraint triangles")
end

"""
    tetgen_mesh(dom::TetDomain; maxvolume::Union{Real,Nothing}=nothing,
                minratio::Union{Real,Nothing}=1.414,
                mindihedral::Union{Real,Nothing}=nothing,
                steiner_on_boundary::Bool=true, quiet::Bool=true,
                switches::Union{String,Nothing}=nothing) -> TetGenResult

Create a constrained tetrahedral mesh of a domain with `TetGen.jl`.

The domain `dom` is described by a `TetDomain`, see for example
`tetgen_domain_box` and `tetgen_domain_add_surface`. The keyword
arguments control the mesh quality and size:
* `maxvolume`: Maximum volume of the tetrahedra. If `nothing`, there is
  no global volume constraint. Regional volume constraints of the domain
  are always applied. Note that `TetGen` treats volume constraints as
  approximate: typically about one percent of the tetrahedra exceed
  them, by up to about fifty percent.
* `minratio`: Maximum ratio of circumradius to shortest edge length of
  the tetrahedra. Smaller values give better shaped tetrahedra, but
  more of them. The default `1.414` is the standard value of `TetGen`.
  If both `minratio` and `mindihedral` are `nothing`, no quality
  refinement is performed.
* `mindihedral`: Minimum dihedral angle of the tetrahedra, in degrees.
  If `nothing`, only `minratio` is used.
* `steiner_on_boundary`: If `false`, no additional points are added on
  the boundary and constraint facets. Note that this also severely
  limits the refinement close to these facets.
* `quiet`: If `true`, the output of `TetGen` is suppressed.
* `switches`: A complete `TetGen` switch string, which overrides all
  other keyword arguments. It should contain `p` if the facets of the
  domain are to be respected, and it should not contain `f`, since
  otherwise interior triangles are reported as boundary triangles
  with marker `0`.

Quality refinement is switched on by default, since flat tetrahedra
can make the construction of multivector fields fragile.

If the facets of the domain intersect, `TetGen` fails, and an
`ArgumentError` is raised. In that case, check the constraint surfaces
and curves, or use `Gmsh.jl` for intersecting geometry. An
`ArgumentError` is also raised if the mesh does not contain any
tetrahedra, which happens if the outer boundary is not closed, or if
a hole point does not lie inside a closed constraint surface.

The resulting `TetGenResult` can be converted into an `EuclideanComplex`
with `tetgen_to_simplicial`.

This function requires `TetGen.jl` to be loaded.

# Example

```julia
using TetGen, ConleyDynamics

dom = tetgen_domain_box([-2,-2,-2], [2,2,2])
V, T = trisurface_sphere([0,0,0], 1.0, 2)
dom = tetgen_domain_add_surface(dom, V, T)
res = tetgen_mesh(dom; maxvolume = 0.01)
ec  = tetgen_to_simplicial(res)
sphere_cells = tetgen_marker_cells(res, ec, 1)
```
"""
function tetgen_mesh end

"""
    create_simplicial_delaunay3d(bmin::AbstractVector{<:Real},
                                 bmax::AbstractVector{<:Real},
                                 pdist::Real, attempts::Int;
                                 p::Int=2, box_corners::Bool=false,
                                 rng::AbstractRNG=Random.default_rng())
                                 -> EuclideanComplex

Create a three-dimensional Delaunay triangulation of a random point
sample inside a box. The complex is over the rationals if `p=0`, and
over `GF(p)` if `p>0`.

The function selects a random sample of points inside the box
`[bmin[1],bmax[1]] x [bmin[2],bmax[2]] x [bmin[3],bmax[3]]`, with
minimum distance `pdist` between the points, see `poisson_points_box`
for the meaning of `attempts`. It then computes the Delaunay
tetrahedralization of these points with `TetGen.jl`, and returns it
as an `EuclideanComplex`. The random number generator can be specified
via `rng` for reproducible results.

As in the planar case `create_simplicial_delaunay`, the complex covers
only the convex hull of the sample, so it does not fill the box
completely near the boundary, and the tetrahedra on the boundary of
the convex hull can be very flat. If `box_corners=true`, the eight
corners of the box are added to the sample, so that the convex hull is
the whole box. For quality meshes of a box, use `tetgen_mesh` with
`tetgen_domain_box` instead.

This function requires `TetGen.jl` to be loaded.
"""
function create_simplicial_delaunay3d end

"""
    create_simplicial_delaunay3d(bmin::AbstractVector{<:Real},
                                 bmax::AbstractVector{<:Real},
                                 npoints::Int;
                                 p::Int=2, box_corners::Bool=false,
                                 rng::AbstractRNG=Random.default_rng())
                                 -> EuclideanComplex

Create a three-dimensional Delaunay triangulation of `npoints` uniformly
distributed random points inside a box. The complex is over the
rationals if `p=0`, and over `GF(p)` if `p>0`.

See the method with arguments `pdist` and `attempts` for details. This
function requires `TetGen.jl` to be loaded.
"""
create_simplicial_delaunay3d(::AbstractVector{<:Real}, ::AbstractVector{<:Real}, ::Int)

"""
    tetgen_to_simplicial(res::TetGenResult; p::Int=2) -> EuclideanComplex

Convert a tetrahedral mesh created by `tetgen_mesh` into an
`EuclideanComplex`. The complex is over the rationals if `p=0`, and
over `GF(p)` if `p>0`.

Vertices of the mesh which do not belong to any tetrahedron are removed.
Vertex labels are zero-padded decimal strings as in
`tetmesh_to_euclidean`.
"""
function tetgen_to_simplicial(res::TetGenResult; p::Int=2)
    #
    # Convert a TetGen mesh into an EuclideanComplex
    #
    return tetmesh_to_euclidean(res.points, res.tetrahedra; p=p, droporphans=true)
end

"""
    tetgen_marker_cells(res::TetGenResult, ec::EuclideanComplex, marker::Int;
                        closure::Bool=true) -> Vector{Int}

Find the cells of a complex which lie on the facets with a given marker.

The complex `ec` has to be created from `res` via `tetgen_to_simplicial`.
The function returns the indices of all triangles of `ec` which lie on
facets of the domain with marker `marker`, as well as all edges on
constraint curves with this marker. If `closure=true` (default), all
faces of these cells are included as well, so the result is a closed
subcomplex. For a constraint surface, this gives the surface as a
two-dimensional subcomplex, and for the marker `0` it gives the outer
boundary of the domain.

The returned cells can be used, for example, to define exit sets or to
restrict a multivector field to a constraint surface.
"""
function tetgen_marker_cells(res::TetGenResult, ec::EuclideanComplex, marker::Int;
                             closure::Bool=true)
    #
    # Find the cells on the facets with a given marker
    #

    labels = _tetgen_point_labels(res)
    cells  = Vector{Int}()

    # Triangles on surface facets

    for (t, m) in zip(res.trifaces, res.trifacemarkers)
        m == marker || continue
        any(isempty, labels[t]) && continue
        k = get(ec.indices, join(sort(labels[t])), 0)
        k > 0 && push!(cells, k)
    end

    # Edges on curve segments. The mesh generator splits every segment
    # into a chain of edges, whose vertices lie on the segment.

    segments = [(collect(a), collect(b)) for (a, b, m) in res.curvesegments if m == marker]
    if !isempty(segments)
        vertices = findall(==(0), ec.dimensions)
        for (a, b) in segments
            d = b .- a
            len2 = sum(abs2, d)
            onseg = Vector{Tuple{Float64,Int}}()
            for v in vertices
                x = ec.coords[v][1]
                t = sum((x .- a) .* d) / len2
                if (-1e-9 <= t <= 1 + 1e-9) &&
                   (sum(abs2, x .- a .- t .* d) <= 1e-18 * len2)
                    push!(onseg, (t, v))
                end
            end
            sort!(onseg)
            for k in 1:length(onseg)-1
                v1, v2 = onseg[k][2], onseg[k+1][2]
                k1 = get(ec.indices, join(sort([ec.labels[v1], ec.labels[v2]])), 0)
                k1 > 0 && push!(cells, k1)
            end
        end
    end

    return closure ? lefschetz_closure(ec, cells) : sort(unique(cells))
end

"""
    tetgen_region_cells(res::TetGenResult, ec::EuclideanComplex, attr::Int;
                        closure::Bool=true) -> Vector{Int}

Find the tetrahedra of a complex which lie in the regions with a given
attribute.

The complex `ec` has to be created from `res` via `tetgen_to_simplicial`,
and the domain has to contain regions, see `tetgen_domain_add_region`
and `tetgen_domain_add_surface`. The function returns the indices of all
tetrahedra with region attribute `attr`. If `closure=true` (default),
all faces of these tetrahedra are included as well, so the result is a
closed subcomplex.

Note that `TetGen` assigns new attributes, larger than all given ones,
to the parts of the domain which do not contain a region point. For
example, if a single region with attribute `1` is defined inside a
closed constraint surface, then the tetrahedra outside of this surface
get attribute `2`. The attributes of all tetrahedra are stored in
`res.tetattributes`.
"""
function tetgen_region_cells(res::TetGenResult, ec::EuclideanComplex, attr::Int;
                             closure::Bool=true)
    #
    # Find the tetrahedra in the regions with a given attribute
    #

    if isempty(res.tetattributes)
        throw(ArgumentError("tetgen_region_cells: the mesh has no region attributes"))
    end

    labels = _tetgen_point_labels(res)
    cells  = Vector{Int}()
    for (t, a) in zip(res.tetrahedra, res.tetattributes)
        a == attr || continue
        k = get(ec.indices, join(sort(labels[t])), 0)
        k > 0 && push!(cells, k)
    end

    return closure ? lefschetz_closure(ec, cells) : sort(unique(cells))
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function _tetgen_point_labels(res::TetGenResult)
    #
    # Vertex labels of the mesh points in the complex created by
    # tetgen_to_simplicial, or an empty string for unused points
    #
    used = falses(length(res.points))
    for t in res.tetrahedra
        used[t] .= true
    end
    newindex = cumsum(used)
    vlabels  = _vertex_labels(newindex[end])
    return [used[k] ? vlabels[newindex[k]] : "" for k in 1:length(res.points)]
end
