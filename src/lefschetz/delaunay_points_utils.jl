export delaunay_points_bnd_rectangle
export delaunay_points_add_segment
export delaunay_points_add_nodes
export delaunay_points_add_split_segs

"""
    delaunay_points_bnd_rectangle(bmin, bmax)
        -> (Vector{Tuple{Float64,Float64}}, Vector{Vector{Int}})

Create the initial point list and closed outer boundary curve for a
rectangular domain, for use as input to `triangulate` from
`DelaunayTriangulation.jl`.

The arguments `bmin = [xmin, ymin]` and `bmax = [xmax, ymax]` give the
lower-left and upper-right corners of the bounding rectangle. The four
corners are stored as `Float64` tuples in counter-clockwise order. The
returned `bndcurve` encodes the closed loop in the format expected by the
`boundary_nodes` keyword of `triangulate`.

Returns a tuple `(points, bndcurve)` where:
- `points` is a `Vector{Tuple{Float64,Float64}}` containing the four
  corner vertices,
- `bndcurve` is a `Vector{Vector{Int}}` that closes the rectangular
  boundary.

After this call, constraint edges and interior points can be appended with
[`delaunay_points_add_segment`](@ref), and the result is then passed to
`triangulate` from `DelaunayTriangulation.jl` followed by
[`delaunay_to_simplicial`](@ref) to obtain a `EuclideanComplex`.

This function is only available when `DelaunayTriangulation.jl` is loaded.

# Example

```julia
using DelaunayTriangulation
using ConleyDynamics

points, bndcurve = delaunay_points_bnd_rectangle([-5, -5], [5, 5])
tri = triangulate(points; boundary_nodes = bndcurve)
sc  = delaunay_to_simplicial(tri)
```
"""
function delaunay_points_bnd_rectangle end

"""
    delaunay_points_add_segment(points, segments, newseg)
        -> (Vector{Tuple{Float64,Float64}}, Set{Tuple{Int,Int}})
    delaunay_points_add_segment(points, newseg)
        -> (Vector{Tuple{Float64,Float64}}, Set{Tuple{Int,Int}})

Append a polygonal curve to the point list and constraint-edge set for a
constrained Delaunay triangulation.

`newseg` is a `Vector{Vector{<:Real}}` where each element `[x, y]` is a
2D point on the curve. If the first and last entries of `newseg` coincide
(Euclidean distance less than `1e-9`) the curve is treated as *closed*: only
one copy of the shared endpoint is stored and the closing edge connects the
last new point back to the first. Otherwise the curve is treated as *open*
and both endpoints are stored separately.

All points in `newseg` are appended to `points` (converted to `Float64`)
and all successive pairs are added as directed constraint edges to
`segments`.

The three-argument form takes an existing `segments::Set{Tuple{Int,Int}}`
and extends it. The two-argument form is a convenience wrapper that
initialises `segments` as an empty set; use it when adding the first
constraint curve after [`delaunay_points_bnd_rectangle`](@ref).

Returns the updated `(points, segments)`.

This function is only available when `DelaunayTriangulation.jl` is loaded.

# Example

```julia
using DelaunayTriangulation
using ConleyDynamics

# Rectangular domain with a circular hole and an open cut
points, bndcurve = delaunay_points_bnd_rectangle([-5, -5], [5, 5])

n = 12
theta  = range(0, 2*pi, length = n+1)
circle = [[4.0*cos(t), 4.0*sin(t)] for t in theta]  # closed curve
cut    = [[-2.0, -2.0], [0.0, 1.0], [2.0, -1.0]]    # open curve

points, segments = delaunay_points_add_segment(points, circle)
points, segments = delaunay_points_add_segment(points, segments, cut)

tri = triangulate(points; boundary_nodes = bndcurve, segments = segments)
sc  = delaunay_to_simplicial(tri)
```
"""
function delaunay_points_add_segment end

"""
    delaunay_points_add_nodes(points, newpoints)
        -> Vector{Tuple{Float64,Float64}}

Append a list of interior nodes to the given points for a Delaunay triangulation.

`newpoints` is a `Vector{Vector{<:Real}}` where each element `[x, y]` is a
2D point which describes an interior node. All points in `newpoints` are appended
to `points` (converted to `Float64`), and the function returns the updated `points`.

This function is only available when `DelaunayTriangulation.jl` is loaded.

# Example

```julia
using DelaunayTriangulation
using ConleyDynamics

# Rectangular domain with a circular hole and an open cut
points, bndcurve = delaunay_points_bnd_rectangle([-5, -5], [5, 5])

newpoints = [-4 .+ 8 .* rand(2) for k in 1:12]   # 12 random points
points    = delaunay_points_add_nodes(points, newpoints)

tri = triangulate(points; boundary_nodes = bndcurve)
sc  = delaunay_to_simplicial(tri)
```
"""
function delaunay_points_add_nodes end

"""
    delaunay_points_add_split_segs(points, newsegments; verbose=false)
        -> (Vector{Tuple{Float64,Float64}}, Set{Tuple{Int,Int}})
    delaunay_points_add_split_segs(points, segments, newsegments; verbose=false)
        -> (Vector{Tuple{Float64,Float64}}, Set{Tuple{Int,Int}})

Compute a planar straight-line arrangement from `newsegments`, add all
resulting vertices and sub-segments to `points` and `segments`, and return
the updated pair. This is the safe alternative to
[`delaunay_points_add_segment`](@ref) when the input curves may intersect
each other or have collinear overlaps.

`newsegments` is a `Vector` of segments, each given as a length-2 vector of
`[x, y]` endpoints, e.g. `[[x1, y1], [x2, y2]]`. The function calls
[`split_segments`](@ref) internally to resolve all crossings and overlaps
into a non-crossing arrangement, then appends the resulting vertices and
edges to the existing `points` and `segments`.

The two-argument form is a convenience wrapper that initialises `segments`
as an empty set. The three-argument form takes an existing
`segments::Set{Tuple{Int,Int}}` and extends it.

If the keyword argument `verbose` is `true`, the function prints a brief
summary: number of input segments, number of intersection points created,
and number of output vertices and edges added to `points`.

Returns the updated `(points, segments)`.

!!! warning "Three-argument form: overlaps with existing segments are not checked"
    When calling the three-argument form, the new segments in `newsegments`
    are split among themselves by [`split_segments`](@ref), but **no check is
    performed for intersections or overlaps with the segments already present
    in `segments`**. It is the caller's responsibility to ensure that the
    newly added sub-segments are compatible with the existing constraint-edge
    set. Incorrect use can produce an invalid segment set that silently causes
    `DelaunayTriangulation.jl` to drop constraints.

This function is only available when `DelaunayTriangulation.jl` is loaded.

# Example

```julia
using DelaunayTriangulation
using ConleyDynamics

points, bndcurve = delaunay_points_bnd_rectangle([-2, -2], [2, 2])

# Two crossing diagonals — split_segments resolves the crossing automatically
segs = [[[-1.0, -1.0], [1.0, 1.0]],
        [[-1.0,  1.0], [1.0, -1.0]]]

points, segments = delaunay_points_add_split_segs(points, segs; verbose=true)

tri = triangulate(points; boundary_nodes = bndcurve, segments = segments)
sc  = delaunay_to_simplicial(tri)
```
"""
function delaunay_points_add_split_segs end

