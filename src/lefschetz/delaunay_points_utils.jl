export delaunay_points_bnd_rectangle
export delaunay_points_add_segment

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
