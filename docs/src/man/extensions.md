# Extensions

[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
exposes optional functionality through *weak dependencies* — packages that
are loaded on demand and extend the library's capabilities without becoming
hard dependencies. This chapter documents the functions that become
available through these extensions. At present, two backends are provided:
`DelaunayTriangulation.jl` for constructing constrained planar triangulations
and converting them to simplicial complexes (described in the section below),
and `Plots.jl` for interactive visualization of Lefschetz complexes and their
dynamics (described in the [Visualization](@ref) chapter).

## DelaunayTriangulation.jl

The
[DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl)
package computes Delaunay triangulations of planar point sets, optionally with
constrained edges and boundary curves.
Loading it alongside `ConleyDynamics.jl` activates three functions:

* [`delaunay_points_bnd_rectangle`](@ref) — initialise a point list and
  outer boundary for a rectangular domain,
* [`delaunay_points_add_segment`](@ref) — append constraint curves (closed
  or open) to that point list, and
* [`delaunay_to_simplicial`](@ref) — convert the finished triangulation into
  an [`EuclideanComplex`](@ref).

The typical workflow is to assemble the geometry with the first two functions,
call `triangulate` (and optionally `refine!`) from `DelaunayTriangulation.jl`,
and then convert the result with `delaunay_to_simplicial`.

### Defining the Bounding Domain

[`delaunay_points_bnd_rectangle`](@ref) creates a four-vertex rectangular
outer boundary. Its two arguments are vectors `bmin = [xmin, ymin]` and
`bmax = [xmax, ymax]` giving the lower-left and upper-right corners. It
returns a point list and a `bndcurve` in the format expected by the
`boundary_nodes` keyword of `triangulate`:

```julia
using DelaunayTriangulation
using ConleyDynamics

points, bndcurve = delaunay_points_bnd_rectangle([-5, -5], [5, 5])
tri = triangulate(points; boundary_nodes = bndcurve)
sc  = delaunay_to_simplicial(tri)
```

This already produces a valid `EuclideanComplex` — a triangulated filled
rectangle without any internal constraints.

### Adding Constraint Curves

[`delaunay_points_add_segment`](@ref) appends a polygonal curve to the
point list and records successive point pairs as constrained edges. A curve
is passed as a `Vector{Vector{<:Real}}` of `[x, y]` coordinates. If the
first and last entries coincide (within a tolerance of `1e-9`), the curve is
treated as *closed* and no duplicate endpoint is stored; otherwise it is
treated as *open*.

The function has two calling forms. The two-argument form initialises a
fresh segment set and is convenient for the first constraint curve; the
three-argument form extends an existing `Set{Tuple{Int,Int}}`:

```julia
using DelaunayTriangulation
using ConleyDynamics

points, bndcurve = delaunay_points_bnd_rectangle([-5, -5], [5, 5])

# Closed circular constraint (n+1 points with repeated first/last)
n     = 16
theta = range(0, 2π, length = n + 1)
circle = [[3.0 * cos(t), 3.0 * sin(t)] for t in theta]

# Open polygonal cut inside the domain
cut = [[-4.0, -4.0], [0.0, 0.0], [4.0, -4.0]]

points, segments = delaunay_points_add_segment(points, circle)          # 2-arg form
points, segments = delaunay_points_add_segment(points, segments, cut)   # 3-arg form

tri = triangulate(points; boundary_nodes = bndcurve, segments = segments)
sc  = delaunay_to_simplicial(tri)
```

Multiple closed curves can be combined freely; for example, an outer
boundary with two circular holes is set up by calling
`delaunay_points_add_segment` twice with circle data.

### Mesh Refinement

After `triangulate`, `DelaunayTriangulation.jl` provides `refine!` to
improve triangle quality. Passing `max_area` as a fraction of the total
domain area and `min_angle` in degrees is a practical starting point:

```julia
triarea = get_area(tri)
refine!(tri; max_area = 0.001 * triarea, min_angle = 25.0)
sc = delaunay_to_simplicial(tri)
```

The refined triangulation is then converted to a `EuclideanComplex` in the
same way as before.

### Complete Example

The following self-contained snippet combines all three steps — bounding
rectangle, two circular constraints, mesh refinement, and conversion:

```julia
using DelaunayTriangulation
using ConleyDynamics

# Bounding box
points, bndcurve = delaunay_points_bnd_rectangle([-5, -5], [5, 5])

# Outer ring (closed, 20 segments)
n1    = 20
th1   = range(0, 2π, length = n1 + 1)
ring  = [[4.0 * cos(t), 4.0 * sin(t)] for t in th1]
points, segments = delaunay_points_add_segment(points, ring)

# Inner disc (closed, 10 segments)
n2    = 10
th2   = range(0, 2π, length = n2 + 1)
disc  = [[1.5 * cos(t), 1.5 * sin(t)] for t in th2]
points, segments = delaunay_points_add_segment(points, segments, disc)

# Triangulate, refine, and convert
tri = triangulate(points; boundary_nodes = bndcurve, segments = segments)
triarea = get_area(tri)
refine!(tri; max_area = 0.005 * triarea, min_angle = 20.0)
sc = delaunay_to_simplicial(tri)
```
