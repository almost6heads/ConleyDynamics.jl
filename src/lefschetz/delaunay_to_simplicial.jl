export delaunay_to_simplicial

"""
    delaunay_to_simplicial(tt; p::Int=2) -> EuclideanComplex

Convert a Delaunay triangulation from `DelaunayTriangulation.jl` into
an `EuclideanComplex`.

The argument `tt` is a triangulation object as returned by `triangulate`
from `DelaunayTriangulation.jl`. Ghost triangles and ghost vertices (the
package's boundary sentinel elements) are automatically excluded via the
`each_solid_*` iterators, so the result contains only the geometric simplices.

The keyword argument `p` specifies the field characteristic of the resulting
complex: `p=0` gives the rationals, `p>0` gives GF(p). The default is `p=2`.

Vertex labels are zero-padded decimal strings (`"001"`, `"002"`, …) whose
width is determined by the total vertex count. The resulting `EuclideanComplex`
carries the embedded 2D coordinates of each simplex.

This function requires `DelaunayTriangulation.jl` to be loaded before
`ConleyDynamics.jl`.

# Example

```julia
using DelaunayTriangulation
using ConleyDynamics

pts = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
tri = triangulate(pts)
sc  = delaunay_to_simplicial(tri)
sc3 = delaunay_to_simplicial(tri; p=3)
```
"""
function delaunay_to_simplicial end
