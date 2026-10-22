# Extension Functions

## DelaunayTriangulation.jl

The following functions are available when `DelaunayTriangulation.jl` is
loaded alongside `ConleyDynamics.jl`.

```@docs
delaunay_points_bnd_rectangle
delaunay_points_add_nodes
delaunay_points_add_segment
delaunay_points_add_split_segs
delaunay_to_simplicial
```

## TetGen.jl

The following functions create tetrahedral meshes when `TetGen.jl` is
loaded alongside `ConleyDynamics.jl`.

```@docs
tetgen_mesh
create_simplicial_delaunay3d
```

The domain descriptions and the mesh results are plain Julia data. The
following types and functions are part of `ConleyDynamics.jl` itself,
and they can be used without loading `TetGen.jl`.

```@docs
TetDomain
TetGenResult
tetgen_domain_box
tetgen_domain_from_surface
tetgen_domain_add_nodes
tetgen_domain_add_surface
tetgen_domain_add_curve
tetgen_domain_add_hole
tetgen_domain_add_region
tetgen_to_simplicial
tetgen_marker_cells
tetgen_region_cells
```

## Gmsh.jl

The following functions are available when `Gmsh.jl` is loaded alongside
`ConleyDynamics.jl`.

```@docs
gmsh_to_simplicial
gmsh_file_to_simplicial
gmsh_marker_cells
```
