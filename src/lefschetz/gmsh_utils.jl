export gmsh_to_simplicial, gmsh_file_to_simplicial, gmsh_marker_cells

"""
    gmsh_to_simplicial(; dim::Int=3, p::Int=2,
                       physical::Union{Tuple{Integer,Integer},Nothing}=nothing)
                       -> EuclideanComplex

Convert the mesh of the current `Gmsh.jl` model into an `EuclideanComplex`.
The complex is over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The function reads the mesh elements of dimension `dim` from the current
model, which has to be meshed already, for example by calling
`gmsh.model.mesh.generate(3)`. For `dim=3`, these are the tetrahedra, and
the result is a tetrahedral mesh. For `dim=2`, the triangles are used,
which gives a surface complex embedded in three dimensions, and for
`dim=1`, the edges are used, which gives a graph. All faces of
these elements are added automatically. If `physical=(pdim, ptag)` is
given, only the elements of the physical group with dimension `pdim` and
tag `ptag` are used, and `pdim` has to equal `dim`.

Only linear simplices are supported. If the mesh contains other elements
of dimension `dim`, such as hexahedra, prisms, or second order
tetrahedra, an `ArgumentError` is raised. Nodes which do not belong to
any of the used elements are dropped. Vertex labels are zero-padded
decimal strings as in `tetmesh_to_euclidean`.

If two distinct nodes of the used elements have the same coordinates,
a warning is printed. This usually means that touching or intersecting
entities were meshed separately, without a `gmsh.model.occ.fragment`
operation, so that the complex is not connected across them.

The function only reads the current model. It neither initializes nor
finalizes `Gmsh`, so the model can be used afterwards, for example to
find the cells of physical groups with `gmsh_marker_cells`.

This function requires `Gmsh.jl` to be loaded.

# Example

```julia
using Gmsh: gmsh
using ConleyDynamics

gmsh.initialize()
gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
gmsh.model.mesh.generate(3)
ec = gmsh_to_simplicial()
gmsh.finalize()
```
"""
function gmsh_to_simplicial end

"""
    gmsh_file_to_simplicial(filename::AbstractString; dim::Int=3, p::Int=2)
                            -> EuclideanComplex

Read a mesh or geometry file with `Gmsh.jl` and convert it into an
`EuclideanComplex`. The complex is over the rationals if `p=0`, and over
`GF(p)` if `p>0`.

The file can be in any format which `Gmsh` can open, such as a mesh file
(`.msh`) or a geometry file (`.geo`). If the file does not contain mesh
elements of dimension `dim`, a mesh of this dimension is generated with
the mesh options specified in the file. The mesh is then converted with
`gmsh_to_simplicial`.

If `Gmsh` is not initialized, the function initializes it and finalizes
it at the end. Otherwise, the file is opened in the running `Gmsh`
session. This creates a new model, which becomes the current model,
while the previous models remain available.

This function requires `Gmsh.jl` to be loaded.
"""
function gmsh_file_to_simplicial end

"""
    gmsh_marker_cells(ec::EuclideanComplex, dim::Integer, tag::Integer;
                      physical::Bool=true, closure::Bool=true) -> Vector{Int}

Find the cells of a complex which belong to a physical group or an
entity of the current `Gmsh.jl` model.

The complex `ec` has to be created from the current model, for example
with `gmsh_to_simplicial`. If `physical=true` (default), the function
returns the indices of the cells of `ec` formed by the mesh elements of
the physical group with dimension `dim` and tag `tag`. If
`physical=false`, the elements of the model entity with dimension `dim`
and tag `tag` are used instead. Depending on `dim`, these are vertices,
edges, triangles, or tetrahedra. Elements which are not cells of `ec`
are ignored. If `closure=true` (default), all faces of these cells are
included as well, so the result is a closed subcomplex.

Physical groups are the most convenient way to keep track of geometry
after boolean operations such as `gmsh.model.occ.fragment`, since these
operations renumber the entities. The cells of `ec` are identified by
the coordinates of their vertices, so the numbering of the nodes does
not matter.

This function requires `Gmsh.jl` to be loaded.
"""
function gmsh_marker_cells end
