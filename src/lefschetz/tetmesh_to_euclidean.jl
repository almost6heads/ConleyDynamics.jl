export tetmesh_to_euclidean

"""
    tetmesh_to_euclidean(points::AbstractVector, cells::AbstractVector;
                         p::Int=2, droporphans::Bool=true) -> EuclideanComplex

Convert a simplicial mesh, given by vertex coordinates and maximal
simplices, into an `EuclideanComplex`. The complex is over the rationals
if `p=0`, and over `GF(p)` if `p>0`.

The vector `points` contains the vertex coordinates, where `points[k]`
can be any vector or tuple of real numbers. Usually these are points in
three dimensions, but any dimension is allowed. The vector `cells`
contains the maximal simplices of the mesh, as vectors of 1-based
indices into `points`. Typically these are the tetrahedra of a
tetrahedral mesh, or the triangles of a triangulated surface, but
simplices of different dimensions can be mixed. All faces of the
given simplices are generated automatically.

If `droporphans=true` (default), vertices which do not belong to any
cell are removed, and the remaining vertices are renumbered in their
original order. Mesh generators frequently leave such unused vertices,
for example geometry points outside of the meshed region. If
`droporphans=false`, unused vertices become isolated vertices of the
complex.

Vertex labels are zero-padded decimal strings (`"001"`, `"002"`, …)
whose width is determined by the number of vertices, and the label of
a higher-dimensional simplex is the concatenation of its sorted vertex
labels. This is the same convention as in `create_simplicial_delaunay`
and `delaunay_to_simplicial`.

This function can be used to import tetrahedral or surface meshes from
any source, such as mesh files or other Julia packages.

# Example

```julia
points = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0],
          [0.0,0.0,1.0], [1.0,1.0,1.0]]
tets   = [[1,2,3,4], [2,3,4,5]]
ec     = tetmesh_to_euclidean(points, tets)
```
"""
function tetmesh_to_euclidean(points::AbstractVector, cells::AbstractVector;
                              p::Int=2, droporphans::Bool=true)
    #
    # Convert a simplicial mesh into an EuclideanComplex
    #

    npts = length(points)
    if npts == 0
        throw(ArgumentError("tetmesh_to_euclidean: the point list is empty"))
    end
    if isempty(cells)
        throw(ArgumentError("tetmesh_to_euclidean: the cell list is empty"))
    end

    # Convert the coordinates and check that they all have the same length

    coords = [Float64.(collect(pt)) for pt in points]
    pdim = length(coords[1])
    if any(length(c) != pdim for c in coords)
        throw(ArgumentError("tetmesh_to_euclidean: all points need to have " *
                            "the same number of coordinates"))
    end
    if !all(all(isfinite, c) for c in coords)
        throw(ArgumentError("tetmesh_to_euclidean: all coordinates need to be finite"))
    end

    # Check the cells

    for cell in cells
        if isempty(cell)
            throw(ArgumentError("tetmesh_to_euclidean: cells cannot be empty"))
        end
        if any(v -> (v < 1) || (v > npts), cell)
            throw(ArgumentError("tetmesh_to_euclidean: cell $(collect(cell)) " *
                                "contains a vertex index outside of 1:$npts"))
        end
        if !allunique(cell)
            throw(ArgumentError("tetmesh_to_euclidean: cell $(collect(cell)) " *
                                "contains a repeated vertex"))
        end
    end

    # Determine the used vertices and renumber them if necessary

    if droporphans
        used = falses(npts)
        for cell in cells
            used[cell] .= true
        end
        newindex = cumsum(used)
        vcoords  = coords[used]
        simplices = [sort([newindex[v] for v in cell]) for cell in cells]
    else
        vcoords   = coords
        simplices = [sort(Int.(collect(cell))) for cell in cells]
    end

    # Create the labels and the complex

    labels = _vertex_labels(length(vcoords))

    return create_simplicial_complex(labels, simplices, vcoords; p=p)
end

"""
    tetmesh_to_euclidean(points::AbstractMatrix{<:Real},
                         cells::AbstractMatrix{<:Integer};
                         p::Int=2, droporphans::Bool=true) -> EuclideanComplex

Convert a simplicial mesh given in matrix form into an `EuclideanComplex`.

Every column of `points` contains the coordinates of one vertex, and
every column of `cells` contains the 1-based vertex indices of one
maximal simplex. For a tetrahedral mesh in three dimensions, `points`
is a `3 x n` and `cells` a `4 x m` matrix. This is the format used by
many mesh generators. See the vector version of `tetmesh_to_euclidean`
for a description of the keyword arguments.
"""
function tetmesh_to_euclidean(points::AbstractMatrix{<:Real},
                              cells::AbstractMatrix{<:Integer};
                              p::Int=2, droporphans::Bool=true)
    #
    # Convert a simplicial mesh in matrix form into an EuclideanComplex
    #
    pointvec = [Vector{Float64}(points[:,k]) for k in 1:size(points,2)]
    cellvec  = [Vector{Int}(cells[:,k]) for k in 1:size(cells,2)]
    return tetmesh_to_euclidean(pointvec, cellvec; p=p, droporphans=droporphans)
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function _vertex_labels(nvert::Int)
    #
    # Create zero-padded decimal vertex labels of equal width
    #
    labwidth = Int(ceil(log(0.5+nvert)/log(10)))
    return [lpad(k, labwidth, '0') for k in 1:nvert]
end
