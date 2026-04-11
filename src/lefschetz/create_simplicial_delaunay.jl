export create_simplicial_delaunay

"""
    create_simplicial_delaunay(boxw::Real, boxh::Real, pdist::Real, attmpt::Int;
                               p::Int=2, euclidean::Bool=false)

Create a planar Delaunay triangulation inside a box. The complex is
over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The function selects a random sample of points inside the rectangular
box `[0,boxw] x [0,boxh]`, while trying to maintain a minimum distance
of `pdist` between the points. The argument `attmpt` specifies the number
of attempts when trying to add points. A standard value is 20, and larger
values tend to fill holes better, but at the expense of runtime. From the
random sample, the function then creates a Delaunay triangulation.

When `euclidean=false` (default), the function returns:
* A simplicial complex `sc::LefschetzComplex`.
* A vector `coords::Vector{Vector{Float64}}` of vertex coordinates.

When `euclidean=true`, it returns a `EuclideanComplex` with embedded
coordinates.

Note that the function does not provide a full triangulation
of the given rectangle. Close to the boundary there will be gaps.
"""
function create_simplicial_delaunay(boxw::Real, boxh::Real, pdist::Real, attmpt::Int;
                                    p::Int=2, euclidean::Bool=false)
    #
    # Create a planar Delaunay triangulation
    #

    # If negative, set number of attemps to 20

    if attmpt < 1
        attmpt = 20
    end

    # Select a random point sample in the box, create point to index dictionary

    vertices = randompointarray(boxw, boxh, pdist, attempts=attmpt)
    nvert = length(vertices)
    vertex2index = Dict{Point,Int}([vertices[k] => k for k in 1:nvert])

    # Create the vector of coordinates

    coords = Vector{Vector{Float64}}()
    for k = 1:nvert
        push!(coords, [vertices[k][1], vertices[k][2]])
    end

    # Create the labels

    labels = Vector{String}()
    labwidth = Int(ceil(log(0.5+nvert)/log(10)))
    labformat = "%0" * string(labwidth) * "d"
    for k = 1:nvert
        clabel = Printf.format(Printf.Format(labformat), k)
        push!(labels, clabel)
    end

    # Create the triangulation

    triangles = polytriangulate(vertices)
    simplices = Vector{Vector{Int}}()
    for tri in triangles
        csimp = Vector{Int}()
        for pts in tri
            push!(csimp, vertex2index[pts])
        end
        push!(simplices, sort(csimp))
    end

    # Create the simplicial complex and return it

    sc = create_simplicial_complex(labels, simplices, p=p)

    if euclidean
        return create_simplicial_complex(labels, simplices, coords; p=p)
    else
        return sc, coords
    end
end

"""
    create_simplicial_delaunay(boxw::Real, boxh::Real, npoints::Int;
                               p::Int=2, euclidean::Bool=false)

Create a planar Delaunay triangulation inside a box. The complex is
over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The function selects a random sample of `npoints` points inside the rectangular
box `[0,boxw] x [0,boxh]`. From the random sample, the function then creates a
Delaunay triangulation.

When `euclidean=false` (default), the function returns:
* A simplicial complex `sc::LefschetzComplex`.
* A vector `coords::Vector{Vector{Float64}}` of vertex coordinates.

When `euclidean=true`, it returns a `EuclideanComplex` with embedded
coordinates.

Note that the function does not provide a full triangulation
of the given rectangle. Close to the boundary there will be gaps.
"""
function create_simplicial_delaunay(boxw::Real, boxh::Real, npoints::Int;
                                    p::Int=2, euclidean::Bool=false)
    #
    # Create a planar Delaunay triangulation
    #

    # Select a random point sample in the box, create point to index dictionary

    vertices = randompointarray(0.0, 0.0, boxw, boxh, npoints)
    nvert = length(vertices)
    vertex2index = Dict{Point,Int}([vertices[k] => k for k in 1:nvert])

    # Create the vector of coordinates

    coords = Vector{Vector{Float64}}()
    for k = 1:nvert
        push!(coords, [vertices[k][1], vertices[k][2]])
    end

    # Create the labels

    labels = Vector{String}()
    labwidth = Int(ceil(log(0.5+nvert)/log(10)))
    labformat = "%0" * string(labwidth) * "d"
    for k = 1:nvert
        clabel = Printf.format(Printf.Format(labformat), k)
        push!(labels, clabel)
    end

    # Create the triangulation

    triangles = polytriangulate(vertices)
    simplices = Vector{Vector{Int}}()
    for tri in triangles
        csimp = Vector{Int}()
        for pts in tri
            push!(csimp, vertex2index[pts])
        end
        push!(simplices, sort(csimp))
    end

    # Create the simplicial complex and return it

    if euclidean
        return create_simplicial_complex(labels, simplices, coords; p=p)
    else
        sc = create_simplicial_complex(labels, simplices, p=p)
        return sc, coords
    end
end

