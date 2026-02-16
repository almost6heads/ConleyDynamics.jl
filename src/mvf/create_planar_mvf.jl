export create_planar_mvf

"""
    create_planar_mvf(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf)

Create a planar multivector field from a regular vector field.

The function expects a planar Lefschetz complex `lc` and a coordinate
vector `coords` of coordinates for all the 0-dimensional cells in the
complex. Moreover, the underlying vector field is specified by the
function `vf(z::Vector{Float64})::Vector{Float64}`, where both the input
and output vectors have length two. The function `create_planar_mvf` returns
a multivector field `mvf` on `lc`, which can then be further analyzed
using for example the function `connection_matrix`.

The input data `lc` and `coords` can be generated using one of the following
methods:

* `create_cubical_rectangle`
* `create_simplicial_rectangle`
* `create_simplicial_delaunay`

In each case, the provided coordinate vector can be transformed
to the correct bounding box values using the conversion function
`convert_planar_coordinates`.

# Example 1

Suppose we define a sample vector field using the commands

```julia
function samplevf(x::Vector{Float64})
    #
    # Sample vector field with nontrivial Morse decomposition
    #
    x1, x2 = x
    y1 = x1 * (1.0 - x1*x1 - 3.0*x2*x2)
    y2 = x2 * (1.0 - 3.0*x1*x1 - x2*x2)
    return [y1, y2]
end
```

One first creates a triangulation of the enclosing box, which in
this case is given by `[-2,2] x [-2,2]` using the commands

```julia
n = 21
lc, coords = create_simplicial_rectangle(n,n);
coordsN = convert_planar_coordinates(coords,[-2.0,-2.0],[2.0,2.0]);
```

The multivector field is then generated using

```julia
mvf = create_planar_mvf(lc,coordsN,samplevf);
```

and the commands

```julia
cm = connection_matrix(lc, mvf);
cm.conley
full_from_sparse(cm.matrix)
```

finally show that this vector field gives rise to a Morse decomposition
with nine Morse sets, and twelve connecting orbits. Using the commands

```julia
fname = "morse_test.pdf"
plot_planar_simplicial_morse(lc, coordsN, fname, cm.morse, pv=true)
```

these Morse sets can be visualized. The image will be saved in `fname`.

# Example 2

An example with periodic orbits can be generated using the
vector field

```julia
function samplevf2(x::Vector{Float64})
    #
    # Sample vector field with nontrivial Morse decomposition
    #
    x1, x2 = x
    c0 = x1*x1 + x2*x2
    c1 = (c0 - 4.0) * (c0 - 1.0)
    y1 = -x2 + x1 * c1
    y2 =  x1 + x2 * c1
    return [-y1, -y2]
end
```

The Morse decomposition can now be computed via

```julia
n2 = 51
lc2, coords2 = create_cubical_rectangle(n2,n2);
coords2N = convert_planar_coordinates(coords2,[-4.0,-4.0],[4.0,4.0]);
mvf2 = create_planar_mvf(lc2,coords2N,samplevf2);
cm2 = connection_matrix(lc2, mvf2);
cm2.conley
cm2.poset
full_from_sparse(cm2.matrix)

fname2 = "morse_test2.pdf"
plot_planar_cubical_morse(lc2, fname2, cm2.morse, pv=true)
```

In this case, one obtains three Morse sets: One is a stable equilibrium,
one is an unstable periodic orbit, and the last is a stable periodic orbit.
"""
function create_planar_mvf(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf)
    #
    # Create a planar multivector field from a regular vector field
    #

    # Find the indices for the vertices and the edges

    vi = findall(isequal(0), lc.dimensions)
    ei = findall(isequal(1), lc.dimensions)

    # Find the vertex to region map

    v2r = [Vector{Int}() for _ in 1:length(vi)]

    Threads.@threads for k = 1:length(vi)
        v2r[k] = planar_mvf_vertex2region(vi[k],lc,coords,vf)
    end

    # Find the edge to region map

    e2r = [Vector{Int}() for _ in 1:length(ei)]

    Threads.@threads for k = 1:length(ei)
        e2r[k] = planar_mvf_edge2region(ei[k],lc,coords,vf)
    end

    # Create the multivector field base

    mvfbase = Vector{Vector{Int}}()

    for k in 1:length(v2r)
        if length(v2r[k]) > 0
            push!(mvfbase, vcat([vi[k]], v2r[k]))
        end
    end

    for k in 1:length(e2r)
        if length(e2r[k]) > 0
            push!(mvfbase, vcat([ei[k]], e2r[k]))
        end
    end

    # Create the multivector field

    mvf = create_mvf_hull(lc, mvfbase)
    return mvf
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function planar_mvf_vertex2region(vindex::Int, lc::LefschetzComplex,
                                  coords::Vector{<:Vector{<:Real}}, vf)
    #
    # Find the 2d regions which can be entered from a vertex
    #

    # Determine the adjacent edges and regions

    openhull = lefschetz_openhull(lc, [vindex])
    openhull_dims = lc.dimensions[openhull]

    adedges = openhull[findall(isequal(1),openhull_dims)]
    adregns = openhull[findall(isequal(2),openhull_dims)]

    # Create a vector of adjacent vertices

    adverts = setdiff(lefschetz_skeleton(lc, adedges, 0), [vindex])

    # Loop through the regions, compute the barycentric
    # coordinates, and add the relevant regions

    targetregions = Vector{Int}()
    c0 = coords[vindex]
    vfc0 = vf(c0)
    tc0 = c0 .+ vfc0

    for k in adregns

        # Determine the vertices of the region

        regionvertices = lefschetz_skeleton(lc, [k], 0)
        adregionvertices = intersect(regionvertices, adverts)

        # Determine the barycentric coordinates

        p1, p2 = adregionvertices
        c1, c2 = planar_barycentric_coords(tc0, coords[p1], coords[p2], c0)

        # Add the face if necessary

        if (c1 >= 0.0) & (c2 >= 0.0)
            push!(targetregions, k)
        end
    end

    # Return the target regions

    return targetregions
end

###############################################################################

function planar_mvf_edge2region(eindex::Int, lc::LefschetzComplex,
                                coords::Vector{<:Vector{<:Real}}, vf)
    #
    # Find the 2d regions which can be entered from an edge
    #

    npts = 10   # Look at npts points along the edge interior to determine the flow

    # Determine the boundary of the edge

    v1, v2 = lefschetz_boundary(lc, eindex)

    # Determine the adjacent region or regions

    adregions = lefschetz_coboundary(lc, eindex)

    # Create normal vector to the edge

    pvec = coords[v2] .- coords[v1]
    nvec = [-pvec[2], pvec[1]]

    # Determine sign of the dot products along the edge

    dp_is_p = false
    dp_is_n = false
    dp_is_z = false

    nptsf = Float64(npts)
    p0 = coords[v1]
    for k = 1:npts
        pc = p0 .+ ((k / (nptsf+1.0)) .* pvec)
        pvf = vf(pc)
        dprod = pvf[1]*nvec[1] + pvf[2]*nvec[2]
        if dprod > 0.0
            dp_is_p = true
        elseif dprod < 0.0
            dp_is_n = true
        else
            dp_is_z = true
        end
    end

    # Decide whether an adjacent region could be a target

    targetregions = Vector{Int}()
    midpoint = 0.5 .* coords[v1] .+ 0.5 .* coords[v2]

    for k in adregions

        # Determine the barycenter of the region

        regionvertices = lefschetz_skeleton(lc, [k], 0)
        nv = 0.0
        bpoint = [0.0, 0.0]
        for m in regionvertices
            bpoint = bpoint .+ coords[m]
            nv = nv + 1.0
        end
        bpoint = bpoint / nv

        # Determine the dot product of the segment from the
        # midpoint to the barycenter with the normal vector

        dprod = (bpoint[1] - midpoint[1]) * nvec[1] +
                (bpoint[2] - midpoint[2]) * nvec[2]

        # Add the region, if necessary

        if (dprod > 0.0) & (dp_is_p | dp_is_z)
            push!(targetregions, k)
        end

        if (dprod < 0.0) & (dp_is_n | dp_is_z)
            push!(targetregions, k)
        end
    end
    
    # Return the target regions

    return targetregions
end

###############################################################################

function planar_barycentric_coords(p::Vector{<:Real},p1::Vector{<:Real},
                                   p2::Vector{<:Real},p3::Vector{<:Real})
    #
    # Compute the barycentric coordinate in the representation
    #
    #    p = a * (p1-p3) + b * (p2-p3) + p3

    x,  y  = p
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) /
        ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3))
    b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) /
        ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3))

    return a, b
end

