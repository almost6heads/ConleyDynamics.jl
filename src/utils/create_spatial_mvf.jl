export create_spatial_mvf

"""
    create_spatial_mvf(lc::LefschetzComplex, coords::Vector{Vector{Float64}}, vf)

Create a spatial multivector field from a regular vector field.

The function expects a three-dimensional Lefschetz complex `lc` and a
coordinate vector `coords` of coordinates for all the 0-dimensional cells
in the complex. Moreover, the underlying vector field is specified by the
function `vf(z::Vector{Float64})::Vector{Float64}`, where both the input
and output vectors have length three. The function `create_spatial_mvf`
returns a multivector field `mvf` on `lc`, which can then be further
analyzed using for example the function `connection_matrix`.

__The following still needs to be adapted to the 3D case!!__

The input data `lc` and `coords` can be generated using one of the following
methods:

* `create_cubical_rectangle`
* `create_simplicial_rectangle`
* `create_simplicial_delaunay`

In each case, the provided coordinate vector can be transformed to the
correct bounding box using `convert_planar_coordinates`.

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
function create_spatial_mvf(lc::LefschetzComplex, coords::Vector{Vector{Float64}}, vf)
    #
    # Create a spatial multivector field from a regular vector field
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

function spatial_mvf_vertex2region(vindex::Int, lc::LefschetzComplex,
                                   coords::Vector{Vector{Float64}}, vf)
    #
    # Find the 3d regions which can be entered from a vertex
    #

    # Determine the adjacent edges and regions

    openhull = lefschetz_openhull(lc, [vindex])
    openhull_dims = lc.dimensions[openhull]

    adedges = openhull[findall(isequal(1),openhull_dims)]
    adregns = openhull[findall(isequal(3),openhull_dims)]

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

        p1, p2, p3 = adregionvertices
        c1, c2, c3 = spatial_barycentric_coords(tc0, coords[p1], coords[p2], coords[p3], c0)

        # Add the region if necessary

        if (c1 >= 0.0) & (c2 >= 0.0) & (c3 >= 0.0)
            push!(targetregions, k)
        end
    end

    # Return the target regions

    return targetregions
end

###############################################################################

function spatial_mvf_edge2region(eindex::Int, lc::LefschetzComplex,
                                 coords::Vector{Vector{Float64}}, vf)
    #
    # Find the 3d regions which can be entered from an edge
    #

    npts = 5   # Look at npts-1 points along the edge to determine the flow

    # Determine the boundary of the edge

    v1, v2 = lefschetz_boundary(lc, eindex)

    # Determine the adjacent faces and regions

    openhull = lefschetz_openhull(lc, [eindex])
    openhull_dims = lc.dimensions[openhull]
    adfaces = openhull[findall(isequal(2),openhull_dims)]
    adregns = openhull[findall(isequal(3),openhull_dims)]

    # Determine the adjacent edges for each edge vertex
    # In each case, remove the edge 'eindex'

    openhull1 = lefschetz_openhull(lc, [v1])
    openhull1_dims = lc.dimensions[openhull1]
    adedges1 = setdiff(openhull1[findall(isequal(1),openhull1_dims)], [eindex])

    openhull2 = lefschetz_openhull(lc, [v2])
    openhull2_dims = lc.dimensions[openhull2]
    adedges2 = setdiff(openhull2[findall(isequal(1),openhull2_dims)], [eindex])

    # Create vectors of adjacent vertices for v1 and v2
    # In each case, remove v2 and v1, respectively

    adverts1 = setdiff(lefschetz_skeleton(lc, adedges1, 0), [v1])
    adverts2 = setdiff(lefschetz_skeleton(lc, adedges2, 0), [v2])

    # Loop through the regions to find all target regions

    targetregions = Vector{Int}()
    cv1 = coords[v1]
    cv2 = coords[v2]

    for r in adregns

        # Determine the faces f1 and f2 of the region which
        # contain the given edge eindex

        f1, f2 = intersect(lefschetz_boundary(lc, r), adfaces)

        # For each of these faces, determine its edges which are
        # adjacent to v1 and v2, but different from eindex

        f1bnd = lefschetz_boundary(lc, f1)
        e11 = first(intersect(f1bnd, adedges1))
        e12 = first(intersect(f1bnd, adedges2))

        f2bnd = lefschetz_boundary(lc, f2)
        e21 = first(intersect(f2bnd, adedges1))
        e22 = first(intersect(f2bnd, adedges2))

        # For each of these edges, determine its boundary vertex
        # which is different from v1 or v2
        # Notation: vXY is vertex on face X opposite to vertex vY

        v11 = first(setdiff(lefschetz_boundary(lc, e11), [v1, v2]))
        v12 = first(setdiff(lefschetz_boundary(lc, e12), [v1, v2]))
        v21 = first(setdiff(lefschetz_boundary(lc, e21), [v1, v2]))
        v22 = first(setdiff(lefschetz_boundary(lc, e22), [v1, v2]))











    end






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

        p1, p2, p3 = adregionvertices
        c1, c2, c3 = spatial_barycentric_coords(tc0, coords[p1], coords[p2], coords[p3], c0)

        # Add the region if necessary

        if (c1 >= 0.0) & (c2 >= 0.0) & (c3 >= 0.0)
            push!(targetregions, k)
        end
    end

    # Return the target regions

    return targetregions






    # Create normal vector to the edge

    pvec = coords[v2] .- coords[v1]
    nvec = [-pvec[2], pvec[1]]

    # Determine sign of the dot products along the edge

    dp_is_p = false
    dp_is_n = false
    dp_is_z = false

    nptsf = Float64(npts)
    p0 = coords[v1]
    for k = 0:npts
        pc = p0 .+ ((k / nptsf) .* pvec)
        pvf = vf(pc)
        dprod = pvf[1]*nvec[1] + pvf[2]*nvec[2]
        if dprod > 0.0
            dp_is_p = true
        elseif dprod < 0.0
            dp_is_n = true
        elseif (k>0) & (k<npts)
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








###############################################################################

function planar_mvf_edge2region(eindex::Int, lc::LefschetzComplex,
                                coords::Vector{Vector{Float64}}, vf)
    #
    # Find the 2d regions which can be entered from an edge
    #

    npts = 4   # Look at npts+1 points along the edge to determine the flow

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
    for k = 0:npts
        pc = p0 .+ ((k / nptsf) .* pvec)
        pvf = vf(pc)
        dprod = pvf[1]*nvec[1] + pvf[2]*nvec[2]
        if dprod > 0.0
            dp_is_p = true
        elseif dprod < 0.0
            dp_is_n = true
        elseif (k>0) & (k<npts)
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

function spatial_barycentric_coords(p::Vector{<:Real},p1::Vector{<:Real},
                                   p2::Vector{<:Real},p3::Vector{<:Real},
                                   p4::Vector{<:Real})
    #
    # Compute the barycentric coordinate in the representation
    #
    #    p = a * (p1-p4) + b * (p2-p4) + c * (p3-p4) + p4

    return [p1-p4;; p2-p4;; p3-p4] \ (p-p4)
end
