export create_spatial_mvf

"""
    create_spatial_mvf(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf)

Create a spatial multivector field from a regular vector field.

The function expects a three-dimensional Lefschetz complex `lc` and a
coordinate vector `coords` of coordinates for all the 0-dimensional cells
in the complex. Moreover, the underlying vector field is specified by the
function `vf(z::Vector{Float64})::Vector{Float64}`, where both the input
and output vectors have length three. The function `create_spatial_mvf`
returns a multivector field `mvf` on `lc`, which can then be further
analyzed using for example the function `connection_matrix`.

The input data `lc` and `coords` has to be of one of the following
two types:

* `lc` is a tetrahedral mesh of a region in three dimensions. In other
  words, the underlying Lefschetz complex is in fact a simplicial 
  complex, and the vector `coords` contains the vertex coordinates.
* `lc` is a three-dimensional cubical complex, i.e., it is the
  closure of a collection of three-dimensional cubes in space.
  The vertex coordinates can be slightly perturbed from the original
  position in the cubical lattice, as long as the overall structure 
  of the complex stays intact. In that case, the faces are interpreted
  as Bezier surfaces with straight edges.

# Example 1

Suppose we define a sample vector field using the commands

```julia
function samplevf(x::Vector{Float64})
    #
    # Sample vector field with nontrivial Morse decomposition
    #
    x1, x2, x3 = x
    y1 = x1 * (1.0 - x1*x1)
    y2 = -x2
    y3 = -x3
    return [y1, y2, y3]
end
```

One first creates a cubical complex covering the interesting dynamics,
say the trapping region `[-1.5,1.5] x [-1,1] x [-1,1]`, using the
commands

```julia
lc, coords = create_cubical_box(3,3,3);
coordsN = convert_spatial_coordinates(coords,[-1.5,-1.0,-1.0],[1.5,1.0,1.0]);
```

The multivector field is then generated using

```julia
mvf = create_spatial_mvf(lc,coordsN,samplevf);
```

and the commands

```julia
cm = connection_matrix(lc, mvf);
cm.conley
full_from_sparse(cm.matrix)
```

finally show that this vector field gives rise to a Morse decomposition
with three Morse sets, and two connecting orbits.
"""
function create_spatial_mvf(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf)
    #
    # Create a spatial multivector field from a regular vector field
    #

    # Find the indices for the vertices, the edges, and the faces

    vi = findall(isequal(0), lc.dimensions)
    ei = findall(isequal(1), lc.dimensions)
    fi = findall(isequal(2), lc.dimensions)

    # Find the vertex to region map

    v2r = [Vector{Int}() for _ in 1:length(vi)]

    Threads.@threads for k = 1:length(vi)
        v2r[k] = spatial_mvf_vertex2region(vi[k],lc,coords,vf)
    end

    # Find the edge to region map

    e2r = [Vector{Int}() for _ in 1:length(ei)]

    Threads.@threads for k = 1:length(ei)
        e2r[k] = spatial_mvf_edge2region(ei[k],lc,coords,vf)
    end

    # Find the face to region map

    f2r = [Vector{Int}() for _ in 1:length(fi)]

    Threads.@threads for k = 1:length(fi)
        f2r[k] = spatial_mvf_face2region(fi[k],lc,coords,vf)
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

    for k in 1:length(f2r)
        if length(f2r[k]) > 0
            push!(mvfbase, vcat([fi[k]], f2r[k]))
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
                                   coords::Vector{<:Vector{<:Real}}, vf)
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

    for k in adregns

        # Determine the vertices of the region

        regionvertices = lefschetz_skeleton(lc, [k], 0)
        adregionvertices = intersect(regionvertices, adverts)

        # Determine the coordinates with respect to the special basis

        p1, p2, p3 = adregionvertices
        c1, c2, c3 = spatial_coords(vfc0, coords[p1].-c0, coords[p2].-c0, coords[p3].-c0)

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
                                 coords::Vector{<:Vector{<:Real}}, vf)
    #
    # Find the 3d regions which can be entered from an edge
    #

    npts = 3   # Look at npts points along the edge interior to determine the flow
    nptsf = Float64(npts)

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

        # Create the coordinate vectors for these points

        cv11 = coords[v11]
        cv12 = coords[v12]
        cv21 = coords[v21]
        cv22 = coords[v22]

        # Move along the edge until first potential flow into region r

        entered_region = false
        edgedir = cv2 .- cv1
        k = 1

        while (k <= npts) & (!entered_region)

            # Find the coordinate of the point along the edge

            t     = k / (nptsf + 1.0)
            ccur  = (1-t) .* cv1 .+ t .* cv2
            vfcur = vf(ccur)

            # Find the tangent vectors to the two faces

            cf1 = (1-t) .* (cv11 .- cv1) .+ t .* (cv12 - cv2)
            cf2 = (1-t) .* (cv21 .- cv1) .+ t .* (cv22 - cv2)

            # Determine the coordinates with respect to the special basis

            c1, c2, c3 = spatial_coords(vfcur, cf1, cf2, edgedir)

            # Determine whether the flow enters

            if (c1 >= 0.0) & (c2 >= 0.0)
                push!(targetregions, r)
                entered_region = true
            end

            # Increment the counter

            k = k + 1
        end
    end

    # Return the target regions

    return targetregions
end

###############################################################################

function spatial_mvf_face2region(findex::Int, lc::LefschetzComplex,
                                 coords::Vector{<:Vector{<:Real}}, vf)
    #
    # Find the 3d regions which can be entered from a face
    #

    npts = 3   # Look at npts points along the edge interior to determine the flow
    nptsf = Float64(npts)
    targetregions = Vector{Int}()

    # Determine the edges and vertices in the boundary

    faceclosure = lefschetz_closure(lc, [findex])
    faceclosure_dims = lc.dimensions[faceclosure]
    facevertices = faceclosure[findall(isequal(0),faceclosure_dims)]
    faceedges    = faceclosure[findall(isequal(1),faceclosure_dims)]

    # Distinguish between triangles and Bezier quadrilaterals
    
    if length(facevertices) == 3

        # Flow through a triangular face
        
        v1, v2, v3 = facevertices
        cv1 = coords[v1]
        cv2 = coords[v2]
        cv3 = coords[v3]

        ce12 = cv2 .- cv1
        ce13 = cv3 .- cv1
        ce23 = cv3 .- cv2
        facenormal = cross(ce12, ce13)

        # Loop through the cofaces

        adregions = lefschetz_coboundary(lc, findex)

        for r in adregions

            # Determine the fourth vertex

            v4 = first(setdiff(lefschetz_skeleton(lc, [r], 0), facevertices))
            cv4 = coords[v4]

            # Determine the normal vector

            if dot(cv4.-cv1, facenormal) > 0
                fn = facenormal
            else
                fn = -facenormal
            end

            # Check for flow directions and add region if it is entered

            entered_region = false
            m1 = 1
            while (m1 <= npts) & (!entered_region)
                t1 = m1 / (1.0 + nptsf)
                m2 = 1
                while (m2 <= m1) & (!entered_region)
                    t2 = m2 / (1.0 + nptsf)
                    cpt = cv1 .+ t1 .* ce12 .+ t2 .* ce23
                    cvf = vf(cpt)
                    if dot(cvf, fn) >= 0
                        push!(targetregions, r)
                        entered_region = true
                    end
                    m2 = m2 + 1
                end
                m1 = m1 + 1
            end
        end

    elseif length(facevertices) == 4

        # Flow through Bezier quadrilateral

        # Determine the vertices so that the edges are v1v2, v2v4, v3v4, v1v3

        v1 = facevertices[1]
        adedges = intersect(lefschetz_coboundary(lc, v1), faceedges)
        v2 = first(setdiff(lefschetz_boundary(lc, adedges[1]), [v1]))
        v3 = first(setdiff(lefschetz_boundary(lc, adedges[2]), [v1]))
        v4 = first(setdiff(facevertices, [v1,v2,v3]))

        cv1 = coords[v1]
        cv2 = coords[v2]
        cv3 = coords[v3]
        cv4 = coords[v4]

        # Precompute the edge vectors

        ce12 = cv2 .- cv1
        ce34 = cv4 .- cv3
        ce13 = cv3 .- cv1
        ce24 = cv4 .- cv2

        # Precompute the face normals at the corners

        fnormal1 = cross(ce12, ce13)
        fnormal2 = cross(ce12, ce24)
        fnormal3 = cross(ce13, ce34)
        fnormal4 = cross(ce24, ce34)

        # Loop through the cofaces

        adregions = lefschetz_coboundary(lc, findex)

        for r in adregions

            # Determine the remaining four vertices
            
            topvertices = setdiff(lefschetz_skeleton(lc, [r], 0), facevertices)
            v1b = first(intersect(lefschetz_skeleton(lc,lefschetz_coboundary(lc,v1),0),topvertices))
            v2b = first(intersect(lefschetz_skeleton(lc,lefschetz_coboundary(lc,v2),0),topvertices))
            v3b = first(intersect(lefschetz_skeleton(lc,lefschetz_coboundary(lc,v3),0),topvertices))
            v4b = first(intersect(lefschetz_skeleton(lc,lefschetz_coboundary(lc,v4),0),topvertices))

            # Find the correct normal vectors at the vertices

            if dot(coords[v1b].-cv1, fnormal1) > 0
                fn1 = fnormal1
            else
                fn1 = -fnormal1
            end
            if dot(coords[v2b].-cv2, fnormal2) > 0
                fn2 = fnormal2
            else
                fn2 = -fnormal2
            end
            if dot(coords[v3b].-cv3, fnormal3) > 0
                fn3 = fnormal3
            else
                fn3 = -fnormal3
            end
            if dot(coords[v4b].-cv4, fnormal4) > 0
                fn4 = fnormal4
            else
                fn4 = -fnormal4
            end

            # Check for flow directions and add region if it is entered

            entered_region = false
            m1 = 1
            while (m1 <= npts) & (!entered_region)
                t1 = m1 / (1.0 + nptsf)
                m2 = 1
                while (m2 <= npts) & (!entered_region)
                    t2 = m2 / (1.0 + nptsf)
                    cpt = (1-t1).*(1-t2).*cv1 .+ t1.*(1-t2).*cv2 .+ (1-t1).*t2.*cv3 .+ t1.*t2.*cv4
                    cfn = (1-t1).*(1-t2).*fn1 .+ t1.*(1-t2).*fn2 .+ (1-t1).*t2.*fn3 .+ t1.*t2.*fn4
                    cvf = vf(cpt)
                    if dot(cvf, cfn) >= 0
                        push!(targetregions, r)
                        entered_region = true
                    end
                    m2 = m2 + 1
                end
                m1 = m1 + 1
            end
        end

    else

        # This case it not covered by the function

        error("Faces have to be triangles or Bezier quadrilaterals!")
    end

    # Return the target regions

    return targetregions
end

###############################################################################

function spatial_coords(p::Vector{<:Real},p1::Vector{<:Real},
                       p2::Vector{<:Real},p3::Vector{<:Real})
    #
    # Compute the spatial coordinates in the representation
    #
    #    p = a * p1 + b * p2 + c * p3

    return [p1;; p2;; p3] \ p
end

