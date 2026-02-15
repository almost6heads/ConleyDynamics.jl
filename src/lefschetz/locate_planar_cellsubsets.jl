export locate_planar_cellsubsets
export cellsubset_location_rectangle
export cellsubset_location_circle
export cellsubset_bounding_box
export signed_distance_rectangle
export signed_distance_circle
export segment_intersects_rectangle
export segment_intersects_circle

"""
    locate_planar_cellsubsets(lc, coords, csubsets, rmin, rmax)

Locate cell subsets relative to a planar rectangle.

For the Lefschetz complex `lc::LefschetzComplex`, whose vertices
have the coordinates given in `coords::Vector{Vector{Real}}`,
and the cell subsets in `csubsets::CellSubsets` this function
extracts the indices of all cell subset closures which lie in the
interior, or intersect the boundary of the rectangle specified
by the minimal and maximal corners `rmin::Vector{Real}`
and `rmax::Vector{Real}`, respectively. The function returns
* `indexI::Vector{Int}`: indices of cell subset closures inside
  the rectangle,
* `indexB::Vector{Int}`: indices of cell subset closures which
  intersect the rectangle boundary.
"""
function locate_planar_cellsubsets(lc::LefschetzComplex,
                                   coords::Vector{Vector{T}},
                                   csubsets::CellSubsets,
                                   rmin::Vector{<:Real},
                                   rmax::Vector{<:Real}) where T<:Real
    #
    # Locate cell subsets relative to a planar rectangle
    #
    Nsub = length(csubsets)
    lindex = fill(0,Nsub)

    Threads.@threads for k in 1:Nsub
        lindex[k] = cellsubset_location_rectangle(lc, coords,
                                       csubsets[k], rmin, rmax)
    end

    indexI = findall(x -> x==1, lindex)
    indexB = findall(x -> x==2, lindex)
    return indexI, indexB
end

"""
    locate_planar_cellsubsets(lc, coords, csubsets, c, r)

Locate cell subsets relative to a planar circle.

For the Lefschetz complex `lc::LefschetzComplex`, whose vertices
have the coordinates given in `coords::Vector{Vector{Real}}`,
and the cell subsets in `csubsets::CellSubsets` this function
extracts the indices of all cell subset closures which lie in the
interior, or intersect the boundary of the circle specified by the
center point `c::Vector{Real}` and radius `r::Real`. The function
returns
* `indexI::Vector{Int}`: indices of cell subset closures inside
  the circle,
* `indexB::Vector{Int}`: indices of cell subset closures which
  intersect the circle.
"""
function locate_planar_cellsubsets(lc::LefschetzComplex,
                                   coords::Vector{Vector{T}},
                                   csubsets::CellSubsets,
                                   c::Vector{<:Real},
                                   r::Real) where T<:Real
    #
    # Locate cell subsets relative to a planar circle
    #
    Nsub = length(csubsets)
    lindex = fill(0,Nsub)

    Threads.@threads for k in 1:Nsub
        lindex[k] = cellsubset_location_circle(lc, coords,
                                       csubsets[k], c, r)
    end

    indexI = findall(x -> x==1, lindex)
    indexB = findall(x -> x==2, lindex)
    return indexI, indexB
end

"""
    cellsubset_location_rectangle(lc, coords, csubset, rmin, rmax)

Determine the location of a cell subset relative to a rectangle.

For the Lefschetz complex `lc::LefschetzComplex`, whose vertices
have the coordinates given in `coords::Vector{Vector{Real}}`,
this function determines the location of the closure of the
cellsubset given in `csubset::Cells` relative to the rectangle
specified by the minimal and maximal corners `rmin::Vector{Real}`
and `rmax::Vector{Real}`, respectively. The function returns
* 1 if the set lies in the interior of the rectangle,
* 2 of the set intersects the rectangle boundary, and
* 3 if the set lies in the exterior of the rectangle.
"""
function cellsubset_location_rectangle(lc::LefschetzComplex,
                                       coords::Vector{Vector{T}},
                                       csubset::Cells,
                                       rmin::Vector{<:Real},
                                       rmax::Vector{<:Real}) where T<:Real
    #
    # Determine the location of a cell subset relative to a rectangle
    #
    @assert length(coords[1])==2 "The complex has to be planar!"

    # Convert to index format if csubset is a string vector
    if typeof(csubset) == Vector{String}
        csubsetI = convert_cells(lc, csubset)
    else
        csubsetI = csubset
    end

    # Extract the vertices in the closure
    vindices = lefschetz_skeleton(lc, csubsetI, 0)

    # Compute the signed distances of the vertices
    sd = [signed_distance_rectangle(coords[k],rmin,rmax) for k in vindices]
    sdmin = minimum(sd)
    sdmax = maximum(sd)
    if sdmax < 0
        return 1
    end
    if (sdmin <= 0) && (sdmax >= 0)
        return 2
    end

    # All points lie outside, but does the bounding box?
    bmin, bmax = cellsubset_bounding_box(lc, coords, csubsetI)
    if (bmin[1] > rmax[1]) || (bmax[1] < rmin[1])
        return 3
    end
    if (bmin[2] > rmax[2]) || (bmax[2] < rmin[2])
        return 3
    end

    # All points lie outside, the bounding box intersects,
    # so I guess we have to check the edges individually
    eindices = lefschetz_skeleton(lc, csubsetI, 1)
    if length(eindices)==0
        return 3
    end

    for eind in eindices
        vindices = lefschetz_boundary(lc, eind)
        p1 = coords[vindices[1]]
        p2 = coords[vindices[2]]
        if segment_intersects_rectangle(p1,p2,rmin,rmax)
            return 2
        end
    end

    return 3
end

"""
    cellsubset_location_circle(lc, coords, csubset, c, r)

Determine the location of a cell subset relative to a circle.

For the Lefschetz complex `lc::LefschetzComplex`, whose vertices
have the coordinates given in `coords::Vector{Vector{Real}}`, this
function determines the location of the closure of the cellsubset
given in `csubset::Cells` relative to the circle specified by the
center point `c::Vector{Real}` and the radius `r::Real`.
The function returns
* 1 if the set lies in the interior of the circle,
* 2 of the set intersects the circle, and
* 3 if the set lies in the exterior of the circle.
"""
function cellsubset_location_circle(lc::LefschetzComplex,
                                    coords::Vector{Vector{T}},
                                    csubset::Cells,
                                    c::Vector{<:Real},
                                    r::Real) where T<:Real
    #
    # Determine the location of a cell subset relative to a circle
    #
    @assert length(coords[1])==2 "The complex has to be planar!"

    # Convert to index format if csubset is a string vector
    if typeof(csubset) == Vector{String}
        csubsetI = convert_cells(lc, csubset)
    else
        csubsetI = csubset
    end

    # Extract the vertices in the closure
    vindices = lefschetz_skeleton(lc, csubsetI, 0)

    # Compute the signed distances of the vertices
    sd = [signed_distance_circle(coords[k],c,r) for k in vindices]
    sdmin = minimum(sd)
    sdmax = maximum(sd)
    if sdmax < 0
        return 1
    end
    if (sdmin <= 0) && (sdmax >= 0)
        return 2
    end

    # All points lie outside, but does the bounding box?
    bmin, bmax = cellsubset_bounding_box(lc, coords, csubsetI)
    if (bmin[1] > c[1]+r) || (bmax[1] < c[1]-r)
        return 3
    end
    if (bmin[2] > c[2]+r) || (bmax[2] < c[2]-r)
        return 3
    end

    # All points lie outside, the bounding box intersects
    # the box determined by the circle, so I guess we have
    # to check the edges individually
    eindices = lefschetz_skeleton(lc, csubsetI, 1)
    if length(eindices)==0
        return 3
    end

    for eind in eindices
        vindices = lefschetz_boundary(lc, eind)
        p1 = coords[vindices[1]]
        p2 = coords[vindices[2]]
        if segment_intersects_circle(p1,p2,c,r)
            return 2
        end
    end

    return 3
end

"""
    cellsubset_bounding_box(lc, coords, csubset)

Compute the bounding box for a cell subset.

For the Lefschetz complex `lc::LefschetzComplex`, whose vertices
have the coordinates given in `coords::Vector{Vector{Float64}}`,
this function computes the smallest enclosing box for the
closure of the cellsubset given in `csubset::Cells`. The function
returns the coordinates `bmin` and `bmax` of the minimal
and maximal corners of the box, respectively.
"""
function cellsubset_bounding_box(lc::LefschetzComplex,
                                 coords::Vector{Vector{T}},
                                 csubset::Cells) where T<:Real
    #
    # Compute the bounding box for a cell subset
    #
    cdim = length(coords[1])
    bmin = fill(0.0, cdim)
    bmax = fill(0.0, cdim)

    # Convert to index format if csubset is a string vector
    if typeof(csubset) == Vector{String}
        csubsetI = convert_cells(lc, csubset)
    else
        csubsetI = csubset
    end

    # Extract the vertices in the closure
    vindices = lefschetz_skeleton(lc, csubsetI, 0)

    # Find the minimal and maximal coordinates
    for k=1:cdim
        ck = [coords[m][k] for m in vindices]
        bmin[k] = minimum(ck)
        bmax[k] = maximum(ck)
    end

    return bmin, bmax
end

"""
    signed_distance_rectangle(p, rmin, rmax)

Compute the signed distance from a point to a rectangle.

The point in question is given as `p::Vector{Real}`.
The rectangle has to be parallel to the coordinate axes and is
given via its lower left corner `rmin::Vector{Real}` and its
upper right corner `rmax::Vector{Real}`. The function returns
the signed distance, which is negative inside the rectangle,
and positive outside.
"""
function signed_distance_rectangle(p::Vector{<:Real},
                                   rmin::Vector{<:Real},
                                   rmax::Vector{<:Real})
    #
    # Compute the signed distance from a point to a rectangle
    #
    xmin, ymin = rmin
    xmax, ymax = rmax
    @assert xmin < xmax "Rectangle points need xmin < xmax!"
    @assert ymin < ymax "Rectangle points need ymin < ymax!"

    # Translate point to rectangle-centered coordinates
    px = p[1] - 0.5 * (xmin + xmax)
    py = p[2] - 0.5 * (ymin + ymax)
    hx = 0.5 * (xmax - xmin)
    hy = 0.5 * (ymax - ymin)

    # Distance to rectangle along each axis
    dx = abs(px) - hx
    dy = abs(py) - hy

    # Outside distance
    outside = hypot(max(dx, 0), max(dy, 0))

    # Inside distance (negative)
    inside = min(max(dx, dy), 0)

    return outside + inside
end

"""
    signed_distance_circle(p, c, r)

Compute the signed distance from a point to a circle.

The point is specified as `p::Vector{Real}`, while the
circle is given by its center point `c::Vector{Real}`
and radius `r::Real`. The function returns the signed 
distance, which is negative inside the circle, and
positive outside.
"""
function signed_distance_circle(p::Vector{<:Real},
                                c::Vector{<:Real},
                                r::Real)
    #
    # Compute the signed distance from a point to a circle
    #
    @assert r >= 0 "The radius has to be nonnegative!"
    dx = p[1] - c[1]
    dy = p[2] - c[2]
    return hypot(dx, dy) - r
end

"""
    segment_intersects_rectangle(p1, p2, rmin, rmax)

Determine whether a line segment intersects a rectangle boundary.

The endpoints of the line segment are specified in the form
`p1,p2::Vector{Real}`. The rectangle is parallel to the coordinate
axes and is given via its lower left corner `rmin::Vector{Real}`
and its upper right corner `rmax::Vector{Real}`. The function
returns `true` if the line segment intersects the rectangle
boundary, otherwise `false`.
"""
function segment_intersects_rectangle(p1::Vector{<:Real},
                                      p2::Vector{<:Real},
                                      rmin::Vector{<:Real},
                                      rmax::Vector{<:Real})
    #
    # Determine whether a line segment intersects a rectangle
    # 
    x1, y1 = p1
    x2, y2 = p2
    xmin, ymin = rmin
    xmax, ymax = rmax

    # Decide answer if either endpoint is inside the rectangle
    p1inside = (xmin ≤ x1 ≤ xmax && ymin ≤ y1 ≤ ymax)
    p2inside = (xmin ≤ x2 ≤ xmax && ymin ≤ y2 ≤ ymax)

    if p1inside && p2inside
        if (min(x1-xmin, xmax-x1, y1-ymin, ymax-y1) == 0) ||
            (min(x2-xmin, xmax-x2, y2-ymin, ymax-y2) == 0)
            return true     # One point lies on the boundary
        else
            return false    # Both points in the interior
        end
    end

    # Now at least one point is outside
    if p1inside || p2inside
        return true
    end

    # Now both points are outside: Apply the
    # Liang–Barsky line clipping algorithm
    dx = x2 - x1
    dy = y2 - y1

    p = [-dx, dx, -dy, dy]
    q = [x1 - xmin, xmax - x1, y1 - ymin, ymax - y1]

    t0, t1 = 0.0, 1.0

    for i in 1:4
        if p[i] == 0
            if q[i] < 0
                return false
            end
        else
            t = q[i] / p[i]
            if p[i] < 0
                t0 = max(t0, t)
            else
                t1 = min(t1, t)
            end
            if t0 > t1
                return false
            end
        end
    end

    return true
end

"""
    segment_intersects_circle(p1, p2, c, r)

Determine whether a line segment intersects a circle.

The endpoints of the line segment are specified in the
form `p1,p2::Vector{Real}`, while the circle is given by
its center point `c::Vector{Real}` and radius `r::Real`. 
The function returns `true` if the line segment intersects
the circle, otherwise it returns `false`.
"""
function segment_intersects_circle(p1::Vector{<:Real},
                                   p2::Vector{<:Real},
                                   c::Vector{<:Real},
                                   r::Real)
    #
    # Determine whether a line segment intersects a circle
    #
    @assert r >= 0 "The radius has to be nonnegative!"
    x1, y1 = p1
    x2, y2 = p2
    cx, cy = c

    # Vector from p1 to p2
    dx = x2 - x1
    dy = y2 - y1

    # Vector from circle center to p1
    fx = x1 - cx
    fy = y1 - cy

    # Quadratic coefficients
    a = dx*dx + dy*dy
    b = 2 * (fx*dx + fy*dy)
    c = fx*fx + fy*fy - r*r

    # Check discriminant
    discriminant = b*b - 4*a*c
    if discriminant < 0
        return false
    end

    # Find intersection points
    discriminant = sqrt(discriminant)
    t1 = (-b - discriminant) / (2*a)
    t2 = (-b + discriminant) / (2*a)

    # Check if either intersection lies on the segment
    return (0.0 ≤ t1 ≤ 1.0) || (0.0 ≤ t2 ≤ 1.0)
end

