
export signed_distance_rectangle
export signed_distance_circle
export segment_intersects_rectangle
export segment_intersects_circle











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
                                   r::Vector{<:Real})
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

