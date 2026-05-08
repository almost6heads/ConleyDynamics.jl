function ConleyDynamics.delaunay_points_bnd_rectangle(bmin, bmax)
    # Convert to Float64: DelaunayTriangulation.jl expects floating-point coordinates.
    xmin, ymin = float.(bmin)
    xmax, ymax = float.(bmax)

    # Four corners in counter-clockwise order starting at lower-left.
    # The repeated index 1 at the end of bndcurve closes the loop, as required
    # by the boundary_nodes convention of DelaunayTriangulation.jl.
    points   = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    bndcurve = [[1, 2, 3, 4, 1]]

    return points, bndcurve
end

function ConleyDynamics.delaunay_points_add_segment(
        pointsIn::Vector{<:Tuple{<:Real,<:Real}},
        segmentsIn::Set{Tuple{Int,Int}},
        newseg::Vector{<:Vector{<:Real}})

    points   = deepcopy(pointsIn)
    segments = deepcopy(segmentsIn)
    lns      = length(newseg)
    @assert lns >= 2 "We need at least two points in the segment!"

    # Floating-point tolerance for coincident-endpoint detection.  We use a
    # distance check rather than exact equality to handle rounding in the
    # caller's parametrisation (e.g. cos/sin roundoff on a closed circle).
    etol = 1.0e-9
    distance = sqrt((newseg[1][1] - newseg[lns][1])^2 +
                    (newseg[1][2] - newseg[lns][2])^2)
    isclosed = distance < etol

    # Append all interior points and their connecting edges.
    # For a closed curve newseg[lns] duplicates newseg[1], so we stop at
    # lns-2 here and handle the closing edge separately below.
    oldlen = length(points)
    for k = 1:lns-2
        x, y = float.(newseg[k])
        push!(points, (x, y))
        push!(segments, (oldlen + k, oldlen + k + 1))
    end

    # Handle the final portion:
    #   closed curve — the last new point connects back to the first new point,
    #                  so no duplicate coordinate is stored.
    #   open curve   — both the penultimate and final points are distinct and
    #                  need to be stored, plus the edge between them.
    k    = lns - 1
    x, y = float.(newseg[k])
    push!(points, (x, y))
    if isclosed
        push!(segments, (oldlen + k, oldlen + 1))
    else
        push!(segments, (oldlen + k, oldlen + k + 1))
        k    = lns
        x, y = float.(newseg[k])
        push!(points, (x, y))
    end

    return points, segments
end

function ConleyDynamics.delaunay_points_add_segment(
        points::Vector{<:Tuple{<:Real,<:Real}},
        newseg::Vector{<:Vector{<:Real}})
    # Convenience wrapper: initialise an empty segment set, then delegate.
    seg0 = Set{Tuple{Int,Int}}()
    return ConleyDynamics.delaunay_points_add_segment(points, seg0, newseg)
end

function ConleyDynamics.delaunay_points_add_nodes(
        pointsIn::Vector{<:Tuple{<:Real,<:Real}},
        newpoints::Vector{<:Vector{<:Real}})
    # Add new interior nodes
    points = deepcopy(pointsIn)
    for newpt in newpoints
        x, y = float.(newpt)
        push!(points, (x,y))
    end
    return points
end

