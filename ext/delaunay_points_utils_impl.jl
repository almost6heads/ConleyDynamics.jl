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

function ConleyDynamics.delaunay_points_add_split_segs(
        pointsIn::Vector{<:Tuple{<:Real,<:Real}},
        segmentsIn::Set{Tuple{Int,Int}},
        newsegments::Vector{<:Vector{<:Vector{<:Real}}};
        verbose::Bool = false)

    # Run the arrangement algorithm on the new segments.
    new_pts, new_edges = ConleyDynamics.split_segments(newsegments)

    n_input    = length(newsegments)
    n_isect    = length(new_pts) - 2 * n_input   # rough count; may be negative
    n_out_pts  = length(new_pts)
    n_out_segs = length(new_edges)

    # Count true intersection points: vertices that are not among the 2*n_input
    # original endpoints (after snap-rounding their duplicates were merged).
    # We use n_out_pts to count cluster centroids.
    n_isect_pts = max(0, n_out_pts - 2 * n_input)

    if verbose
        @printf("split_segments: %d input segment(s)\n", n_input)
        @printf("  intersection points created : %d\n", n_isect_pts)
        @printf("  output vertices             : %d\n", n_out_pts)
        @printf("  output edges                : %d\n", n_out_segs)
    end

    points   = deepcopy(pointsIn)
    segments = deepcopy(segmentsIn)
    offset   = length(points)

    # Append the split-segment vertices (as Float64 tuples).
    for p in new_pts
        push!(points, (Float64(p[1]), Float64(p[2])))
    end

    # Append the split-segment edges, shifted by offset.
    for (i, j) in new_edges
        push!(segments, (offset + i, offset + j))
    end

    return points, segments
end

function ConleyDynamics.delaunay_points_add_split_segs(
        points::Vector{<:Tuple{<:Real,<:Real}},
        newsegments::Vector{<:Vector{<:Vector{<:Real}}};
        verbose::Bool = false)
    seg0 = Set{Tuple{Int,Int}}()
    return ConleyDynamics.delaunay_points_add_split_segs(
        points, seg0, newsegments; verbose=verbose)
end

