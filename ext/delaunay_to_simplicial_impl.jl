function ConleyDynamics.delaunay_to_simplicial(tt; p::Int=2)
    #
    # Convert a DelaunayTriangulation.jl triangulation into an EuclideanComplex.
    #
    # The "solid" iterators exclude ghost simplices — virtual elements used
    # internally by DelaunayTriangulation.jl to represent the unbounded region.

    ttpoints    = get_points(tt)
    ttvertices  = collect(each_solid_vertex(tt))
    tttriangles = collect(each_solid_triangle(tt))

    # Build vertex labels and coordinates.
    # Labels are zero-padded integers so that all labels have the same width,
    # which is required by create_simplicial_complex.

    labwidth = Int(ceil(log(0.5 + length(ttvertices)) / log(10)))
    labformat = "%0" * string(labwidth) * "d"

    labels = Vector{String}()
    coords = Vector{Vector{Float64}}()
    ttv2index = Dict{Int,Int}()   # maps triangulation vertex id -> sequential index

    for k = 1:length(ttvertices)
        clabel = Printf.format(Printf.Format(labformat), k)
        push!(labels, clabel)
        ptcoord = ttpoints[ttvertices[k]]
        push!(coords, [ptcoord[1], ptcoord[2]])
        ttv2index[ttvertices[k]] = k
    end

    # Collect the maximal simplices (triangles). create_simplicial_complex
    # automatically generates all sub-simplices (edges and vertices) from these.

    simplices = Vector{Vector{Int}}()
    for dttri in tttriangles
        v1 = ttv2index[dttri[1]]
        v2 = ttv2index[dttri[2]]
        v3 = ttv2index[dttri[3]]
        push!(simplices, sort([v1, v2, v3]))
    end

    sc = create_simplicial_complex(labels, simplices, coords; p=p)

    return sc
end
