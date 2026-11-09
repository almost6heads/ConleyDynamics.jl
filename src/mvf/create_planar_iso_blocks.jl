export create_planar_dts, close_dts


function create_planar_dts(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf)

    # find vert and edge indices    
    vi = findall(isequal(0), lc.dimensions)
    
    ei_nt = planar_nontransverse_edges(lc, coords, vf)
    ei = findall(isequal(1), lc.dimensions)
    ei = setdiff(ei, ei_nt)

    # Find the vertex to region map
    v2r = [Vector{Int}() for _ in 1:length(vi)]

    Threads.@threads for k in eachindex(vi)
        v2r[k] = planar_mvf_vertex2region(vi[k],lc,coords,vf)
    end

    # Find the edge to region map
    e2r = [Vector{Int}() for _ in 1:length(ei)]

    Threads.@threads for k in eachindex(ei)
        e2r[k] = planar_mvf_edge2region(ei[k],lc,coords,vf)
    end

    # Find each NT edge's toplices
    e2r_nt = [Vector{Int}() for _ in 1:length(ei_nt)]

    Threads.@threads for k in eachindex(ei_nt)
        e2r_nt[k] = lefschetz_coboundary(lc, ei_nt[k])
    end    

    # Create the future dynamics graph
    # the vertex set are the indices of elements in the Lefschetz Complex
    future_graph = SimpleGraph(lc.ncells)

    # map vertex edges
    for k in eachindex(vi)
        if length(v2r[k]) > 0
            for m in eachindex(v2r[k])
                add_edge!(future_graph, vi[k], v2r[k][m])
            end
        end
    end

    # map edge region edges
    for k in eachindex(ei)
        if length(e2r[k]) > 0
            for m in eachindex(e2r[k])
                add_edge!(future_graph, ei[k], e2r[k][m])
            end
        end
    end

    # map nt edges to their neighboring toplices
    for k in eachindex(ei_nt)
        if length(e2r_nt[k]) > 0
            for m in eachindex(e2r_nt[k])
                add_edge!(future_graph, ei_nt[k], e2r_nt[k][m])
            end
        end
    end

    return future_graph, connected_components(future_graph)
end

function create_planar_dts(ec::EuclideanComplex, vf)
    #
    # Build a coords_by_index vector (indexed by cell index, vertex cells filled)
    #
    lc = euclidean_to_lefschetz(ec)
    coords_by_index = [Vector{Float64}() for _ in 1:ec.ncells]
    for k in 1:ec.ncells
        if ec.dimensions[k] == 0
            coords_by_index[k] = ec.coords[k][1]
        end
    end

    return create_planar_dts(lc, coords_by_index, vf)
end


################################################################################

function close_dts(lc::LefschetzComplex, dts::Vector{Vector{Int}})

    merge_flag = true

    while merge_flag
        
        merge_flag = false
        
        for k in eachindex(dts)
            this_lch = lefschetz_lchull(lc, dts[k])
            
            if issetequal(dts[k], this_lch)
                continue
            end
            
            for m in setdiff(eachindex(dts), k)
                if !isempty(intersect(this_lch, dts[m]))
                    temp_union = sort(union(dts[k], dts[m]))
                    deleteat!(dts, k)
                    deleteat!(dts, m)
                    push!(dts, temp_union)
                    
                    merge_flag = true
                    break
                end
            end
        
            if merge_flag
                break
            end
        end
    end

    return dts
end

function close_dts(ec::EuclideanComplex, dts)
    #
    lc = euclidean_to_lefschetz(ec)
    
    return close_dts(lc, dts)
end