export create_dual_set, create_mvf
export graph_merge, dual_graph_merge
export digraph_propagate


################################################################################
function create_dual_set(ec::EuclideanComplex, set::Vector{Int}, vf)

    # Convert EC to LC
    lc = euclidean_to_lefschetz(ec)
    coords = [Vector{Float64}() for _ in 1:ec.ncells]
    for k in 1:ec.ncells
        if ec.dimensions[k] == 0
            coords[k] = ec.coords[k][1]
        end
    end

    # Get boundary elements
    set_bdry = lefschetz_topboundary(lc, set)

    set_bdry_vi = set_bdry[findall(isequal(0), lc.dimensions[set_bdry])]
    set_bdry_ei = set_bdry[findall(isequal(1), lc.dimensions[set_bdry])]


    # Find the dual mouth
    set_dmo = Vector{Int}()
    for idx in set_bdry_vi
        this_v2r = planar_mvf_vertex2region(idx, lc, coords,
                                            vf, vf_rev=true)
        if isempty(intersect(set, this_v2r))
            push!(set_dmo, idx)
        end
    end

    for idx in set_bdry_ei
        this_e2r = planar_mvf_edge2region(idx, lc, coords,
                                          vf, vf_rev=true)
        if isempty(intersect(set, this_e2r))
            push!(set_dmo, idx)
        end
    end

    # Create dual set
    set_intr = lefschetz_interior(lc, set)
    set_dcap = setdiff(set_bdry, set_dmo)

    return union(set_intr, set_dcap)
end


function graph_merge(ec::EuclideanComplex, sg::SimpleGraph)

    # Find future graph conn comp
    sg_cc = connected_components(sg)
    
    # Loop over conn comps and merge as necessary
    for cc in sg_cc
        this_lchull = setdiff(lefschetz_lchull(ec, cc), cc)

        for ele in this_lchull
            this_ohull = lefschetz_openhull(ec, [ele])
            this_ohull_tops = lefschetz_skeleton(ec, this_ohull, ec.dim)
            this_tops = intersect(this_ohull_tops, cc)
            
            for top in this_tops
                add_edge!(sg, ele, top)
            end
        end
    end

    return sg
end


function dual_graph_merge(ec::EuclideanComplex, sg::SimpleGraph, vf)

    # Create dual graph from non-dual conn comp
    sg_dual = SimpleGraph(ec.ncells)

    sg_cc = connected_components(sg)
    for cc in sg_cc
        this_dual = create_dual_set(ec, cc, vf)
        for m = 1:length(this_dual)-1
            add_edge!(sg_dual, this_dual[m], this_dual[m+1])
        end
    end

    # Loop over dual conn comps and merge as necessary
    sg_dual_cc = connected_components(sg_dual)
    for cc in sg_dual_cc
        this_lchull = setdiff(lefschetz_lchull(ec, cc), cc)

        for ele in this_lchull
            this_ohull = lefschetz_openhull(ec, [ele])
            this_ohull_tops = lefschetz_skeleton(ec, this_ohull, ec.dim)
            this_cc_tops = intersect(this_ohull_tops, cc)
            this_tops = setdiff(this_ohull_tops, this_cc_tops)
        
            # Merge in *non-dual* graph
            for top in this_tops
                add_edge!(sg, ele, top)
            end
        end
    end

    # Check if dual graph is locally closed after merging
    dual_is_closed = true
    sg_dual_cc = connected_components(sg_dual)
    for cc in sg_dual_cc
        if !lefschetz_is_locally_closed(ec, cc)
            dual_is_closed = false
            break
        end
    end

    return sg, dual_is_closed
end


function create_mvf(ec, sg, vf; strong_iso::Bool=false)
    counter = 1
    
    merge_again = true
    while merge_again && counter <= 5
        merge_again = false
        

        println("Merge loop ", counter)
        counter += 1
        

        close_loop = true
        while close_loop
            close_loop = false

            graph_cc = connected_components(sg)
            for cc in graph_cc
                if !lefschetz_is_locally_closed(ec, cc)
                    sg = graph_merge(ec, sg)
                    close_loop = true
                    break
                end
            end
        end

        if strong_iso
            sg, dual_is_closed = dual_graph_merge(ec, sg, vf)

            if !dual_is_closed
                merge_again = true
            end
        end
    end

    return sg
end


function digraph_propagate(dg, seeds, depth)
    
    pset = [[seed] for seed in seeds]
    all_parents = [[seed] for seed in seeds]
    visited = Set(seeds)

    idg = [bfs_tree(dg, seed) for seed in seeds]
    

    for _ in 1:depth

        for idx in eachindex(pset)
            this_parents = all_parents[idx]

            union!(visited, Set(this_parents))
            all_kids = Set()

            for parent in this_parents
                kids = outneighbors(idg[idx], parent)
                union!(all_kids, Set(kids))
            end

            eligible_kids = setdiff(all_kids, visited)

            all_parents[idx] = collect(eligible_kids)
            append!(pset[idx], all_parents[idx])
        end
    end

    for idx in eachindex(pset)
        sort!(pset[idx])
    end
    
    return pset
end