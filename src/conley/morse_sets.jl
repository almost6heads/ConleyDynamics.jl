export morse_sets

"""
    morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false)

Find the nontrivial Morse sets of a multivector field on a Lefschetz complex.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The function returns the nontrivial Morse sets as
a `Vector{Vector{Int}}`. If the optional argument `poset=true` is added,
then the function returns both the Morse sets and the adjacency matrix of
the Hasse diagram of the underlying poset.
"""
function morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false)
    #
    # Find the nontrivial Morse sets of a multivector field
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
        return_as_labels = false
    else
        mvfI = convert_cellsubsets(lc, mvf)
        return_as_labels = true
    end

    # Create the cell-to-mv map

    cell2mv = fill(Int(0), lc.ncells)

    for k = 1:length(mvfI)
        for m in mvfI[k]
            if cell2mv[m] == 0
                cell2mv[m] = k
            else
                error(" The multivector field is not a partition!")
            end
        end
    end

    # Create the digraph based on the boundary matrix

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(lc.boundary)

    # Add cycles for the multivectors to the lists c and r

    lenmvf = length(mvfI)
    for k = 1:lenmvf
        mv  = mvfI[k]
        lmv = length(mv)
        for m = 1:lmv-1
            push!(c,mv[m])
            push!(r,mv[m+1])
        end
        push!(c,mv[lmv])
        push!(r,mv[1])
    end

    # Create the edge list and the digraph

    edgelist = Edge.([(c[k],r[k]) for k in 1:length(c)])
    dg = SimpleDiGraph(edgelist)

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    for k=1:length(scc)
        sort!(scc[k])
    end

    # Find all nontrivial Morse sets:
    # This means either nontrivial homology or a
    # Morse set containing at least two multivectors

    isnontrivial = fill(Int(0), length(scc))

    # First check on nontrivial homology in parallel

    Threads.@threads for k = 1:length(scc)
        isnontrivial[k] = sum(conley_index(lc, scc[k]))
    end

    # Then check on number of multivectors

    for k = 1:length(scc)
        if length(unique(cell2mv[scc[k]])) > 1
            isnontrivial[k] = isnontrivial[k] + 1
        end
    end

    # Finally extract the Morse sets

    ntindices = findall(x -> x>0, isnontrivial)
    morsesets = scc[ntindices]

    # Determine the Hasse diagram of the poset

    if poset
        nms = length(morsesets)
        psadj = fill(false,nms,nms)
        for k1 in 1:nms
            vk1 = morsesets[k1][1]
            for k2 in 1:nms
                vk2 = morsesets[k2][1]
                psadj[k1,k2] = has_path(dg, vk1, vk2)
            end
        end
        psdg = SimpleDiGraph(psadj)
        pshasse = transitivereduction(psdg)

        hasse = fill(false,nms,nms)
        for k1 in 1:nms
            for k2 in 1:nms
                hasse[k2,k1] = has_edge(pshasse, k1, k2)
            end
        end
    end

    # Return the results

    if return_as_labels
        morsesets = convert_cellsubsets(lc, morsesets)
    end

    if poset
        return morsesets, hasse
    else
        return morsesets
    end
end

