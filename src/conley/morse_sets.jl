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

    # Create the digraph with one vertex for every cell. Building it
    # from the edge list alone would drop cells with larger indices
    # than any cell incident to an edge.

    dg = SimpleDiGraph(nc)
    for k in 1:length(c)
        add_edge!(dg, c[k], r[k])
    end

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    for k=1:length(scc)
        sort!(scc[k])
    end

    # Find all nontrivial Morse sets:
    # This means either nontrivial homology or a
    # Morse set containing at least two multivectors

    isnontrivial = fill(Int(0), length(scc))

    # A strongly connected component which contains at least two
    # multivectors is always nontrivial. Otherwise, it is a single
    # multivector, which is locally closed, and its Conley index
    # can be determined from the multivector cells alone. The
    # algorithm "auto" uses "rank" for small multivectors and
    # "morse" for large ones.

    Threads.@threads for k = 1:length(scc)
        if length(unique(cell2mv[scc[k]])) > 1
            isnontrivial[k] = 1
        else
            isnontrivial[k] = sum(conley_index(lc, scc[k], algorithm="auto",
                                               validate=false))
        end
    end

    # Finally extract the Morse sets

    ntindices = findall(x -> x>0, isnontrivial)
    morsesets = scc[ntindices]

    # Determine the Hasse diagram of the poset

    if poset
        nms  = length(morsesets)
        nscc = length(scc)

        # Map the cells to the strongly connected components,
        # and the components to the Morse sets (zero if none)

        cell2scc = fill(Int(0), nv(dg))
        for k = 1:nscc
            for v in scc[k]
                cell2scc[v] = k
            end
        end

        scc2ms = fill(Int(0), nscc)
        scc2ms[ntindices] = 1:nms

        # Process the components in reverse topological order
        # of the condensation, i.e., successors come first

        sccorder = reverse(topological_sort(condensation(dg, scc)))

        # For every component, find the Morse sets which can be reached
        # along paths that do not pass through any other Morse set

        firsthit = [Int[] for _ in 1:nscc]
        stamp    = fill(Int(0), nms)

        for k in sccorder
            fh = firsthit[k]
            for v in scc[k]
                for w in outneighbors(dg, v)
                    kw = cell2scc[w]
                    kw == k && continue
                    targets = scc2ms[kw] > 0 ? (scc2ms[kw]:scc2ms[kw]) : firsthit[kw]
                    for m in targets
                        if stamp[m] != k
                            stamp[m] = k
                            push!(fh, m)
                        end
                    end
                end
            end
        end

        # The transitive closure of the first-hit relation is the order
        # between the Morse sets. A first-hit Morse set is in the Hasse
        # diagram unless it can be reached via another first-hit Morse set

        reach   = [falses(nms) for _ in 1:nms]
        covered = falses(nms)
        hasse   = fill(false,nms,nms)

        for k in sccorder
            m1 = scc2ms[k]
            m1 == 0 && continue
            fill!(covered, false)
            for m2 in firsthit[k]
                covered .|= reach[m2]
            end
            for m2 in firsthit[k]
                reach[m1][m2] = true
                if !covered[m2]
                    hasse[m2,m1] = true
                end
            end
            reach[m1] .|= covered
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

