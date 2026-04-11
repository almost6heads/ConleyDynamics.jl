export morse_interval

"""
    morse_interval(lc::AbstractComplex, mvf::CellSubsets,
                   ms::CellSubsets)

Find the isolated invariant set for a Morse set interval.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The collection of Morse sets are contained in`ms`.
All of these sets should be Morse sets in the sense of being strongly
connected components of the flow graph. (Nevertheless, this will be enforced
in the function!) In other words, the sets in `ms` should be determined
using the function `morse_sets`!

The function returns the smallest isolated invariant set which contains
the Morse sets and their connections as a `Vector{Int}`.
"""
function morse_interval(lc::AbstractComplex, mvf::CellSubsets,
                        ms::CellSubsets)
    #
    # Find the isolated invariant set for a Morse set interval
    #

    # Convert the multivector field and the Morse sets to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
    else
        mvfI = convert_cellsubsets(lc, mvf)
    end

    if ms isa Vector{Vector{Int}}
        msI = ms
    else
        msI = convert_cellsubsets(lc, ms)
    end
    
    # Create the digraph based on the boundary matrix

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(lc.boundary)

    # Add cycles for the multivectors
    
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

    # Add cycles for the Morse sets
    
    lenms = length(msI)
    for k = 1:lenms
        cs = msI[k]
        lcs = length(cs)
        for m = 1:lcs-1
            push!(c,cs[m])
            push!(r,cs[m+1])
        end
        push!(c,cs[lcs])
        push!(r,cs[1])
    end

    # Connect the Morse sets to each other

    for k=2:length(msI)
        push!(c,msI[1][1])
        push!(r,msI[k][1])
        push!(c,msI[k][1])
        push!(r,msI[1][1])
    end

    # Create the edge list and the digraph

    edgelist = Edge.([(c[k],r[k]) for k in 1:length(c)])
    dg = SimpleDiGraph(edgelist)

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    # Find the strongly connected component which contains
    # the Morse sets in ms

    invhullindex = first(findall(x -> msI[1][1] in x, scc))
    invhull = sort(scc[invhullindex])

    # Return the results

    return invhull
end

