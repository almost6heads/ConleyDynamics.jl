export mvf_neighborhood

"""
    mvf_neighborhood(lc::AbstractComplex, mvf::CellSubsets,
                     sources::Cells, n::Int)

Determine a multivector neighborhood of a set of cells.

The function returns all multivectors of the multivector field `mvf`
that form a neighborhood of the given cells with a collar width of `n`.
More precisely, for `n = 0` one obtains the multivectors intersecting
the cells in `sources`, while for `n >= 1` the function returns all
multivectors from level `n-1`, together with all multivectors whose
closure intersects the closure of the level `n-1` neighborhood.
Critical cells in this neighborhood will be returned as singletons.
The return type of the multivectors is the same as the type of `mvf`.
"""
function mvf_neighborhood(lc::AbstractComplex, mvf::CellSubsets,
                          sources::Cells, n::Int)
    #
    # Determine a multivector neighborhood of a set of cells
    #
    @assert n>=0 "The neighborhood thickness cannot be negative!"
    @assert length(sources)>0 "Come on, give me at least one cell!"

    # Convert the multivector field and sources to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = deepcopy(mvf)
    else
        mvfI = convert_cellsubsets(lc, mvf)
    end

    if sources isa Vector{Int}
        sourcesI = deepcopy(sources)
    else
        sourcesI = convert_cells(lc, sources)
    end

    # Append the missing singletons to the mvf

    critical = setdiff(collect(1:lc.ncells), reduce(vcat, mvfI; init=Int[]))
    for k in critical
        push!(mvfI, [k])
    end

    # Compute the closures of the multivectors

    mvfIcl = map(t -> lefschetz_closure(lc,t), mvfI)

    # Determine the source multivectors

    cell2mv = cell_to_sets(lc.ncells, mvfI)
    srcmv   = sets_containing(cell2mv, sourcesI, length(mvfI))

    # Determine the edges in the undirected graph:
    # There is an edge whenever the closures of the
    # multivectors intersect, i.e., whenever both
    # closures contain a common cell.

    cell2cl = cell_to_sets(lc.ncells, mvfIcl)

    g = Graph(length(mvfI))
    for c = 1:lc.ncells
        clmv = cell2cl[c]
        for i = 1:length(clmv)
            for j = i+1:length(clmv)
                add_edge!(g, clmv[i], clmv[j])
            end
        end
    end

    # Find the neighborhood

    state = dijkstra_shortest_paths(g, srcmv)
    reachable_vertices = findall(d -> d <= n, state.dists)
    targetsI = mvfI[reachable_vertices]

    # Return the result

    if mvf isa Vector{Vector{Int}}
        return targetsI
    else
        return convert_cellsubsets(lc, targetsI)
    end
end

"""
    mvf_neighborhood(lc::AbstractComplex, mvf::CellSubsets,
                     source::Cell, n::Int)

Determine a multivector neighborhood of a cell.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cell `source` in at most `n` steps
forward or backward time. In other words, it computes a neighborhood
of the multivector containing the cell `source` with a thickness
of `n` layers of multivectors. Critical cells in this neighborhood
will be returned as singletons. The return type of the multivectors
is the same as the type of `mvf`.
"""
function mvf_neighborhood(lc::AbstractComplex, mvf::CellSubsets,
                          source::Cell, n::Int)
    #
    # Determine a multivector neighborhood of a cell
    #
    return mvf_neighborhood(lc, mvf, [source], n)
end

