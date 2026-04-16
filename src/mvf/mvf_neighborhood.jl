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

    critical = setdiff(collect(1:lc.ncells), reduce(vcat, mvfI))
    for k in critical
        push!(mvfI, [k])
    end

    # Compute the closures of the multivectors

    mvfIcl = map(t -> lefschetz_closure(lc,t), mvfI)

    # Determine the source multivectors

    srcmv = Vector{Int}()
    for k in 1:length(mvfI)
        if length(intersect(mvfI[k], sourcesI)) > 0
            push!(srcmv, k)
        end
    end

    # Determine the edges in the undirected graph:
    # There is an edge whenever the closures of the
    # multivectors intersect.

    g = Graph(length(mvfI))
    for k = 1:length(mvfI)
        for m = 1:length(mvfI)
            if (!(k == m)) && length(intersect(mvfIcl[k], mvfIcl[m]))>0
                add_edge!(g, k, m)
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

