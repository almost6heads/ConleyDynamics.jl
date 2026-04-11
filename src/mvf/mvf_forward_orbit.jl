export mvf_forward_orbit
export mvf_backward_orbit

"""
    mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                      sources::Cells)

Determine the forward orbit of a collection of cells.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cells in `sources` in forward time.
Critical cells in the forward orbit will be returned as singletons.
The return type of the multivectors is the same as the type of `mvf`.
"""
function mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                           sources::Cells)
    #
    # Determine the forward orbit of a collection of cells
    #
    return mvf_forward_orbit(lc, mvf, sources, typemax(Int)-1)
end

"""
    mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                      source::Cell)

Determine the forward orbit of a single cell.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cell `source` in forward time. Critical
cells in the forward orbit will be returned as singletons. The return
type of the multivectors is the same as the type of `mvf`.
"""
function mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                           source::Cell)
    #
    # Determine the forward orbit of a cell
    #
    return mvf_forward_orbit(lc, mvf, [source], typemax(Int)-1)
end

"""
    mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                      sources::Cells, nsteps::Int)

Determine a time-restricted forward orbit of a collection of cells.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cells in `sources` in at most `nsteps`
transitions in forward time. Critical cells in the forward orbit will
be returned as singletons. The return type of the multivectors is the
same as the type of `mvf`.
"""
function mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                           sources::Cells, nsteps::Int)
    #
    # Determine the forward orbit of a collection of cells
    #
    @assert nsteps>=0 "The number of steps cannot be negative!"
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

    # Determine the edges in the directed graph

    g = DiGraph(length(mvfI))
    for k = 1:length(mvfI)
        for m = 1:length(mvfI)
            if (!(k == m)) && length(intersect(mvfIcl[k], mvfI[m]))>0
                add_edge!(g, k, m)
            end
        end
    end

    # Find the forward orbit

    state = dijkstra_shortest_paths(g, srcmv)
    reachable_vertices = findall(d -> d <= nsteps, state.dists)
    targetsI = mvfI[reachable_vertices]

    # Return the result

    if mvf isa Vector{Vector{Int}}
        return targetsI
    else
        return convert_cellsubsets(lc, targetsI)
    end
end

"""
    mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                      source::Cell, nsteps::Int)

Determine a time-restricted forward orbit of a single cell.

The function returns all multivectors of the multivector field
`mvf` that can be reached from the cell `source` in at most `nsteps`
transitions in forward time. Critical cells in the forward orbit will
be returned as singletons. The return type of the multivectors is the
same as the type of `mvf`.
"""
function mvf_forward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                           source::Cell, nsteps::Int)
    #
    # Determine the forward orbit of a cell
    #
    return mvf_forward_orbit(lc, mvf, [source], nsteps)
end

######################################################################
######################################################################

"""
    mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                       sources::Cells)

Determine the backward orbit of a collection of cells.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cells in `sources` in backward time.
In other words, it contains all multivectors that flow to `sources`
in forward time. Critical cells in the backward orbit will be returned
as singletons. The return type of the multivectors is the same as
the type of `mvf`.
"""
function mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                            sources::Cells)
    #
    # Determine the backward orbit of a collection of cells
    #
    return mvf_backward_orbit(lc, mvf, sources, typemax(Int)-1)
end

"""
    mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                       source::Cell)

Determine the backward orbit of a single cell.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cell `source` in backward time. In other
words, it contains all multivectors that flow to `source` in forward
time. Critical cells in the backward orbit will be returned as singletons.
The return type of the multivectors is the same as the type of `mvf`.
"""
function mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                            source::Cell)
    #
    # Determine the backward orbit of a cell
    #
    return mvf_backward_orbit(lc, mvf, [source], typemax(Int)-1)
end

"""
    mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                       sources::Cells, nsteps::Int)

Determine a time-restricted backward orbit of a collection of cells.

The function returns all multivectors of the multivector field `mvf`
that can be reached from the cells in `sources` in at most `nsteps`
transitions in backward time. In other words, it contains all
multivectors that flow to `sources` in at most `nsteps` transitions
in forward time. Critical cells in the backward orbit will be returned
as singletons. The return type of the multivectors is the same as the
type of `mvf`.
"""
function mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                            sources::Cells, nsteps::Int)
    #
    # Determine the backward orbit of a collection of cells
    #
    @assert nsteps>=0 "The number of steps cannot be negative!"
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

    # Determine the edges in the directed graph

    g = DiGraph(length(mvfI))
    for k = 1:length(mvfI)
        for m = 1:length(mvfI)
            if (!(k == m)) && length(intersect(mvfIcl[k], mvfI[m]))>0
                add_edge!(g, m, k)   # Reverse direction of time!
            end
        end
    end

    # Find the backward orbit

    state = dijkstra_shortest_paths(g, srcmv)
    reachable_vertices = findall(d -> d <= nsteps, state.dists)
    targetsI = mvfI[reachable_vertices]

    # Return the result

    if mvf isa Vector{Vector{Int}}
        return targetsI
    else
        return convert_cellsubsets(lc, targetsI)
    end
end

"""
    mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                       source::Cell, nsteps::Int)

Determine a time-restricted backward orbit of a single cell.

The function returns all multivectors of the multivector field
`mvf` that can be reached from the cell `source` in at most
`nsteps` transitions in backward time. In other words, it contains
all multivectors that flow to `sources` in at most `nsteps`
transitions in forward time. Critical cells in the backward orbit
will be returned as singletons. The return type of the multivectors
is the same as the type of `mvf`.
"""
function mvf_backward_orbit(lc::AbstractComplex, mvf::CellSubsets,
                            source::Cell, nsteps::Int)
    #
    # Determine the backward orbit of a cell
    #
    return mvf_backward_orbit(lc, mvf, [source], nsteps)
end

