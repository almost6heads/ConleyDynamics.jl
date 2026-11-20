export mvf_digraph, cyclic_partition_bfs, cyclic_partition_angular

"""
    out, cell2mv, mouths = mvf_digraph(lc::AbstractComplex, mvf::Vector{Vector{Int}})

Compute the multivector digraph of a multivector field.

The vertices of the digraph are the indices of the multivectors in `mvf`,
and there is an edge from `k` to `m != k` if the mouth of `mvf[k]`, i.e.,
the set `cl(mvf[k]) ∖ mvf[k]`, intersects `mvf[m]`. The function returns
three vectors:

* `out`: The adjacency lists of the digraph, i.e., `out[k]` contains the
  indices of the targets of all edges starting at `k`, without
  repetitions and in no particular order. Loops are not included.
* `cell2mv`: A vector of length `lc.ncells`, whose entry for a cell is
  the index of the multivector in `mvf` which contains it.
* `mouths`: The vector of the mouths of all multivectors, i.e.,
  `mouths[k]` contains the cells of `cl(mvf[k]) ∖ mvf[k]`.

The multivector field `mvf` has to be complete, i.e., it has to contain
all multivectors, including the singletons, which are usually omitted, see
`mvf_with_singletons`. In other words, every cell of `lc` has to lie in
exactly one set of `mvf`. This is checked, and an `ArgumentError` is
raised if a cell is missing, or if it lies in more than one set. Without
this check, a missing singleton would silently remove all edges into it
from the digraph.
"""
function mvf_digraph(lc::AbstractComplex, mvf::Vector{Vector{Int}})
    #
    # Compute the multivector digraph
    #
    cell2mv = complete_mvf_index(lc, mvf, "mvf_digraph")

    out    = [Vector{Int}() for _ in 1:length(mvf)]
    mouths = [Vector{Int}() for _ in 1:length(mvf)]

    Threads.@threads for k in 1:length(mvf)
        mouth = setdiff(lefschetz_closure(lc, mvf[k]), mvf[k])
        mouths[k] = mouth
        out[k] = unique([cell2mv[c] for c in mouth if (cell2mv[c] != 0) && (cell2mv[c] != k)])
    end

    return out, cell2mv, mouths
end

"""
    labels = cyclic_partition_bfs(out::Vector{Vector{Int}}, C::Vector{Int},
                                  seeds::Vector{Int})

Split a set of multivectors into pieces by a breadth-first search.

The adjacency lists `out` describe a digraph, for example the multivector
digraph from `mvf_digraph`, and `C` is a set of its vertices, usually a
strongly connected component. The search starts simultaneously from all
vertices in `seeds`, which have to belong to `C`, follows only edges
within `C`, and assigns every vertex of `C` to the seed from which it is
reached first. Vertices which are reached from several seeds in the same
step are assigned to the seed which comes first in `seeds`.

The function returns a vector `labels` of the same length as `C`, where
`labels[k]` is the index in `seeds` of the piece containing `C[k]`. Thus,
the pieces are pairwise disjoint and cover `C`. Vertices which cannot be
reached from any seed get the label zero, and a warning is issued. This
cannot happen if `C` is strongly connected.

This function is based on the function `digraph_propagate` by Justin
Thorpe from pull request #10, which could assign a vertex to two seeds
and left vertices beyond a depth limit unassigned.
"""
function cyclic_partition_bfs(out::Vector{Vector{Int}}, C::Vector{Int},
                              seeds::Vector{Int})
    #
    # Split a set of vertices by a simultaneous breadth-first search
    #

    # Map the vertices of C to their positions

    pos = Dict{Int,Int}()
    for k in 1:length(C)
        pos[C[k]] = k
    end
    for s in seeds
        haskey(pos, s) || throw(ArgumentError("cyclic_partition_bfs: all seeds have to belong to C"))
    end

    # Breadth-first search, one layer at a time

    labels = zeros(Int, length(C))
    frontier = Vector{Int}()
    for (i, s) in enumerate(seeds)
        if labels[pos[s]] == 0
            labels[pos[s]] = i
            push!(frontier, s)
        end
    end

    while !isempty(frontier)
        next = Vector{Int}()
        for v in frontier
            for w in out[v]
                k = get(pos, w, 0)
                if (k != 0) && (labels[k] == 0)
                    labels[k] = labels[pos[v]]
                    push!(next, w)
                end
            end
        end
        frontier = next
    end

    if any(isequal(0), labels)
        @warn "cyclic_partition_bfs: $(count(isequal(0), labels)) vertices " *
              "cannot be reached from the seeds"
    end

    return labels
end

"""
    labels = cyclic_partition_angular(lc::AbstractComplex,
                                      coords::Vector{<:Vector{<:Real}},
                                      mvf::Vector{Vector{Int}}, C::Vector{Int},
                                      p::Int; center=[0.0, 0.0],
                                      clockwise::Bool=false)

Split a set of planar multivectors into `p` angular sectors.

For every index `k` in `C`, the barycenter of the vertices of the
multivector `mvf[k]` is computed, where `coords` contains the coordinates
of the vertices as in `create_planar_mvf`. The plane is divided into `p`
sectors of equal angle around `center`, starting at the positive
horizontal direction and numbered counterclockwise, or clockwise if
`clockwise=true`. The function returns a vector `labels` of the same
length as `C`, where `labels[k]` is the number of the sector between 1
and `p` which contains the barycenter of `mvf[C[k]]`.
"""
function cyclic_partition_angular(lc::AbstractComplex,
                                  coords::Vector{<:Vector{<:Real}},
                                  mvf::Vector{Vector{Int}}, C::Vector{Int},
                                  p::Int; center=[0.0, 0.0],
                                  clockwise::Bool=false)
    #
    # Split a set of planar multivectors into angular sectors
    #
    labels = zeros(Int, length(C))
    for k in 1:length(C)
        vertices = lefschetz_skeleton(lc, mvf[C[k]], 0)
        bary = sum(coords[v] for v in vertices) ./ length(vertices) .- center
        theta = (clockwise ? -1 : 1) * atan(bary[2], bary[1])
        labels[k] = 1 + min(p - 1, floor(Int, mod(theta, 2pi) / (2pi / p)))
    end
    return labels
end

"""
    labels = cyclic_partition_angular(ec::EuclideanComplex,
                                      mvf::Vector{Vector{Int}}, C::Vector{Int},
                                      p::Int; center=[0.0, 0.0],
                                      clockwise::Bool=false)

Split a set of planar multivectors on a `EuclideanComplex` into `p`
angular sectors, see the method for Lefschetz complexes.
"""
function cyclic_partition_angular(ec::EuclideanComplex,
                                  mvf::Vector{Vector{Int}}, C::Vector{Int},
                                  p::Int; center=[0.0, 0.0],
                                  clockwise::Bool=false)
    _check_coords2d(ec, "cyclic_partition_angular")
    lc = euclidean_to_lefschetz(ec)
    return cyclic_partition_angular(lc, euclidean_vertex_coords(ec), mvf, C, p,
                                    center=center, clockwise=clockwise)
end

function complete_mvf_index(lc::AbstractComplex, mvf::Vector{Vector{Int}}, fname::String)
    #
    # Check that mvf is a complete multivector field, i.e., that every cell
    # of lc lies in exactly one of its sets, and return for every cell the
    # index of the set containing it
    #
    ncells = lc.ncells
    cell2mv = zeros(Int, ncells)
    multiple = Vector{Int}()
    for k in 1:length(mvf)
        for c in mvf[k]
            if !(1 <= c <= ncells)
                throw(ArgumentError(fname * ": the multivector field contains the " *
                                    "index $c, which is not a cell of the complex"))
            end
            cell2mv[c] == 0 ? (cell2mv[c] = k) : push!(multiple, c)
        end
    end
    if !isempty(multiple)
        multiple = sort(unique(multiple))
        throw(ArgumentError(fname * ": $(length(multiple)) cells lie in more than " *
                            "one multivector, for example the cells " *
                            "$(multiple[1:min(5,end)])"))
    end
    uncovered = findall(iszero, cell2mv)
    if !isempty(uncovered)
        throw(ArgumentError(fname * ": the multivector field is not complete, " *
                            "$(length(uncovered)) cells lie in no multivector, for " *
                            "example the cells $(uncovered[1:min(5,end)]). Use " *
                            "mvf_with_singletons to add the singletons"))
    end
    return cell2mv
end
