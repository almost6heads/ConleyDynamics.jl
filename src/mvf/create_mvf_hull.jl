export create_mvf_hull

"""
    create_mvf_hull(lc::AbstractComplex, mvfbase::Vector{Vector{Int}})

Create the smallest multivector field containing the given sets.

The resulting multivector field has the property that every set of
the form `mvfbase[k]` is contained in a minimal multivector. Notice
that these sets do not have to be disjoint, and that not even their
locally closed hulls have to be disjoint. In the latter case, this
leads to two such sets having to be contained in the same multivector.
If the sets in `mvfbase` are poorly chosen, one might end up with
extremely large multivectors due to the above potential merging of
locally closed hulls. The merging only recomputes the locally closed
hulls of multivectors which have grown, so its cost is governed by the
sizes of the merged multivectors, and not by the size of the complex.
"""
function create_mvf_hull(lc::AbstractComplex, mvfbase::Vector{Vector{Int}})
    #
    # Create the smallest multivector field containing the given sets
    #

    # The multivectors are kept as the blocks of a union-find structure
    # on the cells. Every merge below is forced, since the locally closed
    # hull of a set is contained in every locally closed set containing
    # it. A block is checked again only after it has grown, so the work
    # is proportional to the sizes of the blocks which change, and not to
    # the size of the complex in every merging round. The result is the
    # same as the one of the earlier implementation, which recomputed all
    # locally closed hulls in every round.

    ncells  = lc.ncells
    parent  = collect(1:ncells)
    inblock = falses(ncells)
    members = Vector{Vector{Int}}(undef, ncells)

    # Find the root of the block containing a cell (with path halving)

    function findroot(x::Int)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end

    # Merge two blocks, keeping the larger one as the root

    function mergeblocks!(a::Int, b::Int)
        a = findroot(a)
        b = findroot(b)
        a == b && return a
        if length(members[a]) < length(members[b])
            a, b = b, a
        end
        parent[b] = a
        append!(members[a], members[b])
        members[b] = Vector{Int}()
        return a
    end

    # Put all cells of a set into one block, and return its root

    function absorb!(cellset::Vector{Int})
        root = cellset[1]
        for cell in cellset
            if !inblock[cell]
                inblock[cell] = true
                members[cell] = [cell]
            end
        end
        root = findroot(root)
        for cell in cellset
            root = mergeblocks!(root, cell)
        end
        return root
    end

    # Start with the locally closed hulls of the given sets

    for mvfset in mvfbase
        if length(mvfset) > 1
            absorb!(lefschetz_lchull(lc, mvfset))
        end
    end

    # Check every block once, and again whenever it has grown

    worklist = [k for k in 1:ncells if inblock[k] && (findroot(k) == k)]
    dirty = falses(ncells)
    dirty[worklist] .= true

    while !isempty(worklist)
        root = findroot(pop!(worklist))
        dirty[root] || continue
        dirty[root] = false
        blockhull = lefschetz_lchull(lc, members[root])
        if length(blockhull) > length(members[root])
            newroot = absorb!(blockhull)
            if !dirty[newroot]
                dirty[newroot] = true
                push!(worklist, newroot)
            end
        end
    end

    # Collect the multivectors, each sorted, and ordered by their
    # smallest cell as in the earlier implementation

    mvf = Vector{Vector{Int}}()
    for k in 1:ncells
        if inblock[k] && (findroot(k) == k) && (length(members[k]) > 1)
            push!(mvf, sort(members[k]))
        end
    end
    sort!(mvf, by=first)

    # Return the multivector field

    return mvf
end

"""
    create_mvf_hull(lc::AbstractComplex, mvfbase::Vector{Vector{String}})

Create the smallest multivector field containing the given sets.

The resulting multivector field has the property that every set of
the form `mvfbase[k]` is contained in a minimal multivector. Notice
that these sets do not have to be disjoint, and that not even their
locally closed hulls have to be disjoint. In the latter case, this
leads to two such sets having to be contained in the same multivector.
If the sets in `mvfbase` are poorly chosen, one might end up with
extremely large multivectors due to the above potential merging of
locally closed hulls.
"""
function create_mvf_hull(lc::AbstractComplex, mvfbase::Vector{Vector{String}})
    #
    # Convert a multivector field from label to index form
    #

    mvfbaseI = Vector{Vector{Int}}()

    for k=1:length(mvfbase)
        push!(mvfbaseI,Vector{Int}())
        for m=1:length(mvfbase[k])
            push!(mvfbaseI[k],lc.indices[mvfbase[k][m]])
        end
    end

    mvfI = create_mvf_hull(lc, mvfbaseI)
    mvf  = Vector{Vector{String}}()

    for m in mvfI
        push!(mvf, lc.labels[m])
    end

    return mvf
end

