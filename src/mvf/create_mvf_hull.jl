export create_mvf_hull

"""
    create_mvf_hull(lc::LefschetzComplex, mvfbase::Vector{Vector{Int}})

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
function create_mvf_hull(lc::LefschetzComplex, mvfbase::Vector{Vector{Int}})
    #
    # Create the smallest multivector field containing the given sets
    #

    # Create the locally closed hulls of each set

    lchulls = Vector{Vector{Int}}()

    for mvfset in mvfbase
        mvfset_length = length(mvfset)
        if mvfset_length > 1
            mvfset_lchull = lefschetz_lchull(lc, mvfset)
            push!(lchulls, mvfset_lchull)
        end
    end

    keepgoing = true   # We need to force an initial merge

    # If one of the sets was not locally closed, start merging
    # and then recompute the locally closed hulls

    while keepgoing

        # Create the graph of multivectors for merging

        edgelist = Vector{Tuple{Int,Int}}()
        for k in lchulls
            for m = 1:length(k)-1
                push!(edgelist, (k[m], k[m+1]))
            end
        end
        sgraph = SimpleGraph(lc.ncells)
        add_edges!(sgraph, Edge.(edgelist))

        # Find the connected components

        conncomp = connected_components(sgraph)

        # Create again the locally closed hulls

        lchulls = Vector{Vector{Int}}()
        keepgoing = false

        for mvfset in conncomp
            mvfset_length = length(mvfset)
            if mvfset_length > 1
                mvfset_lchull = lefschetz_lchull(lc, mvfset)
                push!(lchulls, mvfset_lchull)
                if length(mvfset_lchull) > mvfset_length
                    keepgoing = true   # The set was not locally closed
                end
            end
        end
    end

    # Return the multivector field

    return lchulls
end

"""
    create_mvf_hull(lc::LefschetzComplex, mvfbase::Vector{Vector{String}})

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
function create_mvf_hull(lc::LefschetzComplex, mvfbase::Vector{Vector{String}})
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

