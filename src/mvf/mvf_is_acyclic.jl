export mvf_is_acyclic

"""
    mvf_is_acyclic(lc::AbstractComplex, mvf::CellSubsets)

Determine whether a multivector field is acyclic.

The function returns a boolean which indicates whether a 
multivector field forms an acyclic partition of the Lefschetz
complex `lc` or not.
"""
function mvf_is_acyclic(lc::AbstractComplex, mvf::CellSubsets)
    #
    # Determine whether a multivector field is acyclic
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = deepcopy(mvf)
    else
        mvfI = convert_cellsubsets(lc, mvf)
    end

    # Append the missing singletons to the mvf

    critical = setdiff(collect(1:lc.ncells), reduce(vcat, mvfI))
    for k in critical
        push!(mvfI, [k])
    end

    # Compute the closures of the multivectors

    mvfIcl = map(t -> lefschetz_closure(lc,t), mvfI)

    # Determine the edges in the directed graph

    g = DiGraph(length(mvfI))
    for k = 1:length(mvfI)
        for m = 1:length(mvfI)
            if (!(k == m)) && length(intersect(mvfIcl[k], mvfI[m]))>0
                add_edge!(g, k, m)
            end
        end
    end

    # Return the result

    return !is_cyclic(g)
end

