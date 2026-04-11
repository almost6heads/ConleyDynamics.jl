export mvf_critical

"""
    mvf_critical(lc::AbstractComplex, mvf::CellSubsets)

Determines all critical multivectors.

The function returns all critical multivectors of the specified
multivector field `mvf`. This includes both the implied singletons,
as well as the multivectors listed in `mvf` which have nontrivial
homology.

# Example
```jldoctest
julia> lc, mvf= example_forman1d();

julia> mvf
5-element Vector{Vector{String}}:
 ["A", "AD"]
 ["D", "CD"]
 ["C", "AC"]
 ["B", "BE"]
 ["E", "EF"]

julia> mvf_critical(lc, mvf)
3-element Vector{Vector{String}}:
 ["F"]
 ["BF"]
 ["DE"]
```
"""
function mvf_critical(lc::AbstractComplex, mvf::CellSubsets)
    #
    # Find all critical multivectors
    #

    if typeof(mvf) == Vector{Vector{String}}
        mvfI = convert_cellsubsets(lc, mvf)
    else
        mvfI = deepcopy(mvf)
    end

    # Add the singletons to the multivector field

    singletons = setdiff(collect(1:lc.ncells), vcat(mvfI...))
    for cc in singletons
        push!(mvfI, [cc,])
    end

    # Compute the Conley indices and collect the critical ones

    ch = Channel{Vector{Int}}(Inf)

    Threads.@threads for mv in mvfI
        if sum(conley_index(lc, mv)) > 0
            put!(ch, mv)
        end
    end

    # Collect the critical multivectors

    close(ch)
    mvc = sort(collect(ch))

    # Return the results

    if typeof(mvf) == Vector{Vector{String}}
        return convert_cellsubsets(lc, mvc)
    else
        return mvc
    end
end

