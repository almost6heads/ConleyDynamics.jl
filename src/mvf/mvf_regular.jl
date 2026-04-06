export mvf_regular

"""
    mvf_regular(lc::LefschetzComplex, mvf::CellSubsets)

Determines all regular multivectors.

The function returns all regular multivectors of the specified
multivector field `mvf`, i.e., all multivectors whose homology
groups are all trivial.

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

julia> mvf_regular(lc, mvf)
5-element Vector{Vector{String}}:
 ["A", "AD"]
 ["B", "BE"]
 ["C", "AC"]
 ["D", "CD"]
 ["E", "EF"]
```
"""
function mvf_regular(lc::LefschetzComplex, mvf::CellSubsets)
    #
    # Find all regular multivectors
    #

    if typeof(mvf) == Vector{Vector{String}}
        mvfI = convert_cellsubsets(lc, mvf)
    else
        mvfI = deepcopy(mvf)
    end

    # Compute the Conley indices and collect the regular ones

    ch = Channel{Vector{Int}}(Inf)

    Threads.@threads for mv in mvfI
        if sum(conley_index(lc, mv)) == 0
            put!(ch, mv)
        end
    end

    # Collect the regular multivectors

    close(ch)
    mvreg = sort(collect(ch))

    # Return the results

    if typeof(mvf) == Vector{Vector{String}}
        return convert_cellsubsets(lc, mvreg)
    else
        return mvreg
    end
end

