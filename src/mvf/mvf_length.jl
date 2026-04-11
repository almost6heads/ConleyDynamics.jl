export mvf_length

"""
    mvf_length(lc::AbstractComplex, mvf::CellSubsets)

Number of multivectors in a multivector field.

The function returns the number of all multivectors in a given multivector
field. This not only includes the multivectors explicitly listed in the
vector `mvf`, but also the number of all singletons, which are trivially
critical multivectors.

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

julia> mvf_length(lc, mvf)
8
```
"""
function mvf_length(lc::AbstractComplex, mvf::CellSubsets)
    #
    # Find the number of multivectors in a
    # multivector field
    #

    if typeof(mvf) == Vector{Vector{String}}
        mvfI = convert_cellsubsets(lc, mvf)
    else
        mvfI = mvf
    end

    # Determine the number of multivectors and critical cells

    nm = length(mvfI)
    nc = length(setdiff(collect(1:lc.ncells), vcat(mvfI...)))

    return nm + nc
end

