export forman_all_cell_types
export forman_critical_cells

"""
    forman_all_cell_types(lc::AbstractComplex, fvf::CellSubsets)

Find all cell types of a Forman vector field.

This function returns all critical cells, all arrow sources, and all arrow
targets of the Forman vector field `fvf` on the Lefschetz complex `lc`.

# Example
```jldoctest
julia> labels = ["a","b","c","d"];

julia> simplices = [["a","b"], ["b","c"], ["c","d"]];

julia> lc = create_simplicial_complex(labels, simplices, p=5);

julia> fvf = [["ab","b"], ["bc","c"]];

julia> c, s, t = forman_all_cell_types(lc, fvf);

julia> println(c)
["a", "d", "cd"]

julia> println(s)
["b", "c"]

julia> println(t)
["ab", "bc"]
```
"""
function forman_all_cell_types(lc::AbstractComplex, fvf::CellSubsets)
    #
    # Find all cell types of a Forman vector field
    #

    # Perform a quick check that we have a Forman vector field

    if length(fvf)>0
        @assert maximum(length.(fvf))<=2 "Not a Forman vector field!"
    end

    # Deal with the trivial case first

    if length(fvf) == 0
        if typeof(fvf) == Vector{Vector{String}}
            return lc.labels, Vector{Vector{String}}(), Vector{Vector{String}}()
        else
            return Vector{Int}(1:lc.ncells), Vector{Vector{Int}}(), Vector{Vector{Int}}()
        end
    end

    # Convert fvf to index form

    if typeof(fvf) == Vector{Vector{String}}
        fvfI = convert_cellsubsets(lc, fvf)
    else
        fvfI = fvf
    end

    # Collect all cells

    cI = Vector{Int}(1:lc.ncells)
    sI = Vector{Int}()
    tI = Vector{Int}()

    for k = 1:length(fvfI)
        arrow = fvfI[k]
        if length(arrow) == 2
            a1 = arrow[1]
            a2 = arrow[2]
            if lc.dimensions[a1] == 1 + lc.dimensions[a2]
                push!(sI, a2)
                push!(tI, a1)
            elseif 1 + lc.dimensions[a1] == lc.dimensions[a2]
                push!(sI, a1)
                push!(tI, a2)
            else
                error("This is not a Forman vector!")
            end
        end
    end

    cI = sort(unique(setdiff(setdiff(cI, sI), tI)))
    sI = sort(unique(sI))
    tI = sort(unique(tI))

    lsum = length(cI) + length(sI) + length(tI)

    @assert lsum==lc.ncells "Not every cell is critical, source, or target!"

    # Return the results

    if typeof(fvf)==Vector{Vector{String}}
        return convert_cells(lc,cI), convert_cells(lc,sI), convert_cells(lc,tI)
    else
        return cI, sI, tI
    end
end

"""
    forman_critical_cells(lc::AbstractComplex, fvf::CellSubsets)

Find all critical cells of a Forman vector field.

This function returns all critical cells of the Forman vector
field `fvf` on the Lefschetz complex `lc`.

# Example
```jldoctest
julia> labels = ["a","b","c","d"];

julia> simplices = [["a","b"], ["b","c"], ["c","d"]];

julia> lc = create_simplicial_complex(labels, simplices, p=5);

julia> fvf = [["ab","b"], ["bc","c"]];

julia> c = forman_critical_cells(lc, fvf);

julia> println(c)
["a", "d", "cd"]
```
"""
function forman_critical_cells(lc::AbstractComplex, fvf::CellSubsets)
    #
    # Find all critical cells of a Forman vector field
    #
    cc, ss, tt = forman_all_cell_types(lc, fvf)
    return cc
end

