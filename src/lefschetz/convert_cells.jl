export convert_cells
export convert_cellsubsets
export cellsubsets_to_cells

"""
    convert_cells(lc::AbstractComplex, cl::Vector{Int})

Convert cell list `cl` in the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_cells(lc::AbstractComplex, cl::Vector{Int})
    #
    # Convert a cell list from index to label form
    #

    return lc.labels[cl]
end

"""
    convert_cells(lc::AbstractComplex, cl::Vector{String})

Convert cell list `cl` in the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_cells(lc::AbstractComplex, cl::Vector{String})
    #
    # Convert a cell list from label to index form
    #

    return [lc.indices[k] for k in cl]
end

"""
    convert_cellsubsets(lc::AbstractComplex, clsub::Vector{Vector{Int}})

Convert CellSubsets `clsub` in the Lefschetz complex `lc` from 
index form to label form.
"""
function convert_cellsubsets(lc::AbstractComplex, clsub::Vector{Vector{Int}})
    #
    # Convert a CellSubsets from index to label form
    #

    newclsub = Vector{Vector{String}}()

    for k=1:length(clsub)
        push!(newclsub,Vector{String}())
        for m=1:length(clsub[k])
            push!(newclsub[k],lc.labels[clsub[k][m]])
        end
    end

    return newclsub
end

"""
    convert_cellsubsets(lc::AbstractComplex, clsub::Vector{Vector{String}})

Convert CellSubsets `clsub` in the Lefschetz complex `lc` from 
label form to index form.
"""
function convert_cellsubsets(lc::AbstractComplex, clsub::Vector{Vector{String}})
    #
    # Convert a CellSubsets from label to index form
    #

    newclsub = Vector{Vector{Int}}()

    for k=1:length(clsub)
        push!(newclsub,Vector{Int}())
        for m=1:length(clsub[k])
            push!(newclsub[k],lc.indices[clsub[k][m]])
        end
    end

    return newclsub
end

"""
    cellsubsets_to_cells(lc::AbstractComplex, css::CellSubsets)

Convert the cell subsets `css` into a vector of cells, by simply
taking the union of all subsets. The cells are returned in ordered
form, based on their index form in `lc`.
"""
function cellsubsets_to_cells(lc::AbstractComplex, css::CellSubsets)
    #
    # Convert cell subsets to a vector of cells
    #

    # Convert to index form

    if css isa Vector{Vector{Int}}
        cssI = deepcopy(css)
    else
        cssI = convert_cellsubsets(lc, css)
    end

    # Merge the sets

    csI = sort(unique(reduce(vcat, cssI)))

    # Return the result

    if css isa Vector{Vector{Int}}
        return csI
    else
        return convert_cells(lc, csI)
    end
end

