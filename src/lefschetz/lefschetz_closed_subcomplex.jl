export lefschetz_closed_subcomplex

"""
    lefschetz_closed_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})

Extract a closed subcomplex from a Lefschetz complex. The subcomplex is
the closure of the collection of cells given in `subcomp`.
"""
function lefschetz_closed_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    # Determine the closed subcomplex

    clsubcomp = lefschetz_closure(lc, subcomp)

    # Determine the new struct fields

    sub_labels     = lc.labels[clsubcomp]
    sub_dimensions = lc.dimensions[clsubcomp]

    # Construct the new boundary matrix

    bnd = lc.boundary
    sub_boundary = sparse_minor(bnd, clsubcomp, clsubcomp)

    # Create the new Lefschetz complex and return it

    sub_lc = LefschetzComplex(sub_labels, sub_dimensions, sub_boundary; validate=false)
    return sub_lc
end

"""
    lefschetz_closed_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})

Extract a closed subcomplex from a Lefschetz complex. The subcomplex is
the closure of the collection of cells given in `subcomp`.
"""
function lefschetz_closed_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcsub = lefschetz_closed_subcomplex(lc, subcompI)
    return lcsub
end

