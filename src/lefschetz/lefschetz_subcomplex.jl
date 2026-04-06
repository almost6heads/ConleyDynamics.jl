export lefschetz_subcomplex

"""
    lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})

Extract a subcomplex from a Lefschetz complex. The subcomplex has to be
locally closed, and is given by the collection of cells in `subcomp`.
"""
function lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    # Do we have a locally closed set?

    if !lefschetz_is_locally_closed(lc, subcomp)
        error("This is not a locally closed set!")
    end

    # Determine the new struct fields

    subcompsorted  = sort(subcomp)
    sub_labels     = lc.labels[subcompsorted]
    sub_dimensions = lc.dimensions[subcompsorted]

    # Construct the new boundary matrix

    bnd = lc.boundary
    sub_boundary = sparse_minor(bnd, subcompsorted, subcompsorted)

    # Create the new Lefschetz complex and return it

    sub_lc = LefschetzComplex(sub_labels, sub_dimensions, sub_boundary; validate=false)
    return sub_lc
end

"""
    lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})

Extract a subcomplex from a Lefschetz complex. The subcomplex has to be
locally closed, and is given by the collection of cells in `subcomp`.
"""
function lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcsub = lefschetz_subcomplex(lc, subcompI)
    return lcsub
end

