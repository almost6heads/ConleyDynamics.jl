export lefschetz_closed_subcomplex

"""
    lefschetz_closed_subcomplex(lc::AbstractComplex, subcomp::Vector{Int})

Extract a closed subcomplex from a Lefschetz complex. The subcomplex is
the closure of the collection of cells given in `subcomp`.
"""
function lefschetz_closed_subcomplex(lc::AbstractComplex, subcomp::Vector{Int})
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

    # Create the new complex and return it

    if lc isa EuclideanComplex
        sub_coords = lc.coords[clsubcomp]
        return EuclideanComplex(sub_labels, sub_dimensions, sub_boundary,
                                sub_coords; validate=false)
    else
        return LefschetzComplex(sub_labels, sub_dimensions, sub_boundary; validate=false)
    end
end

"""
    lefschetz_closed_subcomplex(lc::AbstractComplex, subcomp::Vector{String})

Extract a closed subcomplex from a Lefschetz complex. The subcomplex is
the closure of the collection of cells given in `subcomp`.

If `lc` is an `EuclideanComplex`, the returned subcomplex is also an
`EuclideanComplex` with the coordinates restricted to the subcomplex cells.
"""
function lefschetz_closed_subcomplex(lc::AbstractComplex, subcomp::Vector{String})
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

