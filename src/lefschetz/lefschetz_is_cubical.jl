export lefschetz_is_cubical
export lefschetz_is_subcubical

"""
    lefschetz_is_cubical(lc::AbstractComplex)

Determine whether a Lefschetz complex is a cubical complex.

This function determines whether a Lefschetz complex is
cubical in the sense of this package. This means that the
labels of the complex have to satisfy the constraints
explained in `create_cubical_complex`.
"""
function lefschetz_is_cubical(lc::AbstractComplex)
    #
    # Determine whether a Lefschetz complex is cubical
    #

    # Deal with some trivial cases first

    if (lc.ncells == 0) || (lc.dim == 0)
        return true
    end

    # Check whether the labels are correct cube labels and
    # whether they are consistent

    dimset = Set{Tuple{Int,Int}}()
    for k in 1:lc.ncells
        if is_cube_label(lc.labels[k])
            push!(dimset, cube_field_size(lc.labels[k]))
        else
            return false
        end
    end

    (length(dimset) > 1) && return false

    # Create the cubical complex based on the labels

    ccfull = create_cubical_complex(lc.labels, p=lc.boundary.char)

    # Return the answer

    if ccfull.ncells == lc.ncells
        return true
    else
        return false
    end
end

"""
    lefschetz_is_subcubical(lc::AbstractComplex)

Determine whether a Lefschetz complex is a locally closed
subset of a cubical complex.

This function determines whether a Lefschetz complex is
a locally cosed subset of a cubical complex in the sense of
this package. This means that the labels of the complex have
to satisfy the constraints explained in `create_cubical_complex`,
and `lc` is locally closed in the complex obtained from the
labels in `lc.labels` via `create_cubical_complex`.
"""
function lefschetz_is_subcubical(lc::AbstractComplex)
    #
    # Determine whether a Lefschetz complex is a
    # locally closed subset of a cubical complex
    #

    # Deal with some trivial cases first

    if (lc.ncells == 0) || (lc.dim == 0)
        return true
    end

    # Check whether the labels are correct cube labels and
    # whether they are consistent

    dimset = Set{Tuple{Int,Int}}()
    for k in 1:lc.ncells
        if is_cube_label(lc.labels[k])
            push!(dimset, cube_field_size(lc.labels[k]))
        else
            return false
        end
    end

    (length(dimset) > 1) && return false

    # Create the cubical complex based on the labels

    ccfull = create_cubical_complex(lc.labels, p=lc.boundary.char)

    # Return the answer

    return lefschetz_is_locally_closed(ccfull, lc.labels)
end

