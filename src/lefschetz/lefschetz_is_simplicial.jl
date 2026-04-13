export lefschetz_is_simplicial

"""
    lefschetz_is_simplicial(lc::AbstractComplex)

Determine whether a Lefschetz complex is a simplicial complex.
"""
function lefschetz_is_simplicial(lc::AbstractComplex)
    #
    # Determine whether a Lefschetz complex is simplicial
    #

    # Deal with some trivial cases first

    if (lc.ncells == 0) || (lc.dim == 0)
        return true
    end

    # Determine the boundary of every cell as sets and make
    # sure every cell has the correct number of faces

    cbnd = [Set(lefschetz_boundary(lc, k)) for k=1:lc.ncells]

    for k = 1:lc.ncells
        kdim    = lc.dimensions[k]
        kbndlen = length(cbnd[k])
        if (kdim == 0) && (kbndlen > 0)
            return false
        end
        if (kdim > 0) && (!(kbndlen == kdim + 1))
            return false
        end
    end

    # Make sure all nontrivial boundaries are distinct

    ind1p = findall(t -> t >= 1, lc.dimensions)

    if !(length(unique(cbnd[ind1p])) == length(ind1p))
        return false
    end

    # Check the size of the 0-skeleton of the closure of a cell

    cskel = [lefschetz_skeleton(lc, [k], 0) for k=1:lc.ncells]

    for k = 1:lc.ncells
        if !(length(cskel[k]) == 1 + lc.dimensions[k])
            return false
        end
    end

    # Once we make it here, it has to be simplicial

    return true
end

