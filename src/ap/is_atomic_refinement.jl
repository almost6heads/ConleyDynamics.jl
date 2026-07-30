export is_atomic_refinement

"""
    is_atomic_refinement(lc::AbstractComplex, mvf1::CellSubsets, mvf2::CellSubsets)

Determine whether `mvf1` is an atomic refinement of `mvf2`, i.e. whether
`mvf1` and `mvf2` are elements of `AP(X)` with `mvf1` obtained from `mvf2`
by splitting exactly one block of `mvf2` into two. Equivalently, `mvf1` and
`mvf2` cover each other in the atomic-refinement order used throughout
`construct_ap_space`/`atomic_distances`/`stratum_adjacency`.
"""
function is_atomic_refinement(lc::AbstractComplex,
                              mvf1::CellSubsets,
                              mvf2::CellSubsets)
    #
    # Check whether mvf1 is an atomic refinement
    # of the multivector field mvf2
    #

    mvset1 = convert_mvf_partition(lc, mvf1)
    mvset2 = convert_mvf_partition(lc, mvf2)
    mvl1 = length(mvset1)
    mvl2 = length(mvset2)

    # Make sure the lengths are ok

    if !(mvl1 == mvl2+1)
        return false
    end

    # Do they agree on enough multivectors?

    if length(intersect(mvset1, mvset2)) == mvl2-1
        return true
    else
        return false
    end
end
