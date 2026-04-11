export relative_homology

"""
    relative_homology(lc::AbstractComplex,subc::Cells)

Compute the relative homology of a Lefschetz complex with
respect to a subcomplex. The computation is performed over
the field associated with the Lefschetz boundary matrix.

The subcomplex is the closure of the cells in `subc`, which
can be given either as indices or labels. The homology is
returned as a vector `betti` of Betti numbers, where `betti[k]`
is the Betti number in dimension `k-1`.
"""
function relative_homology(lc::AbstractComplex,subc::Cells)
    #
    # Compute the homology of a Lefschetz complex
    #

    # First deal with the case of empty subcomplex

    if length(subc) == 0
        return homology(lc)
    end

    # Convert subcomplex list to integers if necessary

    if subc isa Vector{String}
        isubc = [lc.indices[k] for k in subc]
    else
        isubc = subc
    end

    # This is just the Lefschetz homology of the difference set!

    closure = lefschetz_closure(lc, isubc)
    lcdiffset = setdiff(Vector{Int}(1:lc.ncells), closure)

    # Return the results

    if length(lcdiffset) == 0
        return fill(Int(0), lc.dim+1)
    else
        return homology(lc, lcdiffset)
    end
end

function relative_homology(lc::AbstractComplex, subc::Vector{Any})
    if length(subc) == 0
        return relative_homology(lc, Vector{Int}([]))
    else
        error("Unknow subspace type!")
    end
end

"""
    relative_homology(lc::AbstractComplex,subc::Cells,subc0::Cells)

Compute the relative homology of a Lefschetz complex with
respect to a subcomplex. The computation is performed over
the field associated with the Lefschetz boundary matrix.

In this implementation, relative homology of the pair
`cl(subc), cl(subc0))` is computed. An error is raised if
`cl(subc0)` is not a subset of `cl(subc)`. The homology is
returned as a vector `betti` of Betti numbers, where `betti[k]`
is the Betti number in dimension `k-1`.
"""
function relative_homology(lc::AbstractComplex,subc::Cells,subc0::Cells)
    #
    # Compute the homology of a Lefschetz complex
    #

    # First deal with the case of an empty complex subc

    if length(subc) == 0
        return fill(Int(0),lc.dim+1)
    end

    # Convert subcomplex lists to integers if necessary

    if subc isa Vector{String}
        isubc = [lc.indices[k] for k in subc]
    else
        isubc = subc
    end

    if subc0 isa Vector{String}
        isubc0 = [lc.indices[k] for k in subc0]
    else
        isubc0 = subc0
    end

    # Compute the closures of the complexes

    clsubc  = lefschetz_closure(lc, isubc)
    clsubc0 = lefschetz_closure(lc, isubc0)

    # Make sure that cl(subc0) is contained in cl(subc)

    is_contained = true
    for k in clsubc0
        if !(k in clsubc)
            is_contained = false
        end
    end

    if !is_contained
        error("cl(subc0) is not contained in cl(subc)")
    end

    # Compute the relative homology

    lcclsubc = lefschetz_closed_subcomplex(lc, clsubc)
    clsubc0_string = lc.labels[clsubc0]
    bettisub = relative_homology(lcclsubc, clsubc0_string)
    betti = fill(Int(0),lc.dim+1)
    for k=1:length(bettisub)
        betti[k] = bettisub[k]
    end

    # Return the result

    return betti
end

function relative_homology(lc::AbstractComplex,
                           subc::Union{Vector{Int},Vector{String}},
                           subc0::Vector{Any})
    if length(subc0) == 0
        return relative_homology(lc, subc, Vector{Int}([]))
    else
        error("Unknow subspace type!")
    end
end

