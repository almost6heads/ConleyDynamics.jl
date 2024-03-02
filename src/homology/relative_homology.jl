export relative_homology

"""
    relative_homology(lc::LefschetzComplex,
                      subc::Union{Vector{Int},Vector{String}};
                      [p::Int])

Compute the relative homology of a Lefschetz complex with
respect to a subcomplex.

The subcomplex is the closure of the cells in `subc`, which
can be given either as indices or labels. The homology is
computed over the rationals (for `p=0`) or the finite field
`GF(p)` (for prime `p`) and is returned as a vector `betti`
of Betti numbers, where `betti[k]` is the Betti number in
dimension `k-1`. If the Lefschetz complex boundary matrix
already has been specialized to a field, the optional
argument `p` can be omitted.
"""
function relative_homology(lc::LefschetzComplex,
                           subc::Union{Vector{Int},Vector{String}};
                           p::Int=-1)
    #
    # Compute the homology of a Lefschetz complex
    #

    # First deal with the case of empty subcomplex

    if length(subc) == 0
        return homology(lc,p=p)
    end

    # Convert subcomplex list to integers if necessary

    if subc isa Vector{String}
        isubc = [lc.indices[k] for k in subc]
    else
        isubc = subc
    end

    # Create the filtration: 0 for the closed subcomplex,
    # and 1 for the rest

    filtration = fill(Int(1),lc.ncells)
    closure = lefschetz_closure(lc, isubc)

    for k in closure
        filtration[k] = 0
    end

    # Compute the persistent homology

    phs, php = persistent_homology(lc,filtration,p=p)

    # Assemble the Betti numbers:
    # Intervals (1,Inf) in dim p give contribute to dim p
    # Intervals (0,1) in dim p contribute to dim p+1
    #    (there should be none for m=lc.dim+1)

    betti = [Int(0) for _ in 1:lc.dim+1]
    for m=1:lc.dim+1
        for k=1:length(phs[m])
            if phs[m][k] == 1
                betti[m] += 1
            end
        end
    end
    for m=1:lc.dim+1
        for k=1:length(php[m])
            if php[m][k] == (0,1)
                betti[m+1] += 1
            end
        end
    end

    # Return the results

    return betti
end

function relative_homology(lc::LefschetzComplex, subc::Vector{Any}; p::Int=-1)
    if length(subc) == 0
        return relative_homology(lc, Vector{Int}([]), p=p)   
    else
        error("Unknow subspace type!")
    end
end

"""
    relative_homology(lc::LefschetzComplex,
                      subc::Union{Vector{Int},Vector{String}},
                      subc0::Union{Vector{Int},Vector{String}};
                      [p::Int])

Compute the relative homology of a Lefschetz complex with
respect to a subcomplex.

In this implementation, relative homology of the pair
`cl(subc), cl(subc0))` is computed. An error is raised if
`cl(subc0)` is not a subset of `cl(subc)`. The homology is
computed over the rationals (for `p=0`) or the finite field
`GF(p)` (for prime `p`) and is returned as a vector `betti`
of Betti numbers, where `betti[k]` is the Betti number in
dimension `k-1`. If the Lefschetz complex boundary matrix
already has been specialized to a field, the optional
argument `p` can be omitted.
"""
function relative_homology(lc::LefschetzComplex,
                           subc::Union{Vector{Int},Vector{String}},
                           subc0::Union{Vector{Int},Vector{String}};
                           p::Int=-1)
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

    lcclsubc = lefschetz_subcomplex(lc, clsubc)
    clsubc0_string = lc.labels[clsubc0]
    bettisub = relative_homology(lcclsubc, clsubc0_string, p=p)
    betti = fill(Int(0),lc.dim+1)
    for k=1:length(bettisub)
        betti[k] = bettisub[k]
    end

    # Return the result

    return betti
end

function relative_homology(lc::LefschetzComplex,
                           subc::Union{Vector{Int},Vector{String}},
                           subc0::Vector{Any}; p::Int=-1)
    if length(subc0) == 0
        return relative_homology(lc, subc, Vector{Int}([]), p=p)   
    else
        error("Unknow subspace type!")
    end
end





