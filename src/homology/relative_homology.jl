export relative_homology

"""
    betti = relative_homology(lc::LefschetzComplex,
                              subc::Union{Vector{Int},Vector{String}};
                              p::Int=2)

Compute the relative homology of a Lefschetz complex with
respect to a subcomplex.

The subcomplex is the closure of the cells in `subc`, which
can be given either as indices or labels. The homology is
computed over the finite field `GF(p)` and is returned as a
vector `betti` of Betti numbers, where `betti[k]` is the
Betti number in dimension `k-1`.
"""
function relative_homology(lc::LefschetzComplex,
                           subc::Union{Vector{Int},Vector{String}};
                           p::Int=2)
    #
    # Compute the homology of a Lefschetz complex
    #

    # Create the filtration: 0 for the closed subcomplex,
    # and 1 for the rest

    filtration = fill(Int(1),lc.ncells)
    closure = lefschetz_closure(lc, subc)

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
