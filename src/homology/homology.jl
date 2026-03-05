export homology

"""
    homology(lc::LefschetzComplex; [algorithm::String])

Compute the homology of a Lefschetz complex.

The homology is returned as a vector `betti` of Betti numbers,
where `betti[k]` is the Betti number in dimension `k-1`. The
computations are performed over the field associated with the
Lefschetz complex boundary matrix. It is possible to invoke
the computation with one of three different algorithms, by
passing the optional argument `algorithm::String`:

* `algorithm = "matrix"` selects the standard matrix algorithm,
  as for example described in the book by Edelsbrunner, Harer,
* `algorithm = "morse"` selects an algorithm based on Morse
  reductions,
* `algorithm = "pmorse"` selects a parallelized algorithm based
  on Morse reductions.

By default, the function uses the Morse reduction algorithm `morse`.
"""
function homology(lc::LefschetzComplex; algorithm::String="morse")
    #
    # Compute the homology of a Lefschetz complex
    #

    filtration = fill(Int(1), lc.ncells)
    phs, php = persistent_homology(lc, filtration, algorithm=algorithm)

    # Return the results

    betti = [length(phs[k]) for k=1:lc.dim+1]
    return betti
end

"""
    homology(lc::LefschetzComplex, subc::Cells; [algorithm::String])

Compute the homology of a Lefschetz complex, given as a subcomplex.

This method of the `homology` function computes the homology of the
Lefschetz complex given by the locally closed subset `subc` of the
ambient Lefschetz complex `lc`. An error is raised if `subc` is
not locally closed. The homology is returned as a vector `betti`
of Betti numbers, where `betti[k]` is the Betti number in
dimension `k-1`. The computations are performed over the field
associated with the Lefschetz complex boundary matrix. It is
possible to invoke the computation with one of three different
algorithms, by passing the optional argument `algorithm::String`:

* `algorithm = "matrix"` selects the standard matrix algorithm,
  as for example described in the book by Edelsbrunner, Harer,
* `algorithm = "morse"` selects an algorithm based on Morse
  reductions,
* `algorithm = "pmorse"` selects a parallelized algorithm based
  on Morse reductions.

By default, the function uses the Morse reduction algorithm `morse`.
"""
function homology(lc::LefschetzComplex, subc::Cells;
                  algorithm::String="morse")
    #
    # Compute the homology of a Lefschetz complex
    #

    # Make sure the set is locally closed

    @assert lefschetz_is_locally_closed(lc, subc) "Not locally closed!"

    # Create the subcomplex and compute the homology

    lcsub = lefschetz_subcomplex(lc, subc)
    filtration = fill(Int(1), lcsub.ncells)
    phs, php = persistent_homology(lcsub, filtration, algorithm=algorithm)

    # Return the results

    betti = [length(phs[k]) for k=1:lc.dim+1]
    return betti
end

