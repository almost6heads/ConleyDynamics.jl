export persistent_homology

"""
    persistent_homology(lc::AbstractComplex, filtration::Vector{Int};
                        [algorithm::String])

Compute the persistent homology of a Lefschetz complex filtration over
the field associated with the Lefschetz complex boundary matrix.

The function returns the two values
* `phsingles::Vector{Vector{Int}}`
* `phpairs::Vector{Vector{Tuple{Int,Int}}}`
It assumes that the order given by the filtration values is admissible,
i.e., the permuted boundary matrix is strictly upper triangular. The
function returns the starting filtration values for infinite length
persistence intervals in `phsingles`, and the birth- and death-filtration
values for finite length persistence intervals in `phpairs`. It is
possible to invoke the persistent homology computation with one of
three different algorithms, by passing the optional
argument `algorithm::String`:

* `algorithm = "matrix"` selects the standard matrix algorithm,
  as for example described in the book by Edelsbrunner, Harer,
* `algorithm = "morse"` selects an algorithm based on Morse
  reductions,
* `algorithm = "pmorse"` selects a parallelized algorithm based
  on Morse reductions.

By default, the function uses the Morse reduction algorithm `morse`.
"""
function persistent_homology(lc::AbstractComplex, filtration::Vector{Int};
                             algorithm::String="matrix")
    #
    # Compute the persistent homology of a Lefschetz complex filtration
    #

    # Select the method of choice

    if uppercase(algorithm) == "MATRIX"
        phsingles, phpairs = ph_matrix_reduce(lc, filtration)
    elseif uppercase(algorithm) == "MORSE"
        phsingles, phpairs = ph_morse_reduce(lc, filtration)    
    elseif uppercase(algorithm) == "PMORSE"
        phsingles, phpairs = ph_morse_reduce(lc, filtration, parallel=true)
    else
        error("Invalid value of the algorithm flag!")
    end

    # Return the results

    return phsingles, phpairs
end

