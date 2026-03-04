export persistent_homology

"""
    persistent_homology(lc::LefschetzComplex, filtration::Vector{Int};
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
values for finite length persistence intervals in `phpairs`.
"""
function persistent_homology(lc::LefschetzComplex, filtration::Vector{Int};
                             algorithm::String="matrix")
    #
    # Compute the persistent homology of a Lefschetz complex filtration
    #

    # Select the method of choice



    if uppercase(algorithm) == "MATRIX"
        phsingles, phpairs = ph_matrix_reduce(lc, filtration)
    elseif uppercase(algorithm) == "MORSE"
        phsingles, phpairs = ph_morse_reduce(lc, filtration)
    else
        error("Invalid value of the algorithm flag!")
    end

    # Return the results

    return phsingles, phpairs
end

