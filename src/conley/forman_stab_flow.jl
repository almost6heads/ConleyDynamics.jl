export forman_stab_flow

"""
    forman_stab_flow(lc::AbstractComplex, fvf::CellSubsets; maxit::Int=25)

Compute Forman's stabilized combinatorial flow and the associated chain
homotopy which relates it to the identity.

This function returns matrix representations of Forman's stabilized
combinatorial flow `phiI`, and of the associated chain homotopy `gammaI`.
The third return argument is a boolean flag which indicates whether or
not the combinatorial flow stabilized. If it did not, then either the 
underlying Forman vector field is not gradient (this is not checked!),
or the maximal number of iterations has been reached. In the latter case,
one has to pass the optional paramter `maxit` with a larger number of
allowed iterations.

# Example
```jldoctest
julia> labels = ["a","b","c","d"];

julia> simplices = [["a","b"], ["b","c"], ["c","d"]];

julia> lc = create_simplicial_complex(labels, simplices, p=5);

julia> fvf = [["ab","b"], ["bc","c"]];

julia> phiI, gammaI, stabilized = forman_stab_flow(lc, fvf);

julia> stabilized
true

julia> full_from_sparse(gammaI)
7×7 Matrix{Int64}:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  4  4  0  0  0  0
 0  0  4  0  0  0  0
 0  0  0  0  0  0  0

julia> full_from_sparse(phiI)
7×7 Matrix{Int64}:
 1  1  1  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  1  0  0  0
 0  0  0  0  0  0  1
 0  0  0  0  0  0  1
 0  0  0  0  0  0  1
```
"""
function forman_stab_flow(lc::AbstractComplex, fvf::CellSubsets; maxit::Int=25)
    #
    # Compute Forman's stabilized combinatorial flow and the
    # associated chain homotopy which relates it to the identity
    #

    # Compute the combinatorial flow phi, as well as gamma
    
    phi, gamma = forman_comb_flow(lc, fvf)

    # Take powers of phi until they stabilize

    phiI = deepcopy(phi)
    gammaI = deepcopy(gamma)
    k = 1
    stabilized = false

    while (k <= maxit) && (!stabilized)
        phiN = phi * phiI
        gammaN = gamma + phi * gammaI
        if phiN == phiI
            stabilized = true
        else
            k = k + 1
            phiI = deepcopy(phiN)
            gammaI = deepcopy(gammaN)
        end
    end

    # Return the results

    return phiI, gammaI, stabilized
end

