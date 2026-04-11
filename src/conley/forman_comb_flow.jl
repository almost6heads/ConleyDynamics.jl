export forman_comb_flow

"""
    forman_comb_flow(lc::AbstractComplex, fvf::Vector{Vector{String}})

Compute Forman's combinatorial flow and the associated chain homotopy.

This function returns matrix representations of Forman's combinatorial 
flow `phi`, and of the associated chain homotopy `gamma`.

# Example
```jldoctest
julia> labels = ["a","b","c","d"];

julia> simplices = [["a","b"], ["b","c"], ["c","d"]];

julia> lc = create_simplicial_complex(labels, simplices, p=5);

julia> fvf = [["ab","b"], ["bc","c"]];

julia> phi, gamma = forman_comb_flow(lc, fvf);

julia> full_from_sparse(gamma)
7×7 Matrix{Int64}:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  4  0  0  0  0  0
 0  0  4  0  0  0  0
 0  0  0  0  0  0  0

julia> full_from_sparse(phi)
7×7 Matrix{Int64}:
 1  1  0  0  0  0  0
 0  0  1  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  1  0  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
 0  0  0  0  0  0  1

julia> full_from_sparse(phi*phi)
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
function forman_comb_flow(lc::AbstractComplex, fvf::Vector{Vector{Int}})
    #
    # Compute the matrix representation of Forman's combinatorial
    # flow, as well as the associated chain homotopy
    #

    # Create the chain homotopy gamma

    nc = lc.ncells
    pp = lc.boundary.char
    gamma = sparse_zero(nc, nc, p=pp)

    for fvec in fvf
        lenfvec = length(fvec)
        if (lenfvec < 1) || (lenfvec > 2)
            error("This is not a Forman vector field!")
        end
        if lenfvec == 2
            a, b = fvec
            if lc.dimensions[a] == lc.dimensions[b]+1
                gamma[a,b] = scalar_inverse(-lc.boundary[b,a], pp)
            elseif lc.dimensions[b] == lc.dimensions[a]+1
                gamma[b,a] = scalar_inverse(-lc.boundary[a,b], pp)
            else
                error("This is not a Forman vector field!")
            end
        end
    end

    # Create Forman's combinatorial flow

    phi = sparse_identity(nc, p=pp) + lc.boundary * gamma + gamma * lc.boundary

    # Return the results

    return phi, gamma
end

"""
    forman_comb_flow(lc::AbstractComplex, fvf::Vector{Vector{String}})

Compute Forman's combinatorial flow and the associated chain homotopy.

This function returns matrix representations of Forman's combinatorial 
flow `phi`, and of the associated chain homotopy `gamma`.
"""
function forman_comb_flow(lc::AbstractComplex, fvf::Vector{Vector{String}})
    #
    # Compute the matrix representation of Forman's combinatorial
    # flow, as well as the associated chain homotopy
    #
    fvfI = convert_cellsubsets(lc, fvf)
    phi, gamma = forman_comb_flow(lc, fvfI)
    return phi, gamma
end

