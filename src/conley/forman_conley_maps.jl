export forman_conley_maps

"""
    forman_conley_maps(lc::AbstractComplex, fvf::CellSubsets)

Compute the maps associated with the Conley complex of a Forman
gradient field 

This function returns the chain maps and chain homotopy associated with
the Conley complex of a Forman gradient vector field. These maps are
computed via the stabilized combinatorial flow. The function returns
the following variables:
* `pp:SparseMatrix`: The chain equivalence which maps the Lefschetz
  complex to the Conley complex
* `jj:SparseMatrix`: The chain equivalence which maps the Conley
  complex to the Lefschetz complex
* `hh:SparseMatrix`: The chain homotopy between `jj*pp` and the identity
* `cc:Cells`: The list of critical cells which make up the Conley complex

# Example
```jldoctest
julia> labels = ["a","b","c","d"];

julia> simplices = [["a","b"], ["b","c"], ["c","d"]];

julia> lc = create_simplicial_complex(labels, simplices, p=5);

julia> fvf = [["ab","b"], ["bc","c"]];

julia> pp, jj, hh, cc = forman_conley_maps(lc, fvf);

julia> sparse_show(pp, cc, lc.labels)
  ┆  a  b  c  d ab bc cd
--┆---------------------
 a┆  1  1  1  .  .  .  .
 d┆  .  .  .  1  .  .  .
cd┆  .  .  .  .  .  .  1

julia> sparse_show(jj, lc.labels, cc)
  ┆  a  d cd
--┆---------
 a┆  1  .  .
 b┆  .  .  .
 c┆  .  .  .
 d┆  .  1  .
ab┆  .  .  1
bc┆  .  .  1
cd┆  .  .  1

julia> sparse_show(hh, lc.labels, lc.labels)
  ┆  a  b  c  d ab bc cd
--┆---------------------
 a┆  .  .  .  .  .  .  .
 b┆  .  .  .  .  .  .  .
 c┆  .  .  .  .  .  .  .
 d┆  .  .  .  .  .  .  .
ab┆  .  4  4  .  .  .  .
bc┆  .  .  4  .  .  .  .
cd┆  .  .  .  .  .  .  .

julia> sparse_show(jj*pp - lc.boundary*hh - hh*lc.boundary)
 1 . . . . . .
 . 1 . . . . .
 . . 1 . . . .
 . . . 1 . . .
 . . . . 1 . .
 . . . . . 1 .
 . . . . . . 1
```
"""
function forman_conley_maps(lc::AbstractComplex, fvf::CellSubsets)
    #
    # Compute the chain maps and chain homotopy which take a 
    # Lefschetz complex with Forman gradient vector field to
    # its Conley complex, and back. The computations rely on
    # Forman's stabilizing flow
    #

    # Determine the maximal number of iterations
    
    maxits = maximum([length(findall(t -> t==k, lc.dimensions))
                      for k in 0:lc.dim]) + 3

    # Convert Forman field to index form if necessary

    if typeof(fvf) == Vector{Vector{String}}
        fvfi = convert_cellsubsets(lc, fvf)
    else
        fvfi = fvf
    end

    # Compute the stabilized combinatorial flow

    phiI, gammaI, stabilized = forman_stab_flow(lc, fvfi, maxit=maxits)
    @assert stabilized "This Forman vector field does not appear to be gradient!"

    # Determine the critical cells

    ccells = forman_critical_cells(lc, fvfi)

    # Find pp and jj as suitable matrix minors

    allindices = Vector{Int}(1:lc.ncells)
    pp = sparse_minor(phiI, ccells, allindices)
    jj = sparse_minor(phiI, allindices, ccells)

    # Return the results
    
    if typeof(fvf) == Vector{Vector{String}}
        ccr = convert_cells(lc, ccells)
    else
        ccr = ccells
    end

    return pp, jj, gammaI, ccr
end

