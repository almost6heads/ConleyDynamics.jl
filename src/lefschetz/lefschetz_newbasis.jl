export lefschetz_newbasis, lefschetz_newbasis_maps

"""
    lefschetz_newbasis(lc::AbstractComplex, basis::SparseMatrix; maps::Bool=false)

Create a new Lefschetz complex via change of basis.

The new basis has to be specified in the sparse matrix `basis`, whose columns
represent the new basis in terms of the existing one. This matrix has to
respect the grading by dimension, i.e., the cells which are used to form a
new basis chain have to have the same dimensions as the cell which is being
replaced. The function returns the new Lefschetz complex `lcnew`. If the
optional parameter `maps = true` is passed, the the function also returns
the chain maps `pp` and `jj` which are the isomorphisms from `lc` to `lcnew`,
and vice versa, as well as the zero chain homotopy `hh`.
"""
function lefschetz_newbasis(lc::AbstractComplex, basis::SparseMatrix; maps::Bool=false)
    #
    # Create a new Lefschetz complex via change of basis
    #

    # Make sure the basis change is graded

    cellcount, lo, hi = lefschetz_cell_count(lc, bounds=true)
    lored = setdiff(lo, 0)
    hired = setdiff(hi, 0)

    for k=1:length(lored)
        for m=lored[k]:hired[k]
            if length(basis.columns[m]) == 0
                error("The matrix does not represent a basis!")
            end
            cmin = minimum(basis.columns[m])
            cmax = maximum(basis.columns[m])
            if ~((cmin >= lored[k]) && (cmax <= hired[k]))
                error("The matrix is not a graded map!")
            end
        end
    end

    # Determine the two isomorphisms

    pp = sparse_inverse(basis)
    jj = deepcopy(basis)

    # Compute the new boundary map

    bndnew = pp * (lc.boundary * jj)

    # Construct the new Lefschetz complex

    lcnew = LefschetzComplex(lc.labels, lc.dimensions, bndnew; validate=false)

    # Return the results

    if maps == false
        return lcnew
    else
        n  = lcnew.ncells
        p  = lcnew.boundary.char
        hh = sparse_zero(n, n, p=p)
        return lcnew, pp, jj, hh
    end
end

"""
    lefschetz_newbasis_maps(lc::AbstractComplex, basis::SparseMatrix)

Create a new Lefschetz complex via change of basis and return the
associated chain maps.

The new basis has to be specified in the sparse matrix `basis`, whose columns
represent the new basis in terms of the existing one. This matrix has to
respect the grading by dimension, i.e., the cells which are used to form a
new basis chain have to have the same dimensions as the cell which is being
replaced. The function returns the new Lefschetz complex `lcnew`, as well as
the chain maps `pp` and `jj` which are the isomorphisms from `lc` to `lcnew`,
and vice versa.
"""
lefschetz_newbasis_maps(lc::AbstractComplex, basis::SparseMatrix) =
  lefschetz_newbasis(lc, basis, maps=true)

