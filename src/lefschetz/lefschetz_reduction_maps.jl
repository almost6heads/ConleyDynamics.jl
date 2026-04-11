export lefschetz_reduction_maps
export compose_reductions

"""
    lefschetz_reduction_maps(lc::AbstractComplex, redpairs::Vector{Vector{Int}})

Apply a sequence of elementary reductions to a Lefschetz complex and return
the associated chain maps.

The reduction pairs have to be specified in the argument `redpairs`. Each entry
has to be a vector of length two which contains an elementary reduction pair
in index form. In particular, the dimensions of the two cells in the pair have
to differ by one, and once the pair is reached in the reduction sequence, one cell
has to be a face of the other. The function returns a new Lefschetz complex, where
all cells in `redpairs` have been removed, as well as the associated chain maps.

The return values are as follows:
- `lcred`: The first variable contains the reduced Lefschetz complex.
- `pp`: This is a sparse matrix representation of the chain equivalence
  between the original complex and the reduced one.
- `jj`: This is a sparse matrix representation of the chain equivalence
  between the reduced complex and the original one.
- `hh`: This is a sparse matrix representation of the chain homotopy
  which shows that the composition `jj * pp` is chain homotopic to
  the identity.

# Examples
```jldoctest
julia> labels = ["a","b","c","d"];

julia> simplices = [["a","b"], ["b","c"], ["c","d"]];

julia> sc = create_simplicial_complex(labels, simplices, p=0);

julia> redpairs = [["b", "bc"], ["d", "cd"]];

julia> scr, pp, jj, hh = lefschetz_reduction_maps(sc,redpairs);

julia> bnd = deepcopy(sc.boundary);

julia> bndr = deepcopy(scr.boundary);

julia> ii = sparse_identity(sc.ncells, p=0);

julia> sparse_nz_count(pp*bnd - bndr*pp)
0

julia> sparse_nz_count(jj*bndr - bnd*jj)
0

julia> full_from_sparse(pp*jj)
3×3 Matrix{Rational{Int64}}:
 1  0  0
 0  1  0
 0  0  1

julia> full_from_sparse(jj*pp)
7×7 Matrix{Rational{Int64}}:
 1  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  1  1  1  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  1  0  0
 0  0  0  0  1  0  0
 0  0  0  0  0  0  0

julia> sparse_nz_count(ii + bnd*hh + hh*bnd - jj*pp)
0
```
"""
function lefschetz_reduction_maps(lc::AbstractComplex, redpairs::Vector{Vector{Int}})
    #
    # Apply a sequence of elementary reductions to a Lefschetz complex
    # and compute the associated chain maps
    #
    
    # Make sure the reduction pairs are pairs

    for redpair in redpairs
        if !(length(redpair) == 2)
            error("The reduction pairs have to come in pairs!")
        end
    end

    # Make sure the dimensions of the reduction pairs differ by 1

    for redpair in redpairs
        if !(abs(lc.dimensions[redpair[1]] - lc.dimensions[redpair[2]]) == 1)
            error("The reduction pair dimensions have to differ by 1!")
        end
    end

    # Make sure the reduction pair cells are distinct

    if !(length(unique(reduce(vcat,redpairs))) == 2*length(redpairs))
        error("The reduction pair cells have to be distinct!")
    end

    # Convert the reduction pairs to string format. This makes it
    # easier later, since the meaning of the indices will change
    # as we go through the reduction steps. While doing this,
    # make sure that the lower-dimensional cell comes first.

    redpairsSTR = Vector{Vector{String}}([])

    for redpair in redpairs
        a, b = redpair
        if lc.dimensions[a] < lc.dimensions[b]
            push!(redpairsSTR, [lc.labels[a], lc.labels[b]])
        else
            push!(redpairsSTR, [lc.labels[b], lc.labels[a]])
        end
    end

    # Initialize the chain map matrices

    pp = sparse_identity(lc.ncells, p=lc.boundary.char)
    jj = sparse_identity(lc.ncells, p=lc.boundary.char)
    hh = sparse_zero(lc.ncells, lc.ncells, p=lc.boundary.char)

    # Extract the necessary Lefschetz complex information
    
    labels = deepcopy(lc.labels)
    dims   = deepcopy(lc.dimensions)
    bnd    = deepcopy(lc.boundary)

    # Loop through the elementary reductions

    for redpair in redpairsSTR

        # Extract the current number of cells

        n = length(labels)
        present_cells = Vector{Int}(1:n)

        # Extract the two cell indices

        a = findfirst(x -> x==redpair[1], labels)
        b = findfirst(x -> x==redpair[2], labels)
        removed_cells = [a,b]
        new_cells     = setdiff(present_cells, removed_cells)
        
        # Extract the incidence coefficient

        lambda = bnd[a,b]
        if lambda == 0
            error("Reduction pairs have to consist of a face and coface!")
        end
        lambdainv = scalar_inverse(lambda, bnd.char)

        # Modify the columns of the coboundary of a, except for the
        # removed cells. Since we have modified the boundary matrix,
        # we need to get the row of bnd directly.

        bndn = deepcopy(bnd)
        for c in setdiff(sparse_get_nz_row(bnd,a), removed_cells)
            sparse_add_column!(bndn, c, b, -bnd[a,c], lambda)
        end
        bndn = sparse_minor(bndn, new_cells, new_cells)

        # Compute chain map from original complex to reduced one
       
        pn = sparse_identity(n, p=bnd.char)
        for c in sparse_get_nz_column(bnd,b)
            pn[c,a] = pn[c,a] - lambdainv * bnd[c,b]
        end
        pn = sparse_minor(pn, new_cells, present_cells)
        
        # Compute chain map from reduced complex to original one

        jn = sparse_identity(n, p=bnd.char)
        for c in setdiff(sparse_get_nz_row(bnd,a), removed_cells)
            jn[b,c] = jn[b,c] - lambdainv * bnd[a,c]
        end
        jn = sparse_minor(jn, present_cells, new_cells)

        # Compute the chain homotopy defined on the original complex

        hn = sparse_zero(n, n, p=bnd.char)
        hn[b,a] = -lambdainv

        # Compute the new complex information and chain maps

        labelsn = labels[new_cells]
        dimsn   = dims[new_cells]

        ppn, jjn, hhn = compose_reductions(pp,jj,hh,pn,jn,hn)

        # Update the complex information and chain maps

        labels = deepcopy(labelsn)
        dims   = deepcopy(dimsn)
        bnd    = deepcopy(bndn)

        pp = deepcopy(ppn)
        jj = deepcopy(jjn)
        hh = deepcopy(hhn)
    end

    # Create and return the reduced Lefschetz complex

    lcred = LefschetzComplex(labels, dims, bnd; validate=false)
    return lcred, pp, jj, hh
end

"""
    lefschetz_reduction_maps(lc::AbstractComplex, redpairs::Vector{Vector{String}})

Apply a sequence of elementary reductions to a Lefschetz complex and return
the chain maps.

The reduction pairs have to be specified in the argument `redpairs`. Each entry
has to be a vector of length two which contains an elementary reduction pair
in label form. In particular, the dimensions of the two cells in the pair have
to differ by one, and once the pair is reached in the reduction sequence, one cell
has to be a face of the other. The function returns a new Lefschetz complex, where
all cells in `redpairs` have been removed, as well as the involved chain maps.
"""
function lefschetz_reduction_maps(lc::AbstractComplex, redpairs::Vector{Vector{String}})
    #
    # Apply a sequence of elementary reductions to a Lefschetz complex
    # and return the associated chain maps
    #

    redpairsI = Vector{Vector{Int}}()
    for rp in redpairs
        push!(redpairsI, [lc.indices[rp[1]], lc.indices[rp[2]]])
    end

    lcred, pp, jj, hh = lefschetz_reduction_maps(lc, redpairsI)
    return lcred, pp, jj, hh
end

"""
    lefschetz_reduction_maps(lc::AbstractComplex, r1::Int, r2::Int)

Apply a single elementary reduction to a Lefschetz complex and return
the chain maps.

This method expects that the two cells `r1` and `r2` which form the
reduction pair are given in index form. The function returns the reduced
Lefschetz complex, as well as the involved chain maps.
"""
lefschetz_reduction_maps(lc::AbstractComplex, r1::Int, r2::Int) =
   lefschetz_reduction_maps(lc,[[r1,r2]])

"""
    lefschetz_reduction_maps(lc::AbstractComplex, r1::String, r2::String)

Apply a single elementary reduction to a Lefschetz complex and return
the chain maps.

This method expects that the two cells `r1` and `r2` which form the
reduction pair are given in label form. The function returns the reduced
Lefschetz complex, as well as the involved chain maps.
"""
lefschetz_reduction_maps(lc::AbstractComplex, r1::String, r2::String) =
   lefschetz_reduction_maps(lc,[[r1,r2]])

"""
    compose_reductions(pp1::SparseMatrix, jj1::SparseMatrix, hh1::SparseMatrix, 
                       pp2::SparseMatrix, jj2::SparseMatrix, hh2::SparseMatrix)

Combine the chain maps and chain homotopies of two reductions.

The function expects the chain maps and chain homotopies of two Lefschetz
complex reductions, and computes the chain maps and chain homotopy for the
composition of both. The maps `ppX`, `jjX`, and `hhX` can be obtained via
either `lefschetz_reduction_maps` or `lefschetz_newbasis_maps`.
"""
function compose_reductions(pp1::SparseMatrix, jj1::SparseMatrix, hh1::SparseMatrix, 
                            pp2::SparseMatrix, jj2::SparseMatrix, hh2::SparseMatrix)
    #
    # Combine the chain maps and chain homotopies of two reductions
    #

    ppn = pp2 * pp1
    jjn = jj1 * jj2
    hhn = hh1 + jj1 * (hh2 * pp1)

    return ppn, jjn, hhn
end

