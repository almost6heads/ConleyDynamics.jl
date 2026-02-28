export example_small_periodicity

"""
    lc1, lc2, mvf = example_small_periodicity()

Create two representations of the Lefschetz complex and
the multivector field for the example from Figure 4 in
the connection matrix paper by *Mrozek & Wanner*.

The complexes `lc1` and `lc2` are just two representations
of the same complex, but they lead to different connection
matrices. Both Lefschetz complexes are defined over the
finite field GF(2).

The function returns the Lefschetz complexes `lc1` and `lc2`,
as well as the multivector field `mvf`.

# Examples
```jldoctest
julia> lc1, lc2, mvf = example_small_periodicity();

julia> cm1 = connection_matrix(lc1, mvf, algorithm="DLMS24");

julia> cm2 = connection_matrix(lc2, mvf, algorithm="DLMS24");

julia> full_from_sparse(cm1.matrix)
4×4 Matrix{Int64}:
 0  0  0  0
 0  0  0  1
 0  0  0  1
 0  0  0  0

julia> print(cm1.labels)
["A", "a", "b", "alpha"]

julia> full_from_sparse(cm2.matrix)
4×4 Matrix{Int64}:
 0  0  0  0
 0  0  0  0
 0  0  0  1
 0  0  0  0

julia> print(cm2.labels)
["A", "c", "b", "alpha"]
```
"""
function example_small_periodicity()
    #
    # Create two representations of the Lefschetz complex
    # and the multivector field for the example from 
    # Figure 4 in the connection matrix paper
    # by *Mrozek & Wanner*.
    
    nc = 6

    # Create the vector of labels

    labelvec = Vector{String}()
    push!(labelvec,"A","B","a","b","c","alpha")

    # Create the vector of simplex dimensions
    
    sdvec = [0, 0, 1, 1, 1, 2]

    # Create the boundary matrix

    bndmatrix = zeros(Int, nc, nc)

    # Start with the edges

    bndmatrix[[1,2],3] = [1; 1]     # a
    bndmatrix[[1,2],4] = [1; 1]     # b
    bndmatrix[[1,2],5] = [1; 1]     # c

    # Move on to the 2-cell

    bndmatrix[[3,4],6] = [1; 1]     # alpha 

    # Construct the Lefschetz complex struct
    
    lc = LefschetzComplex(labelvec, sdvec,
                          sparse_from_full(bndmatrix, p=2))

    # Create a second version of the Lefschetz complex via permutation

    perm = [1, 2, 5, 3, 4, 6]      # Change a, b, c to the order c, a, b
    labelvec2  = labelvec[perm]
    sdvec2     = sdvec[perm]
    bndmatrix2 = bndmatrix[perm,perm]

    lc2 = LefschetzComplex(labelvec2, sdvec2,
                           sparse_from_full(bndmatrix2, p=2))

    # Create the common part of the combinatorial vector fields
    
    mvf = Vector{Vector{String}}()
    push!(mvf,["A", "a"])       # A - a
    push!(mvf,["B", "c"])       # B - c

    # Return the example data

    return lc, lc2, mvf
end

