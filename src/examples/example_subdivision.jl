export example_subdivision

"""
    lc, mvf = example_subdivision(mvftype)

Create the Lefschetz complex and multivector field for the example
from Figure 7.1 in the connection matrix book by *Mrozek & Wanner*.

Depending on the value of `mvftype`, return the multivector (0=default)
or one of the two combinatorial vector field (1,2) examples.

The function returns the Lefschetz complex `lc` over the rationals
and the multivector field `mvf`.

# Examples
```jldoctest
julia> lc, mvf = example_subdivision(1);

julia> cm = connection_matrix(lc, mvf, algorithm="DHL");

julia> full_from_sparse(cm.matrix)
5×5 Matrix{Rational{Int64}}:
 0  0  -1  -1  -1
 0  0   1   0   0
 0  0   0   0   0
 0  0   0   0   0
 0  0   0   0   0
```
"""
function example_subdivision(mvftype=0)
    #
    # Create the combinatorial vector field information for the example in
    # Figure 11 of the connection matrix paper by Mrozek & Wanner. The function
    # returns the underlying simplicial complex and the combinatorial vector
    # field given by the input variable mvftype. The latter can take the value 0
    # for the left-most example, or 1,2 for the two combinatorial vector field
    # refinements shown in Figure 11.
    
    nc = 9

    # Create the vector of labels

    labelvec = Vector{String}()
    push!(labelvec,"A","B","C")                # 1, 2, 3
    push!(labelvec,"AB","AC","BC","CD","CE")   # 4, 5, 6, 7, 8
    push!(labelvec,"ABC")                      # 9

    # Create the vector of simplex dimensions
    
    sdvec = [0, 0, 0, 1, 1, 1, 1, 1, 2]

    # Create the boundary matrix

    bndmatrix = zeros(Int, nc, nc)

    # Start with the edges

    bndmatrix[[1,2],4] = [-1; 1]     # AB
    bndmatrix[[1,3],5] = [-1; 1]     # AC
    bndmatrix[[2,3],6] = [-1; 1]     # BC
    bndmatrix[3,7] = -1              # CD
    bndmatrix[3,8] = -1              # CE

    # Move on to the triangle

    bndmatrix[[4,5,6],9] = [1; -1; 1]     # ABC 

    # Construct the Lefschetz complex struct
    
    lc = LefschetzComplex(labelvec, sdvec,
                          sparse_from_full(bndmatrix, p=0))

    # Create the common part of the combinatorial vector fields
    
    mvf = Vector{Vector{String}}()
    if mvftype==0
        push!(mvf,["C","AC","BC","ABC"])
    elseif mvftype==1
        push!(mvf,["C","AC"])
        push!(mvf,["BC","ABC"])
    elseif mvftype==2
        push!(mvf,["C","BC"])
        push!(mvf,["AC","ABC"])
    else
        error("mvftype must be 0, 1, or 2!")
    end

    # Return the example data

    return lc, mvf
end

