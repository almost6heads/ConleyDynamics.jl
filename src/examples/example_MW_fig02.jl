export example_MW_fig02

"""
    lc, mvf = example_MW_fig02(mvftype)

Create the simplicial complex and multivector field for the example
from Figure 2 in the connection matrix paper by *Mrozek & Wanner*.

Depending on the value of `mvftype`, return the periodic orbit (0=default)
or one of the three gradient (1,2,3) examples.

# Examples
```jldoctest
julia> lc, mvf = example_MW_fig02(0);

julia> cm = connection_matrix(lc, mvf);

julia> cm.cm
[0   0   0   1   1   0   0   0]
[0   0   0   1   1   0   0   0]
[0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   1   0]
[0   0   0   0   0   0   1   0]
[0   0   0   0   0   0   0   1]
[0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0]

julia> print(cm.labels)
["A", "C", "CE", "AC", "BD", "DF", "ABC", "EFG"]
```
"""
function example_MW_fig02(mvftype=0)
    #
    # Create the combinatorial vector field information for the example in
    # Figure 2 of the connection matrix paper by Mrozek & Wanner. The function
    # returns the underlying simplicial complex and the combinatorial vector
    # field given by the input variable mvftype. The latter can take the value 0
    # for the example with periodic orbit, or 1,2,3 for the three gradient
    # vector field refinements shown in Figure 6.
    
    n0 =  7     # number of vertices
    n1 = 11     # number of edges
    n2 =  4     # number of triangles
    nc = n0 + n1 + n2

    # Create the vector of labels

    labelvec = Vector{String}()
    push!(labelvec,"A","B","C","D","E","F","G")
    push!(labelvec,"AB","AC","BC","BD","CD","CE","DE","DF","EF","EG","FG")
    push!(labelvec,"ABC","BCD","DEF","EFG")

    # Create the label to index dictionary

    indexdict = Dict{String,Int}([(labelvec[k],k) for k in 1:length(labelvec)])

    # Create the vector of Poincare polynomials
    
    PP, t = ZZ["t"]
    ppvec = [t^(length(labelvec[k])-1) for k=1:length(labelvec)]

    # Create the boundary matrix

    bndmatrix = zeros(Int, nc, nc)

    # Start with the edges

    bndmatrix[[1,2], 8] = [-1; 1]     # AB
    bndmatrix[[1,3], 9] = [-1; 1]     # AC
    bndmatrix[[2,3],10] = [-1; 1]     # BC
    bndmatrix[[2,4],11] = [-1; 1]     # BD
    bndmatrix[[3,4],12] = [-1; 1]     # CD
    bndmatrix[[3,5],13] = [-1; 1]     # CE
    bndmatrix[[4,5],14] = [-1; 1]     # DE
    bndmatrix[[4,6],15] = [-1; 1]     # DF
    bndmatrix[[5,6],16] = [-1; 1]     # EF
    bndmatrix[[5,7],17] = [-1; 1]     # EG
    bndmatrix[[6,7],18] = [-1; 1]     # FG

    # Move on to the triangles

    bndmatrix[[ 8, 9,10],19] = [1; -1; 1]     # ABC 
    bndmatrix[[10,11,12],20] = [1; -1; 1]     # BCD
    bndmatrix[[14,15,16],21] = [1; -1; 1]     # DEF
    bndmatrix[[16,17,18],22] = [1; -1; 1]     # EFG

    # Construct the Lefschetz complex struct
    
    lc = LefschetzComplex{typeof(t)}(nc, bndmatrix,
                                     labelvec, indexdict, ppvec)

    # Create the common part of the combinatorial vector fields
    
    mvf = Vector{Vector{Int}}()
    push!(mvf,[2,8])       # B - AB
    push!(mvf,[6,18])      # F - FG
    push!(mvf,[7,17])      # G - EG
    push!(mvf,[10,20])     # BC - BCD
    push!(mvf,[16,21])     # EF - DEF
    
    if mvftype==0
        push!(mvf,[3,12])       # C - CD
        push!(mvf,[4,14])       # D - DE
        push!(mvf,[5,13])       # E - CE
    elseif mvftype==1
        push!(mvf,[4,14])       # D - DE
        push!(mvf,[5,13])       # E - CE
    elseif mvftype==2
        push!(mvf,[3,12])       # C - CD
        push!(mvf,[5,13])       # E - CE
    elseif mvftype==3
        push!(mvf,[3,12])       # C - CD
        push!(mvf,[4,14])       # D - DE
    else
        error("mvftype must be 0, 1, 2, or 3!")
    end

    strmvf = convert_mvf(mvf, lc)

    # Return the example data

    return lc, strmvf
end

