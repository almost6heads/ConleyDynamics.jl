export example_MW_fig03

"""
    lc, mvf = example_MW_fig03()

Create the Lefschetz complex and multivector field for the example
from Figure 3 in the connection matrix paper by *Mrozek & Wanner*.

# Examples
```jldoctest
julia> lc, mvf = example_MW_fig03();

julia> cm = connection_matrix(lc, mvf);

julia> cm.cm
[0   0   0   0]
[0   0   0   0]
[0   0   0   0]
[0   0   0   0]

julia> print(cm.labels)
["BD", "DF", "AC", "CE"]
```
"""
function example_MW_fig03()
    #
    # Create the multivector field information for the example in
    # Figure 3 of the connection matrix paper by Mrozek & Wanner.
    # The function returns the underlying Lefschetz complex and
    # the multivector field.
    
    nc = 10

    # Create the vector of labels

    labelvec = Vector{String}()
    push!(labelvec,"C")                   # 1
    push!(labelvec,"AC","BC","BD","BF")   # 2, 3, 4, 5
    push!(labelvec,"CE","CF","DF")        # 6, 7, 8
    push!(labelvec,"BCF","BDF")           # 9, 10

    # Create the label to index dictionary

    indexdict = Dict{String,Int}([(labelvec[k],k) for k in 1:length(labelvec)])

    # Create the vector of Poincare polynomials
    
    PP, t = ZZ["t"]
    ppvec = [t^0, t, t, t, t, t, t, t, t^2, t^2]

    # Create the boundary matrix

    bndmatrix = zeros(Int, nc, nc)

    # Start with the edges

    bndmatrix[1,2] = 1     # AC
    bndmatrix[1,3] = 1     # BC
    bndmatrix[1,6] = 1     # CE
    bndmatrix[1,7] = 1     # CF

    # Move on to the triangle

    bndmatrix[[3,5,7], 9] = [1; 1; 1]   # BCF
    bndmatrix[[4,5,8],10] = [1; 1; 1]   # BDF

    # Construct the Lefschetz complex struct
    
    lcf = LefschetzComplex{typeof(t)}(nc, bndmatrix,
                                      labelvec, indexdict, ppvec)
    lc = convert_lefschetz_sparse(lcf)

    # Create the common part of the combinatorial vector fields
    
    mvf = Vector{Vector{String}}()
    push!(mvf,["C","BC","CF","BCF"])
    push!(mvf,["BF","BDF"])

    # Return the example data

    return lc, mvf
end
