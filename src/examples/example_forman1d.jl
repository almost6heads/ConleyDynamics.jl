export example_forman1d

"""
    lc, mvf = example_forman1d()

Create the simplicial complex and multivector field
for the example from Figure 1 in the FoCM 2020 paper
by *Batko, Kaczynski, Mrozek, and Wanner*.

The function returns the Lefschetz complex `lc` and the 
multivector field `mvf`. The Lefschetz complex is defined
over the finite field GF(2). If desired for plotting,
the function can be invoked with the optional argument
`euclidean=true`. In that case a `EuclidenComplex`
with embedded coordinates is returned for `lc`.

# Examples
```jldoctest
julia> lc, mvf = example_forman1d();

julia> cm = connection_matrix(lc, mvf, algorithm="DHL");

julia> sparse_show(cm)
  ┆  A CD  F BF DE
--┆---------------
 A┆  .  .  .  .  1
CD┆  .  .  .  .  .
 F┆  .  .  .  .  1
BF┆  .  .  .  .  .
DE┆  .  .  .  .  .

julia> sparse_show(cm.matrix)
 . . . . 1
 . . . . .
 . . . . 1
 . . . . .
 . . . . .

julia> println(cm.labels)
["A", "CD", "F", "BF", "DE"]

julia> lc2, mvf2 = example_forman1d(euclidean=true);

julia> typeof(lc2)
EuclideanComplex

julia> lc2.coords
13-element Vector{Vector{Vector{Float64}}}:
 [[50.0, 100.0]]
 [[250.0, 100.0]]
 [[0.0, 0.0]]
 [[100.0, 0.0]]
 [[200.0, 0.0]]
 [[300.0, 0.0]]
 [[50.0, 100.0], [0.0, 0.0]]
 [[50.0, 100.0], [100.0, 0.0]]
 [[250.0, 100.0], [200.0, 0.0]]
 [[250.0, 100.0], [300.0, 0.0]]
 [[0.0, 0.0], [100.0, 0.0]]
 [[100.0, 0.0], [200.0, 0.0]]
 [[200.0, 0.0], [300.0, 0.0]]
```
"""
function example_forman1d(; euclidean::Bool=false)
    # Create the simplicial complex and multivector field for the example
    # from Figure 1 in the FoCM 2020 paper by Batko, Kaczynski, Mrozek,
    # and Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D", "E", "F"]
    intsimplices = [[1,3],[1,4],[3,4],[2,5],[2,6],[5,6],[4,5]]

    # Create coordinates for plotting

    coords = [[50,100],[250,100],[0,0],[100,0],[200,0],[300,0]]

    # Create the correct complex

    if !euclidean
        lc = create_simplicial_complex(labels, intsimplices, p=2)
    else
        lc = create_simplicial_complex(labels, intsimplices, coords, p=2)
    end

    # Create the multivector field

    mvf = [["A","AD"], ["D","CD"], ["C","AC"], ["B","BE"], ["E","EF"]]

    # Return the example data

    return lc, mvf
end

