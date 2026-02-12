export example_forman1d

"""
    lc, mvf, coords = example_forman1d()

Create the simplicial complex and multivector field
for the example from Figure 1 in the FoCM 2020 paper
by *Batko, Kaczynski, Mrozek, and Wanner*.

The function returns the Lefschetz complex `lc` and the 
multivector field `mvf`. If desired for plotting, the third
return value `coords` gives a vector of coordinates for the
vertices. The Lefschetz complex is defined over the finite
field GF(2).

# Examples
```jldoctest
julia> lc, mvf = example_forman1d();

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show(cm)
  ┆  A AD  F BF DE
--┆---------------
 A┆  0  0  0  0  1
AD┆  0  0  0  0  0
 F┆  0  0  0  0  1
BF┆  0  0  0  0  0
DE┆  0  0  0  0  0

julia> sparse_show(cm.matrix)
 0 0 0 0 1
 0 0 0 0 0
 0 0 0 0 1
 0 0 0 0 0
 0 0 0 0 0

julia> println(cm.labels)
["A", "AD", "F", "BF", "DE"]
```
"""
function example_forman1d()
    # Create the simplicial complex and multivector field for the example
    # from Figure 1 in the FoCM 2020 paper by Batko, Kaczynski, Mrozek,
    # and Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D", "E", "F"]
    intsimplices = [[1,3],[1,4],[3,4],[2,5],[2,6],[5,6],[4,5]]
    lc = create_simplicial_complex(labels, intsimplices, p=2)

    # Create coordinates for plotting

    coords = [[50,100],[250,100],[0,0],[100,0],[200,0],[300,0]]

    # Create the multivector field

    mvf = [["A","AD"], ["D","CD"], ["C","AC"], ["B","BE"], ["E","EF"]]

    # Return the example data

    return lc, mvf, coords
end

