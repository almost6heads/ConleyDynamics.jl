export example_julia_logo

"""
    lc, mvf = example_julia_logo()

Create the simplicial complex and multivector field
for the example from Figure 2.1 in the connection matrix
book by *Mrozek & Wanner*.

The function returns the Lefschetz complex `lc` over GF(2)
and the multivector field `mvf`.

# Examples
```jldoctest
julia> lc, mvf = example_julia_logo();

julia> cm = connection_matrix(lc, mvf, algorithm="DHL");

julia> sparse_show(cm)
   ┆   D  BC ABC
---┆------------
  D┆   .   .   .
 BC┆   .   .   1
ABC┆   .   .   .

julia> print(cm.labels)
["D", "BC", "ABC"]
```
"""
function example_julia_logo()
    # Create the simplicial complex and multivector field for the example
    # from Figure 1 in the connection matrix paper by Mrozek & Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D"]
    intsimplices = [[1,2,3],[2,3,4]]
    lc = create_simplicial_complex(labels, intsimplices, p=2)

    # Create the multivector field

    mvf = [["A","AB"], ["C","AC"], ["B","BC", "BD", "BCD"]]

    # Return the example data

    return lc, mvf
end

