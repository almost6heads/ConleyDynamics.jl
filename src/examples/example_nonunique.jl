export example_nonunique

"""
    lc1, lc2, mvf, coords1, coords2 = example_nonunique()

Create two representations of a simplicial complex and one multivector
field which illustrates nonunique connection matrices.

The two complexes `lc1` and `lc2` represent the same simplicial
complex over GF(2), but differ in the ordering of the labels.

The function returns the Lefschetz complexes `lc1` and `lc2`, as well
as the multivector field `mvf`. If desired for plotting, the fourth
and fifth return values `coords1` and `coords2` give vectors of
coordinates for the vertices of the two complexes.

# Examples
```jldoctest
julia> lc1, lc2, mvf = example_nonunique();

julia> cm1 = connection_matrix(lc1, mvf, algorithm="DLMS24");

julia> cm2 = connection_matrix(lc2, mvf, algorithm="DLMS24");

julia> sparse_show(cm1.matrix)
 0 0 0 1 0 1 0 0 0
 0 0 0 1 0 1 0 0 0
 0 0 0 0 0 0 0 1 1
 0 0 0 0 0 0 1 1 0
 0 0 0 0 0 0 0 1 0
 0 0 0 0 0 0 1 1 0
 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0

julia> print(cm1.labels)
["2", "7", "79", "29", "45", "67", "168", "349", "789"]

julia> sparse_show(cm2.matrix)
 0 0 0 1 0 1 0 0 0
 0 0 0 1 0 1 0 0 0
 0 0 0 0 0 0 1 0 1
 0 0 0 0 0 0 1 1 0
 0 0 0 0 0 0 0 1 0
 0 0 0 0 0 0 1 1 0
 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0

julia> print(cm2.labels)
["2", "8", "78", "29", "45", "67", "168", "349", "789"]
```
"""
function example_nonunique()
    # Create two representations of a simplicial complex and one
    # multivector field which illustrates nonunique connection matrices

    # Initialize the simplicial complex in natural order
    
    labels = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    intsimplices = [[1,2,8],[2,8,9],[2,3,9],
                    [1,6,8],[6,7,8],[7,8,9],[4,7,9],[3,4,9],
                    [5,6,7],[4,5,7]]
    strsimplices = convert_simplices(intsimplices, labels)

    coords1 = [[0,100],[150,0],[300,100],[250,200],[150,300],
               [50,200],[150,200],[100,100],[200,100]]

    # Create the multivector field using string labels

    mvf = [["1","12"], ["3","23"], ["4","34"], ["5","56"], ["6","16"],
           ["7","79"], ["9","89"], ["8","78"],
           ["18","128"], ["28","289"], ["39","239"], ["49","479"],
           ["47","457"], ["57","567"], ["68","678"]]

    # Create the first simplicial complex
    
    labels1 = labels
    lc1 = create_simplicial_complex(labels1, strsimplices, p=2)

    # Create the second simplicial complex, via vertex reordering

    labels2 = ["1", "2", "3", "4", "5", "6", "8", "9", "7"]
    lc2 = create_simplicial_complex(labels2, strsimplices, p=2)

    coords2 = [[0,100],[150,0],[300,100],[250,200],[150,300],
               [50,200],[100,100],[200,100],[150,200]]

    # Return the example data

    return lc1, lc2, mvf, coords1, coords2
end

