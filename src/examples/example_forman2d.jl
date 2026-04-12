export example_forman2d

"""
    lc, mvf = example_forman2d()

Create the simplicial complex and multivector field
for the example from Figure 3 in the FoCM 2020 paper
by *Batko, Kaczynski, Mrozek, and Wanner*.

The function returns the Lefschetz complex `lc` over the
finite field GF(2) and the multivector field `mvf`. If
desired for plotting, the function can be invoked with
the optional argument `euclidean=true`. In that case
a `EuclidenComplex` with embedded coordinates is
returned for `lc`.

# Examples
```jldoctest
julia> lc, mvf = example_forman2d();

julia> cm = connection_matrix(lc, mvf, algorithm="DHL");

julia> sparse_show(cm)
   ┆   D   E   F  IJ  BF  EF  HI ADE FGJ
---┆------------------------------------
  D┆   .   .   .   .   1   .   1   .   .
  E┆   .   .   .   .   .   1   .   .   .
  F┆   .   .   .   .   1   1   1   .   .
 IJ┆   .   .   .   .   .   .   .   .   1
 BF┆   .   .   .   .   .   .   .   1   .
 EF┆   .   .   .   .   .   .   .   .   .
 HI┆   .   .   .   .   .   .   .   1   .
ADE┆   .   .   .   .   .   .   .   .   .
FGJ┆   .   .   .   .   .   .   .   .   .

julia> sparse_show(cm.matrix)
 . . . . 1 . 1 . .
 . . . . . 1 . . .
 . . . . 1 1 1 . .
 . . . . . . . . 1
 . . . . . . . 1 .
 . . . . . . . . .
 . . . . . . . 1 .
 . . . . . . . . .
 . . . . . . . . .

julia> print(cm.labels)
["D", "E", "F", "IJ", "BF", "EF", "HI", "ADE", "FGJ"]

julia> lc2, mvf2 = example_forman2d(euclidean=true);

julia> typeof(lc2)
EuclideanComplex
```
"""
function example_forman2d(; euclidean::Bool=false)
    # Create the simplicial complex and multivector field for the example
    # from Figure 3 in the FoCM 2020 paper by Batko, Kaczynski, Mrozek,
    # and Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    intsimplices = [[1,4,5],[1,2,5],[2,5,6],[2,3,6],[3,6,7],
                    [4,5,8],[5,8,9],[5,6,9],[6,9,10],[6,7,10]]
    strsimplices = convert_simplices(intsimplices, labels)

    # Create the coordinates vector

    coords = [[50,200],[150,200],[250,200],
              [0,100],[100,100],[200,100],[300,100],
              [50,0],[150,0],[250,0]]

    # Create the correct complex

    if !euclidean
        lc = create_simplicial_complex(labels, strsimplices, p=2)
    else
        lc = create_simplicial_complex(labels, strsimplices, coords, p=2)
    end

    # Create the multivector field

    mvf = [["A","AD"], ["B","AB"], ["C","BC"], ["F","FG"], ["G","GJ"],
           ["H","DH"], ["I","FI"], ["J","IJ"],
           ["AE","ABE"], ["BE","BEF"], ["CF","BCF"], ["CG","CFG"],
           ["DE","DEH"], ["EH","EHI"], ["EI","EFI"], ["FJ","FIJ"]]

    # Return the example data

    return lc, mvf
end

