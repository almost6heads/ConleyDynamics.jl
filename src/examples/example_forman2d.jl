export example_forman2d

"""
    lc, mvf, coords = example_forman2d()

Create the simplicial complex and multivector field
for the example from Figure 3 in the FoCM 2020 paper
by *Batko, Kaczynski, Mrozek, and Wanner*.

The function returns the Lefschetz complex `lc` over the
finite field GF(2) and the multivector field `mvf`. If
desired for plotting, the third return value `coords`
gives a vector of coordinates for the vertices.

# Examples
```jldoctest
julia> lc, mvf = example_forman2d();

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show(cm)
   ┆   D   E   F  GJ  BF  EF  HI ADE FGJ
---┆------------------------------------
  D┆   0   0   0   0   1   0   1   0   0
  E┆   0   0   0   0   0   1   0   0   0
  F┆   0   0   0   0   1   1   1   0   0
 GJ┆   0   0   0   0   0   0   0   0   1
 BF┆   0   0   0   0   0   0   0   1   0
 EF┆   0   0   0   0   0   0   0   0   0
 HI┆   0   0   0   0   0   0   0   1   0
ADE┆   0   0   0   0   0   0   0   0   0
FGJ┆   0   0   0   0   0   0   0   0   0

julia> sparse_show(cm.matrix)
 0 0 0 0 1 0 1 0 0
 0 0 0 0 0 1 0 0 0
 0 0 0 0 1 1 1 0 0
 0 0 0 0 0 0 0 0 1
 0 0 0 0 0 0 0 1 0
 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 1 0
 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0

julia> print(cm.labels)
["D", "E", "F", "GJ", "BF", "EF", "HI", "ADE", "FGJ"]
```
"""
function example_forman2d()
    # Create the simplicial complex and multivector field for the example
    # from Figure 3 in the FoCM 2020 paper by Batko, Kaczynski, Mrozek,
    # and Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    intsimplices = [[1,4,5],[1,2,5],[2,5,6],[2,3,6],[3,6,7],
                    [4,5,8],[5,8,9],[5,6,9],[6,9,10],[6,7,10]]
    strsimplices = convert_simplices(intsimplices, labels)
    lc = create_simplicial_complex(labels, strsimplices, p=2)

    # Create the coordinates vector

    coords = [[50,200],[150,200],[250,200],
              [0,100],[100,100],[200,100],[300,100],
              [50,0],[150,0],[250,0]]

    # Create the multivector field

    mvf = [["A","AD"], ["B","AB"], ["C","BC"], ["F","FG"], ["G","GJ"],
           ["H","DH"], ["I","FI"], ["J","IJ"],
           ["AE","ABE"], ["BE","BEF"], ["CF","BCF"], ["CG","CFG"],
           ["DE","DEH"], ["EH","EHI"], ["EI","EFI"], ["FJ","FIJ"]]

    # Return the example data

    return lc, mvf, coords
end

