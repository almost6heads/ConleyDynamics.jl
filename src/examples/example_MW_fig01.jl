export example_MW_fig01

"""
    lcomplex, mvf = example_MW_fig01()

Create the simplicial complex and multivector field
for the example from Figure 1 in the connection matrix
paper by *Mrozek & Wanner*.

# Examples
```jldoctest
julia> lc, mvf = example_MW_fig01();

julia> cm = connection_matrix(lc, mvf);

julia> cm.cm
[0   0   0]
[0   0   1]
[0   0   0]

julia> print(cm.labels)
["D", "AC", "ABC"]
```
"""
function example_MW_fig01()
    # Create the simplicial complex and multivector field for the example
    # from Figure 1 in the connection matrix paper by Mrozek & Wanner.
    
    # Create the simplicial complex

    labels = ["A", "B", "C", "D"]
    intsimplices = [[1,2,3],[2,3,4]]
    lc = create_simplicial_complex(labels, intsimplices)

    # Create the multivector field

    mvf = [["A","AB"], ["C","AC"], ["B","BC", "BD", "BCD"]]

    # Return the example data

    return lc, mvf
end
