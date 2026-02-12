export example_critical_simplex

"""
    lc, mvf = example_critical_simplex(dim)

Create a simplicial complex of dimension `dim` as well as
a multivector field on it in which every cell is critical.

The function returns the Lefschetz complex `lc` over GF(2)
and the multivector field `mvf`.

# Examples
```jldoctest
julia> lc, mvf = example_critical_simplex(2);

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show(cm.matrix)
 0 0 0 1 1 0 0
 0 0 0 1 0 1 0
 0 0 0 0 1 1 0
 0 0 0 0 0 0 1
 0 0 0 0 0 0 1
 0 0 0 0 0 0 1
 0 0 0 0 0 0 0

julia> print(cm.labels)
["A", "B", "C", "AB", "AC", "BC", "ABC"]
```
"""
function example_critical_simplex(dim)
    # Create a simplicial complex of dimension dim as well as
    # a multivector field on it in which every cell is critical.

    # Create the vector of labels and the simplicial complex

    labels = string.(Char.([Int('A') + k for k=0:dim]))
    intsimplices = Vector{Vector{Int}}()
    push!(intsimplices,[k for k=1:dim+1])
    lc = create_simplicial_complex(labels, intsimplices, p=2)

    # Create the multivector field
    
    mvf = Vector{Vector{Int}}()

    # Return the example data

    return lc, mvf
end

