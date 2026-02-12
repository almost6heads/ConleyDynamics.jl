export example_dunce_chaos

"""
    sc, vfG, vfC = example_dunce_chaos()

Create a minimal simplicial complex representation
of the Dunce hat, as well as two Forman vector fields.

The function returns the simplicial representation of
the Dunce hat in `sc` over the finite field `GF(2)`.
The Forman vector field `vfG` is a gradient vector 
field with unique connection matrix. The field `vfC`
is a modification of this field which merges the critical
cells of dimensions 1 and 2 into a Forman arrow. The
resulting Forman vector field is no longer gradient, 
and in fact exhibits Lorez-like chaos.

# Examples
```jldoctest
julia> sc, vfG, vfC = example_dunce_chaos();

julia> homology(sc)
3-element Vector{Int64}:
 1
 0
 0

julia> cmG = connection_matrix(sc, vfG);

julia> sparse_show(cmG.matrix)
 0 0 0
 0 0 1
 0 0 0

julia> print(cmG.labels)
["1", "12", "128"]

julia> cmC = connection_matrix(sc, vfC);

julia> sparse_show(cmC.matrix)
 0

julia> print(cmC.labels)
["1"]

julia> msC, psC = morse_sets(sc, vfC, poset=true);

julia> [conley_index(sc, mset) for mset in msC]
2-element Vector{Vector{Int64}}:
 [1, 0, 0]
 [0, 0, 0]

julia> psC
2×2 Matrix{Bool}:
 0  1
 0  0

julia> msC
2-element Vector{Vector{String}}:
 ["1"]
 ["12", "14", "15", "25", "28", "56", "68", "78", "124", "125", "128", "145", "256", "278", "568", "678"]
```
"""
function example_dunce_chaos()
    # Create a minimal simplicial complex representation
    # of the Dunce hat, as well as two Forman vector fields.
    
    # Create the simplicial complex

    labels = ["1","2","3","4","5","6","7","8"]
    simplices = [["1","3","8"],["1","2","8"],["2","3","4"],["3","4","8"],
                 ["2","7","8"],["2","3","7"],["1","2","4"],["4","5","8"],
                 ["5","6","8"],["6","7","8"],["1","3","7"],["1","4","5"],
                 ["1","6","7"],["1","2","5"],["2","5","6"],["2","3","6"],
                 ["1","3","6"]]

    sc = create_simplicial_complex(labels, simplices, p=2)

    # Create the gradient multivector field

    vfG = [["3","13"],["2","23"],["8","18"],["4","48"],["5","45"],
           ["6","26"],["7","67"],
           ["38","138"],["34","348"],["24","234"],["14","124"],
           ["15","145"],["25","125"],["56","256"],["58","458"],
           ["68","568"],["36","236"],["16","136"],["17","167"],
           ["37","137"],["78","678"],["27","237"],["28","278"]]

    # Create the chaotic multivector field

    vfC = [["3","13"],["2","23"],["8","18"],["4","48"],["5","45"],
           ["6","26"],["7","67"],
           ["38","138"],["34","348"],["24","234"],["14","124"],
           ["15","145"],["25","125"],["56","256"],["58","458"],
           ["68","568"],["36","236"],["16","136"],["17","167"],
           ["37","137"],["78","678"],["27","237"],["28","278"],
           ["12","128"]]

    # Return the example data

    return sc, vfG, vfC
end

