export example_clorenz

"""
    lc, mvf = example_clorenz()

Create the simplicial complex and multivector field
for the example from Figure 3 in the JCD 2016 paper
by Kaczynski, Mrozek, and Wanner.

The function returns the Lefschetz complex `lc` over the
finite field GF(2) and the multivector field `mvf`.

# Examples
```jldoctest
julia> lc, mvf = example_clorenz();

julia> cm = connection_matrix(lc, mvf);

julia> sparse_show(cm.matrix)
 0 0 0 0 1
 0 0 0 0 0
 0 0 0 0 1
 0 0 0 0 0
 0 0 0 0 0

julia> print(cm.labels)
["i", "ip", "g", "gm", "bc"]

julia> ms, ps = morse_sets(lc, mvf, poset=true);

julia> [conley_index(lc, mset) for mset in ms]
4-element Vector{Vector{Int64}}:
 [1, 1, 0]
 [1, 1, 0]
 [0, 1, 0]
 [0, 0, 0]

julia> ps
4×4 Matrix{Bool}:
 0  0  1  0
 0  0  1  0
 0  0  0  1
 0  0  0  0
```
"""
function example_clorenz()
    # Create the simplicial complex and multivector field for the
    # combinatorial Lorenz example from Figure 3 in the JCD 2016
    # paper by Kaczynski, Mrozek, and Wanner.
    
    # Create the simplicial complex

    labels = ["a","b","c","d","e","f","g","h","i","j",
              "k","l","m","n","o","p","q","r","s"]
    simplices = [["l","m","r"],["m","n","r"],["f","l","m"],["f","g","m"],
                 ["h","m","n"],["a","f","g"],["a","g","h"],["a","b","h"],
                 ["o","p","s"],["p","q","s"],["i","o","p"],["j","k","p"],
                 ["k","p","q"],["d","e","i"],["e","i","j"],["e","j","k"],
                 ["b","c","h"],["c","d","i"],["c","h","i"],
                 ["h","i","n"],["h","i","o"]]

    lc = create_simplicial_complex(labels, simplices, p=2)

    # Create the multivector field

    mvf = [["l","lr"],["lm","lmr"],["r","nr"],["mr","mnr"],["f","fl"],
           ["fm","flm"],["fg","fgm"],["g","gm"],["m","hm"],["mn","hmn"],
           ["a","af"],["ag","afg"],["ah","agh"],["h","gh"],["b","ab"],
           ["bh","abh"],["hn","hin"],["n","in"],["o","ho"],["io","hio"],
           ["hi","chi"],["ch","bch"],["ci","cdi"],["c","cd"],["d","de"],
           ["di","dei"],["ei","eij"],["i","ij"],["e","ek"],["ej","ejk"],
           ["j","jp"],["jk","jkp"],["k","kq"],["kp","kpq"],["q","qs"],
           ["pq","pqs"],["s","os"],["ps","ops"],["p","ip"],["op","iop"]]

    # Return the example data

    return lc, mvf
end

