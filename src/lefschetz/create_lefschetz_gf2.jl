export create_lefschetz_gf2

"""
    create_lefschetz_gf2(defcellbnd)

Create a Lefschetz complex over `GF(2)` by specifying its essential
cells and boundaries.

The input argument `defcellbnd` has to be a vector of vectors. Each
entry `defcellbnd[k]` has to be of one of the following two forms:

* `[String, Int, String, String, ...]`: The first `String` contains
   the label for the cell `k`, followed by its dimension in the second
   entry. The remaining entries are for the labels of the cells which
   make up the boundary.
* `[String, Int]`: This shorther form is for cells with empty boundary.
  The first entry denotes the cell label, and the second its dimension.

The cells of the resulting Lefschetz complex correspond to the union of
all occurring labels. Cell labels that only occur in the boundary
specification are assumed to have empty boundary, and they do not have
to be specified separately in the second form above. However, if their
boundary is not empty, they have to be listed via the above first
form as well.

# Examples
```jldoctest
julia> defcellbnd = [["A",0], ["a",1,"B","C"], ["b",1,"B","C"]];

julia> push!(defcellbnd, ["c",1,"B","C"]);

julia> push!(defcellbnd, ["alpha",2,"b","c"]);

julia> lc = create_lefschetz_gf2(defcellbnd);

julia> lc.labels
7-element Vector{String}:
 "A"
 "B"
 "C"
 "a"
 "b"
 "c"
 "alpha"

julia> homology(lc)
3-element Vector{Int64}:
 2
 1
 0
```
"""
function create_lefschetz_gf2(defcellbnd)
    #
    # Create a Lefschetz complex struct over GF(2).
    #

    # Determine the dimension of the complex

    defcelllen = length(defcellbnd)
    defcelldim = [defcellbnd[k][2] for k in 1:defcelllen]

    if minimum(defcelldim) < 0
        error("Dimension of the Lefschetz complex has to be positive!")
    else
        lcdim = maximum(defcelldim)
    end

    if !(typeof(lcdim) == Int)
        error("Dimension of the Lefschetz complex has to be an integer!")
    end

    # Collect the labels for all cells of the complex

    labelsets = [Set{String}() for _ in 0:lcdim]
    
    for k in 1:defcelllen
        kbnd = defcellbnd[k]
        kdim = kbnd[2]
        push!(labelsets[kdim+1], kbnd[1])
        if length(kbnd) > 2
            for m=3:length(kbnd)
                push!(labelsets[kdim], kbnd[m])
            end
        end
    end

    # Order the cell labels and create label to index dictionary

    labelsbydim = [Vector{String}() for _ in 0:lcdim]

    for k=0:lcdim
        append!(labelsbydim[k+1], sort(collect(labelsets[k+1])))
    end

    labelsvec = reduce(vcat,labelsbydim)
    labelindexdict = Dict{String,Int}(labelsvec[j] => j for j=1:length(labelsvec))

    # Create the vector of dimensions for the cells

    ldimvec = Vector{Int}()
    for k=0:lcdim
        for m=1:length(labelsbydim[k+1])
            push!(ldimvec,k)
        end
    end

    # Create the boundary map
    
    ncell  = length(labelsvec)

    Br = Vector{Int}()
    Bc = Vector{Int}()
    Bv = Vector{Int}()
    coeff = Int(1)

    for k in 1:defcelllen
        kbnd = defcellbnd[k]
        if length(kbnd) > 2
            cindex = labelindexdict[kbnd[1]]
            for m in 3:length(kbnd)
                rindex = labelindexdict[kbnd[m]]
                push!(Bc, cindex)   # Column index
                push!(Br, rindex)   # Row index
                push!(Bv, coeff)    # Matrix entry
            end
        end
    end

    B = sparse_from_lists(ncell, ncell, 2, Int(0), Int(1), Br, Bc, Bv)

    # Create and return the Lefschetz complex, validating boundary matrix

    lc = LefschetzComplex(labelsvec, ldimvec, B; validate=true)
    return lc
end

