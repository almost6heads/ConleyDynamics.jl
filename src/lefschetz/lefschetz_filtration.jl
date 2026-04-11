export lefschetz_filtration

"""
    lefschetz_filtration(lc::AbstractComplex, fvalues::Vector{Int})

Compute a filtration on a Lefschetz subset.

The considered Lefschetz complex is given in `lc`. The vector `fvalues`
assigns an integer between 0 and N to every cell in `lc`. For every k
the complex `L_k` is given by the closure of all cells with values
between 1 and k. The function returns the following variables:
* `lcsub`: The subcomplex `L_N`
* `fvalsub`: The filtration on the subcomplex with values 1,...,N

# Example
```jldoctest
julia> labels = ["A","B","C","D","E","F","G"];

julia> simplices = [["A","B","D"],["B","D","E"],["B","C","E"],["C","E","F"],["F","G"]];

julia> sc = create_simplicial_complex(labels,simplices);

julia> filtration = [0,0,0,0,0,0,0,1,1,0,1,2,0,4,2,4,0,5,3,7,6];

julia> lcsub, fvalsub = lefschetz_filtration(sc,filtration);

julia> phinf, phint = persistent_homology(lcsub,fvalsub);

julia> phinf
3-element Vector{Vector{Int64}}:
 [1]
 []
 []

julia> phint
3-element Vector{Vector{Tuple{Int64, Int64}}}:
 []
 [(1, 5), (2, 7), (4, 6)]
 []
```
"""
function lefschetz_filtration(lc::AbstractComplex, fvalues::Vector{Int})
    #
    # Compute a Lefschetz complex filtration
    #

    # Basic consistency check

    maxf = maximum(fvalues)
    fval = sort(unique(fvalues))

    if !((fval == Vector{Int}(0:maxf)) | (fval == Vector{Int}(1:maxf)))
        error("The filtration values have to be consecutive integers starting at 1!")
    end

    # Create an empty filtration array

    newfvalues = fill(Int(0),lc.ncells)

    # Loop through the nonzero filtration values from highest to lowest
    # and fill in the closures

    for k = maxf:-1:1
        indicesk = findall(x -> x==k, fvalues)
        clk = lefschetz_closure(lc, indicesk)
        newfvalues[clk] = fill(k, length(clk))
    end

    # Return the filtration

    if minimum(newfvalues) == 1
        # Every cell appears in at least one filtration complex
        lcsub   = deepcopy(lc)
        fvalsub = deepcopy(newfvalues)
    else
        # Not all cells are used, extract a subcomplex
        indexpos = findall(x -> x>0, newfvalues)
        lcsub = lefschetz_closed_subcomplex(lc, indexpos)
        fvalsub = filter(x -> x>0, newfvalues)
    end
        
    return lcsub, fvalsub
end

"""
    lefschetz_filtration(lc::AbstractComplex, strfilt::Vector{Vector{String}})

Compute a filtration on a Lefschetz subset.

The considered Lefschetz complex is given in `lc`. The vector of string
vectors `strfilt` contains the necessary simplices to build the filtration.
The list `strfilt[k]` contains the simplices that are added at the k-th
step, together with their closures. Thus, for every k the complex `L_k`
is given by the closure of all cells listed in `strfilt[i]` for `i`
between 1 and k. The function returns the following variables:
* `lcsub`: The subcomplex `L_N`, where `N = length(strfilt)`
* `fvalsub`: The filtration on the subcomplex with values 1,...,N

# Example
```jldoctest
julia> labels = ["A","B","C","D","E","F","G"];

julia> simplices = [["A","B","D"],["B","D","E"],["B","C","E"],["C","E","F"],["F","G"]];

julia> sc = create_simplicial_complex(labels,simplices);

julia> strfiltration = [["AB","AD","BD"],["BE","DE"],["BCE"],["CF","EF"],["ABD"],["CEF"],["BDE"]];

julia> lcsub, fvalsub = lefschetz_filtration(sc, strfiltration);

julia> phinf, phint = persistent_homology(lcsub,fvalsub);

julia> phinf
3-element Vector{Vector{Int64}}:
 [1]
 []
 []

julia> phint
3-element Vector{Vector{Tuple{Int64, Int64}}}:
 []
 [(1, 5), (2, 7), (4, 6)]
 []
```
"""
function lefschetz_filtration(lc::AbstractComplex, strfilt::Vector{Vector{String}})
    #
    # Compute a Lefschetz complex filtration
    #

    # Create an empty filtration array

    fvalues = fill(Int(0),lc.ncells)

    # Loop through the string filtration.
    # Start at the last level and go down, to account for
    # simplices that are added more than once. That way they
    # appear first when they should.

    for k = length(strfilt):-1:1
        for label in strfilt[k]
            labelindex = lc.indices[label]
            fvalues[labelindex] = k
        end
    end

    # Invoke the integer filtration list method and return the results

    lcsub, fvalsub = lefschetz_filtration(lc, fvalues)
    return lcsub, fvalsub
end

