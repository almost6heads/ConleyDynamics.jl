export create_simplicial_complex

"""
    create_simplicial_complex(labels::Vector{String},
                              simplices::Vector{Vector{Int}};
                              p::Int=2)

Initialize a Lefschetz complex from a simplicial complex. The
complex is over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The vector `labels` contains a label for every vertex, while
`simplices` contains all the highest-dimensional simplices necessary
to define the simplicial complex. Every simplex is represented as
a vector of `Int`, with entries corresponding to the vertex indices.

!!! warning
    Note that the labels all have to have the same character length!
"""
function create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{Int}};
                                   p::Int=2)
    #
    # Create a Lefschetz complex struct for a simplicial complex.
    #

    # Check whether the labels all have the same length

    sclabels = deepcopy(labels)
    scsimplices = deepcopy(simplices)

    if length(unique(length.(sclabels))) > 1
        error("All labels need to have the same length!")
    end

    # Sort the simplex vertices in ascending order
    
    for k=1:length(scsimplices)
        sort!(scsimplices[k])
    end

    # Find the dimension of the simplicial complex

    scdim = 0
    for k = 1:length(scsimplices)
        scdim = max(scdim, length(scsimplices[k])-1)
    end

    # Create labels for all simplices of dimension at least 1
    # in the simplicial complex and organize them in sets by
    # dimension. In addition, create a dictionary which relates
    # the face labels to the index sequences.
    
    labelsets = [Set{String}() for _ in 1:scdim]
    labelvertexdict = Dict{String,Vector{Int}}()

    for k=1:length(scsimplices)
        csimp = scsimplices[k]
        cdim = length(csimp) - 1
        for m=2:cdim+1
            for faceindices in combinations(csimp,m)
                sortedfaceindices = sort(faceindices)
                facelab = join(sort(sclabels[sortedfaceindices]))
                push!(labelsets[m-1],facelab)
                labelvertexdict[facelab] = sortedfaceindices
            end
        end
    end
    
    # Order the simplex labels and create label to index dictionary

    labelsbydim = [Vector{String}() for _ in 0:scdim]

    append!(labelsbydim[1],sclabels)
    for k=1:scdim
        append!(labelsbydim[k+1],sort(collect(labelsets[k])))
    end

    labelsvec = reduce(vcat,labelsbydim)
    labelindexdict = Dict{String,Int}(labelsvec[j] => j for j=1:length(labelsvec))

    # Create the vector of dimensions for the simplices

    sdimvec = Vector{Int}()
    for k=0:scdim
        for m=1:length(labelsbydim[k+1])
            push!(sdimvec,k)
        end
    end

    # Create the boundary map
    
    nsimp  = length(labelsvec)
    nsimp0 = length(sclabels)

    Br = Vector{Int}()
    Bc = Vector{Int}()
    Bv = Vector{Int}()
    for k=nsimp0+1:nsimp
        csimp = labelvertexdict[labelsvec[k]]
        coeff = Int(1)
        for m=1:length(csimp)
            csimptmp = deepcopy(csimp)
            deleteat!(csimptmp,m)
            bsimp = labelindexdict[join(sort(sclabels[csimptmp]))]
            push!(Br,bsimp)   # Row index
            push!(Bc,k)       # Column index
            push!(Bv,coeff)   # Matrix entry
            coeff = -coeff
        end
    end

    if p > 0
        B = sparse_from_lists(nsimp,nsimp,p,Int(0),Int(1),Br,Bc,Bv)
    else
        tzero = Rational{Int}(0)
        tone  = Rational{Int}(1)
        Bvrational = convert(Vector{Rational{Int}},Bv)
        B = sparse_from_lists(nsimp,nsimp,0,tzero,tone,Br,Bc,Bvrational)
    end
    
    # Create the Lefschetz complex

    lc = LefschetzComplex(labelsvec, sdimvec, B; validate=false)

    # Return the Lefschetz complex

    return lc
end

"""
    create_simplicial_complex(labels::Vector{String},
                              simplices::Vector{Vector{String}};
                              p::Int=2)

Initialize a Lefschetz complex from a simplicial complex. The
complex is over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The vector `labels` contains a label for every vertex, while
`simplices` contains all the highest-dimensional simplices necessary
to define the simplicial complex.
"""
function create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{String}};
                                   p::Int=2)
    #
    # Create a Lefschetz complex struct for a simplicial complex.
    #
    # This is an alternative method for simplices represented
    # in String format.
    #
    newsimplices = convert_simplices(simplices, labels)
    lc = create_simplicial_complex(labels, newsimplices, p=p)
    return lc
end

"""
    create_simplicial_complex(labels::Vector{String},
                              simplices::Vector{Vector{Int}},
                              vertex_coords::Vector{<:Vector{<:Real}};
                              p::Int=2)

Initialize an `EuclideanComplex` from a simplicial complex with embedded
vertex coordinates. The complex is over the rationals if `p=0`, and over
`GF(p)` if `p>0`.

The vector `labels` contains a label for every vertex, `simplices` contains
all the highest-dimensional simplices, and `vertex_coords` contains the
Euclidean coordinates of each vertex (indexed by vertex index, i.e.,
`vertex_coords[k]` gives the coordinate of vertex `k`).
"""
function create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{Int}},
                                   vertex_coords::Vector{<:Vector{<:Real}};
                                   p::Int=2)
    #
    # Create an EuclideanComplex for a simplicial complex with coordinates.
    #
    lc = create_simplicial_complex(labels, simplices; p=p)
    return lefschetz_to_euclidean(lc, vertex_coords)
end

"""
    create_simplicial_complex(labels::Vector{String},
                              simplices::Vector{Vector{String}},
                              vertex_coords::Vector{<:Vector{<:Real}};
                              p::Int=2)

Initialize an `EuclideanComplex` from a simplicial complex with embedded
vertex coordinates (String-simplex variant).
"""
function create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{String}},
                                   vertex_coords::Vector{<:Vector{<:Real}};
                                   p::Int=2)
    #
    # Create an EuclideanComplex for a simplicial complex with coordinates.
    # This is an alternative method for simplices in String format.
    #
    newsimplices = convert_simplices(simplices, labels)
    return create_simplicial_complex(labels, newsimplices, vertex_coords; p=p)
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function convert_simplices(simplices::Vector{Vector{Int}},
                           labels::Vector{String})
    #
    # Convert list of simplices from index form to label form
    #

    newsimplices = Vector{Vector{String}}()

    for k=1:length(simplices)
        push!(newsimplices,Vector{String}())
        for m=1:length(simplices[k])
            push!(newsimplices[k],labels[simplices[k][m]])
        end
    end

    return newsimplices
end

###############################################################################

function convert_simplices(simplices::Vector{Vector{String}},
                           labels::Vector{String})
    #
    # Convert list of simplices from label form to index form
    #

    indices = Dict{String,Int}(labels[k] => k for k in 1:length(labels))
    newsimplices = Vector{Vector{Int}}()

    for k=1:length(simplices)
        push!(newsimplices,Vector{Int}())
        for m=1:length(simplices[k])
            push!(newsimplices[k],indices[simplices[k][m]])
        end
    end

    return newsimplices
end

