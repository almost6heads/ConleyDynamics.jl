export create_cubical_complex, cube_field_size, cube_information, cube_label, get_cubical_coords

"""
    create_cubical_complex(cubes::Vector{String}; p::Int=2)

Initialize a Lefschetz complex from a cubical complex. The complex is
over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The vector `cubes` contains a list of all the highest-dimensional cubes
necessary to define the cubical complex. Every cube is represented as a
string as follows:

* d integers, which correspond to the coordinates of a point in d-dimensional Euclidean space
* a point `.`
* d integers 0 or 1, which give the interval length in the respective dimension

The first d integers all have to occupy the same number of characters. In addition,
if the occupied space is L characters for each coordinate, the coordinates only
can take values from 0 to 10^L - 2. This is due to the fact that the boundary
operator will add one to certain coordinates, and they still need to be 
representable withing the same L digits.

For example, the string `030600.101` corresponds to the point `(3,6,0)` in
three dimensions. The dimensions are 1, 0, and 1, and therefore this string
corresponds to the cube `[3,4] x [6] x [0,1]`. The same cube could have also
been represented by `360.101` or by `003006000.101`.

!!! warning
    Note that the labels all have to have the same format!

# Example
```jldoctest
julia> cubes = ["00.11", "01.01", "02.10", "11.10", "11.01", "22.00"];

julia> lc = create_cubical_complex(cubes);

julia> lc.ncells
17

julia> homology(lc)
3-element Vector{Int64}:
 2
 1
 0
```
"""
function create_cubical_complex(cubes::Vector{String}; p::Int=2)
    #
    # Create a Lefschetz complex struct for a cubical complex.
    #

    # Determine the point dimension and the coordinate field width

    @assert length(cubes)>=1 "We do need at least one cube!"
    pointdim, pointlen = cube_field_size(cubes[1])

    # Create a dictionary with all label to integer data information
    # 
    # If label is any label in cubes, then labelintinfodict[label]
    # is an integer vector with the following entries:
    #
    # 1:pointdim:              coordinates of the anchor point
    # 1+pointdim:2*pointdim:   interval length in each dimension
    # 1+2*pointdim:            dimension of the cube

    labelintinfodict = Dict{String,Vector{Int}}()
    for curlabel in cubes
        labelintinfodict[curlabel] = cube_information(curlabel)
    end

    # Find the dimension of the cubical complex
    
    CCdim = maximum([labelintinfodict[n][1+2*pointdim] for n in cubes])

    # Create labels for all cubes in the complex, organized by
    # dimensions, as well as ordered vectors of labels

    labelsets = [Set{String}() for _ in 0:CCdim]
    labelvect = [Vector{String}() for _ in 0:CCdim]
    
    for curcube in cubes   # Initialize with the given cubes
        curcdim = labelintinfodict[curcube][1+2*pointdim]
        push!(labelsets[curcdim+1], curcube)
    end

    # Work your way from highest to lowest dimension, keep adding faces

    if p==0
        tone  = 1//1
        tzero = 0//1
    else
        tone  = Int(1)
        tzero = Int(0)
    end
    bndchannel = Channel{Tuple{String,String,typeof(tone)}}(Inf)

    for k = CCdim:-1:1
        # Create vector of labels at dimension k

        labelvect[k+1] = sort(collect(labelsets[k+1]))

        # Loop through the cubes and add boundaries

        tmpncubes = length(labelvect[k+1])

        # Create labels and intinfo in parallel

        tmpchannel = Channel{Tuple{String,Vector{Int}}}(Inf)

        Threads.@threads for ck in 1:tmpncubes
            curcube    = labelvect[k+1][ck]
            curintinfo = labelintinfodict[curcube]
            curones = findall(x -> x==1, curintinfo[1+pointdim:2*pointdim])

            bndsign = tone
            for m=1:k
                newintinfoP = copy(curintinfo)
                newintinfoN = copy(curintinfo)
                ione = curones[m]
                newintinfoP[ione] += 1
                newintinfoP[ione+pointdim] = 0
                newintinfoN[ione+pointdim] = 0
                newintinfoP[1+2*pointdim] -= 1
                newintinfoN[1+2*pointdim] -= 1
                newlabelP = cube_label(pointdim, pointlen, newintinfoP)
                newlabelN = cube_label(pointdim, pointlen, newintinfoN)
                put!(tmpchannel, (newlabelP, newintinfoP))
                put!(tmpchannel, (newlabelN, newintinfoN))
                put!(bndchannel, (curcube, newlabelP,  bndsign))
                put!(bndchannel, (curcube, newlabelN, -bndsign))
                bndsign = -bndsign
            end
        end

        # Incorporate data in serial
        
        close(tmpchannel)
        tmppairs = collect(tmpchannel)
        for tmppair in tmppairs
            push!(labelsets[k], tmppair[1])
            labelintinfodict[tmppair[1]] = tmppair[2]
        end
    end

    labelvect[1] = sort(collect(labelsets[1]))

    # After these preparations, create the Lefschetz complex structure!!

    # Organize cube labels in one vector and create label to index dictionary

    CClabelvec = reduce(vcat,labelvect)
    CClabelindexdict = Dict{String,Int}(CClabelvec[j] => j for j=1:length(CClabelvec))
    CCn = length(CClabelvec)

    # Create the vector of dimensions for the cubes

    CCdimvec = Vector{Int}()
    for k=0:CCdim
        for m=1:length(labelvect[k+1])
            push!(CCdimvec,k)
        end
    end

    # Create the boundary map

    close(bndchannel)
    bndtuples = collect(bndchannel)

    Bc = [CClabelindexdict[bndtuples[k][1]] for k = 1:length(bndtuples)]
    Br = [CClabelindexdict[bndtuples[k][2]] for k = 1:length(bndtuples)]
    Bv = [bndtuples[k][3] for k = 1:length(bndtuples)]
    B  = sparse_from_lists(CCn,CCn,p,tzero,tone,Br,Bc,Bv)

    # Create the Lefschetz complex

    lc = LefschetzComplex(CClabelvec, CCdimvec, B; validate=false)

    # Return the Lefschetz complex

    return lc
end

"""
    cube_field_size(cube::String)

Determine the field sizes of a given cube label.

The function returns the dimension of the ambient space in the first output
parameter `pointdim`, and the length of the individual coordinate fields
in the second return variable `pointlen`.

# Example
```jldoctest
julia> cube_field_size("011654003020.0110")
(4, 3)
```
"""
function cube_field_size(cube::String)
    #
    # Determine the field sizes of a given cube label
    #

    pointdim = length(cube) - first(findfirst(".",cube))
    pointlen = Int((length(cube) - 1) / pointdim - 1)

    return pointdim, pointlen
end

"""
    cube_information(cube::String)

Compute a cube's coordinate information.

The function returns an integer vector with the cubes coordinate information.
The return vector `intinfo` contains in its components the following data:

* `1:pointdim`: Coordinates of the anchor point
* `1+pointdim:2*pointdim`: Interval length in each dimension
* `1+2*pointdim`: Dimension of the cube

Note that `pointdim` equals the dimension of the points specifying the cube.

# Example
```jldoctest
julia> cube_information("011654003.011")
7-element Vector{Int64}:
  11
 654
   3
   0
   1
   1
   2
```
"""
function cube_information(cube::String)
    #
    # Create a vector of a cube's integer data information. The
    # vector entries contain the following information:
    #
    # 1:pointdim:              coordinates of the anchor point
    # 1+pointdim:2*pointdim:   interval length in each dimension
    # 1+2*pointdim:            dimension of the cube
    #

    # Get the cube's size information

    pointdim, pointlen = cube_field_size(cube)

    # Extract the anchor point coordinates

    intinfo = Vector{Int}()
    for m = 1:pointdim
        i1 = 1 + (m-1) * pointlen
        i2 = m * pointlen
        push!(intinfo,parse(Int,cube[i1:i2]))
    end

    # Extract the cube's interval length information

    for m = 1:pointdim
        i1 = pointdim * pointlen + 1 + m
        push!(intinfo,parse(Int,cube[i1]))
    end

    # Add the dimension information for convenience

    cdim = sum(intinfo[pointdim+1:2*pointdim])
    push!(intinfo,cdim)

    # Check for admissibility

    if minimum(intinfo[1:pointdim]) < 0
        error("Point coordinates have to be at least 0!")
    end
    if maximum(intinfo[1:pointdim]) > 10^pointlen - 2
        error("Point coordinates overflow! (Cannot be 10^w-1)")
    end

    # Return the integer information

    return intinfo
end

"""
    cube_label(pointdim::Int, pointlen::Int, pointinfo::Vector{Int})

Create a label from a cube's coordinate information.

The dimension of the ambient Eucliden space is `pointdim`, while the field
length for each coordinate is `pointlen`. The vector `pointinfo` has to be of
length at least two times `pointdim`. The first `pointdim` entries contain the
coordinates of the anchor point, while the next `pointdim` entries are either 0
or 1 depending on the size of the interval. For example, if `poindim = 3`
and `pointinfo = [1,2,3,1,0,1]`, then we represent the cube in three-dimensional
space given by `[1,2] x [2] x [3 4]`.

# Example
```jldoctest
julia> cube_label(3,2,[10,23,5,1,1,0])
"102305.110"
```
"""
function cube_label(pointdim::Int, pointlen::Int, pointinfo::Vector{Int})
    #
    # Create the label for a cube
    #

    sp = string.(pointinfo[1:2*pointdim])
    label = ""
    for k = 1:pointdim
        label = label * repeat("0", pointlen-length(sp[k])) * sp[k]
    end
    label = label * "." * join(sp[1+pointdim:2*pointdim])
    
    return label
end

"""
    get_cubical_coords(cc::LefschetzComplex)

Compute the vertex coordinates for a cubical complex.

The variable `cc` has to contain a cubical complex, and the function
returns a vector of coordinates for the vertices of the complex, that
can then be used for plotting.
```
"""
function get_cubical_coords(cc::LefschetzComplex)
    #
    # Compute the vertex coordinates for a cubical complex
    #

    vindices = findall(x -> x==0, cc.dimensions)
    coords = Vector{Vector{Int}}()
    d, L = cube_field_size(cc.labels[1])

    for k in vindices
        cubeinfo  = cube_information(cc.labels[k])
        cubecoord = cubeinfo[1:d]
        push!(coords,cubecoord)
    end

    return coords
end

