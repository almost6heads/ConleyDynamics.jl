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

    pointdim, pointlen = cube_field_size(cubes[1])

    # Create a dictionary with all label to integer data information
    # 
    # If label is any label in cubes, then labelintinfodict[label]
    # is an integer vector with the following entries:
    #
    # 1:pointdim:              coordinates of the anchor point
    # 1+pointdim:2*pointdim:   interval length in each dimension
    # 1+2*pointdim:            dimension of the cube

    @debug "Creating labe dictionary..."
    labelintinfodict = Dict{String,Vector{Int}}()
    for curlabel in cubes
        labelintinfodict[curlabel] = cube_information(curlabel)
    end

    # Find the dimension of the cubical complex
    
    CCdim = maximum([labelintinfodict[n][1+2*pointdim] for n in cubes])

    # Create labels for all cubes in the complex, organized by
    # dimensions, as well as ordered vectors of labels

    @debug "Creating labels..."
    labelsets = [Set{String}() for _ in 0:CCdim]
    labelvect = [Vector{String}() for _ in 0:CCdim]
    
    for curcube in cubes   # Initialize with the given cubes
        curcdim = labelintinfodict[curcube][1+2*pointdim]
        push!(labelsets[curcdim+1], curcube)
    end

    # Work your way from highest to lowest dimension, keep adding faces

    for k = CCdim:-1:1
        @debug "Working on dimension..." k
        # Create vector of labels at dimension k

        labelvect[k+1] = sort(collect(labelsets[k+1]))

        # Loop through the cubes and add boundaries

        tmpncubes   = length(labelvect[k+1])
        tmpnewlabelP = [Vector{String}() for _ in 1:tmpncubes]
        tmpnewlabelN = [Vector{String}() for _ in 1:tmpncubes]
        emptyvecint  = [Vector{Int}() for _ in 1:k]
        tmpnewinfoP  = [deepcopy(emptyvecint) for _ in 1:tmpncubes]
        tmpnewinfoN  = [deepcopy(emptyvecint) for _ in 1:tmpncubes]

        # Create labels and intinfo in parallel

        Threads.@threads for ck in 1:tmpncubes

            curcube    = labelvect[k+1][ck]
            curintinfo = labelintinfodict[curcube]
            curones = findall(x -> x==1, curintinfo[1+pointdim:2*pointdim])

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
                push!(tmpnewlabelP[ck], newlabelP)
                push!(tmpnewlabelN[ck], newlabelN)
                tmpnewinfoP[ck][m] = newintinfoP
                tmpnewinfoN[ck][m] = newintinfoN
            end
        end

        # Incorporate data in serial

        for ck in 1:tmpncubes
            for m in 1:k
                push!(labelsets[k], tmpnewlabelP[ck][m])
                push!(labelsets[k], tmpnewlabelN[ck][m])
                labelintinfodict[tmpnewlabelP[ck][m]] = tmpnewinfoP[ck][m]
                labelintinfodict[tmpnewlabelN[ck][m]] = tmpnewinfoN[ck][m]
            end
        end
    end

    labelvect[1] = sort(collect(labelsets[1]))

    # After these preparations, create the Lefschetz complex structure!!

    # Organize cube labels in one vector and create label to index dictionary

    @debug "Organizing cube labels..."
    CClabelvec = reduce(vcat,labelvect)
    CClabelindexdict = Dict{String,Int}(CClabelvec[j] => j for j=1:length(CClabelvec))
    CCn = length(CClabelvec)

    # Create the vector of dimensions for the cubes

    @debug "creating dimension vector..."
    CCdimvec = Vector{Int}()
    for k=0:CCdim
        for m=1:length(labelvect[k+1])
            push!(CCdimvec,k)
        end
    end

    # Create the boundary map

    @debug "Creating boundary map..."
    pdindices = findall(x -> x>0, CCdimvec)
    pdsize    = length(pdindices)
    tmpBr = [Vector{Int}() for _ in 1:pdsize]
    tmpBc = [Vector{Int}() for _ in 1:pdsize]
    tmpBv = [Vector{Int}() for _ in 1:pdsize]

    Threads.@threads for ki = 1:pdsize
        k = pdindices[ki]
        curcdim = CCdimvec[k]
        curcube = CClabelvec[k]
        curintinfo = labelintinfodict[curcube]
        curones = findall(x -> x==1, curintinfo[1+pointdim:2*pointdim])

        for m=1:curcdim
            if mod(m,2) == 1
                coeff = Int(1)
            else
                coeff = Int(-1)
            end
            newintinfoP = deepcopy(curintinfo)
            newintinfoN = deepcopy(curintinfo)
            ione = curones[m]
            newintinfoP[ione] += 1
            newintinfoP[ione+pointdim] = 0
            newintinfoN[ione+pointdim] = 0
            newintinfoP[1+2*pointdim] -= 1
            newintinfoN[1+2*pointdim] -= 1
            newlabelP = cube_label(pointdim, pointlen, newintinfoP)
            newlabelN = cube_label(pointdim, pointlen, newintinfoN)
            bcellP = CClabelindexdict[newlabelP]
            bcellN = CClabelindexdict[newlabelN]

            push!(tmpBr[ki],bcellP)   # Row index
            push!(tmpBc[ki],k)        # Column index
            push!(tmpBv[ki],coeff)    # Matrix entry
                
            push!(tmpBr[ki],bcellN)   # Row index
            push!(tmpBc[ki],k)        # Column index
            push!(tmpBv[ki],-coeff)   # Matrix entry
        end
    end

    @debug "Creating sparse matrix..."
    Br = reduce(vcat,tmpBr)
    Bc = reduce(vcat,tmpBc)
    Bv = reduce(vcat,tmpBv)

    if p > 0
        B = sparse_from_lists(CCn,CCn,p,Int(0),Int(1),Br,Bc,Bv)
    else
        tzero = Rational{Int}(0)
        tone  = Rational{Int}(1)
        Bvrational = convert(Vector{Rational{Int}},Bv)
        B = sparse_from_lists(CCn,CCn,0,tzero,tone,Br,Bc,Bvrational)
    end

    # Create the Lefschetz complex

    @debug "Creating Lefschetz complex..."
    lc = LefschetzComplex(CClabelvec, CCdimvec, B)

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

    sp = string.(pointinfo)
    label = ""
    for k = 1:pointdim
        label = label * repeat("0", pointlen-length(sp[k])) * sp[k]
    end
    label = label * "." * join(sp[1+pointdim:end])
    
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

