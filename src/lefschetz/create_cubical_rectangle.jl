export create_cubical_rectangle

"""
    create_cubical_rectangle(nx::Int, ny::Int;
                             p::Int=2, randomize::Real=0.0)

Create a cubical complex covering a rectangle in the plane. The
complex is over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The rectangle is given by the subset `[0,nx] x [0,ny]` of the plane,
and each unit square gives a two-dimensional cube in the resulting
cubical complex. The function returns the following objects:

* A cubical complex `cc::LefschetzComplex`
* A vector `coords::Vector{Vector{Float64}}` of vertex coordinates

If the optional parameter `randomize` is assigned a positive real
fraction `r` less that 0.5, then the actual coordinates will be
randomized. They are chosen uniformly from discs of radius `r`
centered at each vertex.
"""
function create_cubical_rectangle(nx::Int, ny::Int;
                                  p::Int=2, randomize::Real=0.0)
    #
    # Create a Lefschetz complex struct for a cubical rectangle.
    #

    # Make sure that we have at least width one in each direction

    @assert min(nx,ny)>=1 "Width and height have to be at least 1!"

    # Determine the field size, and initialize the coordinate labels

    pointdim = 2
    pointlen = Int(ceil(log(max(nx,ny) + 2) / log(10)))

    cl = Vector{String}()
    for k = 0:max(nx,ny)
        push!(cl, repeat("0", pointlen-length(string(k))) * string(k))
    end

    # Create the vertex labels

    cubes  = Vector{String}()
    for k = 0:nx
        for m = 0:ny
            push!(cubes, cl[k+1] * cl[m+1] * ".00")
        end
    end

    labels = sort(cubes)
    dims   = fill(0, length(labels))
    l2i    = Dict(labels[k] => k for k = 1:length(labels))

    # Set up the lists for the boundary matrix creation
    
    if p==0
        tone  = 1//1
        tzero = 0//1
    else
        tone  = Int(1)
        tzero = Int(0)
    end

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(tone)}()

    # Create the edge labels and boundary matrix entries

    cubes   = Vector{String}()
    bndlist = Vector{Tuple{String,Int,typeof(tone)}}()

    for k = 0:nx-1
        for m = 0:ny
            ccube = cl[k+1] * cl[m+1] * ".10"
            push!(cubes, ccube)
            push!(bndlist, (ccube, l2i[cl[k+2] * cl[m+1] * ".00"], tone))
            push!(bndlist, (ccube, l2i[cl[k+1] * cl[m+1] * ".00"], -tone))
        end
    end
    for k = 0:nx
        for m = 0:ny-1
            ccube = cl[k+1] * cl[m+1] * ".01"
            push!(cubes, ccube)
            push!(bndlist, (ccube, l2i[cl[k+1] * cl[m+2] * ".00"], tone))
            push!(bndlist, (ccube, l2i[cl[k+1] * cl[m+1] * ".00"], -tone))
        end
    end

    sort!(cubes)
    ncells = length(labels)
    for k in 1:length(cubes)
        l2i[cubes[k]] = k + ncells
    end

    labels = [labels; cubes]
    dims   = [dims; fill(1, length(cubes))]

    for k in 1:length(bndlist)
        push!(c, l2i[bndlist[k][1]])
        push!(r, bndlist[k][2])
        push!(v, bndlist[k][3])
    end

    # Create the rectangle labels and boundary matrix entries

    cubes   = Vector{String}()
    bndlist = Vector{Tuple{String,Int,typeof(tone)}}()

    for k = 0:nx-1
        for m = 0:ny-1
            ccube = cl[k+1] * cl[m+1] * ".11"
            push!(cubes, ccube)
            push!(bndlist, (ccube, l2i[cl[k+2] * cl[m+1] * ".01"], tone))
            push!(bndlist, (ccube, l2i[cl[k+1] * cl[m+1] * ".01"], -tone))
            push!(bndlist, (ccube, l2i[cl[k+1] * cl[m+2] * ".10"], -tone))
            push!(bndlist, (ccube, l2i[cl[k+1] * cl[m+1] * ".10"], tone))
        end
    end

    sort!(cubes)
    ncells = length(labels)
    for k in 1:length(cubes)
        l2i[cubes[k]] = k + ncells
    end

    labels = [labels; cubes]
    dims   = [dims; fill(2, length(cubes))]
    ncells = length(labels)

    for k in 1:length(bndlist)
        push!(c, l2i[bndlist[k][1]])
        push!(r, bndlist[k][2])
        push!(v, bndlist[k][3])
    end

    # Create the boundary matrix and Lefschetz complex

    B  = sparse_from_lists(ncells,ncells,p,tzero,tone,r,c,v)
    lc = LefschetzComplex(labels, dims, B)

    # Create the coordinate vector

    coords = Vector{Vector{Float64}}()
    for k = 1:lc.ncells
        if lc.dimensions[k] == 0
            vertexk = cube_information(lc.labels[k])
            push!(coords, [Float64(vertexk[1]), Float64(vertexk[2])])
        end
    end

    # Randomize the points if desired

    if randomize >= 0.5
        error("The randomization radius has to be less than 0.5!")
    end

    for k in eachindex(coords)
        rrad = sqrt(rand()) * randomize;
        rtheta = rand() * 2.0 * pi
        coords[k][1] = coords[k][1] + rrad * cos(rtheta);
        coords[k][2] = coords[k][2] + rrad * sin(rtheta);
    end

    # Return the results

    return lc, coords
end

