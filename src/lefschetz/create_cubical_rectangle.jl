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

    if (nx < 1) || (ny < 1)
        error("Width and height have to be at least 1!")
    end

    # Create the vector of two-dimensional cubes

    @debug "Creating vector of cubes..."
    pointdim = 2
    pointlen = Int(ceil(log(maximum([nx,ny]) + 2) / log(10)))
    cubes = Vector{String}()

    for k = 0:nx-1
        for m = 0:ny-1
            pointinfo = [k, m, 1, 1]
            clabel = cube_label(pointdim, pointlen, pointinfo)
            push!(cubes, clabel)
        end
    end

    # Create the Lefschetz complex

    @debug "Creating cubical complex..."
    lc = create_cubical_complex(cubes, p=p)

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

