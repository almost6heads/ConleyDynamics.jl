export create_cubical_box

"""
    create_cubical_box(nx::Int, ny::Int, nz::Int;
                       p::Int=2, randomize::Real=0.0, euclidean::Bool=false)

Create a cubical complex covering a box in space. The complex is
over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The box is given by the subset `[0,nx] x [0,ny] x [0,nz]` of space,
and each unit cube gives a three-dimensional cube in the resulting
cubical complex.

When `euclidean=false` (default), the function returns:
* A cubical complex `cc::LefschetzComplex`
* A vector `coords::Vector{Vector{Float64}}` of vertex coordinates

When `euclidean=true`, it returns a `EuclideanComplex` with embedded
coordinates.

If the optional parameter `randomize` is assigned a positive real
fraction `r` less that 0.5, then the actual coordinates will be
randomized. They are chosen uniformly from balls of radius `r`
centered at each vertex.
"""
function create_cubical_box(nx::Int, ny::Int, nz::Int;
                            p::Int=2, randomize::Real=0.0, euclidean::Bool=false)
    #
    # Create a Lefschetz complex struct for a cubical box
    #

    # Make sure that we have at least width one in each direction

    if (nx < 1) || (ny < 1) || (nz < 1)
        error("Width and height have to be at least 1!")
    end

    # Create the vector of three-dimensional cubes

    pointdim = 3
    pointlen = Int(ceil(log(maximum([nx,ny,nz]) + 2) / log(10)))
    cubes = Vector{String}()

    for kx = 0:nx-1
        for ky = 0:ny-1
            for kz = 0:nz-1
                pointinfo = [kx, ky, kz, 1, 1, 1]
                clabel = cube_label(pointdim, pointlen, pointinfo)
                push!(cubes, clabel)
            end
        end
    end

    # Create the Lefschetz complex

    lc, _ = create_cubical_complex(cubes, p=p)

    # Create the coordinate vector

    coords = Vector{Vector{Float64}}()
    for k = 1:lc.ncells
        if lc.dimensions[k] == 0
            vertexk = cube_information(lc.labels[k])
            push!(coords, [Float64(vertexk[1]), Float64(vertexk[2]), Float64(vertexk[3])])
        end
    end

    # Randomize the points if desired

    if randomize >= 0.5
        error("The randomization radius has to be less than 0.5!")
    end

    for k in eachindex(coords)
        ox = randn()
        oy = randn()
        oz = randn()
        or = randomize * (rand())^(1.0/3.0) / sqrt(ox*ox+oy*oy+oz*oz)
        coords[k][1] = coords[k][1] + or * ox;
        coords[k][2] = coords[k][2] + or * oy;
        coords[k][3] = coords[k][3] + or * oz;
    end

    # Return the results

    if euclidean
        return lefschetz_to_euclidean(lc, coords)
    else
        return lc, coords
    end
end

