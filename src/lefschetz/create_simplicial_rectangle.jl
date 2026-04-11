export create_simplicial_rectangle

"""
    create_simplicial_rectangle(nx::Int, ny::Int; p::Int=2, euclidean::Bool=false)

Create a simplicial complex covering a rectangle in the plane. The
complex is over the rationals if `p=0`, and over `GF(p)` if `p>0`.

The rectangle is given by the subset `[0,nx] x [0,ny]` of the plane.
Each unit square is represented by four triangles, which meet in
the center point of the square. Labels have the following meaning:

- The label `XXXYYYb` corresponds to the point `(XXX, YYY)`.
- The label `XXXYYYc` corresponds to `(XXX + 1/2, YYY + 1/2)`.

The number of characters in `XXX` and `YYY` matches the number
of digits of the larger number of `nx` and `ny`.

When `euclidean=false` (default), the function returns:
* A simplicial complex `sc::LefschetzComplex`.
* A vector `coords::Vector{Vector{Float64}}` of vertex coordinates.

When `euclidean=true`, it returns a `EuclideanComplex` with embedded
coordinates.
"""
function create_simplicial_rectangle(nx::Int, ny::Int; p::Int=2, euclidean::Bool=false)
    #
    # Create a Lefschetz complex struct for a simplicial rectangle.
    #

    # Make sure that we have at least width one in each direction

    if (nx < 1) || (ny < 1)
        error("Width and height have to be at least 1!")
    end

    # Compute the number of digits for the labels

    labwidth = Int(ceil(log(0.5+max(nx,ny))/log(10)))
    labformat = "%0" * string(labwidth) * "d%0" * string(labwidth) * "d%s"

    # Create the vector of labels

    labels = Vector{String}()
    coords = Vector{Vector{Float64}}()
    for ky = 0:ny
        for kx = 0:nx
            vlabel = Printf.format(Printf.Format(labformat), kx, ky, "b")
            push!(labels, vlabel)
            push!(coords, [Float64(kx), Float64(ky)])
        end
    end
    for ky = 0:ny-1
        for kx = 0:nx-1
            vlabel = Printf.format(Printf.Format(labformat), kx, ky, "c")
            push!(labels, vlabel)
            push!(coords, [Float64(kx)+0.5, Float64(ky)+0.5])
        end
    end

    # Create the vector of top-level triangles

    triangles = Vector{Vector{String}}()

    for ky = 0:ny-1
        for kx = 0:nx-1
            vlabel0 = Printf.format(Printf.Format(labformat), kx,   ky,   "b")
            vlabel1 = Printf.format(Printf.Format(labformat), kx+1, ky,   "b")
            vlabel2 = Printf.format(Printf.Format(labformat), kx+1, ky+1, "b")
            vlabel3 = Printf.format(Printf.Format(labformat), kx,   ky+1, "b")
            vlabelc = Printf.format(Printf.Format(labformat), kx, ky, "c")
            push!(triangles, [vlabel0,vlabel1,vlabelc])
            push!(triangles, [vlabel1,vlabel2,vlabelc])
            push!(triangles, [vlabel2,vlabel3,vlabelc])
            push!(triangles, [vlabel0,vlabel3,vlabelc])
        end
    end

    # Create the simplicial complex and return it

    if euclidean
        return create_simplicial_complex(labels, triangles, coords; p=p)
    else
        sc = create_simplicial_complex(labels, triangles, p=p)
        return sc, coords
    end
end

