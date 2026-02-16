export convert_planar_coordinates, convert_spatial_coordinates

"""
    convert_planar_coordinates(coords::Vector{<:Vector{<:Real}},
                               p0::Vector{<:Real},
                               p1::Vector{<:Real})

Convert a given collection of planar coordinates.

The vector `coords` contains pairs of coordinates, which are then
transformed to fit into the box with vertices `p0 = (p0x,p0y)`
and `p1 = (p1x,p1y)`. It is assumed that `p0` denotes the lower
left box corner, while `p1` is the upper right corner. The function
shifts and scales the coordinates in such a way that every side of
the box contains at least one point. Upon completion, it returns
a new coordinate vector `coordsNew`.

More precisely, if the x-coordinates are spanning the interval
`[xmin,xmax]` and the y-coordinates span `[ymin,ymax]`, then the
point `(x,y)` is transformed to `(xn,yn)` with:

* `xn = p0x + (p1x-p0x) * (x-cxmin) / (cxmax-cxmin)`
* `yn = p0y + (p1y-p0y) * (y-cymin) / (cymax-cymin)`
"""
function convert_planar_coordinates(coords::Vector{<:Vector{<:Real}},
                                    p0::Vector{<:Real},
                                    p1::Vector{<:Real})
    #
    # Convert planar coordinates to fit inside a given box
    #

    # Find the minimal and maximal coordinates

    cxmin = minimum([c[1] for c in coords])
    cxmax = maximum([c[1] for c in coords])
    cymin = minimum([c[2] for c in coords])
    cymax = maximum([c[2] for c in coords])

    # Create and return the new coordinate vector

    newcoords = Vector{Vector{Float64}}()

    for c in coords
        cnx = p0[1] + (p1[1] - p0[1]) * (c[1] - cxmin) / (cxmax - cxmin)
        cny = p0[2] + (p1[2] - p0[2]) * (c[2] - cymin) / (cymax - cymin)
        push!(newcoords, [cnx,cny])
    end

    return newcoords
end

"""
    convert_spatial_coordinates(coords::Vector{<:Vector{<:Real}},
                                p0::Vector{<:Real},
                                p1::Vector{<:Real})

Convert a given collection of spatial coordinates.

The vector `coords` contains triples of coordinates, which are then
transformed to fit into the box with vertices `p0 = (p0x,p0y,p0z)`
and `p1 = (p1x,p1y,p1z)`. It is assumed that each coordinate of `p0`
is strictly smaller than the corresponding coordinate of `p1`. The
function shifts and scales the coordinates in such a way that every
side of the box contains at least one point. Upon completion, it
returns a new coordinate vector `coordsNew`.

More precisely, if the x-coordinates are spanning the interval
`[xmin,xmax]`, the y-coordinates span `[ymin,ymax]`, and the
z-coordinates span `[zmin,zmax]`, then the point `(x,y,z)` is
transformed to `(xn,yn,zn)` with:

* `xn = p0x + (p1x-p0x) * (x-cxmin) / (cxmax-cxmin)`
* `yn = p0y + (p1y-p0y) * (y-cymin) / (cymax-cymin)`
* `zn = p0z + (p1z-p0z) * (z-czmin) / (czmax-czmin)`
"""
function convert_spatial_coordinates(coords::Vector{<:Vector{<:Real}},
                                     p0::Vector{<:Real},
                                     p1::Vector{<:Real})
    #
    # Convert spatial coordinates to fit inside a given box
    #

    # Find the minimal and maximal coordinates

    cxmin = minimum([c[1] for c in coords])
    cxmax = maximum([c[1] for c in coords])
    cymin = minimum([c[2] for c in coords])
    cymax = maximum([c[2] for c in coords])
    czmin = minimum([c[3] for c in coords])
    czmax = maximum([c[3] for c in coords])

    # Create and return the new coordinate vector

    newcoords = Vector{Vector{Float64}}()

    for c in coords
        cnx = p0[1] + (p1[1] - p0[1]) * (c[1] - cxmin) / (cxmax - cxmin)
        cny = p0[2] + (p1[2] - p0[2]) * (c[2] - cymin) / (cymax - cymin)
        cnz = p0[3] + (p1[3] - p0[3]) * (c[3] - czmin) / (czmax - czmin)
        push!(newcoords, [cnx,cny,cnz])
    end

    return newcoords
end

