export plot_planar_cubical

"""
    plot_planar_cubical(cc::LefschetzComplex,
                        fname::String;
                        [hfac::Real=1.2,]
                        [vfac::Real=1.2,]
                        [cubefac::Real=0,]
                        [pdim::Vector{Bool}=[true,true,true],]
                        [pv::Bool=false])

Create an image of a planar cubical complex.

The image will be saved in the file with name `fname`, and the
ending determines the image type. Accepted are `.pdf`, `.svg`,
`.png`, and `.eps`. The optional constants `hfac` and `vfac` contain
the horizontal and vertical scale vectors. The optional argument
`cubefac` specifies the side length of an elementary cube for
plotting, and it will be automatically determined otherwise. The
vector `pdim` specifies which cell dimensions should be plotted,
with `pdim[k]` representing dimension `k-1`. Finally if one passes
the argument `pv=true`, then in addition to saving the file
a preview is displayed.

# Examples

Suppose we have created a cubical complex using the commands

```julia
cubes = ["00.11", "01.01", "02.10", "11.10", "11.01", "22.00"]
cc = create_cubical_complex(cubes)
fname = "cc_plot_test.pdf"
```

Then the following code creates an image of the simplicial complex
without labels, but with a preview:

```julia
plot_planar_cubical(cc, fname, pv=true)
```

If one only wants to plot the edges in the complex, but not the
vertices or rectangles, then one can use:

```julia
plot_planar_cubical(cc, fname, pv=true, pdim=[false,true,false])
```
"""
function plot_planar_cubical(cc::LefschetzComplex,
                             fname::String;
                             hfac::Real=1.2,
                             vfac::Real=1.2,
                             cubefac::Real=0,
                             pdim::Vector{Bool}=[true,true,true],
                             pv::Bool=false)
    #
    # Create an image of a planar cubical complex
    #

    # Extract the vertex information

    vertices = lefschetz_skeleton(cc, 0)

    if !(length(vertices) == maximum(vertices))
        error("The vertices are not at the beginning of the cell list!")
    end

    # If desired, determine cubefac automatically

    if iszero(cubefac)
        oxmin = 10^10
        oxmax = 0
        oymin = 10^10
        oymax = 0
        for k = 1:length(vertices)
            intinfo = cube_information(cc.labels[k])
            oxmin = minimum([oxmin, intinfo[1]])
            oxmax = maximum([oxmax, intinfo[1]])
            oymin = minimum([oymin, intinfo[2]])
            oymax = maximum([oymax, intinfo[2]])
        end
        cubefac = 800.0 / maximum([1, oxmax-oxmin, oymax-oymin])
    end

    # Compute the coordinates

    coords = Vector{Vector{Float64}}()
    for k = 1:length(vertices)
        intinfo = cube_information(cc.labels[k])
        push!(coords, [cubefac*intinfo[1], cubefac*intinfo[2]])
    end

    # Create proper coordinates
   
    cx0 = minimum([c[1] for c in coords])
    cx1 = maximum([c[1] for c in coords])
    cy0 = minimum([c[2] for c in coords])
    cy1 = maximum([c[2] for c in coords])
    figw = Int(round((cx1 - cx0) * hfac))
    figh = Int(round((cy1 - cy0) * vfac))
    figdx = (cx1 - cx0) * (hfac-1.0) * 0.5
    figdy = (cy1 - cy0) * (vfac-1.0) * 0.5

    pcoords = [[figdx + (c[1] - cx0) * (figw-2.0*figdx) / (cx1-cx0),
                figdy + (cy1 - c[2]) * (figh-2.0*figdy) / (cy1-cy0)]
               for c in coords]

    # Create a list of vertex points

    points = [Point(pcoords[k][1],pcoords[k][2]) for k in 1:length(pcoords)]

    # Create a list of vertices for each cube

    cellvertices = Vector{Vector{Int}}()
    for k = 1:cc.ncells
        cv = lefschetz_skeleton(cc, [k], 0)
        if !(length(cv) == 2^cc.dimensions[k])
            error("The complex is not a cubical complex!")
        end
        push!(cellvertices, cv)
    end

    # Create the image
   
    Drawing(figw, figh, fname)
    background("white")
    sethue("black")
    
    # Plot the cubical complex

    for k = cc.ncells:-1:1
        cdim = cc.dimensions[k]
        if (cdim == 0) & pdim[1]
            setcolor("royalblue4")
            circle(points[k], 5, action = :fill)
        elseif (cdim == 1) & pdim[2]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            setcolor("royalblue3")
            line(points[k1],points[k2])
            strokepath()
        elseif (cdim == 2) & pdim[3]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][3]
            k3 = cellvertices[k][4]
            k4 = cellvertices[k][2]
            setcolor("steelblue1")
            poly([points[k1],points[k2],points[k3],points[k4]],
                       action = :fill; close=true)
        end
    end

    # Finish the drawing, and preview if desired

    finish()
    if pv
        preview()
    end
end
