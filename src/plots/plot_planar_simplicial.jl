export plot_planar_simplicial

"""
    plot_planar_simplicial(sc::AbstractComplex,
                           coords::Vector{<:Vector{<:Real}},
                           fname::String;
                           [mvf::CellSubsets=Vector{Vector{Int}}([]),]
                           [labeldir::Vector{<:Real}=Vector{Int}([]),]
                           [labeldis::Real=8,]
                           [hfac::Real=1.2,]
                           [vfac::Real=1.2,]
                           [sfac::Real=0,]
                           [pdim::Vector{Bool}=[true,true,true],]
                           [pv::Bool=false])

Create an image of a planar simplicial complex, and if
specified, a Forman vector field on it.

The vector `coords` contains coordinates for every one of the
vertices of the simplicial complex `sc`. The image will be saved
in the file with name `fname`, and the ending determines the image
type. Accepted are `.pdf`, `.svg`, `.png`, and `.eps`.

If the optional `mvf` is specified and is a Forman vector field,
then this Forman vector field is drawn as well. The optional
vector `labeldir` contains directions for the vertex labels,
and `labeldis` the distance from the vertex. The directions
have to be reals between 0 and 4, with 0,1,2,3 corresponding
to E,N,W,S. The optional constants `hfac` and `vfac` contain
the horizontal and vertical scale vectors, while `sfac` describes
a uniform scale. If `sfac=0` the latter is automatically determined.
The vector `pdim` specifies in which dimensions cells are drawn; the
default shows vertices, edges, and triangles. Finally if one passes
the argument `pv=true`, then in addition to saving the file
a preview is displayed.

# Examples

Suppose we have created a simplicial complex using the commands

```julia
sc, coords = create_simplicial_delaunay(300, 300, 30, 20)
fname = "sc_plot_test.pdf"
```

Then the following code creates an image of the simplicial complex
without labels, but with a preview:

```julia
plot_planar_simplicial(sc, coords, fname, pv=true)
```

If we want to see the labels, we can use

```julia
ldir = fill(0.5, sc.ncells);
plot_planar_simplicial(sc, coords, fname, labeldir=ldir, labeldis=10, pv=true)
```

This command puts all labels in the North-East direction at a distance of 10.
"""
function plot_planar_simplicial(sc::AbstractComplex,
                                coords::Vector{<:Vector{<:Real}},
                                fname::String;
                                mvf::CellSubsets=Vector{Vector{Int}}([]),
                                labeldir::Vector{<:Real}=Vector{Int}([]),
                                labeldis::Real=8,
                                hfac::Real=1.2,
                                vfac::Real=1.2,
                                sfac::Real=0,
                                pdim::Vector{Bool}=[true,true,true],
                                pv::Bool=false)
    #
    # Create an image of a planar simplicial complex
    #

    # Check that the filename has a valid extension
    if !(lowercase(fname[end-3:end]) in [".png", ".pdf", ".eps", ".svg"])
        error("The filename must have one of the following extensions: .png, .pdf, .eps, .svg")
    end

    # Create proper coordinates
   
    cx0 = minimum([c[1] for c in coords])
    cx1 = maximum([c[1] for c in coords])
    cy0 = minimum([c[2] for c in coords])
    cy1 = maximum([c[2] for c in coords])

    if iszero(sfac)
        sfac = 800.0 / maximum([1, cx1-cx0, cy1-cy0])
    end

    figw  = Int(round((cx1 - cx0) * hfac * sfac))
    figh  = Int(round((cy1 - cy0) * vfac * sfac))
    figdx = (cx1 - cx0) * (hfac-1.0) * 0.5 * sfac
    figdy = (cy1 - cy0) * (vfac-1.0) * 0.5 * sfac

    pcoords = [[figdx + (c[1] - cx0) * (figw-2.0*figdx) / (cx1-cx0),
                figdy + (cy1 - c[2]) * (figh-2.0*figdy) / (cy1-cy0)]
               for c in coords]

    # Create a list of vertex points

    points = [Point(pcoords[k][1],pcoords[k][2]) for k in 1:length(pcoords)]

    # Create a list of vertices for each simplex

    cellvertices = Vector{Vector{Int}}()
    for k = 1:sc.ncells
        cv = lefschetz_skeleton(sc, [k], 0)
        if !(length(cv) == 1+sc.dimensions[k])
            error("The complex is not a simplicial complex!")
        end
        push!(cellvertices, cv)
    end

    # If necessary prepare for the Forman vector field plot

    if length(mvf) > 0

        # First make sure that the multivector field has type Int

        if typeof(mvf[1][1]) == String
            mvfI = convert_cellsubsets(sc, mvf)
        else
            mvfI = mvf
        end

        # Find the barycenters

        barycs = Vector{Point}()
        for k = 1:sc.ncells
            cdim = sc.dimensions[k]
            if cdim == 0
                push!(barycs, points[k])
            elseif cdim == 1
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                push!(barycs, midpoint(points[k1],points[k2]))
            elseif cdim == 2
                k1 = cellvertices[k][1]
                k2 = cellvertices[k][2]
                k3 = cellvertices[k][3]
                push!(barycs, trianglecenter(points[k1],points[k2],points[k3]))
            end
        end

        # Preprocess the Forman vector field

        invector = fill(false,sc.ncells)
        for k = 1:length(mvfI)
            if !(length(mvfI[k]) == 2)
                error("The multivector field is not a Forman vector field!")
            end
            invector[mvfI[k][1]] = true
            invector[mvfI[k][2]] = true
        end
        
        critical = Vector{Int}()
        for k = 1:sc.ncells
            if invector[k] == false
                push!(critical,k)
            end
        end
    end
    
    # Create the image
   
    Drawing(figw, figh, fname)
    background("white")
    sethue("black")
    
    # Plot the simplicial complex

    for k = sc.ncells:-1:1
        cdim = sc.dimensions[k]
        if cdim == 0
            setcolor("royalblue4")
            circle(points[k], 5, action = :fill)
        elseif cdim == 1
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            setcolor("royalblue3")
            line(points[k1],points[k2])
            strokepath()
        elseif cdim == 2
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            k3 = cellvertices[k][3]
            setcolor("steelblue1")
            poly([points[k1],points[k2],points[k3]],
                       action = :fill; close=true)
        end
    end

    if length(mvf) > 0

        # Plot the critical cells

        for m = 1:length(critical)
            k = critical[m]
            setcolor("red3")
            circle(barycs[k], 5, action = :fill)
        end

        # Plot the arrows

        for m = 1:length(mvf)
            k1 = mvfI[m][1]
            k2 = mvfI[m][2]
            cdim1 = sc.dimensions[k1]
            cdim2 = sc.dimensions[k2]
            setcolor("red3")
            if cdim1 < cdim2
                arrow(barycs[k1],barycs[k2],linewidth=3,arrowheadlength=10)
            elseif cdim1 > cdim2
                arrow(barycs[k2],barycs[k1],linewidth=3,arrowheadlength=10)
            else
                error("The multivector field is not a Forman vector field!")
            end
        end
    end

    # Display the vertex labels

    if length(labeldir) > 0
        setcolor("black")
        for k = 1:length(coords)
            label(sc.labels[k], pi*(4.0-labeldir[k])/2.0, points[k],
                  offset=labeldis)
        end
    end

    # Finish the drawing, and preview if desired

    finish()
    if pv
        preview()
    end
end

"""
    plot_planar_simplicial(ec::EuclideanComplex,
                           fname::String;
                           kwargs...)

Create an image of a planar simplicial complex from a `EuclideanComplex`.

This method extracts the vertex coordinates from the embedded coordinates
stored in `ec` and delegates to the standard plotting method. All keyword
arguments are forwarded unchanged.
"""
function plot_planar_simplicial(ec::EuclideanComplex,
                                fname::String;
                                mvf::CellSubsets=Vector{Vector{Int}}([]),
                                labeldir::Vector{<:Real}=Vector{Int}([]),
                                labeldis::Real=8,
                                hfac::Real=1.2,
                                vfac::Real=1.2,
                                sfac::Real=0,
                                pdim::Vector{Bool}=[true,true,true],
                                pv::Bool=false)
    #
    # Delegate to the LefschetzComplex method after extracting vertex coords
    #
    vertex_coords = [ec.coords[k][1] for k in 1:ec.ncells if ec.dimensions[k] == 0]
    lc = euclidean_to_lefschetz(ec)
    plot_planar_simplicial(lc, vertex_coords, fname;
                           mvf=mvf, labeldir=labeldir, labeldis=labeldis,
                           hfac=hfac, vfac=vfac, sfac=sfac, pdim=pdim, pv=pv)
end
