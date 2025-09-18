export plot_planar_simplicial_morse

"""
    plot_planar_simplicial_morse(sc::LefschetzComplex,
                                 coords::Vector{<:Vector{<:Real}},
                                 fname::String,
                                 morsesets::CellSubsets;
                                 [hfac::Real=1.2,]
                                 [vfac::Real=1.2,]
                                 [sfac::Real=0,]
                                 [pdim::Vector{Bool}=[false,true,true],]
                                 [pv::Bool=false,]
                                 [ci::Bool=false])

Create an image of a planar simplicial complex, together with
Morse sets, or also selected multivectors.

The vector `coords` contains coordinates for every one of the
vertices of the simplicial complex `sc`. The image will be saved
in the file with name `fname`, and the ending determines the image
type. Accepted are `.pdf`, `.svg`, `.png`, and `.eps`.

The vector `morsesets` contains a list of Morse sets, or more general,
subsets of the simplicial complex. For every `k`, the set described
by `morsesets[k]` will be shown in a distinct color.

The optional constants `hfac` and `vfac` contain the horizontal and
vertical scale vectors for the margins, while `sfac` describes a uniform
scale. If `sfac=0` the latter is automatically determined. The vector `pdim`
specifies in which dimensions cells are drawn; the default only shows edges
and triangles. If one passes the argument `pv=true`, then in addition
to saving the file a preview is displayed. Lastly, passing the argument
`ci=true` will color code the Morse sets according to their respective
Conley indices.
"""
function plot_planar_simplicial_morse(sc::LefschetzComplex,
                                       coords::Vector{<:Vector{<:Real}},
                                       fname::String,
                                       morsesets::CellSubsets;
                                       hfac::Real=1.2,
                                       vfac::Real=1.2,
                                       sfac::Real=0,
                                       pdim::Vector{Bool}=[false,true,true],
                                       pv::Bool=false,
                                       ci::Bool=false)
    #
    # Create an image of a planar simplicial complex
    #

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

    # Create the image
    
    Drawing(figw, figh, fname)
    background("white")
    sethue("black")
    
    # Plot the simplicial complex

    for k = sc.ncells:-1:1
        cdim = sc.dimensions[k]
        if (cdim == 0) && pdim[1]
            setcolor("royalblue4")
            setopacity(0.4)
            circle(points[k], 5, action = :fill)
        elseif (cdim == 1) && pdim[2]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            setcolor("royalblue3")
            setopacity(0.4)
            line(points[k1],points[k2])
            strokepath()
        elseif (cdim == 2) && pdim[3]
            k1 = cellvertices[k][1]
            k2 = cellvertices[k][2]
            k3 = cellvertices[k][3]
            setcolor("steelblue1")
            setopacity(0.4)
            poly([points[k1],points[k2],points[k3]],
                 action = :fill; close=true)
        end
    end

    # Plot the Morse sets

    if morsesets isa Vector{Vector{Int}}
        msI = morsesets
    else
        msI = convert_cellsubsets(sc, morsesets)
    end
    
    # Color code Morse sets if 'ci' is flagged
    if ci == true
        col1 = colorant"rgb(51,34,136)"     #color-blind safe 'indigo'
        col2 = colorant"rgb(221,204,119)"   #color-blind safe 'tan' 
        col3 = colorant"rgb(17,119,51)"     #color-blind safe 'moss green'
        col4 = colorant"rgb(204,102,119)"  #color-blind safe 'salmon'
        col5 = colorant"rgb(135,34,85)"     #color-blind safe 'maroon'
        col6 = colorant"rgb(153,153,51)"    #color-blind safe 'asparagus'
        c1 = [1,0,0]  #asym stable
        c2 = [0,1,0]  #1d unstable
        c3 = [0,0,1]  ### unstable
        c4 = [1,1,0]  ##  stable periodic
        c5 = [0,1,1]  # unstable periodic
        
        # Pass the single Morse set `msI[m]`
        for m in eachindex(msI)
            current_ci = conley_index(sc, msI[m])
            
            if (current_ci == c1)
                setcolor(col1)
            elseif (current_ci == c2)     
                setcolor(col2)
            elseif (current_ci == c3)     
                setcolor(col3)
            elseif (current_ci == c4)     
                setcolor(col4)
            elseif (current_ci == c5)     
                setcolor(col5)        
            else 
                setcolor(col6)     
            end
            
            setopacity(0.6)
            for k in msI[m]
                cdim = sc.dimensions[k]
                if (cdim == 0) && pdim[1]
                    circle(points[k], 5, action = :fill)
                elseif (cdim == 1) && pdim[2]
                    k1 = cellvertices[k][1]
                    k2 = cellvertices[k][2]
                    line(points[k1],points[k2])
                    strokepath()
                elseif (cdim == 2) && pdim[3]
                    k1 = cellvertices[k][1]
                    k2 = cellvertices[k][2]
                    k3 = cellvertices[k][3]
                    poly([points[k1],points[k2],points[k3]],
                         action = :fill; close=true)
                end
            end
        end


    else # if ci == false
        col1 = colorant"royalblue4"
        col2 = colorant"royalblue3"
        col3 = colorant"steelblue1"
        cols = distinguishable_colors(length(msI), [col1,col2,col3], dropseed=true)
        
        for m in eachindex(msI)
            setcolor(cols[m])
            setopacity(0.6)
            for k in msI[m]
                cdim = sc.dimensions[k]
                if (cdim == 0) && pdim[1]
                    circle(points[k], 5, action = :fill)
                elseif (cdim == 1) && pdim[2]
                    k1 = cellvertices[k][1]
                    k2 = cellvertices[k][2]
                    line(points[k1],points[k2])
                    strokepath()
                elseif (cdim == 2) && pdim[3]
                    k1 = cellvertices[k][1]
                    k2 = cellvertices[k][2]
                    k3 = cellvertices[k][3]
                    poly([points[k1],points[k2],points[k3]],
                         action = :fill; close=true)
                end
            end
        end
    end

    # Finish the drawing, and preview if desired
    finish()
    if pv
        preview()
    end
end
