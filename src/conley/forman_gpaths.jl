export forman_gpaths

"""
    forman_gpaths(lc::LefschetzComplex, fvf::CellSubsets, x::Cell;
                       extended::Bool=false)

Find all Forman gradient paths starting at a source cell.

For the Forman gradient vector field `fvf` on the Lefschetz complex
`lc` this function finds all maximal gradient paths starting at `x`.
These are solution paths which consist exclusively of Forman vectors,
i.e., they contain an even number of cells whose dimensions alternate
between `dim(x)` and `dim(x) + 1`. Every cell of dimension `dim(x)`
is an arrow source which is followed by its arrow target, while every
cell of dimension `dim(x) + 1` is an arrow target which is succeeded
by a cell in its boundary which is the source of a different arrow,
as long as such a cell exists.

If one passes the optional argument `extended = true`, then the
function determines all maximal gradient paths `p` such that the
subpath `p[1:end-1]` is a gradient path in the above sense,
and `p[end]` is an element in the boundary of `p[end-1]` which
is different from `p[end-2]`. Such paths clearly have an
odd length.
"""
function forman_gpaths(lc::LefschetzComplex, fvf::CellSubsets, x::Cell;
                       extended::Bool=false)
    #
    # Find all Forman gradient paths starting from x
    #

    # Perform a quick check that we have a Forman vector field

    if length(fvf)>0
        @assert maximum(length.(fvf))<=2 "Not a Forman vector field!"
    end

    # Convert x to index form

    if typeof(x)==String
        xi = lc.indices[x]
    else
        xi = x
    end
    xd = lc.dimensions[xi]

    # Convert fvf to index form

    if typeof(fvf)==Vector{Vector{String}}
        fvfi = convert_cellsubsets(lc, fvf)
    else
        fvfi = fvf
    end

    # Create a dictionary of Forman vectors whose source
    # has the same dimension as the cell x

    target = Dict{Int,Int}()
    for vector in fvfi
        if length(vector) == 2
            d1 = lc.dimensions[vector[1]]
            d2 = lc.dimensions[vector[2]]
            if d1 + 1 == d2
                if d1 == xd
                    target[vector[1]] = vector[2]
                end
            elseif d1 == d2 + 1
                if d2 == xd
                    target[vector[2]] = vector[1]
                end
            else
                error("This is not a Forman vector field!")
            end
        end
    end
    sources = collect(keys(target))

    # The cell x is not a source, therefore return the empty path

    if !(xi in sources)
        if typeof(x)==String
            return Vector{Vector{String}}([])
        else
            return Vector{Vector{Int}}([])
        end
    end

    # Find the paths

    paths = Vector{Vector{Int}}([])
    partialpaths = Vector{Vector{Int}}([])
    push!(partialpaths, [xi, target[xi]])

    while length(partialpaths) > 0
        cpath = popfirst!(partialpaths)
        lasttcell = cpath[length(cpath)]
        lastscell = cpath[length(cpath)-1]
        tbnd = setdiff(lefschetz_boundary(lc, lasttcell), [lastscell])
        nextscells = intersect(tbnd, sources)
        
        if length(nextscells) == 0
            push!(paths, cpath)
        else
            for k in nextscells
                if (k in cpath)
                    error("This Forman field is not gradient!")
                else
                    push!(partialpaths, vcat(cpath, [k, target[k]]))
                end
            end
        end
    end

    # If extended==false we are done, return the paths

    if !extended
        if typeof(x)==String
            return convert_cellsubsets(lc, paths)
        else
            return paths
        end
    end

    # If extended==true, we need to find the next elements for each path

    epaths = Vector{Vector{Int}}([])

    while length(paths) > 0
        cpath = popfirst!(paths)
        lasttcell = cpath[length(cpath)]
        lastscell = cpath[length(cpath)-1]
        tbnd = setdiff(lefschetz_boundary(lc, lasttcell), [lastscell])

        # All cells in tbnd are either targets or critical,
        # due to the construction of paths

        if length(tbnd)==0
            push!(epaths, cpath[1:length(cpath)-1])
        else
            for k in tbnd
                push!(epaths, vcat(cpath, [k]))
            end
        end
    end
    
    # Return the results

    if typeof(x)==String
        return convert_cellsubsets(lc, epaths)
    else
        return epaths
    end
end

"""
    forman_gpaths(lc::LefschetzComplex, fvf::CellSubsets,
                       x::Cell, y::Cell)

Find all Forman gradient paths between two cells.

For the Forman gradient vector field `fvf` on the Lefschetz
complex `lc` this function finds all gradient paths between
the cells `x` and `y`. The dimensions of these cells have to
satisfy one of the following three conditions:

* `dim(x) + 1 = dim(y)`: In this case the function
  returns all solution paths between `x` and `y` which
  consist entirely of Forman arrows. All the sources
  have the dimension of `x`, and the targets the
  dimension of `y`.
* `dim(x) = dim(y)`: In this case the function returns
  all solution paths `p` between `x` and `y` for which
  `p[1:end-1]` is a Forman gradient path in the above
  sense, and `p[end]` lies in the boundary of `p[end-1]`
  and is different from the cell `p[end-2]`.
* `dim(x) = dim(y) + 1`: In this case the function
  returns all solution paths `p` between `x` and `y`
  for which `p[2:end]` is a Forman gradient path in
  the sense of the second item, and `p[2]` lies
  in the boundary of `p[1]`.

In all other cases an error is raised.
"""
function forman_gpaths(lc::LefschetzComplex, fvf::CellSubsets,
                       x::Cell, y::Cell)
    #
    # Find all Forman gradient paths starting at x and ending in y
    #

    # Convert x and y to index form, find their dimensions

    if typeof(x)==String
        xi = lc.indices[x]
    else
        xi = x
    end
    xd = lc.dimensions[xi]

    if typeof(y)==String
        yi = lc.indices[y]
    else
        yi = y
    end
    yd = lc.dimensions[yi]

    # Make sure the dimensions are admissible

    @assert abs(xd - yd) <= 1 "Dimensions of x and y are not admissible!"

    # Now find all paths between x and y if dim(y) = dim(x) + 1
    # or if dim(y) = dim(x)

    if (yd == xd + 1) || (yd == xd)

        # Find all paths starting from x

        if yd == xd + 1
            allpaths = forman_gpaths(lc, fvf, xi)
        else
            allpaths = forman_gpaths(lc, fvf, xi, extended=true)
        end

        # Extract all paths that contain y

        paths = Vector{Vector{Int}}([])
        for cpath in allpaths
            yindices = findall(t -> t==yi, cpath)
            @assert length(yindices)<2 "This is not a gradient path!"
            if length(yindices)==1
                push!(paths, cpath[1:yindices[1]])
            end
        end
    end

    # Here we consider the remaining case dim(x) = dim(y) + 1

    if xd == yd + 1

        # Find the boundary of x

        xbnd = lefschetz_boundary(lc, xi)

        # For every boundary point, find all paths to y

        paths = Vector{Vector{Int}}([])
        for cbnd in xbnd
            bndpaths = forman_gpaths(lc, fvf, cbnd, yi)
            for cpath in bndpaths
                push!(paths, vcat([xi], cpath))
            end
        end
    end

    # Finally, return the results

    if typeof(x)==String
        return convert_cellsubsets(lc, unique(paths))
    else
        return unique(paths)
    end
end

