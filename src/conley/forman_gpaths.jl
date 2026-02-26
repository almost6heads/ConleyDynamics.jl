export forman_gpaths
export forman_path_weight

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

    # Find the paths

    paths = Vector{Vector{Int}}([])
    partialpaths = Vector{Vector{Int}}([])

    if xi in sources
        push!(partialpaths, [xi, target[xi]])
    end

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
    push!(epaths, [xi])   # We have to add this trivial path

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
        return convert_cellsubsets(lc, unique(epaths))
    else
        return unique(epaths)
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

In all other cases an empty collection is returned.
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

    # Return the empty collection if |xd - yd| > 1

    if abs(xd - yd) > 1
        if typeof(x)==String
            return Vector{Vector{String}}([])
        else
            return Vector{Vector{String}}([])
        end
    end

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

"""
    forman_path_weight(lc::LefschetzComplex, path::Cells)

Compute the weight of a Forman gradient path.

For the Lefschetz complex `lc` this function computes the 
weight of the Forman gradient path given in `path`.
It is expected that the dimensions of the first and the
last cell in the path differ by at most 1. In case they
have equal dimension, and the path has length larger
than 1, the first cell has to be an arrow source.
"""
function forman_path_weight(lc::LefschetzComplex, path::Cells)
    #
    # Compute the weight of a Forman gradient path
    #

    # Deal with trivial paths

    if length(path) < 2
        return lc.boundary.zero
    end

    # Convert path to index form

    if typeof(path)==Vector{String}
        pathi = convert_cells(lc, path)
    else
        pathi = path
    end
    xd = lc.dimensions[pathi[1]]
    yd = lc.dimensions[pathi[end]]
    p  = lc.boundary.char

    # First consider the case dim(path[1]) = dim(path[end])-1

    if xd == yd - 1
        # Make sure the path length is even

        pathlen = length(pathi)
        @assert mod(pathlen, 2)==0 "The path has to have an even number of cells!"
    
        # Determine the path weight and return it
    
        n      = div(pathlen, 2) - 1     
        kappa1 = (-1)^(n+1)
        kappa2 = lc.boundary[pathi[2*n+1], pathi[2*n+2]]
        weight = scalar_multiply(kappa1, scalar_inverse(kappa2, p), p)

        for k=1:n
            kappa1 = lc.boundary[pathi[2*k+1], pathi[2*k]]
            kappa2 = lc.boundary[pathi[2*k-1], pathi[2*k]]
            weight = scalar_multiply(weight, kappa1, p)
            weight = scalar_multiply(weight, scalar_inverse(kappa2, p), p)
        end
        return weight
    end

    # Now consider the case dim(path[1]) = dim(path[end])

    if xd == yd
        @assert lc.dimensions[pathi[1]]==lc.dimensions[pathi[end-1]]-1 "Dimension mismatch!"
        weight = forman_path_weight(lc, pathi[1:end-1])
        kappa  = lc.boundary[pathi[end], pathi[end-1]]
        weight = scalar_multiply(weight, kappa, p)
        return weight
    end

    # Finally we deal with the case dim(path[1]) = dim(path[end])+1

    if xd == yd + 1
        @assert lc.dimensions[pathi[2]]==lc.dimensions[pathi[end]] "Dimension mismatch!"
        if length(pathi) == 2
            weight = lc.boundary.one
        else
            weight = forman_path_weight(lc, pathi[2:end])
        end
        kappa  = lc.boundary[pathi[2], pathi[1]]
        weight = scalar_multiply(weight, kappa, p)
        return weight
    end
end

"""
    forman_path_weight(lc::LefschetzComplex, paths::CellSubsets)

Accumulated weight of a collection of Forman gradient paths.

For the Lefschetz complex `lc` this function computes the 
sum of the weights of the collection of Forman gradient
paths given in `paths`.
"""
function forman_path_weight(lc::LefschetzComplex, paths::CellSubsets)
    #
    # Compute the accumulated weight of Forman gradient paths
    #
    if length(paths) == 0
        return lc.boundary.zero
    end
    weights = map(x -> forman_path_weight(lc,x), paths)
    addfct = (x,y) -> scalar_add(x, y, lc.boundary.char)
    return reduce(addfct, weights)
end

