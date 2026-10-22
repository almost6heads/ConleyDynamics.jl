# Geometry helpers for the three-dimensional plotting functions of the
# Makie extension. They are part of the core package, so that they can
# be tested without Makie. All geometry is returned as plain vectors:
#
#   tris :: Vector{NTuple{3,Vector{Float64}}}   triangles
#   segs :: Vector{NTuple{2,Vector{Float64}}}   line segments
#   pts  :: Vector{Vector{Float64}}             points

const _Tris3D = Vector{NTuple{3,Vector{Float64}}}
const _Segs3D = Vector{NTuple{2,Vector{Float64}}}
const _Pts3D  = Vector{Vector{Float64}}

function _check_coords3d(ec::EuclideanComplex, fname::String)
    #
    # Check that the complex is embedded in three dimensions
    #
    if (ec.ncells == 0) || (length(ec.coords[1][1]) != 3)
        throw(ArgumentError("$fname: the complex needs to be embedded in " *
                            "three dimensions. For planar complexes, use the " *
                            "Plots.jl functions such as plot_simplicial."))
    end
end

function _check_pdim3d(pdim::AbstractVector, fname::String)
    #
    # Check that pdim has one entry for each cell dimension 0 to 3
    #
    if length(pdim) != 4
        throw(ArgumentError("$fname: pdim needs four entries, for vertices, " *
                            "edges, two-dimensional cells, and three-dimensional " *
                            "cells"))
    end
end

function _cells3d_indices(ec::EuclideanComplex, cells, fname::String)
    #
    # Convert a cell list given by labels or indices into indices
    #
    indices = Vector{Int}()
    for c in cells
        if c isa AbstractString
            k = get(ec.indices, c, 0)
            if k == 0
                throw(ArgumentError("$fname: the complex has no cell with label \"$c\""))
            end
        elseif c isa Integer
            if !(1 <= c <= ec.ncells)
                throw(ArgumentError("$fname: cell index $c outside of 1:$(ec.ncells)"))
            end
            k = Int(c)
        else
            throw(ArgumentError("$fname: cells need to be given by labels or indices"))
        end
        push!(indices, k)
    end
    return indices
end

function _cell_centroid(ec::EuclideanComplex, k::Int)
    #
    # Centroid of the vertices of a cell
    #
    cv = ec.coords[k]
    return sum(cv) ./ length(cv)
end

function _coface_lists(ec::EuclideanComplex)
    #
    # For every cell, the list of cells which contain it in their boundary
    #
    cof = [Int[] for _ in 1:ec.ncells]
    for j in 1:ec.ncells
        for i in ec.boundary.columns[j]
            push!(cof[i], j)
        end
    end
    return cof
end

function _clip_mask(ec::EuclideanComplex, clip, cof::Vector{Vector{Int}})
    #
    # Determine the cells which remain visible after clipping. For
    # clip = (n, offset), a cell without cofaces is removed if its
    # centroid c satisfies n⋅c > offset. Any other cell remains
    # visible if and only if one of its cofaces is visible, so that
    # the visible cells form a closed subcomplex.
    #

    keep = trues(ec.ncells)
    isnothing(clip) && return keep

    nvec, offset = clip
    nv = Float64.(collect(nvec))
    if length(nv) != 3
        throw(ArgumentError("The clip normal needs to have length 3"))
    end

    for k in sortperm(ec.dimensions; rev=true)
        if isempty(cof[k])
            keep[k] = sum(nv .* _cell_centroid(ec, k)) <= offset
        else
            keep[k] = any(keep[c] for c in cof[k])
        end
    end

    return keep
end

function _face_triangles(ec::EuclideanComplex, face::Int,
                         refpoint::Union{Vector{Float64},Nothing}=nothing)
    #
    # Triangulate a two-dimensional cell. Simplicial cells give one
    # triangle, and cubical cells two. The order of the four vertices
    # of a cubical cell is not known, so the vertex farthest from the
    # first vertex is used as the opposite corner. If refpoint is
    # given, the triangles are oriented so that their normals point
    # away from it.
    #

    cv = ec.coords[face]
    if length(cv) == 3
        tris = [(cv[1], cv[2], cv[3])]
    elseif length(cv) == 4
        dists = [sum(abs2, cv[m] .- cv[1]) for m in 2:4]
        opp   = 1 + argmax(dists)
        q, r  = [m for m in 2:4 if m != opp]
        tris  = [(cv[1], cv[q], cv[opp]), (cv[1], cv[opp], cv[r])]
    else
        throw(ArgumentError("Cannot triangulate a cell with $(length(cv)) vertices"))
    end

    if !isnothing(refpoint)
        tris = [_orient_triangle(t, refpoint) for t in tris]
    end

    return tris
end

function _orient_triangle(t::NTuple{3,Vector{Float64}}, refpoint::Vector{Float64})
    #
    # Orient a triangle so that its normal points away from refpoint
    #
    a, b, c = t
    nrm = _cross3(b .- a, c .- a)
    if sum(nrm .* (a .- refpoint)) < 0
        return (a, c, b)
    end
    return t
end

function _cells3d_geometry(ec::EuclideanComplex, cells::Vector{Int},
                           keep::BitVector; hull_only::Bool=true,
                           pdim::Vector{Bool}=[true,true,true,true])
    #
    # Geometry of a cell subset. Three-dimensional cells are shown by
    # their boundary surface (hull_only=true) or by all of their faces.
    # Cells of lower dimension are shown only if they are not faces of
    # higher-dimensional cells of the subset: two-dimensional cells as
    # triangles, edges as segments, and vertices as points. The work is
    # proportional to the size of the subset.
    #

    tris = _Tris3D()
    segs = _Segs3D()
    pts  = _Pts3D()

    allcells = Set{Int}(cells)
    cellset  = Set{Int}(k for k in allcells if keep[k])
    isempty(cellset) && return tris, segs, pts
    bnd = ec.boundary.columns

    # Boundary surface of the three-dimensional cells

    if pdim[4]
        facecount = Dict{Int,Int}()
        facecell  = Dict{Int,Int}()
        for t in cellset
            ec.dimensions[t] == 3 || continue
            for f in bnd[t]
                facecount[f] = get(facecount, f, 0) + 1
                facecell[f]  = t
            end
        end
        for f in sort(collect(keys(facecount)))
            if !hull_only || (facecount[f] == 1)
                append!(tris, _face_triangles(ec, f, _cell_centroid(ec, facecell[f])))
            end
        end
    end

    # Cells which are faces of higher-dimensional cells in the subset,
    # including the ones removed by clipping

    covered = Set{Int}()
    stack = [b for k in allcells for b in bnd[k]]
    while !isempty(stack)
        c = pop!(stack)
        if !(c in covered)
            push!(covered, c)
            append!(stack, bnd[c])
        end
    end

    # Remaining cells of lower dimension

    for k in sort(collect(cellset))
        k in covered && continue
        d = ec.dimensions[k]
        if (d == 2) && pdim[3]
            append!(tris, _face_triangles(ec, k))
        elseif (d == 1) && pdim[2]
            push!(segs, (ec.coords[k][1], ec.coords[k][2]))
        elseif (d == 0) && pdim[1]
            push!(pts, ec.coords[k][1])
        end
    end

    return tris, segs, pts
end

function _background3d_geometry(ec::EuclideanComplex, keep::BitVector,
                                cof::Vector{Vector{Int}}; hull_only::Bool=true,
                                pdim::Vector{Bool}=[true,true,true,true])
    #
    # Geometry of the whole visible complex. For hull_only=true, only
    # the boundary surface of the visible three-dimensional cells is
    # shown, together with its edges and vertices, as well as all
    # visible lower-dimensional cells which are not faces of
    # three-dimensional cells. The vector pdim determines which
    # vertices, edges, free faces, and surfaces of three-dimensional
    # cells are included.
    #

    tris = _Tris3D()
    segs = _Segs3D()
    pts  = _Pts3D()
    bnd  = ec.boundary.columns

    # Faces

    facemask = falses(ec.ncells)
    for f in 1:ec.ncells
        (ec.dimensions[f] == 2) && keep[f] || continue
        cof3 = [c for c in cof[f] if keep[c] && ec.dimensions[c] == 3]
        if isempty(cof3)
            any(ec.dimensions[c] == 3 for c in cof[f]) && continue
            facemask[f] = true
            pdim[3] && append!(tris, _face_triangles(ec, f))
        elseif !hull_only || (length(cof3) == 1)
            facemask[f] = true
            pdim[4] && append!(tris, _face_triangles(ec, f, _cell_centroid(ec, cof3[1])))
        end
    end

    # Edges and vertices

    edgemask = falses(ec.ncells)
    for e in 1:ec.ncells
        (ec.dimensions[e] == 1) && keep[e] || continue
        if !hull_only || any(facemask[c] for c in cof[e]) ||
           !any(ec.dimensions[c] == 2 for c in cof[e])
            edgemask[e] = true
            pdim[2] && push!(segs, (ec.coords[e][1], ec.coords[e][2]))
        end
    end

    for v in 1:ec.ncells
        (ec.dimensions[v] == 0) && keep[v] || continue
        if pdim[1] && (!hull_only || any(edgemask[c] for c in cof[v]) || isempty(cof[v]))
            push!(pts, ec.coords[v][1])
        end
    end

    return tris, segs, pts
end

function _length_scale(ec::EuclideanComplex)
    #
    # Typical edge length of the complex, used for default marker sizes
    #
    lens = Float64[]
    for k in 1:ec.ncells
        if ec.dimensions[k] == 1
            push!(lens, sqrt(sum(abs2, ec.coords[k][1] .- ec.coords[k][2])))
            length(lens) >= 2000 && break
        end
    end
    return isempty(lens) ? 1.0 : sort(lens)[(length(lens)+1) ÷ 2]
end

function _ci_color3d_key(ci::Vector{Int})
    #
    # Pad or truncate a Conley index to four entries for the color lookup
    #
    return [k <= length(ci) ? ci[k] : 0 for k in 1:4]
end
