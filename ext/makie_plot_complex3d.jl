# Plotting of complexes and cell subsets in three dimensions

function _plot_complex3d!(ax, ec::EuclideanComplex, fname::String;
                          pdim::AbstractVector{Bool}=[false,true,true,true],
                          hull_only::Bool=true, clip=nothing, alpha::Real=0.3,
                          facecolor=colorant"steelblue1",
                          edgecolor=colorant"royalblue3",
                          vertexcolor=colorant"royalblue4",
                          linewidth::Real=0.75, linealpha::Real=1.0,
                          markersize=nothing)
    #
    # Draw the visible part of a complex
    #
    CD._check_coords3d(ec, fname)
    CD._check_pdim3d(pdim, fname)
    cof  = CD._coface_lists(ec)
    keep = CD._clip_mask(ec, clip, cof)
    tris, segs, pts = CD._background3d_geometry(ec, keep, cof;
                                                hull_only=hull_only, pdim=pdim)
    msize = isnothing(markersize) ? 0.06 * CD._length_scale(ec) : markersize
    _draw_geometry!(ax, tris, CD._Segs3D(), CD._Pts3D(); color=facecolor, alpha=alpha)
    _draw_geometry!(ax, CD._Tris3D(), segs, CD._Pts3D(); color=edgecolor,
                    linealpha=linealpha, linewidth=linewidth)
    _draw_geometry!(ax, CD._Tris3D(), CD._Segs3D(), pts; color=vertexcolor,
                    markersize=msize)
    return ax
end

function _plot_cells3d!(ax, ec::EuclideanComplex, cells, keep::BitVector;
                        color=colorant"darkorange", alpha::Real=1.0,
                        pdim::AbstractVector{Bool}=[true,true,true,true],
                        hull_only::Bool=true, linewidth::Real=4,
                        markersize::Real=0.05)
    #
    # Draw a cell subset with a precomputed clipping mask
    #
    CD._check_pdim3d(pdim, "plot_cells3d!")
    indices = _cell_indices(ec, cells, "plot_cells3d!")
    tris, segs, pts = CD._cells3d_geometry(ec, indices, keep;
                                           hull_only=hull_only, pdim=pdim)
    _draw_geometry!(ax, tris, segs, pts; color=color, alpha=alpha,
                    linewidth=linewidth, markersize=markersize)
    return ax
end

function ConleyDynamics.plot_cells3d!(ax, ec::EuclideanComplex, cells;
                                      clip=nothing, markersize=nothing, kwargs...)
    CD._check_coords3d(ec, "plot_cells3d!")
    keep  = CD._clip_mask(ec, clip, CD._coface_lists(ec))
    msize = isnothing(markersize) ? 0.12 * CD._length_scale(ec) : markersize
    return _plot_cells3d!(ax, ec, cells, keep; markersize=msize, kwargs...)
end

for (fname, fname!) in ((:plot_simplicial3d, :plot_simplicial3d!),
                        (:plot_cubical3d, :plot_cubical3d!))
    @eval begin
        function ConleyDynamics.$fname!(ax, ec::EuclideanComplex; kwargs...)
            return _plot_complex3d!(ax, ec, $(string(fname!)); kwargs...)
        end
        function ConleyDynamics.$fname(ec::EuclideanComplex; axis::Symbol=:scene,
                                       size=(800,600), kwargs...)
            fig, ax = _new_figure(axis, size)
            _plot_complex3d!(ax, ec, $(string(fname)); kwargs...)
            return fig
        end
    end
end
