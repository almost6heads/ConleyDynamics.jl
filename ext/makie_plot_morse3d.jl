# Plotting of Morse sets and single multivectors in three dimensions

function _plot_subsets3d!(ax, ec::EuclideanComplex, subsets::Vector{Vector{Int}},
                          colors, fname::String;
                          alpha::Real=0.9, bgalpha::Real=0.08,
                          pdim::Vector{Bool}=[false,true,true,true],
                          hull_only::Bool=true, clip=nothing,
                          linewidth::Real=4, markersize=nothing)
    #
    # Draw a faint background complex and colored cell subsets
    #

    CD._check_coords3d(ec, fname)
    cof   = CD._coface_lists(ec)
    keep  = CD._clip_mask(ec, clip, cof)
    msize = isnothing(markersize) ? 0.12 * CD._length_scale(ec) : markersize

    # Background

    tris, segs, pts = CD._background3d_geometry(ec, keep, cof;
                                                hull_only=hull_only, pdim=pdim)
    _draw_geometry!(ax, tris, CD._Segs3D(), CD._Pts3D();
                    color=colorant"steelblue1", alpha=bgalpha)
    _draw_geometry!(ax, CD._Tris3D(), segs, CD._Pts3D(); color=colorant"royalblue3",
                    linealpha=0.25, linewidth=0.5)
    _draw_geometry!(ax, CD._Tris3D(), CD._Segs3D(), pts; color=colorant"royalblue4",
                    markersize=0.4*msize)

    # Cell subsets, combined into one plot call per cell dimension

    tris, tcols = CD._Tris3D(), Any[]
    segs, scols = CD._Segs3D(), Any[]
    pts,  pcols = CD._Pts3D(),  Any[]
    for (s, col) in zip(subsets, colors)
        t, l, p = CD._cells3d_geometry(ec, s, keep; hull_only=hull_only)
        append!(tris, t); append!(tcols, fill(col, length(t)))
        append!(segs, l); append!(scols, fill(col, length(l)))
        append!(pts, p);  append!(pcols, fill(col, length(p)))
    end
    _draw_colored_geometry!(ax, tris, tcols, segs, scols, pts, pcols; alpha=alpha,
                            linewidth=linewidth, markersize=msize)

    return ax
end

function _morse_colors(ec::EuclideanComplex, subsets::Vector{Vector{Int}},
                       ci::Bool, colors)
    #
    # Colors of Morse sets: by Conley index, given, or default
    #
    if ci
        return [_ci_color3d(conley_index(ec, s)) for s in subsets]
    elseif isnothing(colors)
        return _default_colors(length(subsets))
    end
    return [colors[mod1(k, length(colors))] for k in 1:length(subsets)]
end

for (fname, fname!) in ((:plot_simplicial3d_morse, :plot_simplicial3d_morse!),
                        (:plot_cubical3d_morse, :plot_cubical3d_morse!))
    @eval begin
        function ConleyDynamics.$fname!(ax, ec::EuclideanComplex, morsesets;
                                        ci::Bool=false, addcritical::Bool=false,
                                        colors=nothing, kwargs...)
            ms = _cellsubset_indices(ec, morsesets)
            addcritical && (ms = _add_singletons(ec, ms))
            cols = _morse_colors(ec, ms, ci, colors)
            return _plot_subsets3d!(ax, ec, ms, cols, $(string(fname!)); kwargs...)
        end
        function ConleyDynamics.$fname(ec::EuclideanComplex, morsesets;
                                       axis::Symbol=:scene, size=(800,600), kwargs...)
            fig, ax = _new_figure(axis, size)
            ConleyDynamics.$fname!(ax, ec, morsesets; kwargs...)
            return fig
        end
    end
end

for (fname, fname!) in ((:plot_simplicial3d_mv, :plot_simplicial3d_mv!),
                        (:plot_cubical3d_mv, :plot_cubical3d_mv!))
    @eval begin
        function ConleyDynamics.$fname!(ax, ec::EuclideanComplex, mv;
                                        color=colorant"darkorange", kwargs...)
            return _plot_subsets3d!(ax, ec, [_cell_indices(ec, mv)], [color],
                                    $(string(fname!)); kwargs...)
        end
        function ConleyDynamics.$fname(ec::EuclideanComplex, mv;
                                       axis::Symbol=:scene, size=(800,600), kwargs...)
            fig, ax = _new_figure(axis, size)
            ConleyDynamics.$fname!(ax, ec, mv; kwargs...)
            return fig
        end
    end
end
