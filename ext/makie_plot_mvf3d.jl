# Plotting of multivector fields in three dimensions

function _greedy_coloring(ec::EuclideanComplex, mvs::Vector{Vector{Int}}, ncolors::Int)
    #
    # Assign color indices to multivectors such that multivectors with
    # adjacent cells get different colors whenever possible. The
    # candidate colors are tried in a rotating order, so that all
    # colors are used about equally often.
    #
    owner = zeros(Int, ec.ncells)
    for (m, mv) in enumerate(mvs)
        owner[mv] .= m
    end
    cof = CD._coface_lists(ec)
    colidx = zeros(Int, length(mvs))
    for (m, mv) in enumerate(mvs)
        used = falses(ncolors)
        for c in mv, nb in Iterators.flatten((ec.boundary.columns[c], cof[c]))
            o = owner[nb]
            if (o > 0) && (o != m) && (colidx[o] > 0)
                used[colidx[o]] = true
            end
        end
        start = mod1(m, ncolors)
        colidx[m] = start
        for k in 0:ncolors-1
            cand = mod1(start + k, ncolors)
            if !used[cand]
                colidx[m] = cand
                break
            end
        end
    end
    return colidx
end

for (fname, fname!) in ((:plot_simplicial3d_mvf, :plot_simplicial3d_mvf!),
                        (:plot_cubical3d_mvf, :plot_cubical3d_mvf!))
    @eval begin
        function ConleyDynamics.$fname!(ax, ec::EuclideanComplex, mvf;
                                        colors=nothing, ci::Bool=false,
                                        addcritical::Bool=false, alpha::Real=0.8,
                                        kwargs...)
            mvs = filter(mv -> length(mv) > 1, _cellsubset_indices(ec, mvf, $(string(fname!))))
            if addcritical
                mvs = _add_singletons(ec, _cellsubset_indices(ec, mvf, $(string(fname!))))
            end
            if ci
                cols = [_ci_color3d(conley_index(ec, mv)) for mv in mvs]
            else
                palette = isnothing(colors) ? _default_colors(12; light=true) : colors
                cols = palette[_greedy_coloring(ec, mvs, length(palette))]
            end
            return _plot_subsets3d!(ax, ec, mvs, cols, $(string(fname!));
                                    alpha=alpha, kwargs...)
        end
        function ConleyDynamics.$fname(ec::EuclideanComplex, mvf;
                                       axis::Symbol=:scene, size=(800,600), kwargs...)
            fig, ax = _new_figure(axis, size)
            ConleyDynamics.$fname!(ax, ec, mvf; kwargs...)
            return fig
        end
    end
end
