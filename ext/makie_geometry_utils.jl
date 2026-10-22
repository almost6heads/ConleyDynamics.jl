# Conversion of the core geometry helpers into Makie plots

const CD = ConleyDynamics

function _new_figure(axis::Symbol, size)
    #
    # Create a figure with an LScene or an Axis3
    #
    fig = Makie.Figure(size=size)
    if axis == :scene
        ax = Makie.LScene(fig[1,1]; show_axis=false)
    elseif axis == :axis3
        ax = Makie.Axis3(fig[1,1]; aspect=:data)
    else
        throw(ArgumentError("The keyword axis needs to be :scene or :axis3"))
    end
    return fig, ax
end

function _draw_geometry!(ax, tris, segs, pts; color, alpha::Real=1.0,
                         linecolor=color, linealpha::Real=1.0,
                         linewidth::Real=1.0, markersize::Real=0.05)
    #
    # Draw triangles, segments, and points with one plot call each
    #

    if !isempty(tris)
        P = [Makie.Point3f(v...) for t in tris for v in t]
        F = [GB.TriangleFace(3k-2, 3k-1, 3k) for k in 1:length(tris)]
        Makie.mesh!(ax, GB.normal_mesh(P, F); color=(color, alpha),
                    transparency=(alpha < 1), backlight=1.0f0)
    end

    if !isempty(segs)
        S = [Makie.Point3f(v...) for s in segs for v in s]
        Makie.linesegments!(ax, S; color=(linecolor, linealpha), linewidth=linewidth,
                            transparency=(linealpha < 1))
    end

    if !isempty(pts)
        Makie.meshscatter!(ax, [Makie.Point3f(v...) for v in pts];
                           markersize=markersize, color=color)
    end

    return nothing
end

function _draw_colored_geometry!(ax, tris, tricolors, segs, segcolors, pts, ptcolors;
                                 alpha::Real=1.0, linewidth::Real=1.0,
                                 markersize::Real=0.05)
    #
    # Draw triangles, segments, and points with individual colors. All
    # triangles form a single mesh, so that they are depth sorted
    # correctly also by CairoMakie.
    #

    if !isempty(tris)
        P = [Makie.Point3f(v...) for t in tris for v in t]
        F = [GB.TriangleFace(3k-2, 3k-1, 3k) for k in 1:length(tris)]
        C = [Makie.RGBAf(Makie.to_color(c), alpha) for c in tricolors for _ in 1:3]
        Makie.mesh!(ax, GB.normal_mesh(P, F); color=C,
                    transparency=(alpha < 1), backlight=1.0f0)
    end

    if !isempty(segs)
        S = [Makie.Point3f(v...) for s in segs for v in s]
        C = [Makie.to_color(c) for c in segcolors for _ in 1:2]
        Makie.linesegments!(ax, S; color=C, linewidth=linewidth)
    end

    if !isempty(pts)
        Makie.meshscatter!(ax, [Makie.Point3f(v...) for v in pts];
                           markersize=markersize,
                           color=[Makie.to_color(c) for c in ptcolors])
    end

    return nothing
end

function _cell_indices(ec::EuclideanComplex, cells, fname::String)
    #
    # Convert a cell list given by labels or indices into indices
    #
    return CD._cells3d_indices(ec, cells, fname)
end

function _cellsubset_indices(ec::EuclideanComplex, subsets, fname::String)
    #
    # Convert cell subsets given by labels or indices into indices
    #
    return [_cell_indices(ec, s, fname) for s in subsets]
end

function _add_singletons(ec::EuclideanComplex, subsets::Vector{Vector{Int}})
    #
    # Add all cells which do not belong to a subset as singletons
    #
    used = falses(ec.ncells)
    for s in subsets
        used[s] .= true
    end
    return vcat(subsets, [[k] for k in 1:ec.ncells if !used[k]])
end

const _CI3D_COLORS = Dict(
    [0,0,0,0] => colorant"rgb(68,170,153)",
    [1,0,0,0] => colorant"rgb(51,34,136)",
    [0,1,0,0] => colorant"rgb(221,204,119)",
    [0,0,1,0] => colorant"rgb(17,119,51)",
    [0,0,0,1] => colorant"rgb(170,68,153)",
    [1,1,0,0] => colorant"rgb(204,102,119)",
    [0,1,1,0] => colorant"rgb(135,34,85)",
    [0,0,1,1] => colorant"rgb(136,204,238)",
)

function _ci_color3d(ci::Vector{Int})
    #
    # Color of a cell subset with a given Conley index. The colors agree
    # with the planar plots for indices of dimension at most two.
    #
    return get(_CI3D_COLORS, CD._ci_color3d_key(ci), colorant"rgb(153,153,51)")
end

function _default_colors(n::Int; light::Bool=false)
    #
    # The distinguishable colors used for Morse sets in the planar plots.
    # Light colors are used for multivector fields, where many cells of
    # different colors are next to each other.
    #
    seed = [colorant"royalblue4", colorant"royalblue3", colorant"steelblue1"]
    if light
        return distinguishable_colors(n, seed; dropseed=true,
                                      lchoices=range(55, 90; length=8))
    end
    return distinguishable_colors(n, seed; dropseed=true)
end
