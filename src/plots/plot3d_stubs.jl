export plot_simplicial3d, plot_simplicial3d!, plot_cubical3d, plot_cubical3d!
export plot_simplicial3d_morse, plot_simplicial3d_morse!
export plot_cubical3d_morse, plot_cubical3d_morse!
export plot_simplicial3d_mv, plot_simplicial3d_mv!, plot_cubical3d_mv, plot_cubical3d_mv!
export plot_simplicial3d_mvf, plot_simplicial3d_mvf!, plot_cubical3d_mvf, plot_cubical3d_mvf!
export plot_cells3d!

# Stub functions activated when Makie.jl is loaded as a weak dependency,
# usually via GLMakie, CairoMakie, or WGLMakie. Methods are added by
# ConleyDynamicsMakieExt.

"""
    plot_simplicial3d(ec::EuclideanComplex; kwargs...) -> Makie.Figure
    plot_simplicial3d!(ax, ec::EuclideanComplex; kwargs...)

Plot a simplicial complex embedded in three dimensions with Makie.

The first form creates a new figure and returns it, the second form draws
into an existing Makie axis `ax`, such as an `LScene` or an `Axis3`. The
complex can be a tetrahedral mesh, or a two-dimensional complex embedded
in three dimensions, such as a triangulated surface. By default, only the
boundary surface of the tetrahedra is shown, semi-transparent, together
with its edges.

Optional keyword arguments:
- `pdim::Vector{Bool}=[false,true,true,true]`: which cell dimensions
  (0, 1, 2, 3) to draw. Three-dimensional cells are drawn as surfaces.
- `hull_only::Bool=true`: draw only the boundary surface of the
  tetrahedra and its edges, instead of all faces and edges.
- `clip=nothing`: a pair `(n, offset)` with a normal vector `n` and a
  number `offset`. Tetrahedra whose centroid `c` satisfies `n⋅c > offset`
  are removed, which shows the interior of the complex.
- `alpha::Real=0.3`: opacity of the faces.
- `facecolor=colorant"steelblue1"`, `edgecolor=colorant"royalblue3"`,
  `vertexcolor=colorant"royalblue4"`: colors of the cells.
- `linewidth::Real=0.75`: width of the edges.
- `markersize=nothing`: radius of the vertex spheres, in data units. The
  default is a fraction of the typical edge length.

The first form also accepts the following keyword arguments:
- `axis::Symbol=:scene`: `:scene` creates an interactive `LScene`
  without axis decorations, and `:axis3` creates an `Axis3` with ticks
  and equal aspect ratios.
- `size=(800,600)`: the figure size.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`. For interactive rotation and correct transparency, use
`GLMakie`. For figures in files, use `CairoMakie`.
"""
function plot_simplicial3d end

"""
    plot_simplicial3d!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_simplicial3d`.
"""
function plot_simplicial3d! end

"""
    plot_cubical3d(ec::EuclideanComplex; kwargs...) -> Makie.Figure
    plot_cubical3d!(ax, ec::EuclideanComplex; kwargs...)

Plot a cubical complex embedded in three dimensions with Makie.

Squares are drawn as two triangles each, and cubes by their boundary
squares. The keyword arguments are the same as for `plot_simplicial3d`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_cubical3d end

"""
    plot_cubical3d!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_cubical3d`.
"""
function plot_cubical3d! end

"""
    plot_simplicial3d_morse(ec::EuclideanComplex, morsesets::CellSubsets;
                            kwargs...) -> Makie.Figure
    plot_simplicial3d_morse!(ax, ec::EuclideanComplex, morsesets::CellSubsets;
                             kwargs...)

Plot a simplicial complex embedded in three dimensions with Morse sets
highlighted, using Makie.

The complex is drawn faintly in the background, and every Morse set is
drawn in its own color, according to its dimension: tetrahedra by the
boundary surface of their union, triangles which are not faces of
tetrahedra in the Morse set as triangles, edges which are not faces of
higher-dimensional cells as thick lines, and isolated vertices as
spheres.

Optional keyword arguments:
- `ci::Bool=false`: color the Morse sets by their Conley index.
- `addcritical::Bool=false`: cells which do not belong to any Morse set
  are shown as additional singleton Morse sets.
- `colors=nothing`: a vector of colors for the Morse sets. By default,
  the same distinguishable colors as in the planar plots are used.
- `alpha::Real=0.9`: opacity of the Morse sets.
- `bgalpha::Real=0.08`: opacity of the background faces.
- `pdim::Vector{Bool}=[false,true,true,true]`: which cell dimensions to
  draw in the background.
- `hull_only`, `clip`, `linewidth`, `markersize`: as for
  `plot_simplicial3d`. The `clip` keyword applies to the background and
  the Morse sets, and `linewidth=4` applies to the edges of Morse sets.

The first form also accepts the keywords `axis` and `size`, see
`plot_simplicial3d`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_simplicial3d_morse end

"""
    plot_simplicial3d_morse!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_simplicial3d_morse`.
"""
function plot_simplicial3d_morse! end

"""
    plot_cubical3d_morse(ec::EuclideanComplex, morsesets::CellSubsets;
                         kwargs...) -> Makie.Figure
    plot_cubical3d_morse!(ax, ec::EuclideanComplex, morsesets::CellSubsets;
                          kwargs...)

Plot a cubical complex embedded in three dimensions with Morse sets
highlighted, using Makie. The keyword arguments are the same as for
`plot_simplicial3d_morse`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_cubical3d_morse end

"""
    plot_cubical3d_morse!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_cubical3d_morse`.
"""
function plot_cubical3d_morse! end

"""
    plot_simplicial3d_mv(ec::EuclideanComplex, mv; kwargs...) -> Makie.Figure
    plot_simplicial3d_mv!(ax, ec::EuclideanComplex, mv; kwargs...)

Plot a single multivector of a simplicial complex embedded in three
dimensions, using Makie.

The multivector `mv` is a vector of cell indices or labels. It is drawn
in the same way as a Morse set in `plot_simplicial3d_morse`, on top of
the faint background complex. This is useful for checking the
multivectors created by `create_spatial_mvf`.

Optional keyword arguments:
- `color=colorant"darkorange"`: the color of the multivector.
- `alpha`, `bgalpha`, `pdim`, `hull_only`, `clip`, `linewidth`,
  `markersize`: as for `plot_simplicial3d_morse`.

The first form also accepts the keywords `axis` and `size`, see
`plot_simplicial3d`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_simplicial3d_mv end

"""
    plot_simplicial3d_mv!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_simplicial3d_mv`.
"""
function plot_simplicial3d_mv! end

"""
    plot_cubical3d_mv(ec::EuclideanComplex, mv; kwargs...) -> Makie.Figure
    plot_cubical3d_mv!(ax, ec::EuclideanComplex, mv; kwargs...)

Plot a single multivector of a cubical complex embedded in three
dimensions, using Makie. The keyword arguments are the same as for
`plot_simplicial3d_mv`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_cubical3d_mv end

"""
    plot_cubical3d_mv!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_cubical3d_mv`.
"""
function plot_cubical3d_mv! end

"""
    plot_simplicial3d_mvf(ec::EuclideanComplex, mvf::CellSubsets;
                          kwargs...) -> Makie.Figure
    plot_simplicial3d_mvf!(ax, ec::EuclideanComplex, mvf::CellSubsets;
                           kwargs...)

Plot a multivector field on a simplicial complex embedded in three
dimensions, using Makie.

Every multivector with at least two cells is drawn in the same way as a
Morse set in `plot_simplicial3d_morse`, and consecutive multivectors get
different colors. Since multivectors fill the complex, such plots are
only readable for small complexes, or together with the `clip` keyword.

Optional keyword arguments:
- `colors=nothing`: a vector of colors which is cycled through. By
  default, twelve distinguishable colors are used.
- `ci::Bool=false`: color each multivector by its Conley index instead.
- `addcritical::Bool=false`: also draw the cells which do not belong to
  any multivector of `mvf`, which are critical cells.
- `alpha`, `bgalpha`, `pdim`, `hull_only`, `clip`, `linewidth`,
  `markersize`: as for `plot_simplicial3d_morse`.

The first form also accepts the keywords `axis` and `size`, see
`plot_simplicial3d`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_simplicial3d_mvf end

"""
    plot_simplicial3d_mvf!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_simplicial3d_mvf`.
"""
function plot_simplicial3d_mvf! end

"""
    plot_cubical3d_mvf(ec::EuclideanComplex, mvf::CellSubsets;
                       kwargs...) -> Makie.Figure
    plot_cubical3d_mvf!(ax, ec::EuclideanComplex, mvf::CellSubsets;
                        kwargs...)

Plot a multivector field on a cubical complex embedded in three
dimensions, using Makie. The keyword arguments are the same as for
`plot_simplicial3d_mvf`.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_cubical3d_mvf end

"""
    plot_cubical3d_mvf!(ax, args...; kwargs...)

Draw into an existing Makie axis `ax`. See `plot_cubical3d_mvf`.
"""
function plot_cubical3d_mvf! end

"""
    plot_cells3d!(ax, ec::EuclideanComplex, cells; kwargs...)

Draw a subset of cells of a complex embedded in three dimensions into an
existing Makie axis `ax`.

The vector `cells` contains cell indices or labels. The cells are drawn
according to their dimension, in the same way as a Morse set in
`plot_simplicial3d_morse`. This is the basic building block of all
three-dimensional plotting functions, and it can be used to add further
cell subsets to a plot, such as exit sets or constraint surfaces.

Optional keyword arguments:
- `color=colorant"darkorange"`: the color of the cells.
- `alpha::Real=1.0`: opacity of the faces.
- `pdim::Vector{Bool}=[true,true,true,true]`: which cell dimensions to draw.
- `hull_only`, `clip`, `markersize`: as for `plot_simplicial3d`.
- `linewidth::Real=4`: width of the edges.

Requires `Makie.jl`, which is loaded by `GLMakie`, `CairoMakie`, or
`WGLMakie`.
"""
function plot_cells3d! end
