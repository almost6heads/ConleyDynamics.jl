export SimplicialComplexPlot, CubicalComplexPlot,
       SimplicialMorsePlot, CubicalMorsePlot,
       MVFPlot, MVRegionPlot

"""
    SimplicialComplexPlot

Wrapper type for plotting a planar simplicial complex via Plots.jl.
Constructed automatically by `plot_simplicial`; can also be passed
directly to `plot(::SimplicialComplexPlot)`.
"""
struct SimplicialComplexPlot
    complex  :: EuclideanComplex
    mvf      :: CellSubsets
    pdim     :: Vector{Bool}
    labeldir :: Vector{Float64}
    labeldis :: Float64
end

"""
    CubicalComplexPlot

Wrapper type for plotting a planar cubical complex via Plots.jl.
"""
struct CubicalComplexPlot
    complex :: EuclideanComplex
    pdim    :: Vector{Bool}
end

"""
    SimplicialMorsePlot

Wrapper type for plotting a planar simplicial complex with Morse sets via Plots.jl.
"""
struct SimplicialMorsePlot
    complex   :: EuclideanComplex
    morsesets :: CellSubsets
    pdim      :: Vector{Bool}
    ci        :: Bool
end

"""
    CubicalMorsePlot

Wrapper type for plotting a planar cubical complex with Morse sets via Plots.jl.
"""
struct CubicalMorsePlot
    complex   :: EuclideanComplex
    morsesets :: CellSubsets
    pdim      :: Vector{Bool}
end

"""
    MVFPlot

Wrapper type for plotting a multivector field on a planar complex via Plots.jl.
Each multivector is rendered as a stadium (pill) shaped region.
"""
struct MVFPlot
    complex      :: EuclideanComplex
    mvf          :: CellSubsets
    pdim         :: Vector{Bool}
    tubefac      :: Float64
    mvfcolors    :: Vector{String}
    mvfalpha     :: Float64
    addcritical  :: Bool
end

"""
    MVRegionPlot

Wrapper type for plotting a single multivector on a planar complex via Plots.jl.
The multivector is rendered as an inflated convex hull per maximal cell.
"""
struct MVRegionPlot
    complex  :: EuclideanComplex
    mv       :: Vector{Int}
    pdim     :: Vector{Bool}
    tubefac  :: Float64
    mvfcolor :: String
    mvfalpha :: Float64
end

# Stub functions activated when Plots.jl is loaded as a weak dependency.
# Methods are added by ConleyDynamicsPlotsExt.

"""
    plot_simplicial(ec::EuclideanComplex; kwargs...) -> Plots.Plot

Plot a planar simplicial complex. Requires `using Plots` before `using ConleyDynamics`.

Optional keyword arguments:
- `mvf::CellSubsets=[]`: Forman vector field to overlay (arrows + critical cells)
- `labeldir::Vector{<:Real}=[]`: label directions per cell (0=E,1=N,2=W,3=S)
- `labeldis::Real=0.05`: label offset in data-space units
- `pdim::Vector{Bool}=[true,true,true]`: which dimensions (0,1,2) to draw
"""
function plot_simplicial end

"""
    plot_cubical(ec::EuclideanComplex; kwargs...) -> Plots.Plot

Plot a planar cubical complex. Requires `using Plots` before `using ConleyDynamics`.

Optional keyword arguments:
- `pdim::Vector{Bool}=[true,true,true]`: which dimensions (0,1,2) to draw
"""
function plot_cubical end

"""
    plot_simplicial_morse(ec::EuclideanComplex, morsesets::CellSubsets;
                          kwargs...) -> Plots.Plot

Plot a planar simplicial complex with Morse sets highlighted.
Requires `using Plots` before `using ConleyDynamics`.

Optional keyword arguments:
- `pdim::Vector{Bool}=[false,true,true]`: which dimensions (0,1,2) to draw
- `ci::Bool=false`: color Morse sets by their Conley index
"""
function plot_simplicial_morse end

"""
    plot_cubical_morse(ec::EuclideanComplex, morsesets::CellSubsets;
                       kwargs...) -> Plots.Plot

Plot a planar cubical complex with Morse sets highlighted.
Requires `using Plots` before `using ConleyDynamics`.

Optional keyword arguments:
- `pdim::Vector{Bool}=[false,true,true]`: which dimensions (0,1,2) to draw
"""
function plot_cubical_morse end

"""
    plot_simplicial_mvf(ec::EuclideanComplex, mvf::CellSubsets;
                        kwargs...) -> Plots.Plot

Plot a planar simplicial complex with multivector field regions shaded.
Each multivector is rendered as a single colored region, inset from
cell boundaries so adjacent multivectors are visually distinct.
Requires `using Plots` before `using ConleyDynamics`.

Optional keyword arguments:
- `pdim::Vector{Bool}=[true,true,true]`: which dimensions to show in background
- `tubefac::Real=0.05`: tube half-width as fraction of average edge length
- `mvfcolor::Union{String,Vector{String}}="darkorange"`: color(s) for multivector regions;
  a single string applies to all, a vector cycles through multivectors in order
- `mvfalpha::Real=0.3`: fill opacity (0.0 = transparent, 1.0 = opaque)
- `addcritical::Bool=true`: when true, cells absent from the MVF are shown as implicit singletons
"""
function plot_simplicial_mvf end

"""
    plot_cubical_mvf(ec::EuclideanComplex, mvf::CellSubsets;
                     kwargs...) -> Plots.Plot

Plot a planar cubical complex with multivector field regions shaded.
Requires `using Plots` before `using ConleyDynamics`.

Optional keyword arguments:
- `pdim::Vector{Bool}=[true,true,true]`: which dimensions to show in background
- `tubefac::Real=0.05`: tube half-width as fraction of average edge length
- `mvfcolor::Union{String,Vector{String}}="darkorange"`: color(s) for multivector regions
- `mvfalpha::Real=0.3`: fill opacity (0.0 = transparent, 1.0 = opaque)
- `addcritical::Bool=true`: when true, cells absent from the MVF are shown as implicit singletons
"""
function plot_cubical_mvf end

"""
    plot_simplicial_mv(ec::EuclideanComplex, mv; kwargs...) -> Plots.Plot

Plot a single multivector on a planar simplicial complex background.
The multivector is rendered as an inflated convex hull of barycenters,
one hull per maximal cell of the multivector.
Requires `using Plots` before `using ConleyDynamics`.

`mv` may be a `Vector{Int}` (cell indices) or `Vector{String}` (cell labels).

Optional keyword arguments:
- `pdim::Vector{Bool}=[true,true,true]`: which background dimensions to draw
- `tubefac::Real=0.05`: inflation radius as fraction of average edge length
- `mvfcolor::String="darkorange"`: X11 color name for the region
- `mvfalpha::Real=0.3`: fill opacity (0.0 = transparent, 1.0 = opaque)
"""
function plot_simplicial_mv end

"""
    plot_cubical_mv(ec::EuclideanComplex, mv; kwargs...) -> Plots.Plot

Plot a single multivector on a planar cubical complex background.
See `plot_simplicial_mv` for details.
"""
function plot_cubical_mv end

export plot_simplicial, plot_cubical,
       plot_simplicial_morse, plot_cubical_morse,
       plot_simplicial_mvf, plot_cubical_mvf,
       plot_simplicial_mv, plot_cubical_mv
