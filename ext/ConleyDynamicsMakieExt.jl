module ConleyDynamicsMakieExt

using ConleyDynamics
using Colors
using Makie
import Makie.GeometryBasics as GB

include("makie_geometry_utils.jl")
include("makie_plot_complex3d.jl")
include("makie_plot_morse3d.jl")
include("makie_plot_mvf3d.jl")

end
