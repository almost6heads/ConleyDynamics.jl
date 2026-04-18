# Recipe and convenience function for CubicalComplexPlot

@recipe function f(data::CubicalComplexPlot)
    ec = data.complex

    aspect_ratio --> :equal
    legend       --> false
    framestyle   --> :none

    # --- Quads (dim 2) — rendered first ---
    # ec.coords[k] for a cubical 2-cell has 4 entries in order [p1,p2,p3,p4].
    # The correct polygon traversal (matching Luxor convention) is [1,3,4,2].
    if data.pdim[3]
        xs, ys = Float64[], Float64[]
        for k in ec.ncells:-1:1
            ec.dimensions[k] == 2 || continue
            cv = ec.coords[k]
            push!(xs, cv[1][1], cv[3][1], cv[4][1], cv[2][1], cv[1][1], NaN)
            push!(ys, cv[1][2], cv[3][2], cv[4][2], cv[2][2], cv[1][2], NaN)
        end
        if !isempty(xs)
            @series begin
                seriestype := :shape
                fillcolor  --> colorant"steelblue1"
                fillalpha  --> 0.7
                linewidth  --> 0
                xs, ys
            end
        end
    end

    # --- Edges (dim 1) ---
    if data.pdim[2]
        xs, ys = Float64[], Float64[]
        for k in ec.ncells:-1:1
            ec.dimensions[k] == 1 || continue
            cv = ec.coords[k]
            push!(xs, cv[1][1], cv[2][1], NaN)
            push!(ys, cv[1][2], cv[2][2], NaN)
        end
        if !isempty(xs)
            @series begin
                seriestype := :path
                linecolor  --> colorant"royalblue3"
                linewidth  --> 1.5
                xs, ys
            end
        end
    end

    # --- Vertices (dim 0) ---
    if data.pdim[1]
        vxs = [ec.coords[k][1][1] for k in 1:ec.ncells if ec.dimensions[k] == 0]
        vys = [ec.coords[k][1][2] for k in 1:ec.ncells if ec.dimensions[k] == 0]
        if !isempty(vxs)
            @series begin
                seriestype        := :scatter
                markercolor       --> colorant"royalblue4"
                markersize        --> 5
                markerstrokewidth --> 0
                vxs, vys
            end
        end
    end
end

function ConleyDynamics.plot_cubical(ec::EuclideanComplex;
                                     pdim::Vector{Bool}=[true,true,true])
    data = CubicalComplexPlot(ec, pdim)
    return Plots.plot(data)
end
