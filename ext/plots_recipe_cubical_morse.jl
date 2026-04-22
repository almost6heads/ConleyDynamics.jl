# Recipe and convenience function for CubicalMorsePlot

function _ci_color_cubical(ci_val)
    ci_map = Dict(
        [1,0,0] => colorant"rgb(51,34,136)",
        [0,1,0] => colorant"rgb(221,204,119)",
        [0,0,1] => colorant"rgb(17,119,51)",
        [1,1,0] => colorant"rgb(204,102,119)",
        [0,1,1] => colorant"rgb(135,34,85)",
    )
    return get(ci_map, ci_val, colorant"rgb(153,153,51)")
end

@recipe function f(data::CubicalMorsePlot)
    ec = data.complex

    aspect_ratio --> :equal
    legend       --> false
    framestyle   --> :none

    ms_int = data.morsesets isa Vector{Vector{String}} ?
                 convert_cellsubsets(ec, data.morsesets) : data.morsesets

    # Optionally add cells absent from every explicit set as implicit singletons.
    if data.addcritical
        ms_cells = Set{Int}()
        for ms in ms_int; for k in ms; push!(ms_cells, k); end; end
        implicit = [[k] for k in 1:ec.ncells if k ∉ ms_cells]
        ms_int = vcat(ms_int, implicit)
    end

    if data.ci
        cols = [_ci_color_cubical(conley_index(ec, ms_int[m]))
                for m in eachindex(ms_int)]
    else
        seed = [colorant"royalblue4", colorant"royalblue3", colorant"steelblue1"]
        cols = distinguishable_colors(length(ms_int), seed; dropseed=true)
    end

    # --- Background cubical complex (faded) ---
    # ec.coords[k] for dim-2 has 4 entries; polygon traversal [1,3,4,2].
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
                fillalpha  --> 0.4
                linewidth  --> 0
                xs, ys
            end
        end
    end

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
                linealpha  --> 0.4
                linewidth  --> 1.5
                xs, ys
            end
        end
    end

    if data.pdim[1]
        vxs = [ec.coords[k][1][1] for k in 1:ec.ncells if ec.dimensions[k] == 0]
        vys = [ec.coords[k][1][2] for k in 1:ec.ncells if ec.dimensions[k] == 0]
        if !isempty(vxs)
            @series begin
                seriestype        := :scatter
                markercolor       --> colorant"royalblue4"
                markeralpha       --> 0.4
                markersize        --> 5
                markerstrokewidth --> 0
                vxs, vys
            end
        end
    end

    # --- Morse sets ---
    for m in eachindex(ms_int)
        ms  = ms_int[m]
        col = cols[m]

        if data.pdim[3]
            xs, ys = Float64[], Float64[]
            for k in ms
                ec.dimensions[k] == 2 || continue
                cv = ec.coords[k]
                push!(xs, cv[1][1], cv[3][1], cv[4][1], cv[2][1], cv[1][1], NaN)
                push!(ys, cv[1][2], cv[3][2], cv[4][2], cv[2][2], cv[1][2], NaN)
            end
            if !isempty(xs)
                @series begin
                    seriestype := :shape
                    fillcolor  := col
                    fillalpha  --> 0.6
                    linewidth  --> 0
                    xs, ys
                end
            end
        end

        if data.pdim[2]
            xs, ys = Float64[], Float64[]
            for k in ms
                ec.dimensions[k] == 1 || continue
                cv = ec.coords[k]
                push!(xs, cv[1][1], cv[2][1], NaN)
                push!(ys, cv[1][2], cv[2][2], NaN)
            end
            if !isempty(xs)
                @series begin
                    seriestype := :path
                    linecolor  := col
                    linealpha  --> 0.6
                    linewidth  --> 3
                    xs, ys
                end
            end
        end

        if data.pdim[1]
            vxs = [ec.coords[k][1][1] for k in ms if ec.dimensions[k] == 0]
            vys = [ec.coords[k][1][2] for k in ms if ec.dimensions[k] == 0]
            if !isempty(vxs)
                @series begin
                    seriestype        := :scatter
                    markercolor       := col
                    markeralpha       --> 0.6
                    markersize        --> 7
                    markerstrokewidth --> 0
                    vxs, vys
                end
            end
        end
    end
end

function ConleyDynamics.plot_cubical_morse(ec::EuclideanComplex,
                                           morsesets::CellSubsets;
                                           pdim::Vector{Bool}=[false,true,true],
                                           ci::Bool=false,
                                           addcritical::Bool=false)
    data = CubicalMorsePlot(ec, morsesets, pdim, ci, addcritical)
    return Plots.plot(data)
end
