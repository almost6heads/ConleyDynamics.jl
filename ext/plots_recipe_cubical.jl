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

    # --- Optional MVF overlay ---
    if !isempty(data.mvf)
        mvf_int = data.mvf isa Vector{Vector{String}} ?
                      convert_cellsubsets(ec, data.mvf) : data.mvf

        bary_x = [sum(p[1] for p in ec.coords[k]) / length(ec.coords[k])
                  for k in 1:ec.ncells]
        bary_y = [sum(p[2] for p in ec.coords[k]) / length(ec.coords[k])
                  for k in 1:ec.ncells]

        pair_list = Tuple{Int,Int}[]
        critical  = Int[]
        mvf_cells = Set{Int}()
        for mv in mvf_int
            for k in mv; push!(mvf_cells, k); end
            if length(mv) == 1
                push!(critical, mv[1])
            elseif length(mv) == 2
                k1, k2 = mv[1], mv[2]
                if ec.dimensions[k1] < ec.dimensions[k2]
                    push!(pair_list, (k1, k2))
                else
                    push!(pair_list, (k2, k1))
                end
            end
        end
        # Cells absent from the MVF are implicitly critical
        for k in 1:ec.ncells
            k ∉ mvf_cells && push!(critical, k)
        end

        if !isempty(critical)
            @series begin
                seriestype        := :scatter
                markercolor       --> colorant"red3"
                markersize        --> 7
                markerstrokewidth --> 0
                [bary_x[k] for k in critical], [bary_y[k] for k in critical]
            end
        end

        if !isempty(pair_list)
            xs_tail = [bary_x[k1] for (k1, _) in pair_list]
            ys_tail = [bary_y[k1] for (k1, _) in pair_list]
            du      = [bary_x[k2] - bary_x[k1] for (k1, k2) in pair_list]
            dv      = [bary_y[k2] - bary_y[k1] for (k1, k2) in pair_list]
            @series begin
                seriestype := :quiver
                quiver     --> (du, dv)
                linecolor  --> colorant"red3"
                linewidth  --> 2
                xs_tail, ys_tail
            end
        end
    end
end

function ConleyDynamics.plot_cubical(ec::EuclideanComplex;
                                     mvf::CellSubsets=Vector{Vector{Int}}([]),
                                     pdim::Vector{Bool}=[true,true,true])
    data = CubicalComplexPlot(ec, mvf, pdim)
    return Plots.plot(data)
end
