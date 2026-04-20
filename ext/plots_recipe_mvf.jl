# Recipe and convenience functions for MVFPlot

# -----------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------

function _avg_edge_length(ec::EuclideanComplex)
    total = 0.0
    count = 0
    for k in 1:ec.ncells
        ec.dimensions[k] == 1 || continue
        cv = ec.coords[k]
        dx = cv[2][1] - cv[1][1]
        dy = cv[2][2] - cv[1][2]
        total += sqrt(dx*dx + dy*dy)
        count += 1
    end
    return count > 0 ? total / count : 1.0
end

function _barycenter(ec::EuclideanComplex, k::Int)
    cv = ec.coords[k]
    n  = length(cv)
    [sum(p[1] for p in cv) / n, sum(p[2] for p in cv) / n]
end

# All cells strictly in the closure of k (k itself excluded), via BFS on boundary map.
function _closure_cells(ec::EuclideanComplex, k::Int)
    visited = Set{Int}()
    queue   = sparse_get_nz_column(ec.boundary, k)
    while !isempty(queue)
        c = pop!(queue)
        c in visited && continue
        push!(visited, c)
        for bc in sparse_get_nz_column(ec.boundary, c)
            bc ∉ visited && push!(queue, bc)
        end
    end
    return visited
end

# Circle of radius r centered at p, approximated as n-gon, NaN-terminated.
function _circle_polygon(p, r, n=16)
    xs = [p[1] + r * cos(2π * i / n) for i in 0:n]
    ys = [p[2] + r * sin(2π * i / n) for i in 0:n]
    push!(xs, NaN); push!(ys, NaN)
    return xs, ys
end

# Stadium (pill) from p1 to p2 with half-width r, NaN-terminated.
# n_cap points per semicircular cap. Degenerates to a circle when p1 ≈ p2.
function _stadium_polygon(p1, p2, r, n_cap=8)
    dx = p2[1] - p1[1]; dy = p2[2] - p1[2]
    len = sqrt(dx*dx + dy*dy)
    if len < 1e-12
        return _circle_polygon(p1, r)
    end
    d = [dx/len, dy/len]
    n = [-d[2], d[1]]

    xs = Float64[]
    ys = Float64[]

    push!(xs, p1[1] + r*n[1]); push!(ys, p1[2] + r*n[2])
    push!(xs, p2[1] + r*n[1]); push!(ys, p2[2] + r*n[2])

    for i in 1:n_cap-1
        angle = π * i / n_cap
        c = cos(angle); s = sin(angle)
        ex = c*n[1] + s*d[1]
        ey = c*n[2] + s*d[2]
        push!(xs, p2[1] + r*ex); push!(ys, p2[2] + r*ey)
    end
    push!(xs, p2[1] - r*n[1]); push!(ys, p2[2] - r*n[2])

    push!(xs, p1[1] - r*n[1]); push!(ys, p1[2] - r*n[2])

    for i in 1:n_cap-1
        angle = π * i / n_cap
        c = cos(angle); s = sin(angle)
        ex = c*(-n[1]) + s*(-d[1])
        ey = c*(-n[2]) + s*(-d[2])
        push!(xs, p1[1] + r*ex); push!(ys, p1[2] + r*ey)
    end

    push!(xs, p1[1] + r*n[1]); push!(ys, p1[2] + r*n[2])
    push!(xs, NaN); push!(ys, NaN)
    return xs, ys
end

# -----------------------------------------------------------------------
# Recipe
# -----------------------------------------------------------------------

@recipe function f(data::MVFPlot)
    ec = data.complex

    aspect_ratio --> :equal
    legend       --> false
    framestyle   --> :none

    mvf_int = data.mvf isa Vector{Vector{String}} ?
                  convert_cellsubsets(ec, data.mvf) : data.mvf

    w_avg = _avg_edge_length(ec)
    r     = data.tubefac * w_avg

    # --- Background complex (faded) ---
    if data.pdim[3]
        xs, ys = Float64[], Float64[]
        for k in ec.ncells:-1:1
            ec.dimensions[k] == 2 || continue
            cv = ec.coords[k]
            if length(cv) == 3
                append!(xs, [p[1] for p in cv]); push!(xs, cv[1][1], NaN)
                append!(ys, [p[2] for p in cv]); push!(ys, cv[1][2], NaN)
            else
                push!(xs, cv[1][1], cv[3][1], cv[4][1], cv[2][1], cv[1][1], NaN)
                push!(ys, cv[1][2], cv[3][2], cv[4][2], cv[2][2], cv[1][2], NaN)
            end
        end
        if !isempty(xs)
            @series begin
                seriestype := :shape
                fillcolor  --> colorant"steelblue1"
                fillalpha  --> 0.3
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
                linealpha  --> 0.3
                linewidth  --> 1.0
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
                markeralpha       --> 0.3
                markersize        --> 4
                markerstrokewidth --> 0
                vxs, vys
            end
        end
    end

    # --- Multivector regions: barycenter-graph stadiums ---
    col = parse(Colorant, data.mvfcolor)

    # Cells absent from every explicit multivector are implicit singletons.
    mvf_cells = Set{Int}()
    for mv in mvf_int
        for k in mv; push!(mvf_cells, k); end
    end
    implicit_singletons = [[k] for k in 1:ec.ncells if k ∉ mvf_cells]
    all_mvf = vcat(mvf_int, implicit_singletons)

    for mv in all_mvf
        M_set = Set(mv)

        # Stadium for each pair (k, bk) where bk is in the closure of k and in M_set.
        # Only emit when dim(k) > dim(bk) to avoid drawing each pair twice.
        for k in mv
            ec.dimensions[k] == 0 && continue
            for bk in _closure_cells(ec, k)
                bk in M_set || continue
                xs, ys = _stadium_polygon(_barycenter(ec, k), _barycenter(ec, bk), r)
                @series begin
                    seriestype := :shape
                    fillcolor  := col
                    fillalpha  := data.mvfalpha
                    linewidth  --> 0
                    xs, ys
                end
            end
        end

        # Circle at every cell's barycenter (covers isolated cells and endpoints).
        for k in mv
            xs, ys = _circle_polygon(_barycenter(ec, k), r)
            @series begin
                seriestype := :shape
                fillcolor  := col
                fillalpha  := data.mvfalpha
                linewidth  --> 0
                xs, ys
            end
        end
    end
end

# -----------------------------------------------------------------------
# Convenience functions
# -----------------------------------------------------------------------

function ConleyDynamics.plot_simplicial_mvf(ec::EuclideanComplex,
                                            mvf::CellSubsets;
                                            pdim::Vector{Bool}=[true,true,true],
                                            tubefac::Real=0.05,
                                            mvfcolor::String="darkorange",
                                            mvfalpha::Real=0.2)
    data = MVFPlot(ec, mvf, pdim, Float64(tubefac), mvfcolor, Float64(mvfalpha))
    return Plots.plot(data)
end

function ConleyDynamics.plot_cubical_mvf(ec::EuclideanComplex,
                                         mvf::CellSubsets;
                                         pdim::Vector{Bool}=[true,true,true],
                                         tubefac::Real=0.05,
                                         mvfcolor::String="darkorange",
                                         mvfalpha::Real=0.2)
    data = MVFPlot(ec, mvf, pdim, Float64(tubefac), mvfcolor, Float64(mvfalpha))
    return Plots.plot(data)
end
