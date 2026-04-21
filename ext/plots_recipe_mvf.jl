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

# Convex hull of pts (Vector of 2-element Float64 vectors).
# Returns indices into pts in CCW order (Andrew's monotone chain).
function _convex_hull(pts::Vector{Vector{Float64}})
    n = length(pts)
    n <= 1 && return collect(1:n)
    idx = sortperm(pts, by = p -> (p[1], p[2]))
    cross2d(o, a, b) = (a[1]-o[1])*(b[2]-o[2]) - (a[2]-o[2])*(b[1]-o[1])
    lower = Int[]
    for i in idx
        while length(lower) >= 2 && cross2d(pts[lower[end-1]], pts[lower[end]], pts[i]) <= 0
            pop!(lower)
        end
        push!(lower, i)
    end
    upper = Int[]
    for i in reverse(idx)
        while length(upper) >= 2 && cross2d(pts[upper[end-1]], pts[upper[end]], pts[i]) <= 0
            pop!(upper)
        end
        push!(upper, i)
    end
    return vcat(lower[1:end-1], upper[1:end-1])
end

# Inflated convex polygon: Minkowski sum of a CCW convex polygon with a disk of radius r.
# Each vertex becomes a circular arc; each edge becomes a parallel-offset straight side.
# n_cap controls arc resolution. Result is NaN-terminated and has no internal overlaps.
function _inflated_convex_polygon(pts::Vector{Vector{Float64}}, r::Float64, n_cap::Int=6)
    n = length(pts)
    n == 0 && return Float64[], Float64[]
    n == 1 && return _circle_polygon(pts[1], r)
    n == 2 && return _stadium_polygon(pts[1], pts[2], r, n_cap)

    # Outward normal of each edge: right normal of CCW direction = (dy/len, -dx/len)
    nrm = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        p1 = pts[i]; p2 = pts[mod1(i+1, n)]
        dx = p2[1]-p1[1]; dy = p2[2]-p1[2]; len = sqrt(dx^2+dy^2)
        nrm[i] = len < 1e-12 ? [0.0, 0.0] : [dy/len, -dx/len]
    end

    xs = Float64[]; ys = Float64[]
    for i in 1:n
        p      = pts[i]
        p_next = pts[mod1(i+1, n)]
        n_in   = nrm[mod1(i-1, n)]   # outward normal of incoming edge
        n_out  = nrm[i]               # outward normal of outgoing edge

        # Arc at vertex p: sweep CCW from n_in direction to n_out direction
        a0 = atan(n_in[2],  n_in[1])
        a1 = atan(n_out[2], n_out[1])
        while a1 < a0; a1 += 2π; end
        for j in 0:n_cap
            θ = a0 + (a1 - a0) * j / n_cap
            push!(xs, p[1] + r*cos(θ)); push!(ys, p[2] + r*sin(θ))
        end

        # Straight side to next vertex (last arc point already at p + r*n_out)
        push!(xs, p_next[1] + r*n_out[1]); push!(ys, p_next[2] + r*n_out[2])
    end
    push!(xs, xs[1]); push!(ys, ys[1])   # close
    push!(xs, NaN);   push!(ys, NaN)
    return xs, ys
end

# Cells in M_set with empty coboundary within M_set (nothing in M covers them).
function _maximal_cells_in_M(ec::EuclideanComplex, M_set::Set{Int})
    covered = Set{Int}()
    for j in M_set
        for i in sparse_get_nz_column(ec.boundary, j)
            push!(covered, i)
        end
    end
    return [k for k in M_set if k ∉ covered]
end

# NaN-separated polygon data representing one multivector mv.
# One inflated convex hull per maximal cell; reusable from MVFPlot recipe.
function _mvf_region_polygons(ec::EuclideanComplex, mv::Vector{Int}, r::Float64)
    M_set   = Set(mv)
    all_xs  = Float64[]
    all_ys  = Float64[]
    for cmax in _maximal_cells_in_M(ec, M_set)
        local_cells = [cmax]
        for c in _closure_cells(ec, cmax)
            c in M_set && push!(local_cells, c)
        end
        bpts     = [_barycenter(ec, c) for c in local_cells]
        hull_idx = _convex_hull(bpts)
        xs, ys   = _inflated_convex_polygon([bpts[i] for i in hull_idx], r)
        append!(all_xs, xs); append!(all_ys, ys)
    end
    return all_xs, all_ys
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

    # --- Multivector regions: inflated convex hull per maximal cell ---
    col = parse(Colorant, data.mvfcolor)

    # Cells absent from every explicit multivector are implicit singletons.
    mvf_cells = Set{Int}()
    for mv in mvf_int
        for k in mv; push!(mvf_cells, k); end
    end
    implicit_singletons = [[k] for k in 1:ec.ncells if k ∉ mvf_cells]

    for mv in vcat(mvf_int, implicit_singletons)
        xs, ys = _mvf_region_polygons(ec, mv, r)
        isempty(xs) && continue
        @series begin
            seriestype := :shape
            fillcolor  := col
            fillalpha  := data.mvfalpha
            linewidth  --> 0
            xs, ys
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

# -----------------------------------------------------------------------
# MVRegionPlot recipe — single multivector, inflated convex hull geometry
# -----------------------------------------------------------------------

@recipe function f(data::MVRegionPlot)
    ec = data.complex

    aspect_ratio --> :equal
    legend       --> false
    framestyle   --> :none

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

    # --- Multivector region ---
    xs, ys = _mvf_region_polygons(ec, data.mv, r)
    if !isempty(xs)
        @series begin
            seriestype := :shape
            fillcolor  := parse(Colorant, data.mvfcolor)
            fillalpha  := data.mvfalpha
            linewidth  --> 0
            xs, ys
        end
    end
end

function ConleyDynamics.plot_simplicial_mv(ec::EuclideanComplex,
                                           mv::Union{Vector{Int},Vector{String}};
                                           pdim::Vector{Bool}=[true,true,true],
                                           tubefac::Real=0.05,
                                           mvfcolor::String="darkorange",
                                           mvfalpha::Real=0.2)
    mv_int = mv isa Vector{String} ? convert_cellsubsets(ec, [mv])[1] : mv
    data = MVRegionPlot(ec, mv_int, pdim, Float64(tubefac), mvfcolor, Float64(mvfalpha))
    return Plots.plot(data)
end

function ConleyDynamics.plot_cubical_mv(ec::EuclideanComplex,
                                        mv::Union{Vector{Int},Vector{String}};
                                        pdim::Vector{Bool}=[true,true,true],
                                        tubefac::Real=0.05,
                                        mvfcolor::String="darkorange",
                                        mvfalpha::Real=0.2)
    mv_int = mv isa Vector{String} ? convert_cellsubsets(ec, [mv])[1] : mv
    data = MVRegionPlot(ec, mv_int, pdim, Float64(tubefac), mvfcolor, Float64(mvfalpha))
    return Plots.plot(data)
end
