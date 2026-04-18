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

# Returns (xs, ys) polygon for the inset region of cell k.
# w_avg: average edge length (used to set vertex circle radius).
function _inset_polygon(ec::EuclideanComplex, k::Int,
                        insetfac::Float64, edgewidth::Float64, w_avg::Float64)
    cv = ec.coords[k]
    dim = ec.dimensions[k]

    if dim == 0
        # Circle approximated as 16-gon
        r = insetfac * w_avg / 2.0
        p = cv[1]
        n_pts = 16
        xs = [p[1] + r * cos(2π * i / n_pts) for i in 0:n_pts]
        ys = [p[2] + r * sin(2π * i / n_pts) for i in 0:n_pts]
        push!(xs, NaN); push!(ys, NaN)
        return xs, ys

    elseif dim == 1
        p1, p2 = cv[1], cv[2]
        m = [(p1[1]+p2[1])/2, (p1[2]+p2[2])/2]
        a = [p1[1] + insetfac*(m[1]-p1[1]), p1[2] + insetfac*(m[2]-p1[2])]
        b = [p2[1] + insetfac*(m[1]-p2[1]), p2[2] + insetfac*(m[2]-p2[2])]
        dx = p2[1] - p1[1]; dy = p2[2] - p1[2]
        len = sqrt(dx*dx + dy*dy)
        len < 1e-14 && (len = 1.0)
        d = [dx/len, dy/len]
        n = [-d[2], d[1]]
        w = edgewidth
        xs = [a[1]+w*n[1], b[1]+w*n[1], b[1]-w*n[1], a[1]-w*n[1], a[1]+w*n[1], NaN]
        ys = [a[2]+w*n[2], b[2]+w*n[2], b[2]-w*n[2], a[2]-w*n[2], a[2]+w*n[2], NaN]
        return xs, ys

    elseif dim == 2
        nv = length(cv)
        bx = sum(p[1] for p in cv) / nv
        by = sum(p[2] for p in cv) / nv
        # Inset each vertex toward barycenter
        pts = [[p[1] + insetfac*(bx - p[1]), p[2] + insetfac*(by - p[2])] for p in cv]
        if nv == 3
            xs = [pts[1][1], pts[2][1], pts[3][1], pts[1][1], NaN]
            ys = [pts[1][2], pts[2][2], pts[3][2], pts[1][2], NaN]
        else
            # Quad: use [1,3,4,2] traversal (matching Luxor convention)
            xs = [pts[1][1], pts[3][1], pts[4][1], pts[2][1], pts[1][1], NaN]
            ys = [pts[1][2], pts[3][2], pts[4][2], pts[2][2], pts[1][2], NaN]
        end
        return xs, ys
    end
    return Float64[], Float64[]
end

# Returns indices (into cv_high) of vertices shared with cv_low.
# Matching is by coordinate with tolerance.
function _shared_vertex_indices(cv_high, cv_low, tol=1e-8)
    shared = Int[]
    for (i, q) in enumerate(cv_high)
        for p in cv_low
            if (q[1]-p[1])^2 + (q[2]-p[2])^2 < tol^2
                push!(shared, i)
                break
            end
        end
    end
    return shared
end

# Returns (xs, ys) bridge polygon connecting the inset of k_low (face) to k_high.
function _bridge_polygon(ec::EuclideanComplex, k_high::Int, k_low::Int,
                         insetfac::Float64, edgewidth::Float64)
    cv_high = ec.coords[k_high]
    cv_low  = ec.coords[k_low]
    dim_high = ec.dimensions[k_high]
    dim_low  = ec.dimensions[k_low]

    # vertex-edge bridge: no explicit bridge needed (vertex circle and edge
    # rectangle visually overlap when insetfac is small)
    if dim_low == 0 && dim_high == 1
        return Float64[], Float64[]
    end

    # edge-to-2cell bridge (triangle or quad)
    if dim_low == 1 && dim_high == 2
        p1, p2 = cv_low[1], cv_low[2]
        nv_high = length(cv_high)
        bx_high = sum(p[1] for p in cv_high) / nv_high
        by_high = sum(p[2] for p in cv_high) / nv_high

        # Shared vertices from the high cell's perspective
        shared_idx = _shared_vertex_indices(cv_high, cv_low)
        length(shared_idx) < 2 && return Float64[], Float64[]

        s1 = cv_high[shared_idx[1]]
        s2 = cv_high[shared_idx[2]]

        # Match shared vertices to p1, p2 order
        d1_p1 = (s1[1]-p1[1])^2 + (s1[2]-p1[2])^2
        d1_p2 = (s1[1]-p2[1])^2 + (s1[2]-p2[2])^2
        if d1_p1 < d1_p2
            q1_shared, q2_shared = s1, s2  # s1 matches p1, s2 matches p2
        else
            q1_shared, q2_shared = s2, s1  # s1 matches p2, s2 matches p1
        end

        # Triangle-side inset corners (inset toward high-cell barycenter)
        v1_t = [q1_shared[1] + insetfac*(bx_high - q1_shared[1]),
                q1_shared[2] + insetfac*(by_high - q1_shared[2])]
        v2_t = [q2_shared[1] + insetfac*(bx_high - q2_shared[1]),
                q2_shared[2] + insetfac*(by_high - q2_shared[2])]

        # Edge inset endpoints
        mx = (p1[1]+p2[1])/2; my = (p1[2]+p2[2])/2
        a = [p1[1] + insetfac*(mx-p1[1]), p1[2] + insetfac*(my-p1[2])]
        b = [p2[1] + insetfac*(mx-p2[1]), p2[2] + insetfac*(my-p2[2])]

        # Normal of edge pointing toward the high cell's interior
        dx = p2[1]-p1[1]; dy = p2[2]-p1[2]
        len = sqrt(dx*dx + dy*dy); len < 1e-14 && (len = 1.0)
        d = [dx/len, dy/len]
        n = [-d[2], d[1]]
        # Find an opposite vertex to orient the normal
        opp_indices = setdiff(1:nv_high, shared_idx)
        if !isempty(opp_indices)
            opp = cv_high[opp_indices[1]]
            if (n[1]*(opp[1]-p1[1]) + n[2]*(opp[2]-p1[2])) < 0
                n = [-n[1], -n[2]]
            end
        end

        # Edge-side rectangle corners facing the high cell
        v1_e = [a[1] + edgewidth*n[1], a[2] + edgewidth*n[2]]
        v2_e = [b[1] + edgewidth*n[1], b[2] + edgewidth*n[2]]

        xs = [v1_t[1], v2_t[1], v2_e[1], v1_e[1], v1_t[1], NaN]
        ys = [v1_t[2], v2_t[2], v2_e[2], v1_e[2], v1_t[2], NaN]
        return xs, ys
    end

    return Float64[], Float64[]
end

# Returns all (k_high, k_low) incidence pairs within multivector mv.
function _boundary_pairs_in(ec::EuclideanComplex, mv::Vector{Int})
    mv_set = Set(mv)
    pairs  = Tuple{Int,Int}[]
    for k_high in mv
        ec.dimensions[k_high] == 0 && continue
        for k_low in sparse_get_nz_column(ec.boundary, k_high)
            if k_low in mv_set
                push!(pairs, (k_high, k_low))
            end
        end
    end
    return pairs
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

    insetfac  = data.insetfac
    w_avg     = _avg_edge_length(ec)
    edgewidth = iszero(data.edgewidth) ? insetfac * w_avg / 4.0 : data.edgewidth

    # --- Background complex (faded) ---
    # For dim-2, use vertex count to distinguish triangle (3) from quad (4).
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

    # --- Multivector shaded regions ---
    seed = [colorant"royalblue4", colorant"royalblue3", colorant"steelblue1"]
    cols = distinguishable_colors(length(mvf_int), seed; dropseed=true)

    for (m, mv) in enumerate(mvf_int)
        col = cols[m]

        # Inset polygons for each cell in the multivector
        for k in mv
            xs, ys = _inset_polygon(ec, k, insetfac, edgewidth, w_avg)
            isempty(xs) && continue
            @series begin
                seriestype := :shape
                fillcolor  := col
                fillalpha  --> 0.7
                linewidth  --> 0
                xs, ys
            end
        end

        # Bridge polygons connecting adjacent cells within the multivector
        for (k_high, k_low) in _boundary_pairs_in(ec, mv)
            xs, ys = _bridge_polygon(ec, k_high, k_low, insetfac, edgewidth)
            isempty(xs) && continue
            @series begin
                seriestype := :shape
                fillcolor  := col
                fillalpha  --> 0.7
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
                                            insetfac::Real=0.2,
                                            edgewidth::Real=0.0)
    data = MVFPlot(ec, mvf, pdim, Float64(insetfac), Float64(edgewidth))
    return Plots.plot(data)
end

function ConleyDynamics.plot_cubical_mvf(ec::EuclideanComplex,
                                         mvf::CellSubsets;
                                         pdim::Vector{Bool}=[true,true,true],
                                         insetfac::Real=0.2,
                                         edgewidth::Real=0.0)
    data = MVFPlot(ec, mvf, pdim, Float64(insetfac), Float64(edgewidth))
    return Plots.plot(data)
end
