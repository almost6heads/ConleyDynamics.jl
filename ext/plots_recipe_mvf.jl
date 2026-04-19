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

    # Right side: p1+r*n to p2+r*n
    push!(xs, p1[1] + r*n[1]); push!(ys, p1[2] + r*n[2])
    push!(xs, p2[1] + r*n[1]); push!(ys, p2[2] + r*n[2])

    # Right cap at p2: sweep from +n to -n (clockwise, outward from p2)
    for i in 1:n_cap-1
        angle = π * i / n_cap
        c = cos(angle); s = sin(angle)
        ex = c*n[1] + s*d[1]
        ey = c*n[2] + s*d[2]
        push!(xs, p2[1] + r*ex); push!(ys, p2[2] + r*ey)
    end
    push!(xs, p2[1] - r*n[1]); push!(ys, p2[2] - r*n[2])

    # Left side: p2-r*n to p1-r*n
    push!(xs, p1[1] - r*n[1]); push!(ys, p1[2] - r*n[2])

    # Left cap at p1: sweep from -n to +n (outward from p1)
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

# 2D cross product (scalar).
_cross2(a, b) = a[1]*b[2] - a[2]*b[1]

# Intersection of two lines: line1 through q1 in direction d1,
# line2 through q2 in direction d2. Returns the intersection point.
function _line_intersect(q1, d1, q2, d2)
    denom = _cross2(d1, d2)
    if abs(denom) < 1e-12
        return [(q1[1]+q2[1])/2, (q1[2]+q2[2])/2]  # parallel: return midpoint
    end
    t = _cross2([q2[1]-q1[1], q2[2]-q1[2]], d2) / denom
    return [q1[1] + t*d1[1], q1[2] + t*d1[2]]
end

# Inset polygon for a dim-2 face cell k.
# For each boundary edge: if the edge is in M_set, inset = 0 (extend to the edge);
# otherwise inset = r (pull back from that edge toward interior).
# Returns (xs, ys) NaN-terminated polygon, or empty if geometry fails.
function _face_inset_polygon(ec::EuclideanComplex, k::Int,
                              M_set::Set{Int}, r::Float64)
    cv      = ec.coords[k]
    nv      = length(cv)
    bedges  = sparse_get_nz_column(ec.boundary, k)

    # For each boundary edge, find its two vertex coords and the inset amount.
    # Represent the inset line as (anchor, direction, inset_normal * inset_dist).
    edge_lines = []   # (anchor, direction, normal_in * inset)
    for be in bedges
        ecoords = ec.coords[be]
        ea, eb  = ecoords[1], ecoords[2]
        dx = eb[1]-ea[1]; dy = eb[2]-ea[2]
        elen = sqrt(dx*dx + dy*dy)
        elen < 1e-14 && continue
        d = [dx/elen, dy/elen]
        # Normal candidates — pick the one pointing toward interior
        n1 = [-d[2], d[1]]
        # Check by dot with barycenter direction from midpoint
        bary = [sum(p[1] for p in cv)/nv, sum(p[2] for p in cv)/nv]
        mx = (ea[1]+eb[1])/2; my = (ea[2]+eb[2])/2
        dot1 = (bary[1]-mx)*n1[1] + (bary[2]-my)*n1[2]
        n_in = dot1 >= 0 ? n1 : [-n1[1], -n1[2]]
        ins = be in M_set ? 0.0 : r
        # Anchor on the inset line (any point on original edge + ins*n_in)
        anchor = [ea[1] + ins*n_in[1], ea[2] + ins*n_in[2]]
        push!(edge_lines, (anchor, d, ea, eb))
    end

    length(edge_lines) < 3 && return Float64[], Float64[]

    # For a triangle (nv=3): 3 edges, 3 vertices each at the intersection of 2 edges.
    # Determine which edges are adjacent to each face vertex by coordinate matching.
    tol = 1e-8
    new_pts = Vector{Vector{Float64}}(undef, nv)
    for i in 1:nv
        vi = cv[i]
        # Find the (up to) 2 edge_lines that contain vi as an endpoint
        adj = []
        for (anchor, d, ea, eb) in edge_lines
            on_ea = (ea[1]-vi[1])^2 + (ea[2]-vi[2])^2 < tol^2
            on_eb = (eb[1]-vi[1])^2 + (eb[2]-vi[2])^2 < tol^2
            (on_ea || on_eb) && push!(adj, (anchor, d))
        end
        if length(adj) >= 2
            new_pts[i] = _line_intersect(adj[1][1], adj[1][2], adj[2][1], adj[2][2])
        else
            # Fallback: move toward barycenter
            bary = [sum(p[1] for p in cv)/nv, sum(p[2] for p in cv)/nv]
            dvx = bary[1]-vi[1]; dvy = bary[2]-vi[2]
            dv  = sqrt(dvx^2+dvy^2)
            frac = dv > 0 ? min(r/dv, 0.45) : 0.0
            new_pts[i] = [vi[1]+frac*dvx, vi[2]+frac*dvy]
        end
    end

    if nv == 3
        xs = [new_pts[1][1], new_pts[2][1], new_pts[3][1], new_pts[1][1], NaN]
        ys = [new_pts[1][2], new_pts[2][2], new_pts[3][2], new_pts[1][2], NaN]
    else
        # Quad: [1,3,4,2] traversal (matching Luxor convention)
        xs = [new_pts[1][1], new_pts[3][1], new_pts[4][1], new_pts[2][1], new_pts[1][1], NaN]
        ys = [new_pts[1][2], new_pts[3][2], new_pts[4][2], new_pts[2][2], new_pts[1][2], NaN]
    end
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

    # --- Multivector shaded regions ---
    col = parse(Colorant, data.mvfcolor)

    for mv in mvf_int
        M_set = Set(mv)

        for k in mv
            dim = ec.dimensions[k]

            if dim == 0
                xs, ys = _circle_polygon(ec.coords[k][1], r)
                @series begin
                    seriestype := :shape
                    fillcolor  := col
                    fillalpha  --> 1.0
                    linewidth  --> 0
                    xs, ys
                end

            elseif dim == 1
                p1 = ec.coords[k][1]
                p2 = ec.coords[k][2]
                mid = [(p1[1]+p2[1])/2, (p1[2]+p2[2])/2]
                bv = sparse_get_nz_column(ec.boundary, k)
                cap1 = (length(bv) >= 1 && bv[1] in M_set) ? p1 : mid
                cap2 = (length(bv) >= 2 && bv[2] in M_set) ? p2 : mid
                xs, ys = _stadium_polygon(cap1, cap2, r)
                @series begin
                    seriestype := :shape
                    fillcolor  := col
                    fillalpha  --> 1.0
                    linewidth  --> 0
                    xs, ys
                end

            elseif dim == 2
                xs, ys = _face_inset_polygon(ec, k, M_set, r)
                isempty(xs) && continue
                @series begin
                    seriestype := :shape
                    fillcolor  := col
                    fillalpha  --> 1.0
                    linewidth  --> 0
                    xs, ys
                end
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
                                            tubefac::Real=0.15,
                                            mvfcolor::String="darkorange")
    data = MVFPlot(ec, mvf, pdim, Float64(tubefac), mvfcolor)
    return Plots.plot(data)
end

function ConleyDynamics.plot_cubical_mvf(ec::EuclideanComplex,
                                         mvf::CellSubsets;
                                         pdim::Vector{Bool}=[true,true,true],
                                         tubefac::Real=0.15,
                                         mvfcolor::String="darkorange")
    data = MVFPlot(ec, mvf, pdim, Float64(tubefac), mvfcolor)
    return Plots.plot(data)
end
