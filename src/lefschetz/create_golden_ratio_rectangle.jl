export create_golden_ratio_rectangle

# =============================================================================
# Internal helper: recursive golden-ratio subdivision
# =============================================================================

function _subdivide_gr!(
        leaves    :: Vector{NTuple{4,Float64}},
        xmin      :: Float64, xmax :: Float64,
        ymin      :: Float64, ymax :: Float64,
        σ         :: Float64,
        sdfunction,
        sdmin     :: Int,
        sdmax     :: Int,
        depth     :: Int
    )

    if depth > sdmax
        push!(leaves, (xmin, xmax, ymin, ymax))
        return
    end

    # Decide whether to split: forced below sdmin, delegated to sdfunction above.
    should_split = (depth < sdmin) || sdfunction([xmin, ymin], [xmax, ymax])

    if !should_split
        push!(leaves, (xmin, xmax, ymin, ymax))
        return
    end

    # Split along the longer edge (x-direction if tied)
    width  = xmax - xmin
    height = ymax - ymin

    if width >= height
        frac    = rand() < 0.5 ? σ : (1.0 - σ)
        split_x = xmin + frac * width
        _subdivide_gr!(leaves, xmin, split_x, ymin, ymax, σ, sdfunction, sdmin, sdmax, depth + 1)
        _subdivide_gr!(leaves, split_x, xmax, ymin, ymax, σ, sdfunction, sdmin, sdmax, depth + 1)
    else
        frac    = rand() < 0.5 ? σ : (1.0 - σ)
        split_y = ymin + frac * height
        _subdivide_gr!(leaves, xmin, xmax, ymin, split_y, σ, sdfunction, sdmin, sdmax, depth + 1)
        _subdivide_gr!(leaves, xmin, xmax, split_y, ymax, σ, sdfunction, sdmin, sdmax, depth + 1)
    end
end

# =============================================================================
# Public function
# =============================================================================

"""
    create_golden_ratio_rectangle(bmin, bmax;
                                  sigma      = (sqrt(5)-1)/2,
                                  sdfunction = (_,_) -> false,
                                  sdmin      = 0,
                                  sdmax      = 7,
                                  p          = 2)

Create a planar `EuclideanComplex` over a field of characteristic `p` by
recursively subdividing the rectangle `[bmin[1],bmax[1]] × [bmin[2],bmax[2]]`
using the golden-ratio subdivision rule.

# Subdivision rule

At each recursive step the current box is examined:

* If its recursion depth is > `sdmax`, it becomes a leaf (hard cap).
* If its depth is < `sdmin`, it is always split (guaranteed minimum).
* Otherwise `sdfunction([xmin,ymin],[xmax,ymax])` decides: `true` splits,
  `false` makes it a leaf.

The default `sdfunction` always returns `false`, so with default arguments
exactly `sdmin` levels of uniform subdivision are performed.

When a box is split, the longer edge is divided into two segments of lengths
`σ·a` and `(1-σ)·a`; which piece goes to which sub-box is chosen uniformly
at random.  The default `sigma = (√5-1)/2 ≈ 0.618` is the golden ratio
minus 1: starting from a square it produces sub-rectangles whose aspect
ratios stay within the three values {1, φ, φ²} under further subdivision
(see Fig. 3.1 of Cochran et al., SISC 2013).

# Return value

An `EuclideanComplex` whose 2-cells are all rectangles (leaf boxes).  The
boundary of each rectangle is made up of the minimal set of line segments
whose endpoints come from the corners of adjacent leaf boxes, oriented
exactly as in a cubical complex: each bottom or right segment has coefficient
+1 and each top or left segment has coefficient -1 in the boundary of its
enclosing rectangle.
"""
function create_golden_ratio_rectangle(
        bmin :: Vector{<:Real},
        bmax :: Vector{<:Real};
        sigma      :: Real = (sqrt(5.0) - 1.0) / 2.0,
        sdfunction        = (_, _) -> false,
        sdmin      :: Int  = 0,
        sdmax      :: Int  = 7,
        p          :: Int  = 2
    )

    # Validate that sdfunction accepts two Vector{Float64} arguments
    if !applicable(sdfunction, Float64[0.0, 0.0], Float64[1.0, 1.0])
        throw(ArgumentError(
            "create_golden_ratio_rectangle: sdfunction must be callable as " *
            "sdfunction(bmin::Vector{Float64}, bmax::Vector{Float64}); " *
            "got $(typeof(sdfunction))"))
    end

    # -------------------------------------------------------------------------
    # Phase 1: build leaf boxes via recursive subdivision
    # -------------------------------------------------------------------------

    leaves = NTuple{4,Float64}[]
    _subdivide_gr!(leaves,
                   Float64(bmin[1]), Float64(bmax[1]),
                   Float64(bmin[2]), Float64(bmax[2]),
                   Float64(sigma), sdfunction, sdmin, sdmax, 0)

    nr = length(leaves)

    # -------------------------------------------------------------------------
    # Phase 2: collect vertices (corners of leaf boxes)
    # -------------------------------------------------------------------------

    vertex_dict = Dict{Tuple{Float64,Float64}, Int}()
    vertex_xy   = Tuple{Float64,Float64}[]

    for (xmin, xmax, ymin, ymax) in leaves
        for xy in ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax))
            if !haskey(vertex_dict, xy)
                push!(vertex_xy, xy)
                vertex_dict[xy] = length(vertex_xy)
            end
        end
    end

    nv = length(vertex_xy)

    # Build per-coordinate lookup tables for edge-split detection.
    # by_y[y] = vertices at that y-value, sorted by x (ascending).
    # by_x[x] = vertices at that x-value, sorted by y (ascending).

    by_y = Dict{Float64, Vector{Tuple{Float64,Int}}}()
    by_x = Dict{Float64, Vector{Tuple{Float64,Int}}}()

    for (vi, (x, y)) in enumerate(vertex_xy)
        push!(get!(by_y, y, Tuple{Float64,Int}[]), (x, vi))
        push!(get!(by_x, x, Tuple{Float64,Int}[]), (y, vi))
    end
    for vlist in values(by_y); sort!(vlist; by = first); end
    for vlist in values(by_x); sort!(vlist; by = first); end

    # -------------------------------------------------------------------------
    # Phase 3: build edges and rectangle boundary information
    # -------------------------------------------------------------------------
    # Canonical orientation:
    #   horizontal edge — start = geometrically left  (smaller x)
    #   vertical   edge — start = geometrically bottom (smaller y)
    # Edge boundary: +1 on end vertex, -1 on start vertex.
    # Rectangle boundary sign: +1 for bottom / right edges, -1 for top / left.

    edge_start = Int[]
    edge_end   = Int[]
    edge_dict  = Dict{Tuple{Int,Int}, Int}()   # (start_vi, end_vi) → edge index

    function get_or_add_edge_h!(v_left::Int, v_right::Int)
        key = (v_left, v_right)
        if !haskey(edge_dict, key)
            push!(edge_start, v_left)
            push!(edge_end,   v_right)
            edge_dict[key] = length(edge_start)
        end
        return edge_dict[key]
    end

    function get_or_add_edge_v!(v_bot::Int, v_top::Int)
        key = (v_bot, v_top)
        if !haskey(edge_dict, key)
            push!(edge_start, v_bot)
            push!(edge_end,   v_top)
            edge_dict[key] = length(edge_start)
        end
        return edge_dict[key]
    end

    # rect_bnds[i] = list of (edge_index, sign) for the i-th leaf box
    rect_bnds = Vector{Vector{Tuple{Int,Int}}}(undef, nr)

    for (ri, (xmin, xmax, ymin, ymax)) in enumerate(leaves)
        bnd = Tuple{Int,Int}[]

        # Bottom edge: y = ymin, left-to-right, coefficient +1
        bot = filter(t -> xmin <= t[1] <= xmax, by_y[ymin])
        for k in 1:length(bot) - 1
            eidx = get_or_add_edge_h!(bot[k][2], bot[k+1][2])
            push!(bnd, (eidx, 1))
        end

        # Right edge: x = xmax, bottom-to-top, coefficient +1
        rgt = filter(t -> ymin <= t[1] <= ymax, by_x[xmax])
        for k in 1:length(rgt) - 1
            eidx = get_or_add_edge_v!(rgt[k][2], rgt[k+1][2])
            push!(bnd, (eidx, 1))
        end

        # Top edge: y = ymax, canonical left-to-right, coefficient -1
        top = filter(t -> xmin <= t[1] <= xmax, by_y[ymax])
        for k in 1:length(top) - 1
            eidx = get_or_add_edge_h!(top[k][2], top[k+1][2])
            push!(bnd, (eidx, -1))
        end

        # Left edge: x = xmin, canonical bottom-to-top, coefficient -1
        lft = filter(t -> ymin <= t[1] <= ymax, by_x[xmin])
        for k in 1:length(lft) - 1
            eidx = get_or_add_edge_v!(lft[k][2], lft[k+1][2])
            push!(bnd, (eidx, -1))
        end

        rect_bnds[ri] = bnd
    end

    ne     = length(edge_start)
    ntotal = nv + ne + nr

    # -------------------------------------------------------------------------
    # Phase 4: labels, dimensions, boundary matrix, and coords
    # -------------------------------------------------------------------------

    label_width = length(string(max(nv, ne, nr)))

    labels = Vector{String}(undef, ntotal)
    for i in 1:nv
        labels[i]           = "v" * lpad(string(i), label_width, '0')
    end
    for i in 1:ne
        labels[nv + i]      = "e" * lpad(string(i), label_width, '0')
    end
    for i in 1:nr
        labels[nv + ne + i] = "r" * lpad(string(i), label_width, '0')
    end

    dims = [fill(0, nv); fill(1, ne); fill(2, nr)]

    # Field scalars
    if p == 0
        tone  = 1 // 1
        tzero = 0 // 1
    else
        tone  = Int(1)
        tzero = Int(0)
    end

    # Boundary matrix via (row, col, value) triples
    r_idx = Int[]
    c_idx = Int[]
    b_val = Vector{typeof(tone)}()

    # 1-cell boundaries: end_vertex (+1), start_vertex (-1)
    for i in 1:ne
        col = nv + i
        push!(c_idx, col); push!(r_idx, edge_end[i]);   push!(b_val,  tone)
        push!(c_idx, col); push!(r_idx, edge_start[i]); push!(b_val, -tone)
    end

    # 2-cell boundaries
    for i in 1:nr
        col = nv + ne + i
        for (eidx, sgn) in rect_bnds[i]
            push!(c_idx, col)
            push!(r_idx, nv + eidx)
            push!(b_val, sgn == 1 ? tone : -tone)
        end
    end

    B = sparse_from_lists(ntotal, ntotal, p, tzero, tone, r_idx, c_idx, b_val)

    # Per-cell coordinates for the EuclideanComplex
    cell_coords = Vector{Vector{Vector{Float64}}}(undef, ntotal)

    for i in 1:nv
        (x, y) = vertex_xy[i]
        cell_coords[i] = [[x, y]]
    end

    for i in 1:ne
        (xs, ys) = vertex_xy[edge_start[i]]
        (xe, ye) = vertex_xy[edge_end[i]]
        cell_coords[nv + i] = [[xs, ys], [xe, ye]]
    end

    for i in 1:nr
        (xmin, xmax, ymin, ymax) = leaves[i]
        cell_coords[nv + ne + i] = [
            [xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]
        ]
    end

    return EuclideanComplex(labels, dims, B, cell_coords; validate=false)
end
