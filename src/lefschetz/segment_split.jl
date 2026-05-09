export split_segments

using Base.Threads

# =============================================================================
# Geometry primitives
# =============================================================================

@inline _cross2(u, v) = u[1]*v[2] - u[2]*v[1]
@inline _norm2sq(v)   = v[1]^2 + v[2]^2
@inline _norm2(v)     = sqrt(_norm2sq(v))
@inline _distsq(a, b) = (a[1]-b[1])^2 + (a[2]-b[2])^2


#=
    _intersect_point(p1, p2, q1, q2, eps_abs) -> Union{Vector{Float64}, Nothing}

Return the intersection point of segments p1-p2 and q1-q2, or `nothing`.
Endpoint touches (t or s at exactly 0 or 1) are included so that
T-intersections and shared endpoints are recorded as candidate points.
Collinear overlapping segments are handled separately in the main loop.
=#
@inline function _intersect_point(p1, p2, q1, q2, eps_abs)
    d1    = p2 .- p1
    d2    = q2 .- q1
    denom = _cross2(d1, d2)

    # Relative tolerance: scales with segment lengths.
    abs(denom) < eps_abs * (_norm2(d1) * _norm2(d2) + eps_abs) && return nothing

    w = q1 .- p1
    t = _cross2(w, d2) / denom
    s = _cross2(w, d1) / denom

    (t > -eps_abs && t < 1.0 + eps_abs && s > -eps_abs && s < 1.0 + eps_abs) || return nothing
    return p1 .+ clamp(t, 0.0, 1.0) .* d1
end


# =============================================================================
# Union-Find  (path compression + union by rank)
# =============================================================================

struct _UnionFind
    parent :: Vector{Int}
    rank   :: Vector{Int}
end

_UnionFind(n::Int) = _UnionFind(collect(1:n), zeros(Int, n))

function _uf_find!(uf::_UnionFind, x::Int)
    while uf.parent[x] != x
        uf.parent[x] = uf.parent[uf.parent[x]]   # path halving
        x = uf.parent[x]
    end
    return x
end

function _uf_union!(uf::_UnionFind, x::Int, y::Int)
    rx, ry = _uf_find!(uf, x), _uf_find!(uf, y)
    rx == ry && return
    if uf.rank[rx] < uf.rank[ry]; rx, ry = ry, rx; end
    uf.parent[ry] = rx
    uf.rank[rx] == uf.rank[ry] && (uf.rank[rx] += 1)
end


# =============================================================================
# Snap-rounding
# =============================================================================

#=
    _snap_points(pts; eps_snap) -> (reps, label)

Cluster `pts` so that every pair of points within Euclidean distance `eps_snap`
belongs to the same cluster.  Clustering is computed by Union-Find over all
O(m²) pairs with distance <= eps_snap.

Returns
  reps  :: Vector{Vector{Float64}}  — centroid of each cluster, in cluster order
  label :: Vector{Int}              — label[i] is the 1-based cluster index of pts[i]
=#
function _snap_points(pts::Vector{Vector{Float64}}; eps_snap::Float64)
    m      = length(pts)
    uf     = _UnionFind(m)
    eps_snap2 = eps_snap^2

    for i in 1:m
        for j in (i+1):m
            _distsq(pts[i], pts[j]) <= eps_snap2 && _uf_union!(uf, i, j)
        end
    end

    # Assign compact cluster indices.
    root_to_ci = Dict{Int,Int}()
    for i in 1:m
        r = _uf_find!(uf, i)
        haskey(root_to_ci, r) || (root_to_ci[r] = length(root_to_ci) + 1)
    end

    label = [root_to_ci[_uf_find!(uf, i)] for i in 1:m]
    nc    = length(root_to_ci)

    # Centroid of each cluster.
    sums   = [zeros(2) for _ in 1:nc]
    counts = zeros(Int, nc)
    for i in 1:m
        c = label[i]
        sums[c] .+= pts[i]
        counts[c] += 1
    end
    reps = [sums[c] ./ counts[c] for c in 1:nc]

    return reps, label
end


# =============================================================================
# Main function
# =============================================================================

"""
    split_segments(segments; eps_snap=1e-6, eps_abs=1e-10)
        -> (Vector{Vector{Float64}}, Vector{Tuple{Int,Int}})

Compute a planar straight-line arrangement from a collection of line segments:
the output covers the same point-set as the input but no two output segments
share an interior point (no crossings, no overlaps).

# Arguments

- `segments`: a `Vector` of segments, where each segment is a length-2 vector
  of 2D endpoints, e.g. `[[x1, y1], [x2, y2]]`.

# Keyword arguments

- `eps_snap = 1e-6`: distance threshold for snap-rounding. Any two candidate
  points (original endpoints or computed intersection points) within this
  Euclidean distance are identified as the same vertex.
- `eps_abs = 1e-10`: arithmetic guard for degenerate configurations (parallel
  or zero-length segments). Automatically raised to `eps_snap * 1e-4` if
  smaller.

# Returns

A tuple `(points, edges)` where:
- `points :: Vector{Vector{Float64}}` is the deduplicated vertex list
  (centroids of snap-rounding clusters), and
- `edges :: Vector{Tuple{Int,Int}}` lists each output sub-segment as a pair
  `(i, j)` with `i < j`, using 1-based indices into `points`.

# Algorithm

1. **Candidate collection** (parallelised): all `2n` endpoints and every
   pairwise intersection point are collected into a flat candidate array.
   Each thread accumulates into its own private buffer to avoid locking.
2. **Snap-rounding**: all candidate points are clustered with Union-Find so
   that any two within distance `eps_snap` are identified. The representative
   of each cluster is its centroid.
3. **Segment splitting** (parallelised): each input segment is split at every
   snapped candidate point that lies on it (within `eps_snap`), yielding a
   list of consecutive cluster-index pairs.
4. **Output assembly**: duplicate cluster-index pairs are removed and returned
   in canonical `(i, j)` form with `i < j`.

# Threading

The two O(n²) loops are parallelised with `Threads.@threads`. Start Julia
with `--threads auto` or set `JULIA_NUM_THREADS` to use all available cores;
the algorithm is correct with any number of threads, including 1.

# Example

```julia
using ConleyDynamics

# Two crossing diagonals of the unit square
segs = [[[0.0, 0.0], [1.0, 1.0]],
        [[0.0, 1.0], [1.0, 0.0]]]

points, edges = split_segments(segs)
# points has 5 entries: the 4 corners plus the crossing at [0.5, 0.5]
# edges  has 4 entries: one sub-segment on each side of the crossing
```
"""
function split_segments(
        segments :: Vector{Vector{Vector{Float64}}};
        eps_snap    :: Float64 = 1e-6,
        eps_abs  :: Float64 = 1e-10
    )
    n       = length(segments)
    eps_abs = max(eps_abs, eps_snap * 1e-4)

    nt = nthreads()

    # Phase 1: collect all candidate points (parallel)
    #
    # Endpoints: 2n points, one pair per segment.
    endpoint_pts = Vector{Vector{Float64}}(undef, 2n)
    for i in 1:n
        endpoint_pts[2i-1] = copy(segments[i][1])
        endpoint_pts[2i]   = copy(segments[i][2])
    end

    # Intersection points: distribute the upper-triangle pair loop across threads.
    inter_bufs = [Vector{Vector{Float64}}() for _ in 1:nt]

    @threads for i in 1:n
        tid  = threadid()
        p1   = segments[i][1]
        p2   = segments[i][2]
        d1   = p2 .- p1
        len1 = _norm2(d1)

        for j in (i+1):n
            q1 = segments[j][1]
            q2 = segments[j][2]
            d2 = q2 .- q1

            denom = _cross2(d1, d2)

            if abs(denom) < eps_abs * (len1 * _norm2(d2) + eps_abs)
                # Parallel lines: check collinearity.
                w = q1 .- p1
                abs(_cross2(w, d1)) > eps_abs * (len1 + eps_abs) && continue
                # Collinear: record all four endpoints so overlapping intervals
                # are split correctly.
                push!(inter_bufs[tid], copy(q1))
                push!(inter_bufs[tid], copy(q2))
                push!(inter_bufs[tid], copy(p1))
                push!(inter_bufs[tid], copy(p2))
            else
                pt = _intersect_point(p1, p2, q1, q2, eps_abs)
                pt !== nothing && push!(inter_bufs[tid], pt)
            end
        end
    end

    # Concatenate all candidate points: endpoints first, then intersections.
    inter_pts     = reduce(vcat, inter_bufs)
    candidate_pts = vcat(endpoint_pts, inter_pts)
    m             = length(candidate_pts)

    # ep[i]: indices of segment i's endpoints in candidate_pts.
    ep = [(2i-1, 2i) for i in 1:n]

    # Phase 2: snap-rounding
    reps, label = _snap_points(candidate_pts; eps_snap = eps_snap)

    # Phase 3: split each segment at snapped candidate points (parallel)
    #
    # For each segment: project every candidate point onto the segment axis,
    # keep those within eps_snap of the line, sort by projection parameter,
    # deduplicate by cluster label, and emit consecutive cluster-index pairs.
    raw_per_seg = [Tuple{Int,Int}[] for _ in 1:n]

    @threads for i in 1:n
        li, ri   = ep[i]
        p_ci     = label[li]
        q_ci     = label[ri]
        p_snap   = reps[p_ci]
        q_snap   = reps[q_ci]
        seg_vec  = q_snap .- p_snap
        seg_len2 = _norm2sq(seg_vec)

        # Skip segments that collapsed to a point after snapping.
        seg_len2 < eps_snap^2 && continue

        # Collect (parameter t, cluster index) for every candidate on this segment.
        on_seg = Tuple{Float64,Int}[]

        for k in 1:m
            c  = candidate_pts[k]
            dc = c .- p_snap

            t_raw = (dc[1]*seg_vec[1] + dc[2]*seg_vec[2]) / seg_len2

            t_raw < -eps_abs && continue
            t_raw > 1.0 + eps_abs && continue

            dist_to_line = abs(_cross2(dc, seg_vec)) / sqrt(seg_len2)
            dist_to_line > eps_snap + eps_abs && continue

            push!(on_seg, (clamp(t_raw, 0.0, 1.0), label[k]))
        end

        sort!(on_seg; by = first)

        # Deduplicate clusters: for each cluster keep the t closest to the
        # cluster centroid's own projection.
        best = Dict{Int, Float64}()
        for (t, ci) in on_seg
            if !haskey(best, ci)
                best[ci] = t
            else
                rep_dc = reps[ci] .- p_snap
                rep_t  = (rep_dc[1]*seg_vec[1] + rep_dc[2]*seg_vec[2]) / seg_len2
                if abs(t - rep_t) < abs(best[ci] - rep_t)
                    best[ci] = t
                end
            end
        end

        deduped = sort!([(t, ci) for (ci, t) in best]; by = first)

        buf = raw_per_seg[i]
        for k in 1:length(deduped)-1
            ca = deduped[k][2]
            cb = deduped[k+1][2]
            ca != cb && push!(buf, (ca, cb))
        end
    end

    # Phase 4: assemble output — deduplicate edges in canonical (i<j) form.
    edges      = Tuple{Int,Int}[]
    seen_edges = Set{Tuple{Int,Int}}()

    for buf in raw_per_seg
        for (ca, cb) in buf
            key = ca < cb ? (ca, cb) : (cb, ca)
            key ∈ seen_edges && continue
            push!(seen_edges, key)
            push!(edges, key)
        end
    end

    return reps, edges
end

# Convenience wrapper: accept any nested-vector input with numeric element types.
function split_segments(
        segments :: Vector{<:Vector{<:Vector{<:Real}}};
        kwargs...)
    split_segments(
        [[Float64.(p) for p in s] for s in segments]; kwargs...)
end
