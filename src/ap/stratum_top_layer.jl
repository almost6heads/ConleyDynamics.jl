export is_top_layer, stratum_top_layer, construct_ap_stratum_top

"""
    is_top_layer(lc::AbstractComplex, mvf::CellSubsets; connected::Bool=true)

Whether the acyclic partition `mvf` lies in the top layer of its own Morse
stratum, i.e. whether it admits no atomic refinement that preserves its
Morse vector.

By the local characterization of the top layer (every block is either a
hyperbolic singleton, a Forman doubleton, or a genuinely rigid multivector),
this comes down to a per-block test: `mvf` is in the top layer iff none of
its explicit blocks admits an atomic split with `delta = 0`. Implicit
singletons need no check, since they cannot be split further.

If `connected=true` (the default, matching `construct_ap_space`), only
splits into two topologically connected pieces are considered valid atomic
refinements (`connected_split_deltas`, the restricted `AP(X)` sense). If
`connected=false`, any conservative split counts (`block_split_deltas`, the
unrestricted `AP^d(X)` sense). This must match how `mvf` itself was
constructed/enumerated, or the test answers a different question than
intended.

Since `block_split_deltas` brute-forces all subsets of a block, this
inherits its `length(block) <= 20` limitation.

# Example
```jldoctest
julia> lc, mvf= example_forman1d();

julia> is_top_layer(lc, mvf)
true
```
"""
function is_top_layer(lc::AbstractComplex, mvf::CellSubsets; connected::Bool=true)
    #
    # Check whether mvf is in the top layer of its own Morse stratum
    #

    if typeof(mvf) == Vector{Vector{String}}
        mvfI = convert_cellsubsets(lc, mvf)
    else
        mvfI = mvf
    end

    for block in mvfI
        deltas = connected ? connected_split_deltas(lc, block) :
                              block_split_deltas(lc, block)
        if any(r -> iszero(r.delta), deltas)
            return false
        end
    end

    return true
end

"""
    stratum_top_layer(lc::AbstractComplex, ap::Vector{Vector{Vector{Int}}}; connected::Bool=true)

For every Morse-vector stratum represented in `ap`, find the indices of
`ap` belonging to its top layer -- the acyclic partitions of that Morse
vector which cannot be further refined without increasing the Morse
vector.

`ap` is typically the output of `construct_ap_space` or `construct_ap_stratum`.
`connected` must match how `ap` was constructed, since it selects whether a
block's conservative-splittability is tested via `connected_split_deltas`
(restricted `AP(X)`, the default) or `block_split_deltas` (unrestricted
`AP^d(X)`); see `is_top_layer`.

Returns `Dict{Vector{Int},Vector{Int}}`, in the same shape as
`stratum_partition`, but every value is filtered down to just the indices
in the top layer of that stratum (a subset of `stratum_partition(lc,
ap)[M]`).

A block's rigidity is memoized across the whole call (many partitions,
possibly across different strata, reuse identical blocks), so each
distinct block occurring anywhere in `ap` has its split spectrum computed
only once.
"""
function stratum_top_layer(lc::AbstractComplex,
                           ap::Vector{Vector{Vector{Int}}};
                           connected::Bool=true)
    #
    # Filter every stratum of ap down to its top layer
    #

    strata = stratum_partition(lc, ap)

    cache     = Dict{Set{Int},Bool}()
    cachelock = ReentrantLock()

    function block_rigid(block::Vector{Int})
        key = Set(block)
        lock(cachelock) do
            haskey(cache, key) && return cache[key]
            deltas = connected ? connected_split_deltas(lc, block) :
                                  block_split_deltas(lc, block)
            rigid = !any(r -> iszero(r.delta), deltas)
            cache[key] = rigid
            return rigid
        end
    end

    top = Dict{Vector{Int},Vector{Int}}()
    for (M, idx) in strata
        top[M] = filter(i -> all(block_rigid, ap[i]), idx)
    end

    return top
end

"""
    construct_ap_stratum_top(lc::AbstractComplex, target::Vector{Int}; connected::Bool=true)

Directly compute the top layer of the Morse stratum with Morse vector
`target`: the elements of `AP_target(X)` (or `AP^d_target(X)` if
`connected=false`) that admit no atomic refinement preserving `target`.

Equivalent to
```julia
ap = construct_ap_stratum(lc, target; connected=connected)
filter(w -> is_top_layer(lc, w; connected=connected), ap)
```
provided as a single call for convenience. Since the top layer is only
meaningful relative to the stratum it sits inside, this still has to
construct the whole stratum first (via `construct_ap_stratum`) before
filtering it down.

Returns a `Vector{Vector{Vector{Int}}}`, in the same format as
`construct_ap_space`/`construct_ap_stratum`.
"""
function construct_ap_stratum_top(lc::AbstractComplex,
                                  target::Vector{Int};
                                  connected::Bool=true)
    #
    # Construct only the top-layer elements of AP_target(X)/AP^d_target(X)
    #

    ap = construct_ap_stratum(lc, target; connected=connected)
    return filter(w -> is_top_layer(lc, w; connected=connected), ap)
end
