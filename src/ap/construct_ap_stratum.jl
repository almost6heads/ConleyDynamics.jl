export construct_ap_stratum

"""
    construct_ap_stratum(lc::AbstractComplex, target::Vector{Int}; connected::Bool=true)

Construct only the elements of `AP(X)` (or `AP^c(X)` if `connected=true`,
the default) whose Morse vector equals `target`, without enumerating the
rest of the poset.

Uses the same merge-based search as `construct_ap_space`, but prunes a
branch the moment its current Morse vector fails to dominate `target` in
every coordinate: refinement never decreases the Morse vector (`delta =
beta(A)+beta(B)-beta(C) >= 0` for any merge of two blocks `A`, `B` into
`C`, by Lemma 3.3 of doc/morse-vector-strata.md -- "Monotonicity... holds
for all refinements, not only atomic ones"), so once some coordinate has
dropped below `target[k]`, no further merge can ever bring it back and the
branch can safely be dropped. Nodes whose Morse vector matches `target`
exactly are still expanded further, since additional Morse-vector-
preserving (`delta=0`) merges can still produce further distinct
partitions belonging to the same stratum.

Each accepted merge updates the running Morse vector incrementally from
the two merged blocks' already-known Conley indices (`beta(A)+beta(B) -
beta(C)`, memoized in a shared cache keyed by block) rather than
recomputing the whole partition's Morse vector from scratch at every node.
This incremental update is only ever applied to blocks belonging to an
already-`mvf_is_acyclic`-accepted partition -- `conley_index` errors on a
subset that is not locally closed, and a *rejected* merge candidate's
merged block carries no such guarantee, so the acyclicity check always
runs first.

Whether this is actually faster than calling `construct_ap_space` and
filtering by Morse vector afterward is target-dependent, not automatic:
`target = beta(X)` is a *global* lower bound on the Morse vector of every
element of `AP(X)` (the "Extremes" property, Section 3.3), so pruning
against it never triggers and this degenerates to the full enumeration
with extra per-node overhead. For a genuinely intermediate target,
however, most branches do get cut, and the saving can be large -- e.g. for
a 13-cell example where the connected `AP^c(X)` has 40232 elements, fixing
an intermediate stratum of 537 elements takes on the order of a tenth of a
second here versus several seconds for full enumeration plus filtering.
The main intended use is targeting a single stratum of the *unrestricted*
`AP(X)` on complexes where full unrestricted enumeration
(`construct_ap_space(lc; connected=false)`) is no longer tractable at all.

Returns a `Vector{Vector{Vector{Int}}}`, in the same format as
`construct_ap_space`.
"""
function construct_ap_stratum(lc::AbstractComplex, target::Vector{Int}; connected::Bool=true)
    #
    # Construct only the AP(X)/AP^c(X) elements with Morse vector == target
    #

    @assert lc.ncells>2 "You are just testing edge cases, I refuse to play.."
    @assert length(target) == lc.dim+1 "target must have length lc.dim+1 = $(lc.dim+1)"

    # Shared cache of Conley indices, keyed by block (a block's Conley
    # index depends only on its own cell content, not on how it arose, so
    # it can be safely reused across every partition it ever appears in).
    # Guarded by a lock since the search below runs multithreaded.

    betacache = Dict{Set{Int},Vector{Int}}()
    cachelock = ReentrantLock()

    function beta_of(block::Set{Int})
        lock(cachelock) do
            haskey(betacache, block) && return betacache[block]
            b = conley_index(lc, collect(block))
            betacache[block] = b
            return b
        end
    end

    # Create the initial finest multivector field and its Morse vector

    cmvf  = Set(Set.([[k] for k in 1:lc.ncells]))
    total = zeros(Int, lc.dim + 1)
    for c in 1:lc.ncells
        total .+= beta_of(Set([c]))
    end

    hits = Set{Set{Set{Int}}}()
    total == target && push!(hits, cmvf)

    frontier = [(cmvf, total)]

    # Loop through the levels, pruning branches that can no longer reach
    # target and collecting every partition found to match it exactly

    for k in 1:lc.ncells-1

        ch = Channel{Tuple{Set{Set{Int}},Vector{Int}}}(Inf)

        Threads.@threads for (lmvf, ltotal) in frontier

            mvecs = collect(lmvf)

            for m1 in 1:length(mvecs)-1
                mvec1 = mvecs[m1]
                nbhd1 = connected ? Set(lefschetz_neighbors(lc, collect(mvec1))) : nothing
                for m2 in m1+1:length(mvecs)
                    mvec2 = mvecs[m2]

                    if connected && length(intersect(nbhd1, mvec2)) == 0
                        continue
                    end

                    newblock = union(mvec1, mvec2)
                    newmvf   = setdiff(mvecs, [mvec1, mvec2])
                    push!(newmvf, newblock)

                    # Acyclicity first, always -- only an accepted block is
                    # guaranteed locally closed, hence safe for beta_of

                    if !mvf_is_acyclic(lc, collect.(newmvf))
                        continue
                    end

                    beta1 = beta_of(mvec1)
                    beta2 = beta_of(mvec2)
                    betaC = beta_of(newblock)
                    newtotal = ltotal .- beta1 .- beta2 .+ betaC

                    # Monotonicity: refinement never decreases the Morse
                    # vector, so once a coordinate has dropped below
                    # target it can never recover -- prune this branch.

                    any(newtotal .< target) && continue

                    put!(ch, (Set(newmvf), newtotal))
                end
            end
        end

        close(ch)
        frontier = unique(first, collect(ch))

        for (s, t) in frontier
            t == target && push!(hits, s)
        end
    end

    return convert_partition_mvf.(collect(hits))
end
