export construct_ap_space

"""
    construct_ap_space(lc::AbstractComplex; connected::Bool=true)

Construct the space of acyclic partitions of `lc`.

Starting from the partition of `lc` into singletons (Morse vector `f(X)`,
the finest possible acyclic partition), repeatedly merges pairs of blocks
of an already-found acyclic partition into a single block, keeping only
merges that preserve acyclicity (`mvf_is_acyclic`), until every pair of
blocks at every level has been tried. This reaches every element of the
poset `AP(X)` (see `stratum_partition` and `stratum_adjacency` for its
Morse-vector stratification), ordered by atomic refinement
(`is_atomic_refinement`).

If `connected=true` (the default), a merge is only attempted when the two
blocks are topologically adjacent (`lefschetz_neighbors`), so only
partitions built from connected multivectors are reached -- this is
`AP^c(X)`, the sub-poset of practical interest for e.g. Forman-vector-field
style examples, and is dramatically cheaper to enumerate. If
`connected=false`, every pair of blocks is tried regardless of adjacency,
reaching the full, unrestricted `AP(X)`, including partitions with
disconnected multivectors; this is needed whenever a question genuinely
requires disconnected pieces (a weight-1 Morse-vector jump, for instance,
is only *guaranteed* uniform when drawn from the unrestricted `AP(X)` --
see `stratum_adjacency`). The unrestricted enumeration can be substantially
larger than the connected one, and may not be tractable for complexes with
more than a handful of cells; `construct_ap_stratum` targets a single
Morse-vector stratum directly and can remain tractable well past that
point.

Returns a `Vector{Vector{Vector{Int}}}`: each entry is one element of
`AP(X)` (or `AP^c(X)`) as a multivector field in integer form (singletons
omitted, matching the usual `CellSubsets` convention).
"""
function construct_ap_space(lc::AbstractComplex; connected::Bool=true)
    #
    # Construct the space of acyclic partitions
    #

    @assert lc.ncells>2 "You are just testing edge cases, I refuse to play.."

    # Create the initial finest multivector field

    cmvf = Set(Set.([[k] for k in 1:lc.ncells]))

    # Initialize the data structure to hold AP(X)

    ap = Vector{Vector{Set{Set{Int}}}}()
    push!(ap, [cmvf,])

    # Loop through the levels to create all acyclic partitions

    for k in 1:lc.ncells-1

        # Set up a channel to hold all acyclic partitions

        ch = Channel{Set{Set{Int}}}(Inf)

        # Loop through the APs of the previous level

        Threads.@threads for lmvf in ap[k]

            # Extract all multivectors

            mvecs = collect(lmvf)

            # Perform all potential multivector merges

            for m1 in 1:length(mvecs)-1
                mvec1 = mvecs[m1]
                nbhd1 = connected ? Set(lefschetz_neighbors(lc, collect(mvec1))) : nothing
                for m2 in m1+1:length(mvecs)
                    mvec2 = mvecs[m2]

                    # If restricted to connected multivectors, only merge
                    # topologically adjacent blocks

                    if connected && length(intersect(nbhd1, mvec2)) == 0
                        continue
                    end

                    newmvf = setdiff(mvecs, [mvec1, mvec2])
                    push!(newmvf, union(mvec1, mvec2))

                    # We are only interested in acyclic partitions
                    if mvf_is_acyclic(lc, collect.(newmvf))
                        put!(ch, Set(newmvf))
                    end
                end
            end
        end

        # Collect the multivector fields from the channel

        close(ch)
        push!(ap, collect(unique(collect(ch))))
    end

    # Return the space of acyclic partitions

    return convert_partition_mvf.(vcat(ap...))
end
