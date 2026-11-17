export mvf_dual_set, mvf_dual

"""
    mvf_dual_set(lc::AbstractComplex, mv::Vector{Int}, bwd::Vector{Vector{Int}})

Compute the dual of a multivector with respect to given backward transitions.

The backward transitions `bwd` have the format of the output of
`planar_transitions` or `spatial_transitions`, applied to the reversed
vector field, i.e., every entry has the form `[c; targets]`, where
`targets` approximates the immediate past `Ipt(c)` of the cell `c`. Cells
without an entry have an empty immediate past, which is the case for the
entrance cells through which the flow enters the domain. The immediate
past of a toplex is the toplex itself.

The function returns the sorted vector of all cells `c` whose immediate
past intersects `mv`, i.e., the set `mv* = { c : Ipt(c) ∩ mv ≠ ∅ }`. If
`mv` belongs to a strongly accumulating multivector field, then `mv*` is
the multivector of the dual multivector field with the same toplexes, and
the mouth of `mv*` describes the entry set of the strong isolating block
given by the closure of `mv`.

Since the transitions list only toplexes as targets, the function in fact
returns the cells whose immediate past contains a toplex of `mv`. If `mv`
is a multivector of an accumulating multivector field, this is the same
set: If a cell `ρ` of `mv` which is not a toplex lies in `Ipt(c)`, then
`mv` contains a toplex of `Ift(ρ)`, which lies in the star of `ρ`, and
therefore also in the open set `Ipt(c)`. For other sets `mv`, however,
cells whose immediate past meets `mv` only in cells of lower dimension
are not included.

The dual was introduced by Justin Thorpe in his dissertation *Periodic
Orbits in the Setting of Combinatorial Dynamics* (George Mason University,
2026) as the union of the interior of `mv` with the cells on its boundary
whose immediate past meets `mv`, and it was implemented in pull request
#10 as `create_dual_set`. For cells with nonempty immediate past, both
descriptions agree. The present version excludes the entrance cells, and it
uses precomputed backward transitions.
"""
function mvf_dual_set(lc::AbstractComplex, mv::Vector{Int}, bwd::Vector{Vector{Int}})
    #
    # Compute the dual of a single multivector
    #
    return mvf_dual(lc, [mv], bwd)[1]
end

"""
    mvf_dual(lc::AbstractComplex, mvf::Vector{Vector{Int}}, bwd::Vector{Vector{Int}})

Compute the duals of all multivectors in `mvf` with respect to the backward
transitions `bwd`.

The result is a vector whose `k`-th entry is `mvf_dual_set(lc, mvf[k], bwd)`,
but all duals are computed in a single pass over the transitions. The sets
in `mvf` have to be disjoint. If the multivector field `mvf` is strongly
accumulating with respect to the given transitions, then the duals are
disjoint, and they form its dual multivector field. Conversely, a cell
which lies in two duals shows that the corresponding multivectors have to
be merged in every strongly accumulating multivector field. Toplexes which
do not belong to any set in `mvf` are singletons, and their duals, which
are again singletons, are not computed.
"""
function mvf_dual(lc::AbstractComplex, mvf::Vector{Vector{Int}}, bwd::Vector{Vector{Int}})
    #
    # Compute the duals of all multivectors in a single pass
    #

    # Find the multivector containing every cell

    cell2mv = zeros(Int, lc.ncells)
    for k in 1:length(mvf)
        for c in mvf[k]
            if cell2mv[c] != 0
                throw(ArgumentError("mvf_dual: the multivectors have to be disjoint"))
            end
            cell2mv[c] = k
        end
    end

    # Every toplex is its own immediate past

    duals = [Vector{Int}() for _ in 1:length(mvf)]
    for k in 1:length(mvf)
        for c in mvf[k]
            if lc.dimensions[c] == lc.dim
                push!(duals[k], c)
            end
        end
    end

    # Add every cell to the duals of the multivectors which
    # its immediate past meets

    for t in bwd
        c = t[1]
        for m in 2:length(t)
            k = cell2mv[t[m]]
            if (k != 0) && (isempty(duals[k]) || (last(duals[k]) != c))
                push!(duals[k], c)
            end
        end
    end

    # Return the sorted duals

    for k in 1:length(duals)
        unique!(sort!(duals[k]))
    end
    return duals
end

###############################################################################
#
# The dual merging procedure of pull request #10, as a cross-check
#
###############################################################################

"""
    strong_mvf_dual_merging(lc::AbstractComplex, fwd::Vector{Vector{Int}},
                            bwd::Vector{Vector{Int}}; maxrounds::Int=1000,
                            exitcells=nothing, entrancecells=nothing)

Merge multivectors until all duals are locally closed and disjoint.

This is the merging procedure of Justin Thorpe from pull request #10 (the
former function `create_mvf(...; strong_iso=true)`), formulated for
precomputed forward and backward transitions `fwd` and `bwd`, and with
three changes. Multivectors whose duals overlap are merged, the duals do
not contain entrance cells, and a warning is issued if the procedure has
not converged after `maxrounds` rounds. The exit and entrance cells are
determined as in `create_strong_mvf`, unless they are given, and the
forward transitions of exit cells and the backward transitions of entrance
cells are ignored.

Every round forms the smallest multivector field containing the forward
transitions and all merges so far, and computes the duals of its
multivectors with `mvf_dual`. If a cell lies in several duals, then the
corresponding multivectors are merged. If a dual `D` is not locally closed,
then every cell of the remainder `lch(D) ∖ D` is joined to the toplexes in
its star which are not in `D`. For planar complexes, every merge of this
kind is forced, i.e., it occurs in every strongly accumulating multivector
field. Thus, if the procedure converges, its result is the finest strongly
accumulating multivector field for the given transitions. In three
dimensions, the second rule may merge more than necessary.

The function returns the multivector field in the format of
`create_mvf_hull`. It is not exported, and is meant as an independent
check of the planar case.
"""
function strong_mvf_dual_merging(lc::AbstractComplex, fwd::Vector{Vector{Int}},
                                 bwd::Vector{Vector{Int}}; maxrounds::Int=1000,
                                 exitcells=nothing, entrancecells=nothing)
    #
    # The dual merging procedure of pull request #10
    #

    # Remove the transitions of the boundary cells

    if (exitcells === nothing) && (entrancecells === nothing)
        exitcells, entrancecells = strong_mvf_boundary_cells(lc, fwd, bwd)
    end
    isexit = falses(lc.ncells)
    isexit[exitcells] .= true
    isentr = falses(lc.ncells)
    isentr[entrancecells] .= true
    bwd = [t for t in bwd if !isentr[t[1]]]

    base = [t for t in fwd if !isexit[t[1]]]
    toplexes = findall(isequal(lc.dim), lc.dimensions)
    mvf = Vector{Vector{Int}}()
    converged = false
    stalled = false

    for round in 1:maxrounds

        # Form the multivector field, including the singleton toplexes

        mvf = create_mvf_hull(lc, base)
        inmv = falses(lc.ncells)
        for mv in mvf
            inmv[mv] .= true
        end
        allmvf = vcat(mvf, [[t] for t in toplexes if !inmv[t]])
        duals = mvf_dual(lc, allmvf, bwd)

        # Merge multivectors with overlapping duals

        nbase = length(base)
        clean = true
        firstdual = zeros(Int, lc.ncells)
        for k in 1:length(duals)
            for c in duals[k]
                if firstdual[c] == 0
                    firstdual[c] = k
                elseif firstdual[c] != k
                    clean = false
                    t1 = first(c2 for c2 in allmvf[firstdual[c]] if lc.dimensions[c2] == lc.dim)
                    t2 = first(c2 for c2 in allmvf[k] if lc.dimensions[c2] == lc.dim)
                    push!(base, [t1, t2])
                end
            end
        end

        # Join the remainder cells of duals which are not locally closed
        # to the toplexes in their star outside of the dual

        for k in 1:length(duals)
            dual = duals[k]
            remainder = setdiff(lefschetz_lchull(lc, dual), dual)
            isempty(remainder) && continue
            clean = false
            indual = Set(dual)
            for c in remainder
                tops = [t for t in lefschetz_openhull(lc, [c])
                        if (lc.dimensions[t] == lc.dim) && !(t in indual)]
                if !isempty(tops)
                    push!(base, vcat([c], tops))
                end
            end
        end

        # Stop if all duals are disjoint and locally closed, or if
        # nothing could be merged

        if clean
            converged = true
            break
        elseif length(base) == nbase
            stalled = true
            break
        end
    end

    if stalled
        @warn "strong_mvf_dual_merging: a dual is not locally closed, " *
              "but no merge applies"
    elseif !converged
        @warn "strong_mvf_dual_merging: no convergence after $maxrounds rounds"
    end

    # Return the multivector field

    return mvf
end
