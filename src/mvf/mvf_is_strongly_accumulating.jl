export mvf_is_strongly_accumulating, strong_isoblock_sets

"""
    mvf_is_strongly_accumulating(lc::AbstractComplex, mvf::Vector{Vector{Int}},
                                 dualmvf::Vector{Vector{Int}},
                                 fwd::Vector{Vector{Int}}, bwd::Vector{Vector{Int}};
                                 exitcells=nothing, entrancecells=nothing,
                                 merge=Vector{Vector{Int}}(), details::Bool=false)

Check whether a multivector field is strongly accumulating with dual
`dualmvf`, with respect to the forward transitions `fwd` and the backward
transitions `bwd`.

The arguments have the format of the input and output of
`create_strong_mvf`, and the boundary cells are inferred from the
transitions in the same way, unless both `exitcells` and `entrancecells`
are given. Cells which are not listed in `mvf` or `dualmvf` are
singletons. The following conditions are checked:

* The sets in `mvf` are disjoint, locally closed, and contain no exit
  cells, and the sets in `dualmvf` are disjoint, locally closed, and
  contain no entrance cells.
* Every forward transition of a cell which is not an exit cell lies in a
  single multivector of `mvf`, every backward transition of a cell which
  is not an entrance cell lies in a single multivector of `dualmvf`, and
  these cells all have transitions.
* Both multivector fields induce the same partition of the toplexes.
* Every dual multivector `V*` consists of the cells whose backward
  transitions meet the corresponding multivector `V`, and `V` consists of
  the cells whose forward transitions meet `V*`. These sets are computed
  independently with `mvf_dual`.
* The closures of `V` and `V*` agree.
* Every set in `merge` contains no exit cells and lies in a single
  multivector of `mvf`, see `create_strong_mvf`.

The function returns `true` if all conditions hold. If `details=true`, it
returns a tuple `(ok, report)`, where `report` is a dictionary which lists
the violations of every condition. Every entry is a vector, which is empty
if and only if the corresponding condition is satisfied, and `ok` is `true`
if and only if all entries except the last two are empty. Cells and
multivectors are given by their indices, where multivector indices refer
to the positions in `mvf` or `dualmvf`. The dictionary has the following
keys:

* `"overlapping cells"`, `"dual overlapping cells"`: The cells which lie
  in more than one multivector of `mvf` or `dualmvf`, respectively.
* `"not locally closed"`, `"dual not locally closed"`: The multivectors of
  `mvf` or `dualmvf`, respectively, which are not locally closed.
* `"exit cells in mvf"`: The exit cells which lie in a multivector of
  `mvf`.
* `"entrance cells in dual"`: The entrance cells which lie in a
  multivector of `dualmvf`.
* `"split forward transitions"`: The cells which are not exit cells, and
  whose forward transition does not lie in a single multivector of `mvf`.
  This includes cells which are singletons of `mvf`, but which have a
  nontrivial forward transition.
* `"split backward transitions"`: The analogous cells for the backward
  transitions, the entrance cells, and `dualmvf`.
* `"missing forward transitions"`: The cells which are neither toplexes
  nor exit cells, and which have no forward transition.
* `"missing backward transitions"`: The cells which are neither toplexes
  nor entrance cells, and which have no backward transition.
* `"split merge sets"`: The indices of the sets in `merge` which contain
  an exit cell, or which do not lie in a single multivector of `mvf`.
* `"toplexes with different partitions"`: The toplexes at which the
  partitions of the toplexes induced by `mvf` and `dualmvf` differ.
* `"duals not given by immediate pasts"`: The multivectors `V` for which
  the dual computed with `mvf_dual` differs from the multivector of
  `dualmvf` which contains the same toplexes. Here the indices refer to
  the list which consists of `mvf`, followed by the singletons of `mvf`
  which are toplexes, in increasing order.
* `"multivectors not given by immediate futures"`: The analogous check in
  the reverse direction. The indices refer to the list which consists of
  `dualmvf`, followed by the singletons of `dualmvf` which are toplexes,
  in increasing order. If the multivectors overlap, then these two
  entries are not computed, and they contain the string `"not checked,
  since the multivectors overlap"` instead.
* `"closures differ"`: The multivectors of `mvf` which contain a toplex,
  and whose closure differs from the closure of the corresponding
  multivector of `dualmvf`.
* `"disconnected"`, `"dual disconnected"`: The multivectors of `mvf` or
  `dualmvf`, respectively, which are not connected. The theory guarantees
  connectedness only if the transitions over-approximate the immediate
  futures and pasts, so a disconnected multivector indicates that the
  sampling missed some transitions. This diagnostic does not affect the
  result `ok`.
"""
function mvf_is_strongly_accumulating(lc::AbstractComplex, mvf::Vector{Vector{Int}},
                                      dualmvf::Vector{Vector{Int}},
                                      fwd::Vector{Vector{Int}}, bwd::Vector{Vector{Int}};
                                      exitcells=nothing, entrancecells=nothing,
                                      merge::Vector{Vector{Int}}=Vector{Vector{Int}}(),
                                      details::Bool=false)
    #
    # Check whether a multivector field is strongly accumulating
    #

    ncells = lc.ncells
    istop(c) = lc.dimensions[c] == lc.dim
    report = Dict{String,Any}()

    # Determine the boundary cells

    if (exitcells === nothing) && (entrancecells === nothing)
        exitcells, entrancecells = strong_mvf_boundary_cells(lc, fwd, bwd)
    elseif (exitcells === nothing) || (entrancecells === nothing)
        throw(ArgumentError("mvf_is_strongly_accumulating: specify both exitcells " *
                            "and entrancecells, or neither"))
    end
    isexit = falses(ncells)
    isexit[exitcells] .= true
    isentr = falses(ncells)
    isentr[entrancecells] .= true

    # The multivectors: disjointness, local closedness, boundary cells

    function celltable(sets)
        table = zeros(Int, ncells)
        overlaps = Vector{Int}()
        for k in 1:length(sets)
            for c in sets[k]
                table[c] == 0 ? (table[c] = k) : push!(overlaps, c)
            end
        end
        return table, sort(unique(overlaps))
    end
    c2m, report["overlapping cells"] = celltable(mvf)
    c2d, report["dual overlapping cells"] = celltable(dualmvf)
    report["not locally closed"] = [k for k in 1:length(mvf)
                                    if !lefschetz_is_locally_closed(lc, mvf[k])]
    report["dual not locally closed"] = [k for k in 1:length(dualmvf)
                                         if !lefschetz_is_locally_closed(lc, dualmvf[k])]
    report["exit cells in mvf"] = [c for c in 1:ncells if isexit[c] && (c2m[c] != 0)]
    report["entrance cells in dual"] = [c for c in 1:ncells if isentr[c] && (c2d[c] != 0)]

    # Every transition lies in a single multivector

    function splittransitions(transitions, table, excluded)
        split = Vector{Int}()
        hastransition = falses(ncells)
        for t in transitions
            excluded[t[1]] && continue
            hastransition[t[1]] = true
            k = table[t[1]]
            if (k == 0) || any(c -> table[c] != k, t)
                push!(split, t[1])
            end
        end
        missing = [c for c in 1:ncells if !istop(c) && !excluded[c] && !hastransition[c]]
        return sort(unique(split)), missing
    end
    report["split forward transitions"], report["missing forward transitions"] =
        splittransitions(fwd, c2m, isexit)
    report["split backward transitions"], report["missing backward transitions"] =
        splittransitions(bwd, c2d, isentr)
    report["split merge sets"] = [k for k in 1:length(merge)
                                  if any(c -> isexit[c], merge[k]) ||
                                     (length(merge[k]) > 1 &&
                                      (c2m[merge[k][1]] == 0 ||
                                       any(c -> c2m[c] != c2m[merge[k][1]], merge[k])))]

    # Both fields induce the same partition of the toplexes. Singletons
    # get the negative toplex as their label.

    toplexes = [c for c in 1:ncells if istop(c)]
    flabel(t) = c2m[t] != 0 ? c2m[t] : -t
    dlabel(t) = c2d[t] != 0 ? c2d[t] : -t
    f2d = Dict{Int,Int}()
    d2f = Dict{Int,Int}()
    mismatched = Vector{Int}()
    for t in toplexes
        a, b = flabel(t), dlabel(t)
        if (get!(f2d, a, b) != b) || (get!(d2f, b, a) != a)
            push!(mismatched, t)
        end
    end
    report["toplexes with different partitions"] = mismatched

    # The duals are given by the immediate pasts, and vice versa. All
    # multivectors containing toplexes are compared, including singletons.

    allmvf = vcat(mvf, [[t] for t in toplexes if c2m[t] == 0])
    alldual = vcat(dualmvf, [[t] for t in toplexes if c2d[t] == 0])
    function blocks_with_toplex(sets, table)
        result = Dict{Int,Vector{Int}}()
        for t in toplexes
            k = table[t]
            result[t] = k != 0 ? sort(sets[k]) : [t]
        end
        return result
    end
    mvfoftop = blocks_with_toplex(mvf, c2m)
    dualoftop = blocks_with_toplex(dualmvf, c2d)
    fwdX = [t for t in fwd if !isexit[t[1]]]
    bwdX = [t for t in bwd if !isentr[t[1]]]
    if isempty(report["overlapping cells"]) && isempty(report["dual overlapping cells"])
        computeddual = mvf_dual(lc, allmvf, bwdX)
        computedmvf = mvf_dual(lc, alldual, fwdX)
        report["duals not given by immediate pasts"] =
            [k for k in 1:length(allmvf)
             if computeddual[k] != dualoftop[first(c for c in allmvf[k] if istop(c))]]
        report["multivectors not given by immediate futures"] =
            [k for k in 1:length(alldual)
             if computedmvf[k] != mvfoftop[first(c for c in alldual[k] if istop(c))]]
    else
        notchecked = "not checked, since the multivectors overlap"
        report["duals not given by immediate pasts"] = notchecked
        report["multivectors not given by immediate futures"] = notchecked
    end

    # The closures of corresponding multivectors agree

    report["closures differ"] =
        [k for k in 1:length(mvf)
         if any(istop, mvf[k]) &&
            sort(lefschetz_closure(lc, mvf[k])) !=
            sort(lefschetz_closure(lc, dualoftop[first(c for c in mvf[k] if istop(c))]))]

    # Collect the result

    ok = all(isempty, values(report))

    # Diagnostic: connectedness of the multivectors

    report["disconnected"] = [k for k in 1:length(mvf) if !is_connected_block(lc, mvf[k])]
    report["dual disconnected"] = [k for k in 1:length(dualmvf)
                                   if !is_connected_block(lc, dualmvf[k])]

    return details ? (ok, report) : ok
end

"""
    mouth, dualset, dualmouth, islc = strong_isoblock_sets(lc::AbstractComplex,
                                          U::Vector{Int}, mvf::Vector{Vector{Int}},
                                          dualmvf::Vector{Vector{Int}})

Determine the sets which describe the isolating block given by a union of
multivectors.

Let `mvf` be a strongly accumulating multivector field with dual
`dualmvf`, as returned by `create_strong_mvf`, and let `U` be a locally
closed union of multivectors of `mvf`, where singletons are allowed. Then
the closure of `U` is an isolating block whose exit set is given by the
mouth of `U`, i.e., by `cl(U) ∖ U`. The *dual* of `U` is the union of the
duals of the multivectors in `U`. It has the same closure as `U`. If
the dual is locally closed, then the isolating block is a *strong*
isolating block, whose entry set is given by the mouth of the dual.

The function returns the mouth of `U`, the dual `dualset` of `U`, its mouth
`dualmouth`, and a boolean `islc` which indicates whether the dual is
locally closed. All sets are sorted. The function raises an error if `U` is
not a locally closed union of multivectors. For a single multivector, the
dual is always locally closed.
"""
function strong_isoblock_sets(lc::AbstractComplex, U::Vector{Int},
                              mvf::Vector{Vector{Int}}, dualmvf::Vector{Vector{Int}})
    #
    # Determine the mouth and the dual of a union of multivectors
    #

    ncells = lc.ncells
    c2m = zeros(Int, ncells)
    for k in 1:length(mvf)
        for c in mvf[k]
            c2m[c] = k
        end
    end
    c2d = zeros(Int, ncells)
    for k in 1:length(dualmvf)
        for c in dualmvf[k]
            c2d[c] = k
        end
    end

    # Check that U is a locally closed union of multivectors

    inU = falses(ncells)
    inU[U] .= true
    for c in U
        if (c2m[c] != 0) && !all(inU[mvf[c2m[c]]])
            throw(ArgumentError("strong_isoblock_sets: the set is not a union " *
                                "of multivectors"))
        end
    end
    lefschetz_is_locally_closed(lc, U) ||
        throw(ArgumentError("strong_isoblock_sets: the set is not locally closed"))

    # The dual is the union of the duals of the multivectors, which are
    # determined by their toplexes

    indual = falses(ncells)
    for c in U
        if lc.dimensions[c] == lc.dim
            if c2d[c] != 0
                indual[dualmvf[c2d[c]]] .= true
            else
                indual[c] = true
            end
        end
    end
    dualset = findall(indual)

    # The mouths and the local closedness of the dual

    Usorted = sort(unique(U))
    mouth = sort(setdiff(lefschetz_closure(lc, Usorted), Usorted))
    dualmouth = isempty(dualset) ? Vector{Int}() :
                sort(setdiff(lefschetz_closure(lc, dualset), dualset))
    islc = isempty(dualset) || lefschetz_is_locally_closed(lc, dualset)

    return mouth, dualset, dualmouth, islc
end
