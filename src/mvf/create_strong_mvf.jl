export create_strong_mvf, create_planar_strong_mvf, create_spatial_strong_mvf

"""
    mvf, dualmvf, info = create_strong_mvf(lc::AbstractComplex,
                                           fwd::Vector{Vector{Int}},
                                           bwd::Vector{Vector{Int}};
                                           exitcells=nothing,
                                           entrancecells=nothing,
                                           merge=Vector{Vector{Int}}())

Create the finest strongly accumulating multivector field for given
forward and backward transitions.

The forward transitions `fwd` and the backward transitions `bwd` have the
format of the output of `planar_transitions` or `spatial_transitions`,
i.e., every entry has the form `[c; targets]`, where `targets` approximates
the immediate future `Ift(c)`, or the immediate past `Ipt(c)`, of the cell
`c`. The backward transitions are usually obtained by applying the same
function to the reversed vector field. Cells without an entry have no
transitions, and the transitions of the toplexes are trivial.

A multivector field `mvf` is *accumulating* if every forward transition is
contained in a single multivector. It is *strongly accumulating* if, in
addition, there is a *dual* multivector field which contains every
backward transition in a single multivector, and whose multivectors
contain the same toplexes as the multivectors of `mvf`. In this case, the
closure of every multivector is a strong isolating block, i.e., both its
exit set and its entry set are closed, and the entry set is given by the
mouth of the dual multivector. The function returns the finest strongly
accumulating multivector field `mvf` together with its dual `dualmvf`,
which are uniquely determined. Both are given in the format of
`create_mvf_hull`, i.e., only multivectors with at least two cells are
listed. The dictionary `info` contains the following entries:

* `"exit cells"`, `"entrance cells"`: The boundary cells, see below.
* `"mvf to dual"`: For every multivector in `mvf`, the index of the
  multivector in `dualmvf` with the same toplexes, or zero if the latter
  is a single toplex.
* `"dual to mvf"`: The analogous index vector for `dualmvf`.
* `"forward merges"`, `"backward merges"`, `"couplings"`, `"time"`:
  Statistics of the computation.

The optional keyword `merge` is a list of sets of cells, each of which has
to lie in a single multivector of `mvf`. The function then returns the
finest strongly accumulating multivector field with this additional
property, together with its dual. The sets must not contain exit cells,
and they are only imposed on `mvf`, i.e., the dual field contains the
toplexes of each set in a single multivector, but not necessarily its
other cells. If every set is connected, then all multivectors of the
result are connected (up to the sampling of the transitions). This can be
used to remove combinatorial artifacts of a mesh, for example a small
cycle in the multivector digraph (see `mvf_digraph`) which does not
reflect the dynamics, by merging the union of its multivectors.

The boundary of the domain is described by the *exit cells*, through
which the flow leaves the domain, and the *entrance cells*, through which
it enters. Cells at which the flow touches the boundary from outside belong
to both families. If `exitcells` and `entrancecells` are not given, they
are inferred from the transitions of the boundary faces of codimension
one: Such a face is an exit cell if only its backward transition reaches
its toplex, an entrance cell if only its forward transition does, and both
if neither does. The two families are then the closures of these faces. A
boundary face whose transitions reach its toplex in both directions
violates the assumptions, and leads to an error. The forward field is
defined on all cells except the exit cells, and the dual field on all cells
except the entrance cells. Thus, the exit cells do not belong to any
multivector of `mvf`, and if there are exit cells, then Morse sets and
Conley indices have to be computed on the complex without them, see for
example `remove_exit_set`. The function raises an error if one of the
families is not closed, if a cell outside the exit cells has no forward
transition, or if a cell outside the entrance cells has no backward
transition. The last two conditions are sampled versions of the
assumption that every cell can be left in forward and in backward time.

The construction, including the keyword `merge`, is Theorem "Finest
strongly accumulating multivector field" of the paper *Global dynamics via multivector fields on Lefschetz
complexes*. It formalizes the merging procedure of Justin Thorpe in his
dissertation *Periodic Orbits in the Setting of Combinatorial Dynamics*
(George Mason University, 2026), which was implemented in pull request #10.
The implementation keeps one union-find structure for each of the two
multivector fields. Whenever two blocks containing toplexes are merged in
one field, the corresponding toplexes are merged in the other field, and
local closedness is restored incrementally, by checking only the closure
and the star of the cells of the smaller block after each merge. Every
merge is forced, so the result is the finest multivector field with the
required properties.
"""
function create_strong_mvf(lc::AbstractComplex, fwd::Vector{Vector{Int}},
                           bwd::Vector{Vector{Int}};
                           exitcells=nothing, entrancecells=nothing,
                           merge::Vector{Vector{Int}}=Vector{Vector{Int}}())
    #
    # Create the finest strongly accumulating multivector field
    #

    tstart = time()
    ncells = lc.ncells
    dims = lc.dimensions
    faces = lc.boundary.columns
    cofaces = lc.boundary.rows

    # Check homogeneity: every cell without cofaces is a toplex

    for c in 1:ncells
        if isempty(cofaces[c]) && (dims[c] != lc.dim)
            throw(ArgumentError("create_strong_mvf: the complex is not homogeneous, " *
                                "cell $c of dimension $(dims[c]) has no cofaces"))
        end
    end

    # Determine the boundary cells

    if (exitcells === nothing) && (entrancecells === nothing)
        exitcells, entrancecells = strong_mvf_boundary_cells(lc, fwd, bwd)
    elseif (exitcells === nothing) || (entrancecells === nothing)
        throw(ArgumentError("create_strong_mvf: specify both exitcells and " *
                            "entrancecells, or neither"))
    else
        exitcells = sort(unique(Vector{Int}(exitcells)))
        entrancecells = sort(unique(Vector{Int}(entrancecells)))
        lefschetz_is_closed(lc, exitcells) ||
            throw(ArgumentError("create_strong_mvf: the exit cells are not closed"))
        lefschetz_is_closed(lc, entrancecells) ||
            throw(ArgumentError("create_strong_mvf: the entrance cells are not closed"))
    end
    isexit = falses(ncells)
    isexit[exitcells] .= true
    isentr = falses(ncells)
    isentr[entrancecells] .= true
    if any(c -> (dims[c] == lc.dim) && (isexit[c] || isentr[c]), 1:ncells)
        throw(ArgumentError("create_strong_mvf: toplexes cannot be boundary cells"))
    end

    # The prescribed merge sets must not contain exit cells

    for S in merge
        all(c -> 1 <= c <= ncells, S) ||
            throw(ArgumentError("create_strong_mvf: a merge set contains an invalid cell index"))
        if any(c -> isexit[c], S)
            throw(ArgumentError("create_strong_mvf: the merge set starting with cell " *
                                "$(S[1]) contains exit cells"))
        end
    end

    # Check that every cell can be left in forward and in backward time

    hasfwd = falses(ncells)
    hasbwd = falses(ncells)
    for t in fwd
        length(t) > 1 && (hasfwd[t[1]] = true)
    end
    for t in bwd
        length(t) > 1 && (hasbwd[t[1]] = true)
    end
    nofwd = [c for c in 1:ncells if (dims[c] < lc.dim) && !isexit[c] && !hasfwd[c]]
    nobwd = [c for c in 1:ncells if (dims[c] < lc.dim) && !isentr[c] && !hasbwd[c]]
    if !isempty(nofwd)
        throw(ArgumentError("create_strong_mvf: $(length(nofwd)) cells which are " *
                            "not exit cells have no forward transition, for example " *
                            "the cells $(nofwd[1:min(5,end)])"))
    end
    if !isempty(nobwd)
        throw(ArgumentError("create_strong_mvf: $(length(nobwd)) cells which are " *
                            "not entrance cells have no backward transition, for example " *
                            "the cells $(nobwd[1:min(5,end)])"))
    end

    # Union-find structures for the forward (1) and backward (2) fields,
    # with member lists and a representative toplex for every block

    parent  = (collect(1:ncells), collect(1:ncells))
    members = (Vector{Vector{Int}}(undef, ncells), Vector{Vector{Int}}(undef, ncells))
    hasmem  = (falses(ncells), falses(ncells))
    toprep  = ([dims[c] == lc.dim ? c : 0 for c in 1:ncells],
               [dims[c] == lc.dim ? c : 0 for c in 1:ncells])
    pending = (Vector{Int}(), Vector{Int}())
    queue   = Vector{Tuple{Int,Int,Int}}()
    stats   = zeros(Int, 3)

    function findroot(f::Int, x::Int)
        p = parent[f]
        while p[x] != x
            p[x] = p[p[x]]
            x = p[x]
        end
        return x
    end

    function getmembers(f::Int, r::Int)
        if !hasmem[f][r]
            members[f][r] = [r]
            hasmem[f][r] = true
        end
        return members[f][r]
    end

    # Merge two blocks. The cells of the smaller block are checked for
    # local closedness later, and toplexes are coupled to the other field

    function unite!(f::Int, a::Int, b::Int)
        ra = findroot(f, a)
        rb = findroot(f, b)
        ra == rb && return
        ma = getmembers(f, ra)
        mb = getmembers(f, rb)
        if length(ma) < length(mb)
            ra, rb = rb, ra
            ma, mb = mb, ma
        end
        parent[f][rb] = ra
        append!(ma, mb)
        append!(pending[f], mb)
        members[f][rb] = Vector{Int}()
        stats[f] += 1
        ta = toprep[f][ra]
        tb = toprep[f][rb]
        if (ta != 0) && (tb != 0)
            push!(queue, (3 - f, ta, tb))
            stats[3] += 1
        elseif ta == 0
            toprep[f][ra] = tb
        end
        return
    end

    # Restore local closedness at a cell c of a block B: Add every cell y
    # outside of B with c <= y <= z for some z in B, and every cell y
    # outside of B with x <= y <= c for some x in B. The stars and closures
    # are traversed in order of dimension, so that it suffices to look at
    # the immediate cofaces and faces.

    mark  = zeros(Int, ncells)
    above = zeros(Int, ncells)
    stamp = Ref(0)

    function localcells(c::Int, nbrs)
        stamp[] += 1
        s = stamp[]
        cells = [c]
        mark[c] = s
        i = 1
        while i <= length(cells)
            for w in nbrs[cells[i]]
                if mark[w] != s
                    mark[w] = s
                    push!(cells, w)
                end
            end
            i += 1
        end
        return cells, s
    end

    function closeat!(f::Int, c::Int)
        for (nbrs, upward) in ((cofaces, true), (faces, false))
            cells, s = localcells(c, nbrs)
            sort!(cells, by = y -> dims[y], rev = upward)
            r = findroot(f, c)
            toadd = Vector{Int}()
            for y in cells
                inblock = (findroot(f, y) == r)
                hit = inblock
                if !hit
                    for w in nbrs[y]
                        if above[w] == s
                            hit = true
                            break
                        end
                    end
                end
                if hit
                    above[y] = s
                    inblock || push!(toadd, y)
                end
            end
            for y in toadd
                unite!(f, c, y)
            end
        end
        return
    end

    # Enter the transitions, without the exit cells in the forward
    # field and without the entrance cells in the backward field

    for (f, transitions, excluded, name) in ((1, fwd, isexit, "forward"),
                                             (2, bwd, isentr, "backward"))
        for t in transitions
            excluded[t[1]] && continue
            for m in 2:length(t)
                if excluded[t[m]]
                    throw(ArgumentError("create_strong_mvf: the $name transition " *
                                        "of cell $(t[1]) reaches the boundary cell $(t[m])"))
                end
                push!(queue, (f, t[1], t[m]))
            end
        end
    end

    # Enter the prescribed merge sets into the forward field

    for S in merge
        for m in 2:length(S)
            push!(queue, (1, S[1], S[m]))
        end
    end

    # Process merges and local closedness checks until nothing changes

    while true
        if !isempty(queue)
            f, a, b = pop!(queue)
            unite!(f, a, b)
        elseif !isempty(pending[1])
            closeat!(1, pop!(pending[1]))
        elseif !isempty(pending[2])
            closeat!(2, pop!(pending[2]))
        else
            break
        end
    end

    # Collect the multivectors of both fields

    function collectblocks(f::Int, excluded)
        blocks = Vector{Vector{Int}}()
        for c in 1:ncells
            if !excluded[c] && (findroot(f, c) == c) && hasmem[f][c] &&
               (length(members[f][c]) > 1)
                push!(blocks, sort(members[f][c]))
            end
        end
        sort!(blocks, by=first)
        return blocks
    end

    mvf = collectblocks(1, isexit)
    dualmvf = collectblocks(2, isentr)

    # Match the multivectors via their toplexes

    dualindex = zeros(Int, ncells)
    for k in 1:length(dualmvf)
        dualindex[findroot(2, dualmvf[k][1])] = k
    end
    mvfindex = zeros(Int, ncells)
    for k in 1:length(mvf)
        mvfindex[findroot(1, mvf[k][1])] = k
    end
    mvf2dual = zeros(Int, length(mvf))
    for k in 1:length(mvf)
        t = toprep[1][findroot(1, mvf[k][1])]
        t == 0 || (mvf2dual[k] = dualindex[findroot(2, t)])
    end
    dual2mvf = zeros(Int, length(dualmvf))
    for k in 1:length(dualmvf)
        t = toprep[2][findroot(2, dualmvf[k][1])]
        t == 0 || (dual2mvf[k] = mvfindex[findroot(1, t)])
    end

    # Return the multivector fields and the information

    info = Dict{String,Any}()
    info["exit cells"] = exitcells
    info["entrance cells"] = entrancecells
    info["mvf to dual"] = mvf2dual
    info["dual to mvf"] = dual2mvf
    info["forward merges"] = stats[1]
    info["backward merges"] = stats[2]
    info["couplings"] = stats[3]
    info["time"] = time() - tstart

    return mvf, dualmvf, info
end

"""
    mvf, dualmvf, info = create_planar_strong_mvf(lc::AbstractComplex,
                             coords::Vector{<:Vector{<:Real}}, vf;
                             npts::Int=10, nontransverse::Bool=false,
                             exitcells=nothing, entrancecells=nothing,
                             merge=Vector{Vector{Int}}())

Create the finest strongly accumulating multivector field for a planar
vector field.

The arguments `lc`, `coords`, and `vf` are as in `create_planar_mvf`. The
forward and backward transitions are computed by `planar_transitions` for
the vector field `vf` and for the reversed field `x -> -vf(x)`, with the
same keywords `npts` and `nontransverse`, and they are passed to
`create_strong_mvf`, together with the optional boundary cells and merge sets. See
`create_strong_mvf` for a description of the result.
"""
function create_planar_strong_mvf(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}}, vf;
                                  npts::Int=10, nontransverse::Bool=false,
                                  exitcells=nothing, entrancecells=nothing,
                                  merge::Vector{Vector{Int}}=Vector{Vector{Int}}())
    fwd = planar_transitions(lc, coords, vf, npts=npts, nontransverse=nontransverse)
    bwd = planar_transitions(lc, coords, x -> -vf(x), npts=npts,
                             nontransverse=nontransverse)
    return create_strong_mvf(lc, fwd, bwd, exitcells=exitcells,
                             entrancecells=entrancecells, merge=merge)
end

"""
    mvf, dualmvf, info = create_planar_strong_mvf(ec::EuclideanComplex, vf;
                             npts::Int=10, nontransverse::Bool=false,
                             exitcells=nothing, entrancecells=nothing,
                             merge=Vector{Vector{Int}}())

Create the finest strongly accumulating multivector field for a planar
vector field on a `EuclideanComplex`, see the method for Lefschetz
complexes.
"""
function create_planar_strong_mvf(ec::EuclideanComplex, vf; npts::Int=10,
                                  nontransverse::Bool=false,
                                  exitcells=nothing, entrancecells=nothing,
                                  merge::Vector{Vector{Int}}=Vector{Vector{Int}}())
    _check_coords2d(ec, "create_planar_strong_mvf",
                    "For complexes in three dimensions, use create_spatial_strong_mvf.")
    lc = euclidean_to_lefschetz(ec)
    return create_planar_strong_mvf(lc, euclidean_vertex_coords(ec), vf, npts=npts,
                                    nontransverse=nontransverse, exitcells=exitcells,
                                    entrancecells=entrancecells, merge=merge)
end

"""
    mvf, dualmvf, info = create_spatial_strong_mvf(lc::AbstractComplex,
                             coords::Vector{<:Vector{<:Real}}, vf;
                             npts::Int=3, exitcells=nothing,
                             entrancecells=nothing, merge=Vector{Vector{Int}}())

Create the finest strongly accumulating multivector field for a vector
field in three dimensions.

The arguments `lc`, `coords`, and `vf` are as in `create_spatial_mvf`. The
forward and backward transitions are computed by `spatial_transitions` for
the vector field `vf` and for the reversed field `x -> -vf(x)`, with the
same keyword `npts`, and they are passed to `create_strong_mvf`, together
with the optional boundary cells. See `create_strong_mvf` for a
description of the result.
"""
function create_spatial_strong_mvf(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}}, vf;
                                   npts::Int=3, exitcells=nothing, entrancecells=nothing,
                                   merge::Vector{Vector{Int}}=Vector{Vector{Int}}())
    fwd = spatial_transitions(lc, coords, vf, npts=npts)
    bwd = spatial_transitions(lc, coords, x -> -vf(x), npts=npts)
    return create_strong_mvf(lc, fwd, bwd, exitcells=exitcells,
                             entrancecells=entrancecells, merge=merge)
end

"""
    mvf, dualmvf, info = create_spatial_strong_mvf(ec::EuclideanComplex, vf;
                             npts::Int=3, exitcells=nothing,
                             entrancecells=nothing, merge=Vector{Vector{Int}}())

Create the finest strongly accumulating multivector field for a vector
field in three dimensions on a `EuclideanComplex`, see the method for
Lefschetz complexes.
"""
function create_spatial_strong_mvf(ec::EuclideanComplex, vf; npts::Int=3,
                                   exitcells=nothing, entrancecells=nothing,
                                   merge::Vector{Vector{Int}}=Vector{Vector{Int}}())
    _check_spatial_complex(ec, "create_spatial_strong_mvf")
    lc = euclidean_to_lefschetz(ec)
    return create_spatial_strong_mvf(lc, euclidean_vertex_coords(ec), vf, npts=npts,
                                     exitcells=exitcells, entrancecells=entrancecells,
                                     merge=merge)
end

###############################################################################
#
# Internal functions
#
###############################################################################

function strong_mvf_boundary_cells(lc::AbstractComplex, fwd::Vector{Vector{Int}},
                                   bwd::Vector{Vector{Int}})
    #
    # Infer the exit and entrance cells from the transitions of the
    # boundary faces of codimension one, i.e., of the faces of exactly
    # one toplex. Returns the closures of the exit and entrance faces,
    # where faces with no transitions count as both (external tangency).
    #
    ncells = lc.ncells
    cofaces = lc.boundary.rows
    reachfwd = falses(ncells)
    reachbwd = falses(ncells)
    for (transitions, reach) in ((fwd, reachfwd), (bwd, reachbwd))
        for t in transitions
            c = t[1]
            if (lc.dimensions[c] == lc.dim - 1) && (length(cofaces[c]) == 1)
                if cofaces[c][1] in view(t, 2:length(t))
                    reach[c] = true
                end
            end
        end
    end

    exitfaces = Vector{Int}()
    entrfaces = Vector{Int}()
    mixed = Vector{Int}()
    for c in 1:ncells
        if (lc.dimensions[c] == lc.dim - 1) && (length(cofaces[c]) == 1)
            if reachfwd[c] && reachbwd[c]
                push!(mixed, c)
            elseif reachfwd[c]
                push!(entrfaces, c)
            elseif reachbwd[c]
                push!(exitfaces, c)
            else
                push!(exitfaces, c)
                push!(entrfaces, c)
            end
        end
    end

    if !isempty(mixed)
        throw(ArgumentError("create_strong_mvf: the flow enters and leaves the " *
                            "domain through $(length(mixed)) boundary faces, for " *
                            "example the cells $(mixed[1:min(5,end)]). Refine the " *
                            "mesh near these faces, or specify exitcells and " *
                            "entrancecells"))
    end

    exitcells = isempty(exitfaces) ? Vector{Int}() : lefschetz_closure(lc, exitfaces)
    entrancecells = isempty(entrfaces) ? Vector{Int}() : lefschetz_closure(lc, entrfaces)
    return sort(exitcells), sort(entrancecells)
end

function strong_mvf_reference(lc::AbstractComplex, fwd::Vector{Vector{Int}},
                              bwd::Vector{Vector{Int}};
                              exitcells=nothing, entrancecells=nothing,
                              merge::Vector{Vector{Int}}=Vector{Vector{Int}}())
    #
    # Reference implementation of the finest strongly accumulating
    # multivector field, which follows the iteration in the proof of
    # the theorem literally: Starting with the partition of the toplexes
    # into singletons, compute the minimal forward and backward fields
    # containing the transitions and the current partition, and replace
    # the partition by the one generated by the toplexes of both fields,
    # until it is stable. Returns the two fields and the number of rounds.
    #
    if (exitcells === nothing) && (entrancecells === nothing)
        exitcells, entrancecells = strong_mvf_boundary_cells(lc, fwd, bwd)
    end
    isexit = falses(lc.ncells)
    isexit[exitcells] .= true
    isentr = falses(lc.ncells)
    isentr[entrancecells] .= true
    fwdX = vcat([t for t in fwd if !isexit[t[1]]], merge)
    bwdX = [t for t in bwd if !isentr[t[1]]]
    toplexes = findall(isequal(lc.dim), lc.dimensions)

    partition = Vector{Vector{Int}}()
    rounds = 0
    while true
        mvf = create_mvf_hull(lc, vcat(fwdX, partition))
        dualmvf = create_mvf_hull(lc, vcat(bwdX, partition))

        # The partition of the toplexes generated by both fields

        label = collect(1:lc.ncells)
        function root(x)
            while label[x] != x
                x = label[x]
            end
            return x
        end
        for mv in vcat(mvf, dualmvf)
            tops = filter(c -> lc.dimensions[c] == lc.dim, mv)
            for t in tops[2:end]
                r1, r2 = root(tops[1]), root(t)
                r1 == r2 || (label[r2] = r1)
            end
        end
        groups = Dict{Int,Vector{Int}}()
        for t in toplexes
            push!(get!(groups, root(t), Vector{Int}()), t)
        end
        newpartition = sort([sort(g) for g in values(groups) if length(g) > 1], by=first)

        if newpartition == partition
            return mvf, dualmvf, rounds
        end
        partition = newpartition
        rounds += 1
    end
end
