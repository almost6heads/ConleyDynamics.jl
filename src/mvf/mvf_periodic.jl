export mvf_with_singletons, check_cyclic_conditions, check_brick_lyapunov
export cyclic_local_cycles

"""
    fullmvf = mvf_with_singletons(lc::AbstractComplex, mvf::Vector{Vector{Int}})

Add the singletons to a multivector field.

Multivector fields are usually stored without their singletons, see for
example `create_mvf_hull`. This function returns the multivector field
`fullmvf`, which consists of the sets of `mvf`, followed by the singletons
of all cells of `lc` which do not belong to a set of `mvf`. Since the
order of the sets in `mvf` is preserved, the indices of these sets do not
change. The complete multivector field is needed for the multivector
digraph of `mvf_digraph`, and for the functions `check_cyclic_conditions`
and `check_brick_lyapunov`.
"""
function mvf_with_singletons(lc::AbstractComplex, mvf::Vector{Vector{Int}})
    #
    # Add the singletons to a multivector field
    #
    inmvf = falses(lc.ncells)
    for V in mvf
        inmvf[V] .= true
    end
    return vcat(mvf, [[c] for c in 1:lc.ncells if !inmvf[c]])
end

"""
    result = check_cyclic_conditions(lc::AbstractComplex, mvf::Vector{Vector{Int}},
                                     C::Vector{Int}, labels::Vector{Int}, p::Int;
                                     out=nothing, mouths=nothing)

Check the combinatorial hypotheses for the existence of a periodic orbit.

The multivector field `mvf` has to contain all multivectors, including
the singletons, see `mvf_with_singletons`. As in `mvf_digraph`, this is
checked, and an `ArgumentError` is raised if a cell lies in no set or in
more than one set of `mvf`. The vector `C` contains the
indices of the multivectors of a locally closed and `mvf`-compatible set
`U`, for example a Morse set, and `labels[k]` is the number between `1`
and `p` of the subcollection `C_i` which contains the multivector
`mvf[C[k]]`, for example as computed by `cyclic_partition_angular`. The
subcollections are numbered cyclically, i.e., `C_{p+1} = C_1`. The
following hypotheses of the theorem "Periodic orbits via multivector
fields" in the paper *Global dynamics via multivector fields on Lefschetz
complexes* are checked:

* `(C2)`: The closures of the unions `U_i` of consecutive subcollections
  intersect.
* `(C3)`: The closures of non-consecutive unions are disjoint, and for
  `p = 3` the closures of all three unions have no common cell.
* `(C4)`: For every `i`, every multivector `V_0` in `C_i`, and every
  maximal directed path `(V_0, V_1, ...)` in the multivector digraph
  restricted to `C`, there is a `k >= 0` such that either the mouth of
  `V_k` lies in the mouth of `U`, or `V_k` lies in `C_{i+1}` and none of
  `V_1, ..., V_{k-1}` lies in `C_{i+2}`.

Here a path is maximal if it is infinite, or if its last vertex has no
outgoing edge within `C`. Hypothesis `(C1)`, i.e., `p >= 3`, is required.
The function returns a named tuple with the boolean entries `C2`, `C3`,
and `C4`, the vector `fails`, which contains for every subcollection
`C_i` the number of start multivectors for which `(C4)` fails, and the
vector `sizes` of the numbers of multivectors in the subcollections. The
multivector digraph and the mouths of the multivectors can be passed as
`out` and `mouths` if they have been computed with `mvf_digraph` already.

The condition `(C4)` is decided by a backward fixed-point iteration on
the product of the digraph with a flag, which records whether the path
has already visited `C_{i+2}`. Its cost is linear in the size of the
digraph for every `i`.
"""
function check_cyclic_conditions(lc::AbstractComplex, mvf::Vector{Vector{Int}},
                                 C::Vector{Int}, labels::Vector{Int}, p::Int;
                                 out=nothing, mouths=nothing)
    #
    # Check the hypotheses (C2) to (C4)
    #
    p >= 3 || throw(ArgumentError("check_cyclic_conditions: p has to be at least 3"))
    length(labels) == length(C) ||
        throw(ArgumentError("check_cyclic_conditions: labels and C must have the same length"))
    all(l -> 1 <= l <= p, labels) ||
        throw(ArgumentError("check_cyclic_conditions: the labels have to lie between 1 and p"))
    if (out === nothing) || (mouths === nothing)
        out, _, mouths = mvf_digraph(lc, mvf)
    else
        complete_mvf_index(lc, mvf, "check_cyclic_conditions")
    end
    nC = length(C)

    # (C2) and (C3): closures of the unions of the subcollections

    Ucells = reduce(vcat, mvf[C])
    moU = Set(setdiff(lefschetz_closure(lc, Ucells), Ucells))
    cls = [Set(lefschetz_closure(lc, reduce(vcat, mvf[C[labels .== i]]; init=Int[])))
           for i in 1:p]
    C2 = all(!isempty(intersect(cls[i], cls[mod1(i + 1, p)])) for i in 1:p)
    C3 = all(isempty(intersect(cls[i], cls[j])) for i in 1:p, j in 1:p
             if !(mod(j - i, p) in (0, 1, p - 1)))
    if p == 3
        C3 = C3 && isempty(intersect(cls[1], cls[2], cls[3]))
    end

    # (C4): A state (k, f) consists of a multivector C[k] and a flag f,
    # which is 1 if the path has visited C_{i+2} after its start. A state
    # is safe if every maximal path from it reaches a good state. Safe
    # states are found by a backward fixed-point iteration.

    pos = Dict{Int,Int}()
    for k in 1:nC
        pos[C[k]] = k
    end
    moinU = [all(c -> c in moU, mouths[v]) for v in C]
    succ = [[pos[w] for w in out[C[k]] if haskey(pos, w)] for k in 1:nC]
    pred = [Vector{Int}() for _ in 1:nC]
    for k in 1:nC
        for l in succ[k]
            push!(pred[l], k)
        end
    end

    fails = zeros(Int, p)
    for i in 1:p
        inext = mod1(i + 1, p)
        iskip = mod1(i + 2, p)
        good(k, f) = moinU[k] || ((labels[k] == inext) && (f == 0))
        newflag(l, f) = ((f == 1) || (labels[l] == iskip)) ? 1 : 0
        idx(k, f) = k + nC * f
        safe = falses(2 * nC)
        remaining = zeros(Int, 2 * nC)
        queue = Vector{Int}()
        for k in 1:nC, f in 0:1
            if good(k, f)
                safe[idx(k, f)] = true
                push!(queue, idx(k, f))
            else
                remaining[idx(k, f)] = length(succ[k])
            end
        end
        while !isempty(queue)
            q = pop!(queue)
            l = (q - 1) % nC + 1
            fl = (q - 1) ÷ nC
            for k in pred[l], f in 0:1
                newflag(l, f) == fl || continue
                s = idx(k, f)
                safe[s] && continue
                remaining[s] -= 1
                if remaining[s] == 0
                    safe[s] = true
                    push!(queue, s)
                end
            end
        end
        fails[i] = count(k -> (labels[k] == i) && !safe[idx(k, 0)], 1:nC)
    end

    return (C2=C2, C3=C3, C4=all(isequal(0), fails), fails=fails,
            sizes=[count(isequal(i), labels) for i in 1:p])
end

"""
    values, diameters = check_brick_lyapunov(lc::AbstractComplex,
                            coords::Vector{<:Vector{<:Real}},
                            mvf::Vector{Vector{Int}}, C::Vector{Int}, vf;
                            order::Int=2, method::Symbol=:optimal,
                            iters::Int=2000)

Sampled test for the hypothesis that solutions leave the closures of
multivectors.

For every index `k` in `C`, the vector field `vf` is evaluated at the
points of a barycentric grid of order `order` on every cell of the
multivector `mvf[C[k]]`, where `coords` contains the vertex coordinates
as in `create_planar_mvf` and `create_spatial_mvf`. The samples are
normalized, and for a unit vector `v` the value `values[k]` is the
smallest inner product of `v` with a normalized sample. If it is
positive, then the function `x -> v ⋅ x` increases along all solutions in
the closure `B_V` of the multivector, up to the sampling, and therefore
no solution stays in `B_V` for all times. This is hypothesis (D) of the
theorem "Periodic orbits via multivector fields" in the paper *Global
dynamics via multivector fields on Lefschetz complexes*. For
`method=:mean`, the vector `v` is the normalized mean of the normalized
samples. For `method=:optimal`, it is the direction of the point of
smallest norm in the convex hull of the normalized samples, which is
approximated by the Frank-Wolfe (Gilbert) iteration with at most `iters`
steps. This direction maximizes the smallest inner product, and the value
is `-1` if the convex hull contains the origin. The function returns two
vectors of the same length as `C`. The entry `values[k]` is the smallest
inner product described above for the multivector `mvf[C[k]]`, and
`diameters[k]` is the Euclidean diameter of the set of vertices of its
closure, i.e., of the set `B_V`.

The test is not rigorous, since the vector field is only evaluated at
finitely many points.
"""
function check_brick_lyapunov(lc::AbstractComplex, coords::Vector{<:Vector{<:Real}},
                              mvf::Vector{Vector{Int}}, C::Vector{Int}, vf;
                              order::Int=2, method::Symbol=:optimal, iters::Int=2000)
    #
    # Sampled direction test for the closures of multivectors
    #
    method in (:mean, :optimal) ||
        throw(ArgumentError("check_brick_lyapunov: method has to be :mean or :optimal"))
    values = fill(-1.0, length(C))
    diameters = zeros(length(C))

    Threads.@threads for k in 1:length(C)
        G = Vector{Vector{Float64}}()
        X = Vector{Vector{Float64}}()
        for c in mvf[C[k]]
            verts = lefschetz_skeleton(lc, [c], 0)
            Y = [Vector{Float64}(coords[v]) for v in verts]
            append!(X, Y)
            d = length(Y)
            for w in Iterators.product(ntuple(_ -> 0:order, d)...)
                sum(w) == order || continue
                g = vf(sum(w[i] / order * Y[i] for i in 1:d))
                ng = norm(g)
                ng > 0 && push!(G, g / ng)
            end
        end
        diameters[k] = isempty(X) ? 0.0 : maximum(norm(a - b) for a in X, b in X)
        isempty(G) && continue

        if method == :mean
            s = sum(G)
            norm(s) > 0 && (values[k] = minimum(dot(s / norm(s), g) for g in G))
        else
            x = sum(G) / length(G)
            for _ in 1:iters
                j = argmin([dot(x, g) for g in G])
                dot(x, x) - dot(x, G[j]) < 1e-14 && break
                dd = x - G[j]
                t = clamp(dot(x, dd) / dot(dd, dd), 0, 1)
                x = x - t * dd
            end
            if norm(x) >= 1e-10
                v = x / norm(x)
                values[k] = minimum(dot(v, g) for g in G)
            end
        end
    end

    return values, diameters
end

"""
    values, diameters = check_brick_lyapunov(ec::EuclideanComplex,
                            mvf::Vector{Vector{Int}}, C::Vector{Int}, vf;
                            order::Int=2, method::Symbol=:optimal,
                            iters::Int=2000)

Sampled test for the hypothesis that solutions leave the closures of
multivectors on a `EuclideanComplex`, see the method for Lefschetz
complexes.
"""
function check_brick_lyapunov(ec::EuclideanComplex, mvf::Vector{Vector{Int}},
                              C::Vector{Int}, vf; order::Int=2,
                              method::Symbol=:optimal, iters::Int=2000)
    return check_brick_lyapunov(ec, euclidean_vertex_coords(ec), mvf, C, vf,
                                order=order, method=method, iters=iters)
end

"""
    cycles = cyclic_local_cycles(out::Vector{Vector{Int}}, C::Vector{Int},
                                 param::Vector{<:Real}; width::Real=0.02,
                                 nwindows::Int=100)

Find the cycles of a digraph which stay within short windows of a cyclic
parameter.

The adjacency lists `out` describe a digraph, for example the multivector
digraph from `mvf_digraph`, `C` is a set of its vertices, and `param[k]`
is a parameter in `[0,1)` of the vertex `C[k]`, for example the relative
position along a periodic orbit. For each of the `nwindows` windows of
length `width`, which start at the points `w / nwindows`, the function
computes the strongly connected components of the subgraph induced by
the vertices of `C` with parameter in the window. It returns the vector
`cycles` of the components with more than one vertex, without
repetitions. Every component is a sorted vector of vertices of the
digraph, i.e., of elements of `C` such as multivector indices, and not of
positions in `C`. The components are sorted by their smallest vertex, and
the result is empty if there are no such cycles. Such a component is a cycle in the digraph along
which the parameter does not advance, and it prevents hypothesis (C4) of
`check_cyclic_conditions`. If it is an artifact of the mesh, the union
of its multivectors can be passed to `create_strong_mvf` as a merge set.
"""
function cyclic_local_cycles(out::Vector{Vector{Int}}, C::Vector{Int},
                             param::Vector{<:Real}; width::Real=0.02,
                             nwindows::Int=100)
    #
    # Strongly connected components within windows of the parameter
    #
    length(param) == length(C) ||
        throw(ArgumentError("cyclic_local_cycles: param and C must have the same length"))
    pos = Dict{Int,Int}()
    for k in 1:length(C)
        pos[C[k]] = k
    end
    found = Set{Vector{Int}}()
    for w in 0:(nwindows - 1)
        S = [k for k in 1:length(C) if mod(param[k] - w / nwindows, 1.0) < width]
        isempty(S) && continue
        loc = Dict{Int,Int}()
        for (j, k) in enumerate(S)
            loc[k] = j
        end
        g = Graphs.SimpleDiGraph(length(S))
        for k in S, x in out[C[k]]
            l = get(pos, x, 0)
            if (l != 0) && haskey(loc, l)
                Graphs.add_edge!(g, loc[k], loc[l])
            end
        end
        for comp in Graphs.strongly_connected_components(g)
            length(comp) > 1 && push!(found, sort(C[S[comp]]))
        end
    end
    return sort(collect(found), by=first)
end
