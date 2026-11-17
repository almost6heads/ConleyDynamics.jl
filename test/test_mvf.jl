@testset "MVF acyclicity" begin
    # Single gradient pair on a one-edge complex: a→ab, b is critical → acyclic
    sc_edge = create_simplicial_complex(["a","b"], [["a","b"]])
    @test mvf_is_acyclic(sc_edge, [["a","ab"]])

    # Gradient variant from example_three_cm(1): no periodic orbits → acyclic
    lc_grad, mvf_grad = example_three_cm(1)
    @test mvf_is_acyclic(lc_grad, mvf_grad)

    # Periodic orbit example from example_three_cm(0) is NOT acyclic
    lc_three, mvf_periodic = example_three_cm(0)
    @test !mvf_is_acyclic(lc_three, mvf_periodic)

    # Empty mvf (all cells critical singletons) must not crash and is acyclic
    sc_singletons = create_simplicial_complex(["a","b","c"], [["a"],["b"],["c"]])
    @test mvf_is_acyclic(sc_singletons, Vector{Int}[])
end

@testset "MVF graphs agree with brute force" begin
    # Brute-force references based on pairwise comparisons of multivectors
    function full_mvf(lc, mvfI)
        covered = isempty(mvfI) ? Int[] : reduce(vcat, mvfI)
        return vcat(mvfI, [[k] for k in setdiff(1:lc.ncells, covered)])
    end
    function ref_adjacency(lc, mvf, undirected)
        cl = [lefschetz_closure(lc, mv) for mv in mvf]
        other = undirected ? cl : mvf
        return [[m for m in 1:length(mvf) if m != k && !isempty(intersect(cl[k], other[m]))]
                for k in 1:length(mvf)]
    end
    function ref_acyclic(adj)
        state = fill(0, length(adj))          # 0 new, 1 active, 2 done
        function visit(k)
            state[k] = 1
            for m in adj[k]
                state[m] == 1 && return false
                state[m] == 0 && !visit(m) && return false
            end
            state[k] = 2
            return true
        end
        return all(k -> state[k] != 0 || visit(k), 1:length(adj))
    end
    function ref_reach(adj, mvf, sources, n)
        dist = fill(typemax(Int), length(adj))
        front = [k for k in 1:length(mvf) if !isempty(intersect(mvf[k], sources))]
        dist[front] .= 0
        for d in 1:n
            front = unique([m for k in front for m in adj[k] if dist[m] > d])
            dist[front] .= d
        end
        return mvf[findall(d -> d <= n, dist)]
    end

    function samplevf(x::Vector{Float64})
        x1, x2 = x
        return [x1 * (1.0 - x1*x1 - 3.0*x2*x2), x2 * (1.0 - 3.0*x1*x1 - x2*x2)]
    end
    lc6, coords6 = create_simplicial_rectangle(6, 6)
    mvf6 = create_planar_mvf(lc6, convert_planar_coordinates(coords6, [-2.0,-2.0], [2.0,2.0]), samplevf)

    for (lc, mvf) in (example_three_cm(), example_three_cm(1), example_clorenz(),
                      example_torsion_chaos(3, 3)[[1, 3]], (lc6, mvf6))
        mvfI = mvf isa Vector{Vector{Int}} ? mvf : convert_cellsubsets(lc, mvf)
        mvfF = full_mvf(lc, mvfI)
        fwd = ref_adjacency(lc, mvfF, false)
        bwd = [[k for k in 1:length(mvfF) if m in fwd[k]] for m in 1:length(mvfF)]
        und = ref_adjacency(lc, mvfF, true)
        @test mvf_is_acyclic(lc, mvfI) == ref_acyclic(fwd)
        for sources in ([1], [lc.ncells], collect(1:5:lc.ncells)), n in 0:3
            @test mvf_forward_orbit(lc, mvfI, sources, n) == ref_reach(fwd, mvfF, sources, n)
            @test mvf_backward_orbit(lc, mvfI, sources, n) == ref_reach(bwd, mvfF, sources, n)
            @test mvf_neighborhood(lc, mvfI, sources, n) == ref_reach(und, mvfF, sources, n)
        end
    end
end

@testset "MVF information" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)
    formanvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]

    info = mvf_information(sc, formanvf)

    # Required keys exist
    @test haskey(info, "N mv")
    @test haskey(info, "N regular")
    @test haskey(info, "N critical")
    @test haskey(info, "Lengths regular")
    @test haskey(info, "Lengths critical")

    # Total multivectors = regular + critical
    @test info["N mv"] == info["N regular"] + info["N critical"]
end

@testset "Morse sets" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)
    formanvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]

    ms = morse_sets(sc, formanvf)

    # The gradient field has critical cells: F, DF, EF, DEF, and the connected component {A,B,C,...}
    # Number of nontrivial Morse sets equals number in cm.morse from Tutorial testset
    cm = connection_matrix(sc, formanvf)
    @test length(ms) == length(cm.morse)

    # Periodic orbit example: one Morse set containing the orbit
    lc_three, mvf_periodic = example_three_cm(0)
    ms_per = morse_sets(lc_three, mvf_periodic)
    @test length(ms_per) >= 1
end

@testset "MVF hull agrees with brute force" begin
    # Brute-force least fixed point: merge every base set into one block,
    # then merge the locally closed hull of every block until all blocks
    # are locally closed
    function hull_bruteforce(lc, base)
        lab = collect(1:lc.ncells)
        function unify!(cells)
            t = lab[cells[1]]
            for c in cells
                o = lab[c]
                o == t || (lab[lab .== o] .= t)
            end
        end
        for s in base
            length(s) > 1 && unify!(s)
        end
        changed = true
        while changed
            changed = false
            for l in unique(lab)
                B = findall(==(l), lab)
                length(B) > 1 || continue
                H = lefschetz_lchull(lc, B)
                if length(H) > length(B)
                    unify!(H)
                    changed = true
                end
            end
        end
        mvf = [findall(==(l), lab) for l in unique(lab)]
        return sort!(filter(m -> length(m) > 1, mvf), by=first)
    end

    # Deterministic pseudo-random base sets: a cell together with
    # some cells from its closure or open hull
    state = Ref(UInt64(20261105))
    nextint(n) = (state[] = state[] * 6364136223846793005 + 1442695040888963407;
                  Int((state[] >> 33) % UInt64(n)) + 1)
    function random_base(lc, nsets)
        base = Vector{Vector{Int}}()
        for _ in 1:nsets
            c = nextint(lc.ncells)
            nb = setdiff(union(lefschetz_closure(lc, [c]), lefschetz_openhull(lc, [c])), [c])
            s = [c]
            if !isempty(nb)
                for _ in 1:nextint(3)
                    push!(s, nb[nextint(length(nb))])
                end
            end
            push!(base, isodd(nextint(2)) ? s : reverse(s))
        end
        return base
    end

    complexes = [example_moebius(2)[1], example_three_cm(0)[1],
                 create_cubical_rectangle(5, 4)[1], create_simplicial_rectangle(4, 4)[1],
                 create_cubical_box(3, 3, 2)[1]]
    for lc in complexes
        for trial in 1:20
            nsets = [1, 3, lc.ncells ÷ 20 + 1, lc.ncells ÷ 8 + 1, lc.ncells ÷ 3][1 + trial % 5]
            base = random_base(lc, nsets)
            mvf = create_mvf_hull(lc, base)
            @test mvf == hull_bruteforce(lc, base)
        end
        # Edge cases: empty base, only singletons, idempotence, label input
        @test isempty(create_mvf_hull(lc, Vector{Vector{Int}}()))
        @test isempty(create_mvf_hull(lc, [[k] for k in 1:lc.ncells]))
        mvf = create_mvf_hull(lc, random_base(lc, lc.ncells ÷ 8 + 1))
        @test create_mvf_hull(lc, mvf) == mvf
        base = random_base(lc, lc.ncells ÷ 8 + 1)
        @test create_mvf_hull(lc, [lc.labels[s] for s in base]) ==
              [lc.labels[m] for m in create_mvf_hull(lc, base)]
    end
end

@testset "Exit set removal agrees with brute force" begin
    # The exit set consists of the critical cells on the manifold boundary
    for vf in (x -> [x[1], -x[2], -x[3]], x -> -x, x -> x)   # saddle, sink, source
        ec = create_cubical_box(4, 4, 4; euclidean=true)
        svf(x) = vf(x .- 2.0)
        mvf = create_spatial_mvf(ec, svf)
        bnd = manifold_boundary(ec)
        crit = [c for c in bnd if !any(mv -> length(mv) > 1 && c in mv, mvf)]
        exitset = lefschetz_closure(ec, crit)
        lcr, mvfr = remove_exit_set(ec, mvf)
        @test lcr.ncells == ec.ncells - length(exitset)
        @test sort(lcr.labels) == sort(ec.labels[setdiff(1:ec.ncells, exitset)])
    end
end

@testset "Flow transitions" begin
    #
    # The transitions reproduce create_planar_mvf and create_spatial_mvf
    #
    function pvf(x::Vector{Float64})
        x1, x2 = x
        return [x1 * (1.0 - x1*x1 - 3.0*x2*x2), x2 * (1.0 - 3.0*x1*x1 - x2*x2)]
    end
    lc, coords = create_simplicial_rectangle(8, 8)
    coordsN = convert_planar_coordinates(coords, [-2.0,-2.0], [2.0,2.0])
    fwd = planar_transitions(lc, coordsN, pvf)
    @test create_mvf_hull(lc, fwd) == create_planar_mvf(lc, coordsN, pvf)
    @test all(t -> lc.dimensions[t[1]] < 2 && all(lc.dimensions[t[2:end]] .== 2), fwd)

    fwd5 = planar_transitions(lc, coordsN, pvf, npts=5)
    @test create_mvf_hull(lc, fwd5) == create_planar_mvf(lc, coordsN, pvf, npts=5)

    # The stricter test only enlarges the transitions, and it flags
    # the same edges for the reversed vector field

    fwdnt = planar_transitions(lc, coordsN, pvf, nontransverse=true)
    @test length(fwdnt) == length(fwd)
    @test all(issubset(fwd[k], fwdnt[k]) for k in eachindex(fwd))
    @test sort(planar_nontransverse_edges(lc, coordsN, pvf)) ==
          sort(planar_nontransverse_edges(lc, coordsN, x -> -pvf(x)))

    # The EuclideanComplex method agrees with the Lefschetz complex method

    ec = create_simplicial_rectangle(8, 8, euclidean=true)
    lce = euclidean_to_lefschetz(ec)
    evf(x) = pvf(x ./ 2.0 .- 2.0)
    @test create_mvf_hull(lce, planar_transitions(ec, evf)) == create_planar_mvf(ec, evf)

    # The vertex test uses the direction of the vector field directly, so
    # that a small transverse component is not lost for large coordinates

    lcv, coordsv = create_simplicial_rectangle(2, 2)
    coordsL = [isempty(c) ? c : c .+ 1.0e6 for c in coordsv]
    tvf(x) = [1.0, 1.0e-12]
    v = findfirst(k -> lcv.dimensions[k] == 0 && coordsv[k] == [1.0, 1.0], 1:lcv.ncells)
    tv = only(filter(t -> t[1] == v, planar_transitions(lcv, coordsL, tvf)))
    @test length(tv) == 2

    # Spatial transitions

    function svf(x::Vector{Float64})
        x1, x2, x3 = x
        return [x1 * (1.0 - x1*x1), -x2, -x3]
    end
    lc3, coords3 = create_cubical_box(3, 3, 3)
    coords3N = convert_spatial_coordinates(coords3, [-1.5,-1.0,-1.0], [1.5,1.0,1.0])
    fwd3 = spatial_transitions(lc3, coords3N, svf)
    @test create_mvf_hull(lc3, fwd3) == create_spatial_mvf(lc3, coords3N, svf)
    @test create_mvf_hull(lc3, spatial_transitions(lc3, coords3N, svf, npts=5)) ==
          create_spatial_mvf(lc3, coords3N, svf, npts=5)
end

@testset "Dual multivectors" begin
    #
    # Duals for the Ouroboros example, compared with the definition
    # in the dissertation of Justin Thorpe (interior and dual cap)
    #
    ec, vf = example_ouroboros()
    lc = euclidean_to_lefschetz(ec)
    @test lc.ncells == 259
    fwd = planar_transitions(ec, vf)
    bwd = planar_transitions(ec, x -> -vf(x))
    mvf = create_mvf_hull(lc, fwd)
    duals = mvf_dual(lc, mvf, bwd)
    @test duals[3] == mvf_dual_set(lc, mvf[3], bwd)

    ipt = [Vector{Int}() for _ in 1:lc.ncells]
    for t in bwd
        ipt[t[1]] = t[2:end]
    end
    for k in 1:lc.ncells
        if lc.dimensions[k] == 2
            ipt[k] = [k]
        end
    end
    for k in eachindex(mvf)
        V = mvf[k]
        interior = lefschetz_interior(lc, V)
        boundary = lefschetz_topboundary(lc, V)
        thesisdual = union(interior, [c for c in boundary if !isempty(intersect(ipt[c], V))])
        @test duals[k] == sort(filter(c -> !isempty(ipt[c]), thesisdual))
    end

    # The flow enters through all 18 boundary edges, and these
    # entrance cells belong to no dual

    bedges = [c for c in 1:lc.ncells if lc.dimensions[c] == 1 && length(lefschetz_coboundary(lc, c)) == 1]
    @test length(bedges) == 18
    @test all(isempty(ipt[c]) for c in bedges)
    @test all(!(c in D) for D in duals, c in bedges)
end

@testset "Dual merging cross-check" begin
    #
    # The merging procedure of pull request #10 on the Ouroboros example
    # merges exactly one pair of multivectors, and its result has disjoint
    # and locally closed duals with the same toplexes
    #
    ec, vf = example_ouroboros()
    lc = euclidean_to_lefschetz(ec)
    fwd = planar_transitions(ec, vf)
    bwd = planar_transitions(ec, x -> -vf(x))
    mvf0 = create_mvf_hull(lc, fwd)
    mvf = ConleyDynamics.strong_mvf_dual_merging(lc, fwd, bwd)
    @test length(mvf0) == 65
    @test length(mvf) == 64
    @test all(any(issubset(V0, V) for V in mvf) for V0 in mvf0)

    toplexes = findall(isequal(2), lc.dimensions)
    inmv = falses(lc.ncells)
    for V in mvf
        inmv[V] .= true
    end
    allmvf = vcat(mvf, [[t] for t in toplexes if !inmv[t]])
    duals = mvf_dual(lc, allmvf, bwd)
    count = zeros(Int, lc.ncells)
    for D in duals, c in D
        count[c] += 1
    end
    @test maximum(count) == 1
    @test all(D -> lefschetz_is_locally_closed(lc, D), duals)
    @test all(filter(c -> lc.dimensions[c] == 2, duals[k]) ==
              filter(c -> lc.dimensions[c] == 2, allmvf[k]) for k in eachindex(duals))
end

@testset "Multivector digraphs and cyclic partitions" begin
    #
    # The breadth-first split always returns a partition
    #
    out = [[3], [3], Int[]]
    @test cyclic_partition_bfs(out, [1, 2, 3], [1, 2]) == [1, 2, 1]
    cycle = [[2], [3], [4], [1]]
    @test cyclic_partition_bfs(cycle, [1, 2, 3, 4], [1]) == [1, 1, 1, 1]
    @test cyclic_partition_bfs(cycle, [1, 2, 3, 4], [1, 3]) == [1, 1, 2, 2]
    @test cyclic_partition_bfs(cycle, [2, 3, 4], [2]) == [1, 1, 1]
    @test_logs (:warn,) cyclic_partition_bfs(cycle, [1, 2, 3], [2])

    # The digraph agrees with the definition via mouths

    ec, vf = example_ouroboros()
    lc = euclidean_to_lefschetz(ec)
    mvf = create_mvf_hull(lc, planar_transitions(ec, vf))
    inmv = falses(lc.ncells)
    for V in mvf
        inmv[V] .= true
    end
    allmvf = vcat(mvf, [[c] for c in 1:lc.ncells if !inmv[c]])
    dg, cell2mv, mouths = mvf_digraph(lc, allmvf)
    @test all(cell2mv .> 0)
    for k in eachindex(allmvf)
        mo = setdiff(lefschetz_closure(lc, allmvf[k]), allmvf[k])
        @test sort(mouths[k]) == sort(mo)
        @test sort(dg[k]) == [m for m in eachindex(allmvf) if m != k && !isempty(intersect(mo, allmvf[m]))]
    end

    # The digraph requires a complete multivector field: missing
    # singletons and overlapping sets lead to an error

    @test mvf_digraph(lc, mvf_with_singletons(lc, mvf)) == (dg, cell2mv, mouths)
    incomplete = allmvf[2:end]
    @test_throws ArgumentError mvf_digraph(lc, incomplete)
    msg = try mvf_digraph(lc, incomplete); "" catch err; err.msg end
    @test occursin("not complete", msg) && occursin("mvf_with_singletons", msg)
    overlapping = vcat(allmvf, [[allmvf[1][1]]])
    @test_throws ArgumentError mvf_digraph(lc, overlapping)
    msg = try mvf_digraph(lc, overlapping); "" catch err; err.msg end
    @test occursin("more than one multivector", msg)
    @test_throws ArgumentError mvf_digraph(lc, vcat(allmvf, [[lc.ncells + 1]]))

    # The same check in check_cyclic_conditions, also if the digraph is given

    Call = collect(eachindex(allmvf))
    lab3 = [mod1(k, 3) for k in Call]
    @test_throws ArgumentError check_cyclic_conditions(lc, incomplete, [1, 2, 3],
                                                       [1, 2, 3], 3)
    @test_throws ArgumentError check_cyclic_conditions(lc, incomplete, [1, 2, 3],
                                                       [1, 2, 3], 3; out=dg, mouths=mouths)
    @test_throws ArgumentError check_cyclic_conditions(lc, overlapping, Call, lab3, 3;
                                                       out=dg, mouths=mouths)
    @test check_cyclic_conditions(lc, allmvf, Call, lab3, 3; out=dg, mouths=mouths) isa NamedTuple

    # Angular sectors around the center of the rectangle

    labels = cyclic_partition_angular(ec, allmvf, collect(eachindex(allmvf)), 4, center=[2.0, 2.5])
    @test sort(unique(labels)) == [1, 2, 3, 4]
    vertex(p) = findfirst(k -> lc.dimensions[k] == 0 && ec.coords[k][1] == p, 1:lc.ncells)
    corners = [[vertex([4.0, 5.0])], [vertex([0.0, 0.0])], [vertex([4.0, 0.0]), vertex([4.0, 1.0])]]
    @test cyclic_partition_angular(ec, corners, [1, 2, 3], 4, center=[2.0, 2.5]) == [1, 3, 4]
    @test cyclic_partition_angular(ec, corners, [1, 2, 3], 4, center=[2.0, 2.5],
                                   clockwise=true) == [4, 2, 1]
end

@testset "Finest strongly accumulating multivector fields" begin
    #
    # The worklist implementation agrees with the iteration in the proof
    # of the theorem, and in the planar case with the dual merging of
    # pull request #10. The result is strongly accumulating.
    #
    canon(m) = Set(sort.(m))

    function check_strong(lc, fwd, bwd)
        mvf, dual, info = create_strong_mvf(lc, fwd, bwd, details=true)
        rmvf, rdual, rounds = ConleyDynamics.strong_mvf_reference(lc, fwd, bwd)
        @test canon(mvf) == canon(rmvf)
        @test canon(dual) == canon(rdual)

        # Every forward and backward transition lies in one multivector

        isexit = falses(lc.ncells)
        isexit[info["exit cells"]] .= true
        isentr = falses(lc.ncells)
        isentr[info["entrance cells"]] .= true
        c2m = zeros(Int, lc.ncells)
        for (k, V) in enumerate(mvf), c in V
            c2m[c] = k
        end
        c2d = zeros(Int, lc.ncells)
        for (k, V) in enumerate(dual), c in V
            c2d[c] = k
        end
        @test all(t -> isexit[t[1]] || all(c -> c == t[1] || c2m[c] == c2m[t[1]], t) &&
                       (c2m[t[1]] != 0), fwd)
        @test all(t -> isentr[t[1]] || all(c -> c == t[1] || c2d[c] == c2d[t[1]], t) &&
                       (c2d[t[1]] != 0), bwd)
        @test all(V -> lefschetz_is_locally_closed(lc, V), mvf)
        @test all(V -> lefschetz_is_locally_closed(lc, V), dual)
        @test all(V -> !any(isexit[V]), mvf)
        @test all(V -> !any(isentr[V]), dual)

        # The multivectors and their duals have the same toplexes, and the
        # duals agree with the characterization via immediate pasts

        tops(V) = filter(c -> lc.dimensions[c] == lc.dim, V)
        for k in eachindex(mvf)
            j = info["mvf to dual"][k]
            if j == 0
                @test length(tops(mvf[k])) == 1
            else
                @test tops(mvf[k]) == tops(dual[j])
                @test info["dual to mvf"][j] == k
            end
        end
        toplexes = findall(isequal(lc.dim), lc.dimensions)
        inmv = falses(lc.ncells)
        for V in mvf
            inmv[V] .= true
        end
        allmvf = vcat(mvf, [[t] for t in toplexes if !inmv[t]])
        @test canon(filter(D -> length(D) > 1, mvf_dual(lc, allmvf, bwd))) == canon(dual)

        # The result is coarser than the accumulating multivector field

        mvf0 = create_mvf_hull(lc, [t for t in fwd if !isexit[t[1]]])
        @test all(any(issubset(V0, V) for V in mvf) for V0 in mvf0)
        return mvf, dual, info, rounds
    end

    # The Ouroboros example (attractor, all boundary cells are entrance cells)

    ec, vf = example_ouroboros()
    lc = euclidean_to_lefschetz(ec)
    fwd = planar_transitions(ec, vf)
    bwd = planar_transitions(ec, x -> -vf(x))
    mvf, dual, info, rounds = check_strong(lc, fwd, bwd)
    @test length(mvf) == 64
    @test rounds == 1
    @test isempty(info["exit cells"])
    @test length(info["entrance cells"]) == 36
    @test canon(mvf) == canon(ConleyDynamics.strong_mvf_dual_merging(lc, fwd, bwd))
    mvfw, dualw = create_planar_strong_mvf(ec, vf)
    @test mvfw == mvf
    @test dualw == dual

    # The dictionary info is only returned with details=true

    res = create_strong_mvf(lc, fwd, bwd)
    @test length(res) == 2
    @test res == (mvf, dual)
    res = create_strong_mvf(lc, fwd, bwd, details=true)
    @test length(res) == 3
    @test res[1] == mvf && res[2] == dual && res[3] isa Dict
    resw = create_planar_strong_mvf(ec, vf, details=true)
    @test length(resw) == 3
    @test resw[1] == mvf && resw[2] == dual
    @test resw[3]["entrance cells"] == res[3]["entrance cells"]

    # A tilted flow with entrance, exit, and tangency cells

    ec = create_simplicial_rectangle(6, 4, euclidean=true)
    lc = euclidean_to_lefschetz(ec)
    tvf(x) = [1.0, 0.3]
    fwd = planar_transitions(ec, tvf)
    bwd = planar_transitions(ec, x -> -tvf(x))
    mvf, dual, info, rounds = check_strong(lc, fwd, bwd)
    @test canon(mvf) == canon(ConleyDynamics.strong_mvf_dual_merging(lc, fwd, bwd))
    vertex(p) = findfirst(k -> lc.dimensions[k] == 0 && ec.coords[k][1] == p, 1:lc.ncells)
    tangency = intersect(info["exit cells"], info["entrance cells"])
    @test sort(tangency) == sort([vertex([6.0, 0.0]), vertex([0.0, 4.0])])
    @test lefschetz_is_closed(lc, info["exit cells"])
    @test lefschetz_is_closed(lc, info["entrance cells"])

    # A saddle at the center of the rectangle

    ec = create_simplicial_rectangle(8, 8, euclidean=true)
    lc = euclidean_to_lefschetz(ec)
    svf(x) = [x[1] - 4.0, 4.0 - x[2]]
    fwd = planar_transitions(ec, svf)
    bwd = planar_transitions(ec, x -> -svf(x))
    mvf, dual, info, rounds = check_strong(lc, fwd, bwd)
    @test canon(mvf) == canon(ConleyDynamics.strong_mvf_dual_merging(lc, fwd, bwd))

    # A cubical box in three dimensions

    lc3, coords3 = create_cubical_box(4, 4, 4)
    coords3N = convert_spatial_coordinates(coords3, [-1.5,-1.0,-1.0], [1.5,1.0,1.0])
    v3(x) = [x[1] * (1.0 - x[1]^2) - 0.3*x[2], -x[2] + 0.3*x[1], -x[3]]
    fwd3 = spatial_transitions(lc3, coords3N, v3)
    bwd3 = spatial_transitions(lc3, coords3N, x -> -v3(x))
    mvf3, dual3, info3, rounds3 = check_strong(lc3, fwd3, bwd3)
    mvf3w, dual3w = create_spatial_strong_mvf(lc3, coords3N, v3)
    @test mvf3w == mvf3

    # Errors: a flow along the boundary, boundary cells which are not
    # closed, and cells without transitions

    ec = create_simplicial_rectangle(6, 4, euclidean=true)
    lc = euclidean_to_lefschetz(ec)
    hvf(x) = [1.0, 0.0]
    @test_throws ArgumentError create_planar_strong_mvf(ec, hvf)
    fwd = planar_transitions(ec, tvf)
    bwd = planar_transitions(ec, x -> -tvf(x))
    exitcells, entrancecells = ConleyDynamics.strong_mvf_boundary_cells(lc, fwd, bwd)
    edge = first(c for c in exitcells if lc.dimensions[c] == 1)
    @test_throws ArgumentError create_strong_mvf(lc, fwd, bwd, exitcells=[edge],
                                                 entrancecells=entrancecells)
    @test_throws ArgumentError create_strong_mvf(lc, fwd, bwd, exitcells=Int[],
                                                 entrancecells=entrancecells)
    @test_throws ArgumentError create_strong_mvf(lc, fwd[2:end], bwd)
end

@testset "Checking strongly accumulating multivector fields" begin
    #
    # The result of create_strong_mvf passes all checks, while the
    # accumulating field of the Ouroboros example fails at its dual
    #
    ec, vf = example_ouroboros()
    lc = euclidean_to_lefschetz(ec)
    fwd = planar_transitions(ec, vf)
    bwd = planar_transitions(ec, x -> -vf(x))
    mvf, dual, info = create_strong_mvf(lc, fwd, bwd, details=true)
    ok, report = mvf_is_strongly_accumulating(lc, mvf, dual, fwd, bwd, details=true)
    @test ok
    @test all(isempty, values(report))
    @test mvf_is_strongly_accumulating(lc, mvf, dual, fwd, bwd)

    toplexes = findall(isequal(2), lc.dimensions)
    mvf0 = create_mvf_hull(lc, fwd)
    inmv = falses(lc.ncells)
    for V in mvf0
        inmv[V] .= true
    end
    dual0 = filter(D -> length(D) > 1, mvf_dual(lc, vcat(mvf0, [[t] for t in toplexes if !inmv[t]]), bwd))
    ok0, report0 = mvf_is_strongly_accumulating(lc, mvf0, dual0, fwd, bwd, details=true)
    @test !ok0
    @test length(report0["dual not locally closed"]) == 1
    @test isempty(report0["toplexes with different partitions"])

    # Merging two adjacent multivectors without changing the dual is detected

    m = findfirst(j -> j != 1 && !isempty(intersect(lefschetz_closure(lc, mvf[j]), mvf[1])), eachindex(mvf))
    mvf2 = vcat([sort(vcat(mvf[1], mvf[m]))], mvf[setdiff(eachindex(mvf), [1, m])])
    ok2, report2 = mvf_is_strongly_accumulating(lc, mvf2, dual, fwd, bwd, details=true)
    @test !ok2
    @test !isempty(report2["toplexes with different partitions"])

    # A multivector which does not contain its forward transitions

    mvf3 = deepcopy(mvf)
    k = findfirst(V -> any(c -> lc.dimensions[c] == 1, V), mvf3)
    e = first(c for c in mvf3[k] if lc.dimensions[c] == 1)
    mvf3[k] = setdiff(mvf3[k], [e])
    ok3, report3 = mvf_is_strongly_accumulating(lc, mvf3, dual, fwd, bwd, details=true)
    @test !ok3
    @test e in report3["split forward transitions"]

    # Mixed boundary types and three dimensions

    ecT = create_simplicial_rectangle(6, 4, euclidean=true)
    lcT = euclidean_to_lefschetz(ecT)
    tvf(x) = [1.0, 0.3]
    fwdT = planar_transitions(ecT, tvf)
    bwdT = planar_transitions(ecT, x -> -tvf(x))
    mvfT, dualT, infoT = create_strong_mvf(lcT, fwdT, bwdT, details=true)
    @test mvf_is_strongly_accumulating(lcT, mvfT, dualT, fwdT, bwdT)
    @test mvf_is_strongly_accumulating(lcT, mvfT, dualT, fwdT, bwdT,
                                       exitcells=infoT["exit cells"],
                                       entrancecells=infoT["entrance cells"])

    lc3, coords3 = create_cubical_box(4, 4, 4)
    coords3N = convert_spatial_coordinates(coords3, [-1.5,-1.0,-1.0], [1.5,1.0,1.0])
    v3(x) = [x[1] * (1.0 - x[1]^2) - 0.3*x[2], -x[2] + 0.3*x[1], -x[3]]
    fwd3 = spatial_transitions(lc3, coords3N, v3)
    bwd3 = spatial_transitions(lc3, coords3N, x -> -v3(x))
    mvf3d, dual3d = create_strong_mvf(lc3, fwd3, bwd3)
    ok3d, report3d = mvf_is_strongly_accumulating(lc3, mvf3d, dual3d, fwd3, bwd3, details=true)
    @test ok3d
    @test isempty(report3d["disconnected"])

    # Isolating blocks: a single multivector, and the union of all of them

    for V in mvf
        mouth, dualset, dualmouth, islc = strong_isoblock_sets(lc, V, mvf, dual)
        @test islc
        @test mouth == sort(setdiff(lefschetz_closure(lc, V), V))
        @test sort(lefschetz_closure(lc, dualset)) == sort(lefschetz_closure(lc, V))
        @test sort(union(dualset, dualmouth)) == sort(union(V, mouth))
    end
    U = collect(1:lc.ncells)
    mouth, dualset, dualmouth, islc = strong_isoblock_sets(lc, U, mvf, dual)
    @test isempty(mouth)
    @test islc
    @test sort(dualmouth) == sort(info["entrance cells"])
    @test_throws ArgumentError strong_isoblock_sets(lc, mvf[1][1:1], mvf, dual)
end

@testset "Strongly accumulating fields: minimality and invariants" begin
    #
    # A small deterministic random number generator (splitmix64), so that
    # the test instances do not depend on the Julia version
    #
    rngstate = Ref(UInt64(2026))
    function rand01()
        rngstate[] += 0x9e3779b97f4a7c15
        z = rngstate[]
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        z = z ⊻ (z >> 31)
        return (z >> 11) / 2.0^53
    end
    randelement(v) = v[1 + floor(Int, rand01() * length(v))]

    # Flow-like transitions: every toplex gets a random height, every cell
    # flows forward to the highest toplex in its star and backward to the
    # lowest one, and with probability pextra a further toplex is added

    function height_transitions(lc; pextra=0.15)
        h = [rand01() for _ in 1:lc.ncells]
        fwd = Vector{Vector{Int}}()
        bwd = Vector{Vector{Int}}()
        for c in 1:lc.ncells
            lc.dimensions[c] == lc.dim && continue
            tops = filter(x -> lc.dimensions[x] == lc.dim, lefschetz_openhull(lc, [c]))
            f = [tops[argmax(h[tops])]]
            b = [tops[argmin(h[tops])]]
            rand01() < pextra && push!(f, randelement(tops))
            rand01() < pextra && push!(b, randelement(tops))
            push!(fwd, vcat([c], unique(f)))
            push!(bwd, vcat([c], unique(b)))
        end
        return fwd, bwd
    end

    # Random boundary types: every boundary face of codimension one becomes
    # an exit face (no forward transition), an entrance face (no backward
    # transition), or a tangency face (neither)

    function random_boundary!(lc, fwd, bwd)
        bfaces = [c for c in 1:lc.ncells if lc.dimensions[c] == lc.dim - 1 &&
                                            length(lefschetz_coboundary(lc, c)) == 1]
        nofwd = Set{Int}()
        nobwd = Set{Int}()
        for c in bfaces
            r = rand01()
            r < 0.4 && push!(nofwd, c)
            (r >= 0.4) && (r < 0.8) && push!(nobwd, c)
            r >= 0.8 && (push!(nofwd, c); push!(nobwd, c))
        end
        filter!(t -> !(t[1] in nofwd), fwd)
        filter!(t -> !(t[1] in nobwd), bwd)
        return fwd, bwd
    end

    # Partitions of the toplexes, as canonical label vectors

    canonical(lab) = (d = Dict{Int,Int}(); [get!(d, x, length(d)) for x in lab])
    function toplabels(lc, mvf, tops)
        c2m = zeros(Int, lc.ncells)
        for (k, V) in enumerate(mvf), c in V
            c2m[c] = k
        end
        return canonical([c2m[t] == 0 ? -t : c2m[t] for t in tops])
    end
    refines(a, b) = all(a[i] != a[j] || b[i] == b[j] for i in eachindex(a), j in eachindex(a))
    canon(m) = Set(sort.(m))

    # Brute force: a partition T of the toplexes belongs to a strongly
    # accumulating multivector field if and only if the minimal forward and
    # backward fields containing the transitions and T both induce T. The
    # output has to be the unique finest such partition.

    function all_valid_partitions(lc, fwd, bwd)
        tops = findall(isequal(lc.dim), lc.dimensions)
        n = length(tops)
        valid = Vector{Vector{Int}}()
        a = zeros(Int, n)
        function visit(i, maxlabel)
            if i > n
                groups = Dict{Int,Vector{Int}}()
                for k in 1:n
                    push!(get!(groups, a[k], Int[]), tops[k])
                end
                T = [g for g in values(groups) if length(g) > 1]
                V = create_mvf_hull(lc, vcat(fwd, T))
                W = create_mvf_hull(lc, vcat(bwd, T))
                if toplabels(lc, V, tops) == a && toplabels(lc, W, tops) == a
                    push!(valid, copy(a))
                end
                return
            end
            for v in 0:(maxlabel + 1)
                a[i] = v
                visit(i + 1, max(maxlabel, v))
            end
        end
        a[1] = 0
        visit(2, 0)
        return valid
    end

    nontrivial = 0
    for lc in (create_simplicial_rectangle(1, 2)[1], create_cubical_rectangle(2, 4)[1],
               create_cubical_box(2, 2, 1)[1])
        tops = findall(isequal(lc.dim), lc.dimensions)
        for trial in 1:4
            fwd, bwd = height_transitions(lc)
            mvf, dual = create_strong_mvf(lc, fwd, bwd, exitcells=Int[], entrancecells=Int[])
            valid = all_valid_partitions(lc, fwd, bwd)
            lab = toplabels(lc, mvf, tops)
            @test lab in valid
            @test all(refines(lab, v) for v in valid)
            @test toplabels(lc, dual, tops) == lab
            nontrivial += (length(valid) > 1)
        end
    end
    @test nontrivial >= 6

    # General transitions, which are open hulls of random subsets of the
    # star and may contain lower-dimensional cells

    lc = create_simplicial_rectangle(1, 1)[1]
    tops = findall(isequal(2), lc.dimensions)
    for trial in 1:4
        transitions = Vector{Vector{Vector{Int}}}()
        for direction in 1:2
            tr = Vector{Vector{Int}}()
            for c in 1:lc.ncells
                lc.dimensions[c] == 2 && continue
                star = lefschetz_openhull(lc, [c])
                S = filter(_ -> rand01() < 0.3, star)
                isempty(S) && (S = [randelement(star)])
                targets = setdiff(lefschetz_openhull(lc, S), [c])
                push!(tr, vcat([c], targets))
            end
            push!(transitions, tr)
        end
        fwd, bwd = transitions
        mvf, dual = create_strong_mvf(lc, fwd, bwd, exitcells=Int[], entrancecells=Int[])
        valid = all_valid_partitions(lc, fwd, bwd)
        lab = toplabels(lc, mvf, tops)
        @test lab in valid
        @test all(refines(lab, v) for v in valid)
        @test mvf_is_strongly_accumulating(lc, mvf, dual, fwd, bwd, exitcells=Int[],
                                           entrancecells=Int[])
    end

    # Larger complexes with random boundary types: agreement with the
    # reference implementation (and in the plane with the dual merging of
    # pull request #10), all conditions, idempotence, and monotonicity

    for lc in (create_simplicial_rectangle(4, 4)[1], create_cubical_rectangle(7, 6)[1],
               create_cubical_box(3, 3, 2)[1])
        for trial in 1:3
            fwd, bwd = random_boundary!(lc, height_transitions(lc)...)
            mvf, dual, info = create_strong_mvf(lc, fwd, bwd, details=true)
            rmvf, rdual, rounds = ConleyDynamics.strong_mvf_reference(lc, fwd, bwd)
            @test canon(mvf) == canon(rmvf)
            @test canon(dual) == canon(rdual)
            if lc.dim == 2
                @test canon(mvf) == canon(ConleyDynamics.strong_mvf_dual_merging(lc, fwd, bwd))
            end
            @test mvf_is_strongly_accumulating(lc, mvf, dual, fwd, bwd)

            # Idempotence: adding the toplex partition of the result as
            # transitions does not change the result

            groups = [filter(c -> lc.dimensions[c] == lc.dim, V) for V in mvf]
            seeded = vcat(fwd, filter(g -> length(g) > 1, groups))
            mvfs, duals = create_strong_mvf(lc, seeded, bwd)
            @test canon(mvfs) == canon(mvf)
            @test canon(duals) == canon(dual)

            # Monotonicity: additional transitions lead to a coarser field

            extra = deepcopy(fwd)
            for t in extra
                if rand01() < 0.1
                    tops = filter(x -> lc.dimensions[x] == lc.dim, lefschetz_openhull(lc, [t[1]]))
                    push!(t, randelement(tops))
                    unique!(t)
                end
            end
            mvfm, dualm = create_strong_mvf(lc, extra, bwd,
                                                   exitcells=info["exit cells"],
                                                   entrancecells=info["entrance cells"])
            @test all(any(issubset(V, W) for W in mvfm) for V in mvf)
        end
    end

    # Prescribed merge sets (part (b) of the theorem): brute-force
    # minimality on small complexes, agreement with the reference
    # implementation, the checking function, and connectedness for
    # connected merge sets

    function all_valid_partitions_merge(lc, fwd, bwd, merge)
        tops = findall(isequal(lc.dim), lc.dimensions)
        n = length(tops)
        valid = Vector{Vector{Int}}()
        a = zeros(Int, n)
        function visit(i, maxlabel)
            if i > n
                groups = Dict{Int,Vector{Int}}()
                for k in 1:n
                    push!(get!(groups, a[k], Int[]), tops[k])
                end
                T = [g for g in values(groups) if length(g) > 1]
                V = create_mvf_hull(lc, vcat(fwd, merge, T))
                W = create_mvf_hull(lc, vcat(bwd, T))
                if toplabels(lc, V, tops) == a && toplabels(lc, W, tops) == a
                    push!(valid, copy(a))
                end
                return
            end
            for v in 0:(maxlabel + 1)
                a[i] = v
                visit(i + 1, max(maxlabel, v))
            end
        end
        a[1] = 0
        visit(2, 0)
        return valid
    end

    # A connected merge set: a toplex, a facet, and the other toplex at it

    function random_pair(lc)
        facets = [c for c in 1:lc.ncells if lc.dimensions[c] == lc.dim - 1 &&
                                            length(lefschetz_coboundary(lc, c)) == 2]
        f = randelement(facets)
        return sort(vcat([f], lefschetz_coboundary(lc, f)))
    end

    for lc in (create_simplicial_rectangle(1, 2)[1], create_cubical_rectangle(2, 4)[1])
        tops = findall(isequal(lc.dim), lc.dimensions)
        for trial in 1:3
            fwd, bwd = height_transitions(lc)
            merge = [random_pair(lc)]
            mvf, dual = create_strong_mvf(lc, fwd, bwd, exitcells=Int[],
                                                entrancecells=Int[], merge=merge)
            valid = all_valid_partitions_merge(lc, fwd, bwd, merge)
            lab = toplabels(lc, mvf, tops)
            @test lab in valid
            @test all(refines(lab, v) for v in valid)
            @test any(issubset(merge[1], V) for V in mvf)
        end
    end

    for lc in (create_simplicial_rectangle(4, 4)[1], create_cubical_box(3, 3, 2)[1])
        for trial in 1:3
            fwd, bwd = random_boundary!(lc, height_transitions(lc)...)
            mvf0, dual0, info0 = create_strong_mvf(lc, fwd, bwd, details=true)
            bc = (exitcells=info0["exit cells"], entrancecells=info0["entrance cells"])
            isexit = falses(lc.ncells)
            isexit[bc.exitcells] .= true
            merge = filter(S -> !any(c -> isexit[c], S), [random_pair(lc) for _ in 1:3])
            mvf, dual = create_strong_mvf(lc, fwd, bwd; bc..., merge=merge)
            rmvf, rdual, rounds = ConleyDynamics.strong_mvf_reference(lc, fwd, bwd; bc...,
                                                                      merge=merge)
            @test canon(mvf) == canon(rmvf)
            @test canon(dual) == canon(rdual)
            @test all(S -> any(issubset(S, V) for V in mvf), merge)
            @test all(any(issubset(V, W) for W in mvf) for V in mvf0)
            @test mvf_is_strongly_accumulating(lc, mvf, dual, fwd, bwd; bc..., merge=merge)
            @test !mvf_is_strongly_accumulating(lc, mvf0, dual0, fwd, bwd; bc...,
                                                merge=merge) || all(S -> any(issubset(S, V)
                                                for V in mvf0), merge)
        end
    end

    lc = create_cubical_rectangle(3, 3)[1]
    fwd, bwd = random_boundary!(lc, height_transitions(lc)...)
    ex = ConleyDynamics.strong_mvf_boundary_cells(lc, fwd, bwd)[1]
    if !isempty(ex)
        @test_throws ArgumentError create_strong_mvf(lc, fwd, bwd, merge=[[ex[1]]])
    end

    # Local connectedness at dual cells (Lemma dualstar of the paper): For
    # c in V*, all toplexes of V in the star of c lie in one component of V.
    # The lemma requires exact transitions, but it holds in these examples.

    function components_of(lc, V)
        inV = falses(lc.ncells)
        inV[V] .= true
        comp = zeros(Int, lc.ncells)
        ncomp = 0
        for s in V
            comp[s] != 0 && continue
            ncomp += 1
            comp[s] = ncomp
            stack = [s]
            while !isempty(stack)
                x = pop!(stack)
                for y in vcat(lefschetz_boundary(lc, x), lefschetz_coboundary(lc, x))
                    if inV[y] && comp[y] == 0
                        comp[y] = ncomp
                        push!(stack, y)
                    end
                end
            end
        end
        return comp
    end
    function circle_vf(v::Vector{Float64})
        x, y = v / 5 - [2.0, 2.0]
        return [-y + x * (x^2 + y^2 - 1.0), x + y * (x^2 + y^2 - 1.0)]
    end
    examples = [example_ouroboros(),
                (create_simplicial_rectangle(8, 8, euclidean=true), x -> [x[1] - 4.0, 4.0 - x[2]]),
                (create_simplicial_rectangle(20, 20, euclidean=true), circle_vf)]

    # Connected merge sets lead to connected multivectors (part (b) of the
    # theorem), for transitions which are computed from vector fields

    for (ec, vf) in examples[2:3]
        lc = euclidean_to_lefschetz(ec)
        mvf0, dual0, info0 = create_planar_strong_mvf(ec, vf, details=true)
        isexit = falses(lc.ncells)
        isexit[info0["exit cells"]] .= true
        merge = filter(S -> !any(c -> isexit[c], S), [random_pair(lc) for _ in 1:6])
        mvf, dual = create_planar_strong_mvf(ec, vf, merge=merge)
        fwd = planar_transitions(ec, vf)
        bwd = planar_transitions(ec, x -> -vf(x))
        ok, rep = mvf_is_strongly_accumulating(lc, mvf, dual, fwd, bwd, merge=merge,
                                               details=true)
        @test ok
        @test isempty(rep["disconnected"])
        @test all(S -> any(issubset(S, V) for V in mvf), merge)
    end
    for (ec, vf) in examples
        lc = euclidean_to_lefschetz(ec)
        mvf, dual, info = create_planar_strong_mvf(ec, vf, details=true)
        for (j, D) in enumerate(dual)
            k = info["dual to mvf"][j]
            k == 0 && continue
            comp = components_of(lc, mvf[k])
            inV = falses(lc.ncells)
            inV[mvf[k]] .= true
            @test all(length(unique(comp[[t for t in lefschetz_openhull(lc, [c])
                                          if lc.dimensions[t] == 2 && inV[t]]])) <= 1 for c in D)
        end
    end

    # Conley indices over GF(2) and over the rationals: a saddle, and the
    # circle field with an attracting equilibrium and a repelling periodic
    # orbit, where the exit set has to be removed first

    for p in (2, 0)
        ec = create_simplicial_rectangle(2, 2, euclidean=true, p=p)
        lc = euclidean_to_lefschetz(ec)
        mvf, dual = create_planar_strong_mvf(ec, x -> [x[1] - 1.0, 1.0 - x[2]])
        center = findfirst(k -> lc.dimensions[k] == 0 && ec.coords[k][1] == [1.0, 1.0], 1:lc.ncells)
        V = mvf[findfirst(V -> center in V, mvf)]
        @test conley_index(lc, V) == [0, 1, 0]

        ec = create_simplicial_rectangle(20, 20, euclidean=true, p=p)
        mvf, dual, info = create_planar_strong_mvf(ec, circle_vf, details=true)
        @test length(info["exit cells"]) == 160
        @test isempty(info["entrance cells"])
        lcr, mvfr = remove_exit_set(ec, mvf)
        indices = sort([conley_index(lcr, M) for M in morse_sets(lcr, mvfr)])
        @test [1, 0, 0] in indices
        @test [0, 1, 1] in indices
        lcr0, mvfr0 = remove_exit_set(ec, create_planar_mvf(ec, circle_vf))
        @test indices == sort([conley_index(lcr0, M) for M in morse_sets(lcr0, mvfr0)])
    end
end

@testset "Periodic-orbit checks" begin
    #
    # An attracting periodic orbit, the circle of radius 1 around the
    # center of the square [-2,2]^2, which is mapped to [0,n]^2
    #
    function cycle_vf(v::Vector{Float64}, n)
        x, y = 4 * v / n - [2.0, 2.0]
        r2 = x^2 + y^2
        return [-y - x * (r2 - 1.0), x - y * (r2 - 1.0)]
    end

    # (C4) by brute force: a start state is bad if a maximal path through
    # non-good states exists, i.e., if the non-good states reachable from
    # it contain a dead end or a cycle

    function c4_bruteforce(out, mvf, C, labels, p, moU)
        nC = length(C)
        pos = Dict(v => k for (k, v) in enumerate(C))
        succ = [[pos[w] for w in out[C[k]] if haskey(pos, w)] for k in 1:nC]
        moinU = [all(c -> c in moU, setdiff(lefschetz_closure(ec, mvf[v]), mvf[v])) for v in C]
        fails = zeros(Int, p)
        for i in 1:p
            good(k, f) = moinU[k] || (labels[k] == mod1(i + 1, p) && f == 0)
            nxt(k, f) = [(l, (f == 1 || labels[l] == mod1(i + 2, p)) ? 1 : 0) for l in succ[k]]
            for k0 in 1:nC
                labels[k0] == i || continue
                good(k0, 0) && continue
                seen = Set([(k0, 0)])
                stack = [(k0, 0)]
                bad = false
                while !isempty(stack)
                    s = pop!(stack)
                    isempty(succ[s[1]]) && (bad = true)
                    for t in nxt(s...)
                        good(t...) && continue
                        t in seen || (push!(seen, t); push!(stack, t))
                    end
                end
                if !bad
                    # cycle in the non-good states reachable from the start
                    S = collect(seen)
                    idx = Dict(s => j for (j, s) in enumerate(S))
                    indeg = zeros(Int, length(S))
                    for s in S, t in nxt(s...)
                        haskey(idx, t) && (indeg[idx[t]] += 1)
                    end
                    queue = [j for j in 1:length(S) if indeg[j] == 0]
                    removed = 0
                    while !isempty(queue)
                        j = pop!(queue); removed += 1
                        for t in nxt(S[j]...)
                            haskey(idx, t) || continue
                            indeg[idx[t]] -= 1
                            indeg[idx[t]] == 0 && push!(queue, idx[t])
                        end
                    end
                    bad = removed < length(S)
                end
                bad && (fails[i] += 1)
            end
        end
        return fails
    end

    n = 20
    ec = create_simplicial_rectangle(n, n, euclidean=true)
    vf = v -> cycle_vf(v, n)
    mvf, dual = create_planar_strong_mvf(ec, vf)
    ms = morse_sets(ec, mvf)
    ci = [conley_index(ec, M) for M in ms]
    @test [1, 1, 0] in ci
    @test [0, 0, 1] in ci
    full = mvf_with_singletons(ec, mvf)
    @test full[1:length(mvf)] == mvf
    @test sort(reduce(vcat, full)) == collect(1:ec.ncells)
    out, c2m, mouths = mvf_digraph(ec, full)
    C = sort(unique(c2m[ms[findfirst(isequal([1, 1, 0]), ci)]]))
    Ucells = reduce(vcat, full[C])
    moU = Set(setdiff(lefschetz_closure(ec, Ucells), Ucells))

    for p in (3, 4, 6)
        lab = cyclic_partition_angular(ec, full, C, p; center=[n / 2, n / 2])
        r = check_cyclic_conditions(ec, full, C, lab, p; out=out, mouths=mouths)
        @test r.C2 && r.C3 && r.C4
        @test r == check_cyclic_conditions(ec, full, C, lab, p)
        @test sum(r.sizes) == length(C)
        labr = cyclic_partition_angular(ec, full, C, p; center=[n / 2, n / 2], clockwise=true)
        rr = check_cyclic_conditions(ec, full, C, labr, p; out=out, mouths=mouths)
        @test !rr.C4
        @test rr.fails == c4_bruteforce(out, full, C, labr, p, moU)
    end

    # Random labels: agreement of (C4) with the brute-force check

    rngstate = Ref(UInt64(17))
    function rand01()
        rngstate[] += 0x9e3779b97f4a7c15
        z = rngstate[]
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        z = z ⊻ (z >> 31)
        return (z >> 11) / 2.0^53
    end
    for p in (3, 4, 5)
        base = cyclic_partition_angular(ec, full, C, p; center=[n / 2, n / 2])
        for trial in 1:3
            lab = [rand01() < 0.1 ? 1 + floor(Int, rand01() * p) : l for l in base]
            r = check_cyclic_conditions(ec, full, C, lab, p; out=out, mouths=mouths)
            @test r.fails == c4_bruteforce(out, full, C, lab, p, moU)
        end
    end
    @test_throws ArgumentError check_cyclic_conditions(ec, full, C, ones(Int, length(C)), 2)

    # The direction test: positive on the periodic orbit, -1 at the
    # repelling focus, and exactly 1 for a constant vector field

    vals, diam = check_brick_lyapunov(ec, full, C, vf)
    valm, _ = check_brick_lyapunov(ec, full, C, vf; method=:mean)
    @test minimum(vals) > 0.9
    @test all(vals .>= valm .- 1e-12)
    @test maximum(diam) <= sqrt(2) + 1e-12
    Ceq = sort(unique(c2m[ms[findfirst(isequal([0, 0, 1]), ci)]]))
    @test any(isequal(-1.0), check_brick_lyapunov(ec, full, Ceq, vf)[1])
    valc, _ = check_brick_lyapunov(ec, full, C, v -> [1.0, 0.0]; method=:mean)
    @test all(x -> isapprox(x, 1.0), valc)

    # Local cycles: a short cycle is found, a cycle around the circle is not

    outs = [[2], [3], [1], [5], [6], [4]]
    param = [0.10, 0.11, 0.12, 0.0, 0.33, 0.66]
    @test cyclic_local_cycles(outs, collect(1:6), param; width=0.05) == [[1, 2, 3]]
    @test isempty(cyclic_local_cycles(outs, collect(4:6), param[4:6]; width=0.05))
end
