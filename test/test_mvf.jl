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
        mvf, dual, info = create_strong_mvf(lc, fwd, bwd)
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
    mvfw, dualw, infow = create_planar_strong_mvf(ec, vf)
    @test mvfw == mvf
    @test dualw == dual

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
    mvf3w, dual3w, info3w = create_spatial_strong_mvf(lc3, coords3N, v3)
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
