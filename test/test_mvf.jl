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
