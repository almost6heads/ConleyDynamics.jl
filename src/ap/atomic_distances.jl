export atomic_distances

"""
    atomic_distances(lc::AbstractComplex, ap::Vector{Vector{Vector{Int}}}; p::Int=99991)

Compute the atomic-refinement adjacency matrix of `ap` (typically the
output of `construct_ap_space`).

Returns a sparse matrix `A` (in the `ConleyDynamics` sparse format, over
the prime field `p`) with `A[i,j] == 1` iff `ap[i]` is an atomic refinement
of `ap[j]` (`is_atomic_refinement(lc, ap[i], ap[j])`). This is the
expensive step of the pipeline (checking every pair of elements at
adjacent lengths), and is what `stratum_adjacency` uses to determine
coverage between Morse-vector strata.
"""
function atomic_distances(lc::AbstractComplex,
                          ap::Vector{Vector{Vector{Int}}};
                          p::Int=99991)
    #
    # Compute the matrix of atomic distances
    #
    nap   = length(ap)
    aplen = map(t -> mvf_length(lc,t), ap)
    Drow  = Vector{Int}()
    Dcol  = Vector{Int}()
    Dval  = Vector{Int}()

    aplenval = sort(unique(aplen))
    ch = Channel{Vector{Int}}(Inf)

    for k = 1:length(aplenval)-1
        aplen1 = aplenval[k+1]
        aplen2 = aplenval[k]
        ind1 = findall(t -> t == aplen1, aplen)
        ind2 = findall(t -> t == aplen2, aplen)
        Threads.@threads for m1 in ind1
            for m2 in ind2
                if is_atomic_refinement(lc, ap[m1], ap[m2])
                    put!(ch, [m1, m2])
                end
            end
        end
        @info "Finished level " k
    end

    close(ch)
    chpairs = collect(ch)
    for pair in chpairs
        push!(Dcol, pair[1])
        push!(Drow, pair[2])
        push!(Dval, 1)
    end

    return sparse_from_lists(nap,nap,p,Int(0),Int(1),Drow,Dcol,Dval)
end
