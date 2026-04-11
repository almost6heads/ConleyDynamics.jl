export remove_exit_set

"""
    remove_exit_set(lc::AbstractComplex, mvf::CellSubsets)

Exit set removal for a multivector field on a Lefschetz subcomplex.

It is assumed that the Lefschetz complex `lc` is a topological manifold
and that `mvf` contains a multivector field that is created via either
`create_planar_mvf` or `create_spatial_mvf`. The function identifies
cells on the boundary at which the flows exits the region covered by
the Lefschetz complex. If this exit set is closed, we have found an 
isolated invariant set and the function returns a Lefschetz complex `lcr`
restricted to it, as well as the restricted multivector field `mvfr`.
If the exit set is not closed, a warning is displayed and the function
returns the restricted Lefschetz complex and multivector field obtained
by removing the closure of the exit set. In the latter case, unexpected
results might be obtained.
"""
function remove_exit_set(lc::AbstractComplex, mvf::CellSubsets)
    #
    # Exit set removal for a multivector field on a Lefschetz subcomplex
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
        return_as_labels = false
    else
        mvfI = convert_cellsubsets(lc, mvf)
        return_as_labels = true
    end

    # Find the manifold boundary

    lcbnd = manifold_boundary(lc)

    # Find all critical cells on the boundary

    mvfT = deepcopy(mvfI)
    Threads.@threads for k in eachindex(mvfT)
        mv = mvfT[k]
        if length(mv) > 1
            intersect!(mvfT[k], lcbnd)
        else
            mvfT[k] = Vector{Int}()
        end
    end
    criticalcells = setdiff(lcbnd, unique(reduce(vcat,mvfT)))

    # Below is the older and slower code for finding
    # the critical cells. This should be removed eventually.
    # 
    # iscritical = fill(Int(1), length(lcbnd)) 
    # for mv in mvfI
    #     if length(mv) > 1
    #         mvbnd = intersect(mv, lcbnd)
    #         if length(mvbnd) > 0
    #             for k in mvbnd
    #                 ki = findfirst(x -> x==k, lcbnd)
    #                 iscritical[ki] = 0
    #             end
    #         end
    #     end
    # end   
    # criticalindices = findall(x -> x==1, iscritical)
    # criticalcells = lcbnd[criticalindices]

    # Check whether the critical boundary cells form a closed set

    if !lefschetz_is_closed(lc, criticalcells)
        closurecc = lefschetz_closure(lc, criticalcells)
        println(" ")
        println(" Warning! The exit set is not closed!")
        println(" Automatically using the closure..")
        println(" ")
    else
        closurecc = criticalcells
    end

    lcsub = setdiff(1:lc.ncells, closurecc)
    lcr, mvfr = restrict_dynamics(lc, mvf, lcsub)

    # Return the results

    return lcr, mvfr
end

