export sparse_minor

"""
    smp = sparse_minor(sm::SparseMatrix, rvec::Vector{Int}, cvec::Vector{Int})

Create sparse submatrix by specifying the desired row and column indices.
"""
function sparse_minor(sm::SparseMatrix, rvec::Vector{Int}, cvec::Vector{Int})
    #
    # Create sparse submatrix by specifying the desired row and column indices
    #

    # Create the new dimensions

    newnr = length(rvec)
    newnc = length(cvec)

    # Determine the new nonzero entry lists

    r = Vector{Int}([])
    c = Vector{Int}([])
    vals = Vector{typeof(sm.zero)}([])

    for m=1:newnc
        colm = sparse_get_nz_column(sm, cvec[m])
        for centry in colm
            rowind = findall(x -> x==centry, rvec)
            for k in rowind
                push!(r,k)
                push!(c,m)
                push!(vals,sm[rvec[k],cvec[m]])
            end
        end
    end

    # Create and return new sparse matrix

    msm = sparse_from_lists(newnr, newnc, sm.char, sm.zero, sm.one, r, c, vals)
    return msm
end

