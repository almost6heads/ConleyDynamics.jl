export sparse_minor

"""
    smp = sparse_minor(sm::SparseMatrix, rvec::Vector{Int}, cvec::Vector{Int})

Create sparse submatrix by specifying the desired row and column indices.
"""
function sparse_minor(sm::SparseMatrix, rvec::Vector{Int}, cvec::Vector{Int})
    #
    # Create sparse submatrix by specifying the desired row and column indices
    #
    @assert length(rvec)==length(unique(rvec)) "Rows cannot repeat!"
    @assert length(cvec)==length(unique(cvec)) "Columns cannot repeat!"

    rdict = Dict{Int,Int}([rvec[k] => k for k in 1:length(rvec)])
    cdict = Dict{Int,Int}([cvec[k] => k for k in 1:length(cvec)])
    rset  = Set(rvec)
    cset  = Set(cvec)

    # Extract lists from the sparse matrix

    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(sm)

    # Filter out the correct part

    rcv    = [(r[k],c[k],v[k]) for k in 1:length(r)]
    newrcv = filter(v -> (v[1] in rset) && (v[2] in cset), rcv)
    rs     = getindex.(newrcv, 1)
    cs     = getindex.(newrcv, 2)
    newv   = getindex.(newrcv, 3)
    newr   = [rdict[m] for m in rs]
    newc   = [cdict[m] for m in cs]

    # Create the new dimensions

    newnr = length(rvec)
    newnc = length(cvec)

    # Create and return new sparse matrix

    msm = sparse_from_lists(newnr, newnc, sm.char, sm.zero, sm.one, newr, newc, newv)
    return msm
end

