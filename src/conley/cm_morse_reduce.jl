export cm_reduce_hms21
export cm_reduce_pmorse26

"""
    cm_reduce_hms21(matrix::SparseMatrix, psetvec::Vector{Int})

Compute the connection matrix.

This function uses the algorithm from the paper Harker, Mischaikow,
Spendlove (Journal of Applied and Computational Topology, 2021).
Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
"""
function cm_reduce_hms21(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Find the initial matching

    current_matrix  = matrix
    current_psetvec = psetvec
    current_cvector = Vector{Int}(1:matrix.nrow)

    current_m = matching(current_matrix, current_psetvec)
    
    # If this matching produced fewer critical cells than the
    # number we started out with, we can determine the smaller
    # boundary matrix

    while length(current_m.critical) < current_matrix.nrow

        # Determine the new matrix, whose columns are indexed
        # by the critical cells

        new_matrix, new_cvector = morse_delta(current_matrix, current_m)

        current_matrix  = new_matrix
        current_psetvec = current_psetvec[new_cvector]
        current_cvector = current_cvector[new_cvector]

        # Try to find another matching

        current_m = matching(current_matrix, current_psetvec)
    end

    # Return the result

    return current_matrix, current_cvector
end

"""
    cm_reduce_pmorse26(matrix::SparseMatrix, psetvec::Vector{Int})

Compute the connection matrix.

This function uses a parallelized Morse matching algorithm.
Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
"""
function cm_reduce_pmorse26(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Find the initial matching

    current_matrix  = matrix
    current_psetvec = psetvec
    current_cvector = Vector{Int}(1:matrix.nrow)

    current_m = parallel_matching(current_matrix, current_psetvec)
    
    # If this matching produced fewer critical cells than the
    # number we started out with, we can determine the smaller
    # boundary matrix

    while length(current_m.critical) < current_matrix.nrow

        # Determine the new matrix, whose columns are indexed
        # by the critical cells

        new_matrix, new_cvector = parallel_morse_delta(current_matrix, current_m)

        current_matrix  = new_matrix
        current_psetvec = current_psetvec[new_cvector]
        current_cvector = current_cvector[new_cvector]

        # Try to find another matching

        current_m = parallel_matching(current_matrix, current_psetvec)
    end

    # Return the result

    return current_matrix, current_cvector
end

######################################################################################

struct Matching
    critical::Set{Int}
    sources::Set{Int}
    targets::Set{Int}
    arrows::Dict{Int,Int}
end


function matching(matrix::SparseMatrix, filter::Vector{Int})
    #
    # Find a matching for a given boundary matrix, based
    # on a provided filter for the cells.
    #

    @assert filter==sort(filter) "The filter has to be increasing!"
    @assert length(filter)==matrix.nrow "Wrong filter size!"
    @assert matrix.ncol==matrix.nrow "Wrong matrix size!"

    # Initialize the return variables

    critical = Set{Int}()
    sources  = Set{Int}()
    targets  = Set{Int}()
    arrows   = Dict{Int,Int}()

    # Find the filter values

    fvalues = sort(unique(filter))

    # Loop through the values and create the matching

    for fval in fvalues
        # Find the subcomplex and a partial matching

        findices = findall(t -> t==fval, filter)
        cp, ap = partial_matching(matrix, findices)
        
        # Update the return variables

        sp = Vector{Int}()
        tp = Vector{Int}()
        for (s,t) in ap
            push!(sp, s)
            push!(tp, t)
        end

        union!(critical, cp)
        union!(sources, sp)
        union!(targets, tp)
        merge!(arrows, ap)
    end

    return Matching(critical, sources, targets, arrows)
end


function parallel_matching(matrix::SparseMatrix, filter::Vector{Int})
    #
    # Find a matching for a given boundary matrix, based
    # on a provided filter for the cells. This function is
    # parallelized.
    #

    @assert filter==sort(filter) "The filter has to be increasing!"
    @assert length(filter)==matrix.nrow "Wrong filter size!"
    @assert matrix.ncol==matrix.nrow "Wrong matrix size!"

    # Initialize the return variables

    critical = Set{Int}()
    sources  = Set{Int}()
    targets  = Set{Int}()
    arrows   = Dict{Int,Int}()

    # Find the filter values

    fvalues = sort(unique(filter))

    # Loop through the values and create the matching

    critical_channel = Channel{Vector{Int}}(Inf)
    arrow_channel = Channel{Dict{Int,Int}}(Inf)

    Threads.@threads for fval in fvalues

        # Find the subcomplex and a partial matching

        findices = findall(t -> t==fval, filter)
        cp, ap = partial_matching(matrix, findices)
        put!(critical_channel, cp)
        put!(arrow_channel, ap)
    end

    close(critical_channel)
    close(arrow_channel)
    critical_chvec = collect(critical_channel)
    arrow_chvec = collect(arrow_channel)

    # Update the return variables

    for cp in critical_chvec
        union!(critical, cp)
    end

    for ap in arrow_chvec
        sp = Vector{Int}()
        tp = Vector{Int}()
        for (s,t) in ap
            push!(sp, s)
            push!(tp, t)
        end
        union!(sources, sp)
        union!(targets, tp)
        merge!(arrows, ap)
    end

    # Return the results

    return Matching(critical, sources, targets, arrows)
end


function partial_matching(matrix::SparseMatrix, mind::Vector{Int})
    #
    #  Find a partial matching in the minor `matrix` determined
    #  by the columns in `mind`
    #

    @assert length(mind)>0 "I need something to work with!"
    critical = Vector{Int}([])
    arrows   = Dict{Int,Int}()

    # Create a boundary dictionary for faster manipulation
    # from the given sparse matrix. The dictionary only
    # contains the cells in `mind`, and for every cell it
    # points to its boundary within `mind`.

    mindset = Set{Int}(mind)
    bnddict = Dict{Int, Set{Int}}(c => Set{Int}() for c in mindset)
    for k in mind
        if length(matrix.columns[k]) > 0
            bndset = intersect(Set(matrix.columns[k]), mindset)
            if length(bndset) > 0
                bnddict[k] = bndset
            end
        end
    end
    for k in keys(bnddict)
        if length(bnddict[k]) == 0
            delete!(bnddict, k)
        end
    end

    # Create a second dictionary which records the coboundaries
    # of cells
    
    cbddict = Dict{Int, Set{Int}}(c => Set{Int}() for c in mindset)
    for (k, bndk) in bnddict
        for m in bndk
            push!(cbddict[m], k)
        end
    end
    for k in keys(cbddict)
        if length(cbddict[k]) == 0
            delete!(cbddict, k)
        end
    end

    tomatch = mindset
    removed = Set{Int}()

    # Collect all cells with boundary of size 1

    queue = Vector{Int}()
    for k in keys(bnddict)
        if length(bnddict[k]) == 1
            push!(queue, k)
        end
    end

    # Find the matching

    while length(tomatch) > 0

        # Find all coreduction pairs

        while !isempty(queue)
            t = popfirst!(queue)
            
            # Skip if t was already removed, as target of someone else
            (t in removed || length(get(bnddict, t, Set())) != 1) && continue
    
            # Identify the source
            s = first(bnddict[t])
            arrows[s] = t
            
            # Now we have to remove s and t
            for c in (t, s)
                c in removed && error("s error")
                push!(removed, c)
                delete!(tomatch, c)
                
                # Update the boundary of every cell which has c in
                # its boundary
                if haskey(cbddict, c)
                    for cc in cbddict[c]
                        cc in removed && continue

                        # Remove c from the boundary of bc
                        delete!(bnddict[cc], c)
                        
                        # If cc now has a boundary of size 1,
                        # add it to the queue
                        if length(bnddict[cc]) == 1
                            push!(queue, cc)
                        end
                    end
                end
                
                # Clean up the boundary dictionary entry for
                # the removed cell
                delete!(bnddict, c)
            end
        end

        # Find one element with zero boundary

        if length(tomatch) > 0
            foundit = false
            for k in tomatch
                if length(intersect(get(bnddict, k, Set{Int}()), tomatch)) == 0

                    # We found a critical cell!
                    push!(critical, k)
                    foundit = true
                    
                    # Update the boundary of every cell which has k in
                    # its boundary
                    if haskey(cbddict, k)
                        for ck in cbddict[k]
                            ck in removed && continue

                            # Remove k from the boundary of ck
                            delete!(bnddict[ck], k)

                            # If ck now has a boundary of size 1,
                            # add it to the queue
                            if length(bnddict[ck]) == 1
                                push!(queue, ck)
                            end
                        end
                    end

                    # Update the cell lists and coboundary dictionary
                    push!(removed, k)
                    delete!(tomatch, k)
                    delete!(bnddict, k)
                    break
                end
            end
            @assert foundit "There should have been a critical cell!"
        end
    end

    return critical, arrows
end


function morse_gamma(bnd::SparseMatrix, m::Matching, xin::SparseMatrix)
    #
    # This function implements the Gamma algorithm from the paper
    # by Harker, Mischaikow, Mrozek, Nanda (FOCM, 2014). Note that the
    # version in the paper by Harker, Mischaikow, Spendlove (JACT, 2021)
    # is not correct for general fields.
    #
    @assert bnd.nrow==bnd.ncol "Wrong matrix dimensions!"
    @assert xin.nrow==bnd.ncol "Matrix-vector product not defined!"
    @assert xin.ncol==1 "The second argument is not a vector!"

    # Initialization

    p    = xin.char
    n    = xin.nrow
    x    = fill(xin.zero, n)
    c    = fill(xin.zero, n)
    qset = sort(collect(m.sources))

    for j in 1:length(xin.columns[1])
        x[xin.columns[1][j]] = xin.entries[1][j]
    end

    # Main loop

    while true
        qi = findlast(t -> !iszero(t), x[qset])
        qi == nothing && break

        Q     = qset[qi]
        K     = m.arrows[Q]
        omega = scalar_multiply(-x[Q], scalar_inverse(bnd[Q,K],p), p)
        c[K]  = scalar_add(c[K], omega, p)

        for j in 1:length(bnd.columns[K])
            jr = bnd.columns[K][j]
            vr = bnd.entries[K][j]
            x[jr] = scalar_add(x[jr], scalar_multiply(omega, vr, p), p)
        end
    end

    # Construct the result as sparse matrix

    rl = Vector{Int}()
    cl = Vector{Int}()
    vl = Vector{typeof(bnd.zero)}()
 
    for j = 1:n
        if !iszero(c[j])
            push!(rl, j)
            push!(cl, 1)
            push!(vl, c[j])
        end
    end

    cs = sparse_from_lists(n,1,p,bnd.zero,bnd.one,rl,cl,vl)
    return cs
end


function morse_delta(bnd::SparseMatrix, match::Matching)
    #
    # This function computes the reduced boundary matrix
    # induced by the specified matching. In essence, we
    # compute for each unit vector v corresponding to a
    # critical cell the vector
    #
    #    psi(delta(phi(x))),
    #
    # as defined in (5) on page 166 in Harker, Mischaikow,
    # Mrozek, Nanda (FOCM, 2014). Notice that their formula
    # can be simplified using the identity
    #
    #    gamma delta gamma = -gamma.
    #
    # Thus, we actually compute
    #
    #    delta v + delta(gamma(delta v)),
    #
    # where delta denotes the boundary operator.
    #
    @assert bnd.nrow==bnd.ncol "Wrong matrix dimensions!"

    # Initialization

    p = bnd.char
    cproj = sort(collect(match.critical))
    lproj = length(cproj)
    brows = Vector{Int}(1:bnd.nrow)
    rlist = Vector{Int}()
    clist = Vector{Int}()
    vlist = Vector{typeof(bnd.zero)}()

    for ell = 1:length(cproj)

        # Construct the column vector

        k    = cproj[ell]
        bndk = sparse_minor(bnd, brows, [k])
        psi  = bndk + bnd * morse_gamma(bnd, match, bndk)

        # Extract the nonzero entries in the critical part

        for j = 1:lproj
            if !(psi[cproj[j],1] == 0)
                push!(rlist, j)
                push!(clist, ell)
                push!(vlist, psi[cproj[j],1])
            end
        end
    end

    # Construct the new boundary matrix

    D = sparse_from_lists(lproj, lproj, bnd.char,
               bnd.zero, bnd.one, rlist, clist, vlist)
    return D, cproj
end


function parallel_morse_delta(bnd::SparseMatrix, match::Matching)
    #
    # This function is just a parallelized version of the
    # one right above.
    #
    @assert bnd.nrow==bnd.ncol "Wrong matrix dimensions!"

    # Initialization

    p = bnd.char
    cproj = sort(collect(match.critical))
    lproj = length(cproj)
    brows = Vector{Int}(1:bnd.nrow)
    rlist = Vector{Int}()
    clist = Vector{Int}()
    vlist = Vector{typeof(bnd.zero)}()

    column_channel = Channel{Tuple{Int,SparseMatrix}}(Inf)

    Threads.@threads for ell = 1:length(cproj)

        # Construct the column vector

        k    = cproj[ell]
        bndk = sparse_minor(bnd, brows, [k])
        psi  = bndk + bnd * morse_gamma(bnd, match, bndk)
        put!(column_channel, (ell, psi))
    end

    close(column_channel)
    column_chvec = collect(column_channel)

    # Extract the nonzero entries in the critical part

    for (ell, psi) in column_chvec
        for j = 1:lproj
            if !(psi[cproj[j],1] == 0)
                push!(rlist, j)
                push!(clist, ell)
                push!(vlist, psi[cproj[j],1])
            end
        end
    end

    # Construct the new boundary matrix

    D = sparse_from_lists(lproj, lproj, bnd.char,
               bnd.zero, bnd.one, rlist, clist, vlist)
    return D, cproj
end

