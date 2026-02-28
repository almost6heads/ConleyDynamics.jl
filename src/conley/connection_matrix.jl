export connection_matrix

"""
    connection_matrix(lc::LefschetzComplex, mvf::CellSubsets;
                      [algorithm::String],
                      [returnbasis::Bool])

Compute a connection matrix for the multivector field `mvf` on the
Lefschetz complex `lc` over the field associated with the Lefschetz
complex boundary matrix.

The function returns an object of type `ConleyMorseCM`. If the optional
argument `returnbasis::Bool=true` is given, then the function also returns
a dictionary which gives the basis for the connection matrix columns in
terms of the original labels. Finally, it is possible to invoke the
connection matrix computation with one of two different algorithms,
by passing the optional argument `algorithm::String`:

* `algorithm = "DLMS24"`
* `algorithm = "DHL26"`

If the flag `returnbasis::Bool=true` is given the function automatically
chooses `algorithm = "DLMS24"`.
"""
function connection_matrix(lc::LefschetzComplex, mvfarg::CellSubsets;
                           algorithm::String="DLMS24",
                           returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Convert the multivector field to Int if necessary

    if mvfarg isa Vector{Vector{String}}
        mvf = convert_cellsubsets(lc, mvfarg)
    else
        mvf = mvfarg
    end

    # Copy boundary and label data
    
    lcdim     = lc.dim
    bndmatrix = lc.boundary
    labels    = lc.labels
    celldims  = lc.dimensions

    # Find an admissible order

    adorder, adbnd, psetvec, scc = admissible_order(bndmatrix, mvf)

    # Compute the connection matrix in reordered form

    cmRedP = deepcopy(adbnd)
    if returnbasis
        cmMatrP, cmColsP, cmBasisP = cm_reduce_dlms24!(cmRedP, psetvec; returnbasis=true)
    elseif algorithm == "DHL26"
        cmMatrP, cmColsP = cm_reduce_dhl26!(cmRedP, psetvec)
    elseif algorithm == "DLMS24"
        cmMatrP, cmColsP = cm_reduce_dlms24!(cmRedP, psetvec)
    else
        error("Invalid value of the algorithm flag!")
    end

    # Compute first few return variables
    
    cmMatr      = cmMatrP
    cmCols      = adorder[cmColsP]
    cmPoset     = psetvec[cmColsP]
    cmLabels    = labels[cmCols]

    # If desired, determine the basis dictionary

    if returnbasis
        basisdict = Dict{String,Vector{String}}()
        for k=1:length(cmBasisP)
            basisStr = labels[cmCols[k]]
            basisVec = Vector{String}()
            for m=1:length(cmBasisP[k])
                push!(basisVec,labels[adorder[cmBasisP[k][m]]])
            end
            basisdict[basisStr] = basisVec
        end
    end

    # Determine the Morse sets and their Poincare polynomials
    # using the original labels

    poset_indices = sort(union(cmPoset))
    cmMorseSets = Vector{Vector{String}}()
    for k = 1:length(poset_indices)
        push!(cmMorseSets,labels[scc[poset_indices[k]]])
    end

    cmPoincare  = Vector{Vector{Int}}()
    for k=1:length(poset_indices)
        ppolytemp = fill(Int(0),lcdim+1)
        for j=1:length(cmPoset)
            if cmPoset[j] == poset_indices[k]
                ppolytemp[celldims[cmCols[j]]+1] += 1
            end
        end
        push!(cmPoincare,ppolytemp)
    end

    # Renumber the poset indices in consecutive order

    renumber_poset!(cmPoset)

    # Construct the Conley complex as a Lefschetz complex

    cc_ncells  = length(cmLabels)
    cc_bnd     = cmMatr
    cc_labels  = cmLabels
    cc_dims    = Vector{Int}([lc.dimensions[lc.indices[cmLabels[k]]]
                              for k in 1:cc_ncells])

    # Make sure the cells are ordered by dimension

    cdorder    = sortperm(cc_dims)
    cc_labelsO = cc_labels[cdorder]
    cc_dimsO   = cc_dims[cdorder]
    cc_bndO    = sparse_permute(cc_bnd, cdorder, cdorder)

    # Create the new Lefschetz complex and return it

    cmCC = LefschetzComplex(cc_labelsO, cc_dimsO, cc_bndO)
 
    # Return the connection matrix information
    
    cm = ConleyMorseCM{typeof(cmMatr.zero)}(
                cmMatr, cmCols, cmPoset, cmLabels, cmMorseSets, cmPoincare, cmCC)

    if returnbasis
        return cm, basisdict
    else
        return cm
    end
end

###############################################################################

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

###############################################################################

## export admissible_order
##
## """
##     admissible_order(bndmatrix::SparseMatrix, mvf::Vector{Vector{Int}})
##
## Find an admissible order based on the boundary matrix and the multivector field.
##
## The vector `sccnumber` contains the strongly connected component for each column,
## and the matrix `admibnd` is the reordered boundary matrix corresponding to the
## order given by `admiorder`.
##
## # Return values:
## * `admiorder`: Admissible order
## * `admibnd`: Reordered boundary matrix
## * `sccnumber`: Strongly connected component index of column
## * `scc`: Strongly connected components
## """
function admissible_order(bndmatrix::SparseMatrix, mvf::Vector{Vector{Int}})
    #
    # Find an admissible order based on the boundary matrix and the multivector
    # field
    #
    
    # Create the digraph based on the boundary matrix

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(bndmatrix)

    # Add cycles for the multivectors to the lists c and r

    lenmvf = length(mvf)
    for k = 1:lenmvf
        mv  = mvf[k]
        lmv = length(mv)
        for m = 1:lmv-1
            push!(c,mv[m])
            push!(r,mv[m+1])
        end
        push!(c,mv[lmv])
        push!(r,mv[1])
    end

    # Create the edge list and the digraph

    edgelist = Edge.([(c[k],r[k]) for k in 1:length(c)])
    dg = SimpleDiGraph(edgelist)

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    for k=1:length(scc)
        sort!(scc[k])
    end

    # Create the admissible order and the vector of scc numbers

    admiorder = Vector{Int}()
    sccnumber = Vector{Int}()

    scccount = 0

    for k=1:length(scc)
        scccount += 1
        append!(admiorder,scc[k])
        append!(sccnumber,scc[k] .* 0 .+ scccount)
    end

    # Create the reordered boundary matrix

    admibnd = sparse_permute(bndmatrix, admiorder, admiorder)

    # Return the results

    return admiorder, admibnd, sccnumber, scc
end

###############################################################################

## export renumber_poset!
##
## """
##     renumber_poset!(poset::Vector{Int})
##
## Renumber the poset given by the increasing integer vector `poset`.
## """
function renumber_poset!(poset::Vector{Int})
    #
    # Renumber the poset vector
    #
    
    # Check whether the vector is increasing

    len_poset = length(poset)
    is_increasing = true

    for i = 1:len_poset-1
        if poset[i] > poset[i+1]
            is_increasing = false
        end
    end

    if !is_increasing
        error("Poset must be increasing!")
    end

    # Renumber the poset

    current_value = poset[1]
    new_value = 1

    for i = 1:len_poset
        if poset[i] == current_value
            poset[i] = new_value
        else
            current_value = poset[i]
            new_value += 1
            poset[i] = new_value
        end
    end
end

###############################################################################

