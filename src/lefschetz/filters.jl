export create_random_filter
export filter_shallow_pairs
export filter_induced_mvf

"""
    create_random_filter(lc::AbstractComplex)

Create a random injective filter on a Lefschetz complex.

The function creates a random injective filter on a Lefschetz
complex. The filter is created by assigning integers to cell
groups, increasing with dimension. Within each dimension the
assignment is random, but all filter values of cells of dimension
`k` are less than all filter values of cells with dimension `k+1`.
The function returns the filter as `Vector{Int}`, with indices
corresponding to the cell indices in the Lefschetz complex.
"""
function create_random_filter(lc::AbstractComplex)
    #
    # Create a random injective filter on a Lefschetz complex
    #

    # Initialize the filter and the dimension

    lcdim = lc.dim
    phi = Vector{Int}(undef, lc.ncells)
    voffset = 0

    # Loop through the dimensions

    for kdim in 0:lcdim
        kindices = findall(x -> x==kdim, lc.dimensions)
        klen = length(kindices)
        if klen > 0
            phi[kindices] = randperm(klen) .+ voffset
        end
        voffset = voffset + klen
    end

    # Return the filter

    return phi
end

"""
    filter_shallow_pairs(lc::AbstractComplex, phi)

Find all shallow pairs for a filter.

This function finds all shallow pairs for the filter `phi`. These
are face-coface pairs `(x,y)` whose dimensions differ by one, and
such that `y` has the smallest filter value on the coboundary of `x`,
and `x` has the largest filter value on the boundary of `y`.

If the filter is injective, these pairs give rise to a Forman
vector field on the underlying Lefschetz complex. For noninjective
filters this is not true in general.
"""
function filter_shallow_pairs(lc::AbstractComplex, phi)
    #
    # Find all shallow pairs for a filter
    #
    phish = Vector{Vector{Int}}([])

    for x=1:lc.ncells
        phicbdx = phi_cbd(lc, phi, x)
        for y in phicbdx
            if x in phi_bd(lc, phi, y)
                push!(phish, [x,y])
            end
        end
    end

    return phish
end

"""
    filter_induced_mvf(lc::AbstractComplex, phi)

Compute the multivector field induced by a filter.

This function returns the smallest multivector field which has the
property that every shallow pair is contained in a multivector.
For injective filters this is a Forman vector field, but in the
noninjective case it can be a general multivector field.
"""
function filter_induced_mvf(lc::AbstractComplex, phi)
    #
    # Compute the multivector field induced by a filter
    #
    phish = filter_shallow_pairs(lc, phi)
    phimvf = create_mvf_hull(lc, phish)
    return phimvf
end

###############################
#                             #
#     Auxiliary functions     #
#                             #
###############################

function phi_bd(lc::AbstractComplex, phi, y::Int)
    #
    # Compute the phi boundary of a cell y
    #
    bdy = lefschetz_boundary(lc, y)
    phi_bd_cells = Vector{Int}([])
    
    if length(bdy) > 0
        maxphival = maximum(phi[bdy])
        for x in bdy
            if phi[x] == maxphival
                push!(phi_bd_cells, x)
            end
        end
    end

    return phi_bd_cells
end

function phi_cbd(lc::AbstractComplex, phi, x::Int)
    #
    # Compute the phi coboundary of a cell x
    #
    cbdx = lefschetz_coboundary(lc, x)
    phi_cbd_cells = Vector{Int}([])
    
    if length(cbdx) > 0
        minphival = minimum(phi[cbdx])
        for y in cbdx
            if phi[y] == minphival
                push!(phi_cbd_cells, y)
            end
        end
    end

    return phi_cbd_cells
end

