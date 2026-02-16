export planar_nontransverse_edges

"""
    planar_nontransverse_edges(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf;
                               npts::Int=100)

Find all edges of a planar Lefschetz complex which are not flow transverse.

The Lefschetz complex is given in `lc`, the coordinates of all vertices of the complex
in `coords`, and the vector field is specified in `vf`. The optional parameter `npts`
determines how many points along an edge are evaluated for the transversality check.
The function returns a list of nontransverse edges as `Vector{Int}`, which contains
the edge indices.
"""
function planar_nontransverse_edges(lc::LefschetzComplex, coords::Vector{<:Vector{<:Real}}, vf;
                                    npts::Int=100)
    #
    # Find all edges of a planar Lefschetz complex which are not flow transverse
    #

    # Find a list of all edges

    edges = findall(isequal(1), lc.dimensions)
    nedges = length(edges)

    # Check transversality for every edge
    
    is_transverse = fill(false,nedges)
    Threads.@threads for k in eachindex(edges)
        bnd = lefschetz_boundary(lc, edges[k])
        p0 = coords[bnd[1]]
        p1 = coords[bnd[2]]
        is_transverse[k] = planar_edge_is_transverse(p0, p1, vf, np=npts)
    end

    # Return the non-transverse edges

    ntindices = findall(isequal(false), is_transverse)
    return edges[ntindices]
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function planar_edge_is_transverse(p0::Vector{<:Real}, p1::Vector{<:Real}, vf;
                                   np::Int=100)
    #
    # Test whether an edge is transversal to the flow
    #

    # Create normal vector to the edge

    pvec = p1 .- p0
    nvec = [-pvec[2], pvec[1]]

    # Determine the inner product of flow and normal vector at p0

    p0vf = vf(p0)
    dprod = p0vf[1]*nvec[1] + p0vf[2]*nvec[2]

    sfac = -1.0
    if iszero(dprod)
        return false
    elseif dprod > 0.0
        sfac = 1.0
    end

    # Determine sign of the dot products along the edge

    dp_is_n = false
    npf = Float64(np)
    k = 0

    while (k < np) && (dp_is_n == false)
        k = k + 1
        pc = p0 .+ ((k / npf) .* pvec)
        pvf = vf(pc)
        dprod = pvf[1]*nvec[1] + pvf[2]*nvec[2]
        if (dprod * sfac <= 0.0)
            dp_is_n = true   # There was a sign change
        end
    end

    # If there was a sign change, return false

    if dp_is_n
        return false
    else
        return true
    end
end

###############################################################################

