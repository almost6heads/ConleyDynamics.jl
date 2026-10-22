function ConleyDynamics.create_simplicial_delaunay3d(bmin::AbstractVector{<:Real},
                                                     bmax::AbstractVector{<:Real},
                                                     pdist::Real, attempts::Int;
                                                     p::Int=2, box_corners::Bool=false,
                                                     rng::AbstractRNG=Random.default_rng())
    #
    # Delaunay tetrahedralization of a Poisson-disk sample in a box
    #
    points = poisson_points_box(bmin, bmax, pdist; attempts=attempts, rng=rng)
    return _delaunay3d_from_points(points, bmin, bmax, p, box_corners)
end

function ConleyDynamics.create_simplicial_delaunay3d(bmin::AbstractVector{<:Real},
                                                     bmax::AbstractVector{<:Real},
                                                     npoints::Int;
                                                     p::Int=2, box_corners::Bool=false,
                                                     rng::AbstractRNG=Random.default_rng())
    #
    # Delaunay tetrahedralization of uniformly distributed points in a box
    #
    if (length(bmin) != 3) || (length(bmax) != 3) || any(bmax .<= bmin)
        throw(ArgumentError("create_simplicial_delaunay3d: need boxes in three " *
                            "dimensions with bmin < bmax"))
    end
    if npoints < 4
        throw(ArgumentError("create_simplicial_delaunay3d: need at least 4 points"))
    end
    lo = Float64.(bmin)
    hi = Float64.(bmax)
    points = [lo .+ (hi .- lo) .* rand(rng, 3) for _ in 1:npoints]
    return _delaunay3d_from_points(points, bmin, bmax, p, box_corners)
end

function _delaunay3d_from_points(points::Vector{Vector{Float64}}, bmin, bmax,
                                 p::Int, box_corners::Bool)
    #
    # Unconstrained Delaunay tetrahedralization of a point set
    #

    if length(bmin) != 3
        throw(ArgumentError("create_simplicial_delaunay3d: need boxes in three dimensions"))
    end

    if box_corners
        corners = [[x, y, z] for z in (bmin[3], bmax[3]),
                                 y in (bmin[2], bmax[2]),
                                 x in (bmin[1], bmax[1])]
        points = vcat([Float64.(c) for c in vec(corners)], points)
    end

    io = RawTetGenIO{Cdouble}()
    io.pointlist = reduce(hcat, points)
    out = tetrahedralize(io, "Q")

    opoints = [Vector{Float64}(out.pointlist[:,k]) for k in 1:size(out.pointlist,2)]
    tets    = [Vector{Int}(out.tetrahedronlist[1:4,k]) for k in 1:size(out.tetrahedronlist,2)]

    return tetmesh_to_euclidean(opoints, tets; p=p)
end
