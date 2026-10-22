export poisson_points_box

"""
    poisson_points_box(bmin::AbstractVector{<:Real}, bmax::AbstractVector{<:Real},
                       pdist::Real; attempts::Int=20,
                       rng::AbstractRNG=Random.default_rng())
                       -> Vector{Vector{Float64}}

Create a random sample of points in a box with a minimum distance
between the points (Poisson-disk sampling).

The box is given by `[bmin[1],bmax[1]] x ... x [bmin[d],bmax[d]]`, and
it can have any dimension `d`. The function selects random points inside
this box such that any two points have distance at least `pdist`, while
trying to fill the box as densely as possible. The argument `attempts`
specifies how many candidate points are tried around each accepted point
before giving up on it. The standard value is 20, and larger values tend
to fill holes better, but at the expense of runtime. The function uses
Bridson's algorithm, and its runtime is proportional to the number of
generated points. Every candidate point is compared with the points in
`(2*ceil(sqrt(d))+1)^d` neighboring grid cells, so the function is only
practical in low dimensions, up to about `d = 5`.

The random number generator can be specified via `rng`. Passing a seeded
generator, such as `Random.Xoshiro(1234)`, makes the result reproducible.

The points are returned as a vector of coordinate vectors. They can be
used, for example, as nodes for a Delaunay triangulation.

# Example

```julia
using Random
pts = poisson_points_box([0,0,0], [2,1,1], 0.1; rng=Xoshiro(42))
```
"""
function poisson_points_box(bmin::AbstractVector{<:Real}, bmax::AbstractVector{<:Real},
                            pdist::Real; attempts::Int=20,
                            rng::AbstractRNG=Random.default_rng())
    #
    # Poisson-disk sampling in a box via Bridson's algorithm
    #

    # Check the arguments

    d = length(bmin)
    if (d == 0) || (length(bmax) != d)
        throw(ArgumentError("poisson_points_box: bmin and bmax need to have " *
                            "the same positive length"))
    end
    if any(bmax .<= bmin)
        throw(ArgumentError("poisson_points_box: bmax needs to be larger " *
                            "than bmin in every coordinate"))
    end
    if pdist <= 0
        throw(ArgumentError("poisson_points_box: pdist needs to be positive"))
    end
    if attempts < 1
        throw(ArgumentError("poisson_points_box: attempts needs to be positive"))
    end

    lo = Float64.(collect(bmin))
    hi = Float64.(collect(bmax))
    r  = Float64(pdist)

    # Set up the background grid. Its cells have diameter r, so every
    # grid cell contains at most one point. The grid is stored as a flat
    # vector of point indices, with 0 denoting an empty cell.

    cellsize = r / sqrt(d)
    gdims    = [max(1, Int(ceil((hi[k]-lo[k]) / cellsize))) for k in 1:d]
    gstrides = cumprod(vcat(1, gdims[1:end-1]))
    grid     = zeros(Int, prod(gdims))
    reach    = Int(ceil(sqrt(d)))
    offsets  = [collect(o) for o in Iterators.product(ntuple(_ -> -reach:reach, d)...)]
    gp       = zeros(Int, d)

    function gridpos!(x)
        for k in 1:d
            gp[k] = min(gdims[k]-1, Int(floor((x[k]-lo[k]) / cellsize)))
        end
    end

    # Check whether a candidate point is far enough from all others

    points = Vector{Vector{Float64}}()

    function isvalid(x)
        gridpos!(x)
        for offset in offsets
            idx = 1
            inside = true
            for k in 1:d
                nk = gp[k] + offset[k]
                if (nk < 0) || (nk >= gdims[k])
                    inside = false
                    break
                end
                idx += nk * gstrides[k]
            end
            inside || continue
            q = grid[idx]
            if q > 0
                pq = points[q]
                dist2 = 0.0
                for k in 1:d
                    dist2 += (pq[k] - x[k])^2
                end
                dist2 < r*r && return false
            end
        end
        return true
    end

    function addpoint!(x, active)
        push!(points, x)
        gridpos!(x)
        grid[1 + sum(gp .* gstrides)] = length(points)
        push!(active, length(points))
    end

    # Start with a uniformly distributed random point

    active = Vector{Int}()
    addpoint!(lo .+ (hi .- lo) .* rand(rng, d), active)

    # Process the active list

    while !isempty(active)
        ka = rand(rng, 1:length(active))
        center = points[active[ka]]
        found = false
        for _ in 1:attempts
            # Random point in the spherical shell with radii r and 2r
            dir = randn(rng, d)
            nrm = sqrt(sum(abs2, dir))
            nrm == 0 && continue
            x = center .+ (r * (1 + rand(rng)) / nrm) .* dir
            if all(lo .<= x .<= hi) && isvalid(x)
                addpoint!(x, active)
                found = true
                break
            end
        end
        if !found
            active[ka] = active[end]
            pop!(active)
        end
    end

    return points
end
