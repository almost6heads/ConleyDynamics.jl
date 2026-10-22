export trisurface_sphere, trisurface_ellipsoid, trisurface_cylinder
export trisurface_torus, trisurface_plane_patch, trisurface_parametric

"""
    trisurface_sphere(center::AbstractVector{<:Real}, r::Real, n::Int)
                      -> (V::Vector{Vector{Float64}}, T::Vector{Vector{Int}})

Create a triangulated sphere with center `center` and radius `r`.

The triangulation is an icosphere: an icosahedron whose triangles are
subdivided `n` times, with all vertices projected onto the sphere. The
result has `10*4^n + 2` vertices and `20*4^n` triangles, which are all of
roughly the same size and shape.

The function returns the vertex coordinates `V` and the triangles `T`,
where every triangle is a vector of three indices into `V`. The triangles
are consistently oriented, with normals pointing outward.

The returned surface can be used as a constraint or boundary surface for
mesh generation, or it can be turned into a two-dimensional simplicial
complex in three dimensions via `tetmesh_to_euclidean(V, T)`.
"""
function trisurface_sphere(center::AbstractVector{<:Real}, r::Real, n::Int)
    #
    # Create a triangulated sphere via icosahedron subdivision
    #
    return trisurface_ellipsoid(center, [r, r, r], n)
end

"""
    trisurface_ellipsoid(center::AbstractVector{<:Real},
                         radii::AbstractVector{<:Real}, n::Int)
                         -> (V::Vector{Vector{Float64}}, T::Vector{Vector{Int}})

Create a triangulated ellipsoid with center `center` and semi-axes
`radii` along the three coordinate directions.

The triangulation is obtained by scaling an icosphere with `n`
subdivisions, see `trisurface_sphere`. The function returns the vertex
coordinates `V` and the triangles `T`, consistently oriented with
outward normals.
"""
function trisurface_ellipsoid(center::AbstractVector{<:Real},
                              radii::AbstractVector{<:Real}, n::Int)
    #
    # Create a triangulated ellipsoid via icosahedron subdivision
    #

    if (length(center) != 3) || (length(radii) != 3)
        throw(ArgumentError("trisurface_ellipsoid: center and radii need " *
                            "to have length 3"))
    end
    if any(radii .<= 0)
        throw(ArgumentError("trisurface_ellipsoid: the radii need to be positive"))
    end
    if n < 0
        throw(ArgumentError("trisurface_ellipsoid: n needs to be nonnegative"))
    end

    # Start with an icosahedron

    t = (1.0 + sqrt(5.0)) / 2.0
    V = [[-1.0, t, 0.0], [1.0, t, 0.0], [-1.0, -t, 0.0], [1.0, -t, 0.0],
         [0.0, -1.0, t], [0.0, 1.0, t], [0.0, -1.0, -t], [0.0, 1.0, -t],
         [t, 0.0, -1.0], [t, 0.0, 1.0], [-t, 0.0, -1.0], [-t, 0.0, 1.0]]
    T = [[1,12,6], [1,6,2], [1,2,8], [1,8,11], [1,11,12],
         [2,6,10], [6,12,5], [12,11,3], [11,8,7], [8,2,9],
         [4,10,5], [4,5,3], [4,3,7], [4,7,9], [4,9,10],
         [5,10,6], [3,5,12], [7,3,11], [9,7,8], [10,9,2]]
    V = [v ./ sqrt(sum(abs2, v)) for v in V]

    # Subdivide every triangle into four, projecting the new
    # vertices onto the unit sphere

    for _ in 1:n
        midpoints = Dict{Tuple{Int,Int},Int}()
        function midpoint!(a, b)
            key = (min(a,b), max(a,b))
            if !haskey(midpoints, key)
                m = V[a] .+ V[b]
                push!(V, m ./ sqrt(sum(abs2, m)))
                midpoints[key] = length(V)
            end
            return midpoints[key]
        end
        newT = Vector{Vector{Int}}()
        for (a, b, c) in T
            ab = midpoint!(a, b)
            bc = midpoint!(b, c)
            ca = midpoint!(c, a)
            push!(newT, [a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca])
        end
        T = newT
    end

    # Scale and translate

    cv = Float64.(collect(center))
    rv = Float64.(collect(radii))
    V  = [cv .+ rv .* v for v in V]

    return V, T
end

"""
    trisurface_cylinder(p0::AbstractVector{<:Real}, p1::AbstractVector{<:Real},
                        r::Real, ntheta::Int, nz::Int; capped::Bool=true)
                        -> (V::Vector{Vector{Float64}}, T::Vector{Vector{Int}})

Create a triangulated cylinder of radius `r` around the axis from `p0`
to `p1`.

The lateral surface is subdivided into `ntheta` intervals around the
axis and `nz` intervals along the axis. If `capped=true` (default), the
two ends are closed by triangle fans around the points `p0` and `p1`,
so that the result is a closed surface. Otherwise, the result is an open
tube.

The function returns the vertex coordinates `V` and the triangles `T`,
consistently oriented with outward normals.
"""
function trisurface_cylinder(p0::AbstractVector{<:Real}, p1::AbstractVector{<:Real},
                             r::Real, ntheta::Int, nz::Int; capped::Bool=true)
    #
    # Create a triangulated cylinder
    #

    if (length(p0) != 3) || (length(p1) != 3)
        throw(ArgumentError("trisurface_cylinder: p0 and p1 need to have length 3"))
    end
    if r <= 0
        throw(ArgumentError("trisurface_cylinder: the radius needs to be positive"))
    end
    if (ntheta < 3) || (nz < 1)
        throw(ArgumentError("trisurface_cylinder: need ntheta >= 3 and nz >= 1"))
    end

    q0 = Float64.(collect(p0))
    ax = Float64.(collect(p1)) .- q0
    e1, e2, _ = _orthonormal_frame(ax)

    # Lateral surface

    f(u, v) = q0 .+ v .* ax .+ r .* (cos(u) .* e1 .+ sin(u) .* e2)
    V, T = _parametric_grid(f, (0.0, 2.0*pi), (0.0, 1.0), ntheta, nz, true, false)

    # Caps. The ring vertices at the bottom are 1:ntheta, and the
    # ring vertices at the top are nz*ntheta+1:(nz+1)*ntheta.

    if capped
        push!(V, q0)
        c0 = length(V)
        push!(V, q0 .+ ax)
        c1 = length(V)
        for i in 1:ntheta
            inext = mod1(i+1, ntheta)
            push!(T, [c0, inext, i])
            push!(T, [c1, nz*ntheta + i, nz*ntheta + inext])
        end
    end

    return V, T
end

"""
    trisurface_torus(center::AbstractVector{<:Real}, R::Real, r::Real,
                     nu::Int, nv::Int; axis::AbstractVector{<:Real}=[0,0,1])
                     -> (V::Vector{Vector{Float64}}, T::Vector{Vector{Int}})

Create a triangulated torus with center `center`, major radius `R`,
and minor radius `r`.

The torus is rotationally symmetric around the line through `center`
in the direction `axis`. It is subdivided into `nu` intervals around
this axis, and into `nv` intervals around the tube. The function
returns the vertex coordinates `V` and the triangles `T`, consistently
oriented with outward normals. The result is a closed surface of
genus one.
"""
function trisurface_torus(center::AbstractVector{<:Real}, R::Real, r::Real,
                          nu::Int, nv::Int; axis::AbstractVector{<:Real}=[0,0,1])
    #
    # Create a triangulated torus
    #

    if (length(center) != 3) || (length(axis) != 3)
        throw(ArgumentError("trisurface_torus: center and axis need to have length 3"))
    end
    if (r <= 0) || (R <= r)
        throw(ArgumentError("trisurface_torus: need 0 < r < R"))
    end
    if (nu < 3) || (nv < 3)
        throw(ArgumentError("trisurface_torus: need nu >= 3 and nv >= 3"))
    end

    cv = Float64.(collect(center))
    e1, e2, e3 = _orthonormal_frame(Float64.(collect(axis)))

    function f(u, v)
        rho = R + r * cos(v)
        return cv .+ (rho * cos(u)) .* e1 .+ (rho * sin(u)) .* e2 .+ (r * sin(v)) .* e3
    end

    return _parametric_grid(f, (0.0, 2.0*pi), (0.0, 2.0*pi), nu, nv, true, true)
end

"""
    trisurface_plane_patch(origin::AbstractVector{<:Real}, u::AbstractVector{<:Real},
                           v::AbstractVector{<:Real}, nu::Int, nv::Int)
                           -> (V::Vector{Vector{Float64}}, T::Vector{Vector{Int}})

Create a triangulated planar parallelogram.

The parallelogram consists of the points `origin + s*u + t*v` for
`s` and `t` in `[0,1]`, and it is subdivided into `nu` intervals in the
direction `u` and `nv` intervals in the direction `v`. The function
returns the vertex coordinates `V` and the triangles `T`, consistently
oriented with normals in the direction of the cross product of `u`
and `v`. The result is an open surface.
"""
function trisurface_plane_patch(origin::AbstractVector{<:Real}, u::AbstractVector{<:Real},
                                v::AbstractVector{<:Real}, nu::Int, nv::Int)
    #
    # Create a triangulated planar parallelogram
    #

    if any(length.((origin, u, v)) .!= 3)
        throw(ArgumentError("trisurface_plane_patch: origin, u, and v need " *
                            "to have length 3"))
    end
    if (nu < 1) || (nv < 1)
        throw(ArgumentError("trisurface_plane_patch: need nu >= 1 and nv >= 1"))
    end

    ov = Float64.(collect(origin))
    uv = Float64.(collect(u))
    vv = Float64.(collect(v))
    f(s, t) = ov .+ s .* uv .+ t .* vv

    return _parametric_grid(f, (0.0, 1.0), (0.0, 1.0), nu, nv, false, false)
end

"""
    trisurface_parametric(f, urange, vrange, nu::Int, nv::Int;
                          uperiodic::Bool=false, vperiodic::Bool=false,
                          tol::Real=1e-10)
                          -> (V::Vector{Vector{Float64}}, T::Vector{Vector{Int}})

Create a triangulated parametric surface.

The surface is the image of the function `f(u,v)`, which returns a
coordinate vector for parameters `u` in `urange = (umin, umax)` and `v`
in `vrange = (vmin, vmax)`. The parameter rectangle is subdivided into
`nu` intervals in `u` and `nv` intervals in `v`, and every subrectangle
is split into two triangles. The triangles are oriented so that their
normals point in the direction of the cross product of the partial
derivatives of `f` with respect to `u` and `v`.

If `uperiodic=true`, then `f(umin,v)` and `f(umax,v)` are assumed to be
the same point, and the two seams are glued together. The keyword
`vperiodic` does the same in the `v` direction. In addition, vertices
which are closer than `tol` are merged, and triangles which then become
degenerate are removed. This takes care of collapsed edges, such as the
poles of a sphere in spherical coordinates. The tolerance `tol` has to
be positive, and since it is an absolute distance, it has to be scaled
with the size of the surface.

# Example

```julia
fsph(u,v) = [cos(u)*sin(v), sin(u)*sin(v), cos(v)]
V, T = trisurface_parametric(fsph, (0,2*pi), (pi,0), 24, 12; uperiodic=true)
```

In this example, the poles are merged automatically, and the range
`(pi,0)` for `v` gives outward normals.
"""
function trisurface_parametric(f, urange, vrange, nu::Int, nv::Int;
                               uperiodic::Bool=false, vperiodic::Bool=false,
                               tol::Real=1e-10)
    #
    # Create a triangulated parametric surface
    #

    if (length(urange) != 2) || (length(vrange) != 2)
        throw(ArgumentError("trisurface_parametric: urange and vrange need " *
                            "to have length 2"))
    end
    if (nu < 1) || (nv < 1) || (uperiodic && nu < 3) || (vperiodic && nv < 3)
        throw(ArgumentError("trisurface_parametric: need nu, nv >= 1, " *
                            "and at least 3 intervals in periodic directions"))
    end
    if !(tol > 0)
        throw(ArgumentError("trisurface_parametric: tol needs to be positive"))
    end

    g(u, v) = Float64.(collect(f(u, v)))
    ur = (Float64(urange[1]), Float64(urange[2]))
    vr = (Float64(vrange[1]), Float64(vrange[2]))
    V, T = _parametric_grid(g, ur, vr, nu, nv, uperiodic, vperiodic)

    return _merge_close_vertices(V, T, Float64(tol))
end

############################
#                          #
#   Auxilliary functions   #
#                          #
############################

function _parametric_grid(f, urange::Tuple{Float64,Float64},
                          vrange::Tuple{Float64,Float64},
                          nu::Int, nv::Int, uper::Bool, vper::Bool)
    #
    # Triangulate a parameter grid. The vertex with grid position (i,j)
    # has index 1 + i + j*nui, where nui is the number of distinct
    # vertices in the u direction. Periodic directions wrap around.
    #

    nui = uper ? nu : nu + 1
    nvi = vper ? nv : nv + 1

    V = Vector{Vector{Float64}}()
    for j in 0:nvi-1, i in 0:nui-1
        u = urange[1] + i * (urange[2] - urange[1]) / nu
        v = vrange[1] + j * (vrange[2] - vrange[1]) / nv
        push!(V, f(u, v))
    end

    id(i, j) = 1 + mod(i, nui) + mod(j, nvi) * nui

    T = Vector{Vector{Int}}()
    for j in 0:nv-1, i in 0:nu-1
        push!(T, [id(i,j), id(i+1,j), id(i+1,j+1)])
        push!(T, [id(i,j), id(i+1,j+1), id(i,j+1)])
    end

    return V, T
end

function _merge_close_vertices(V::Vector{Vector{Float64}},
                               T::Vector{Vector{Int}}, tol::Float64)
    #
    # Merge vertices closer than tol, remove degenerate triangles,
    # and remove vertices which are no longer used
    #

    # The bucket keys are floating point numbers, so that large
    # coordinates cannot overflow. Adding 0.0 turns -0.0 into 0.0.

    d = length(V[1])
    buckets = Dict{Vector{Float64},Vector{Int}}()
    rep = zeros(Int, length(V))
    offsets = [Float64.(collect(o)) for o in Iterators.product(ntuple(_ -> -1:1, d)...)]

    for k in 1:length(V)
        key = floor.(V[k] ./ tol) .+ 0.0
        for o in offsets
            for m in get(buckets, key .+ o, Int[])
                if sum(abs2, V[m] .- V[k]) <= tol*tol
                    rep[k] = m
                    break
                end
            end
            rep[k] > 0 && break
        end
        if rep[k] == 0
            rep[k] = k
            push!(get!(buckets, key, Int[]), k)
        end
    end

    newT = [rep[t] for t in T if allunique(rep[t])]

    used = falses(length(V))
    for t in newT
        used[t] .= true
    end
    newindex = cumsum(used)

    return V[used], [newindex[t] for t in newT]
end

function _orthonormal_frame(a::Vector{Float64})
    #
    # Return a right-handed orthonormal frame (e1, e2, e3) with e3
    # in the direction of a
    #
    na = sqrt(sum(abs2, a))
    if na == 0
        throw(ArgumentError("The axis vector cannot be zero"))
    end
    e3 = a ./ na
    h  = abs(e3[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    e1 = h .- sum(h .* e3) .* e3
    e1 = e1 ./ sqrt(sum(abs2, e1))
    e2 = _cross3(e3, e1)
    return e1, e2, e3
end

_cross3(a, b) = [a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]]

function _trisurface_check(T::AbstractVector)
    #
    # Check whether a triangulated surface is closed (every edge lies
    # in exactly two triangles), and whether it is consistently
    # oriented (every directed edge occurs at most once)
    #
    edgecount = Dict{Tuple{Int,Int},Int}()
    directed  = Set{Tuple{Int,Int}}()
    oriented  = true
    for t in T
        for (a, b) in ((t[1],t[2]), (t[2],t[3]), (t[3],t[1]))
            key = (min(a,b), max(a,b))
            edgecount[key] = get(edgecount, key, 0) + 1
            if (a, b) in directed
                oriented = false
            end
            push!(directed, (a, b))
        end
    end
    closed = all(c == 2 for c in values(edgecount))
    return closed, oriented
end

function _trisurface_volume(V::AbstractVector, T::AbstractVector)
    #
    # Signed volume enclosed by a closed, oriented triangulated surface.
    # It is positive if the normals point outward.
    #
    vol = 0.0
    for t in T
        a, b, c = V[t[1]], V[t[2]], V[t[3]]
        vol += sum(a .* _cross3(b, c)) / 6.0
    end
    return vol
end

function _trisurface_interior_point(V::AbstractVector, T::AbstractVector,
                                    Tall::AbstractVector=T)
    #
    # Find a point strictly inside a closed, consistently oriented
    # triangulated surface T. Start at the centroid of the largest
    # triangle, follow the inward normal to the next intersection with
    # the triangles in Tall (which contains T, and possibly further
    # components), and return the midpoint of this segment.
    #

    Vf = [Float64.(collect(v)) for v in V]
    sgn = _trisurface_volume(Vf, T) >= 0 ? 1.0 : -1.0

    areas = [sqrt(sum(abs2, _cross3(Vf[t[2]] .- Vf[t[1]], Vf[t[3]] .- Vf[t[1]])))
             for t in T]
    k0 = argmax(areas)
    a, b, c = Vf[T[k0][1]], Vf[T[k0][2]], Vf[T[k0][3]]
    center = (a .+ b .+ c) ./ 3.0
    nrm = _cross3(b .- a, c .- a)
    dir = (-sgn / sqrt(sum(abs2, nrm))) .* nrm

    # Moller-Trumbore ray-triangle intersection

    tmin = Inf
    t0 = T[k0]
    for t in Tall
        t === t0 && continue
        p0, p1, p2 = Vf[t[1]], Vf[t[2]], Vf[t[3]]
        e1 = p1 .- p0
        e2 = p2 .- p0
        pv = _cross3(dir, e2)
        det = sum(e1 .* pv)
        abs(det) < 1e-14 && continue
        tv = center .- p0
        u = sum(tv .* pv) / det
        (u < 0 || u > 1) && continue
        qv = _cross3(tv, e1)
        v = sum(dir .* qv) / det
        (v < 0 || u + v > 1) && continue
        s = sum(e2 .* qv) / det
        if (s > 1e-12) && (s < tmin)
            tmin = s
        end
    end

    if !isfinite(tmin)
        throw(ArgumentError("Could not find an interior point; " *
                            "is the surface closed?"))
    end

    return center .+ (tmin / 2.0) .* dir
end

function _trisurface_components(T::AbstractVector)
    #
    # Split a triangulated surface into its connected components, where
    # triangles are connected if they share a vertex. Returns vectors of
    # triangle indices.
    #
    parent = Dict{Int,Int}()
    function root(v)
        while parent[v] != v
            parent[v] = parent[parent[v]]
            v = parent[v]
        end
        return v
    end
    for t in T, v in t
        haskey(parent, v) || (parent[v] = v)
    end
    for t in T
        r1 = root(t[1])
        for v in t[2:end]
            r2 = root(v)
            r1 != r2 && (parent[r2] = r1)
        end
    end
    comps = Dict{Int,Vector{Int}}()
    for (k, t) in enumerate(T)
        push!(get!(comps, root(t[1]), Int[]), k)
    end
    return sort(collect(values(comps)); by=first)
end
