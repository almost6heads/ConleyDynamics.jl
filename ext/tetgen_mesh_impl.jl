function ConleyDynamics.tetgen_mesh(dom::TetDomain;
                                    maxvolume::Union{Real,Nothing}=nothing,
                                    minratio::Union{Real,Nothing}=1.414,
                                    mindihedral::Union{Real,Nothing}=nothing,
                                    steiner_on_boundary::Bool=true,
                                    quiet::Bool=true,
                                    switches::Union{String,Nothing}=nothing)
    #
    # Create a constrained tetrahedral mesh of a TetDomain with TetGen
    #

    if isempty(dom.facets)
        throw(ArgumentError("tetgen_mesh: the domain has no facets"))
    end
    if !isnothing(minratio) && (minratio < 1)
        throw(ArgumentError("tetgen_mesh: minratio needs to be at least 1, since " *
                            "TetGen does not terminate for smaller values"))
    end
    if !isnothing(mindihedral) && !(0 < mindihedral < 90)
        throw(ArgumentError("tetgen_mesh: mindihedral needs to be between 0 and 90"))
    end

    # Fill the TetGen input structure. Indices are 1-based on the
    # Julia side.

    io = RawTetGenIO{Cdouble}()
    io.pointlist = Float64[pt[k] for k in 1:3, pt in dom.points]
    facetlist!(io, dom.facets)
    io.facetmarkerlist = Cint.(dom.markers)
    if !isempty(dom.holes)
        io.holelist = Float64[pt[k] for k in 1:3, pt in dom.holes]
    end
    if !isempty(dom.regions)
        io.regionlist = Float64[k <= 3 ? r[1][k] : (k == 4 ? r[2] : r[3])
                                for k in 1:5, r in dom.regions]
    end

    # Build the switch string

    if isnothing(switches)
        sw = "p"
        if !isnothing(minratio) || !isnothing(mindihedral)
            sw *= "q" * string(Float64(something(minratio, 2.0)))
            if !isnothing(mindihedral)
                sw *= "/" * string(Float64(mindihedral))
            end
        end
        if !isnothing(maxvolume)
            if maxvolume <= 0
                throw(ArgumentError("tetgen_mesh: maxvolume needs to be positive"))
            end
            sw *= "a" * string(Float64(maxvolume))
        end
        if !isempty(dom.regions)
            sw *= "A"
            if any(r[3] > 0 for r in dom.regions)
                sw *= "a"
            end
        end
        if !steiner_on_boundary
            sw *= "Y"
        end
        if quiet
            sw *= "Q"
        end
    else
        sw = switches
    end

    # Call TetGen and translate its errors

    out = _tetgen_run(io, sw, "tetgen_mesh",
                      "Possible causes are constraint surfaces or curves which " *
                      "intersect each other or the boundary (for intersecting " *
                      "geometry, use Gmsh.jl instead), or quality bounds which " *
                      "are too strict, such as minratio close to 1 or a large " *
                      "mindihedral.")

    ntet = size(out.tetrahedronlist, 2)
    if ntet == 0
        throw(ArgumentError("tetgen_mesh: the mesh contains no tetrahedra. " *
                            "Is the outer boundary closed, and does every hole " *
                            "point lie inside a closed constraint surface?"))
    end

    # Collect the results

    points     = [Vector{Float64}(out.pointlist[:,k]) for k in 1:size(out.pointlist,2)]
    tetrahedra = [Vector{Int}(out.tetrahedronlist[1:4,k]) for k in 1:ntet]
    trifaces   = [Vector{Int}(out.trifacelist[:,k]) for k in 1:size(out.trifacelist,2)]
    trimarkers = Vector{Int}(out.trifacemarkerlist)
    if length(trimarkers) != length(trifaces)
        trimarkers = zeros(Int, length(trifaces))
    end
    tetattrs = size(out.tetrahedronattributelist, 1) > 0 ?
                   round.(Int, vec(out.tetrahedronattributelist[1,:])) : Int[]
    segments = [(dom.points[f[1]], dom.points[f[2]], m)
                for (f, m) in zip(dom.facets, dom.markers) if length(f) == 2]

    return TetGenResult(points, tetrahedra, trifaces, trimarkers, tetattrs, segments)
end

function _tetgen_tmpfiles()
    #
    # Debugging files which TetGen writes into the current directory
    # when it detects intersecting facets
    #
    return Set(f for f in readdir(pwd()) if startswith(f, "tetgen-tmpfile"))
end

function _tetgen_run(io::RawTetGenIO, switches::String, fname::String, hint::String)
    #
    # Call TetGen, turn its errors into an ArgumentError, and remove
    # the debugging files it created in the current directory
    #
    before = _tetgen_tmpfiles()
    try
        return tetrahedralize(io, switches)
    catch e
        if e isa TetGen.TetGenError
            throw(ArgumentError("$fname: TetGen failed (" *
                                strip(sprint(showerror, e)) * "). " * hint))
        end
        rethrow()
    finally
        for f in setdiff(_tetgen_tmpfiles(), before)
            rm(joinpath(pwd(), f); force=true)
        end
    end
end
