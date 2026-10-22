# Gmsh element types of linear simplices, indexed by dimension + 1
const GMSH_SIMPLEX_TYPES = (15, 1, 2, 4)

function _gmsh_check_initialized(fname::String)
    #
    # Make sure that Gmsh is running
    #
    if gmsh.isInitialized() != 1
        throw(ArgumentError("$fname: Gmsh is not initialized"))
    end
end

function _gmsh_elements(dim::Integer, entities, fname::String)
    #
    # Collect the node tags of the linear simplices of dimension dim on
    # the given entities, and reject all other element types
    #
    stype = GMSH_SIMPLEX_TYPES[dim+1]
    elements = Vector{Vector{UInt64}}()
    for ent in entities
        etypes, _, enodes = gmsh.model.mesh.getElements(dim, ent)
        for (et, nodes) in zip(etypes, enodes)
            if et != stype
                ename = gmsh.model.mesh.getElementProperties(et)[1]
                throw(ArgumentError("$fname: the mesh contains elements of type " *
                                    "\"$ename\", but only linear simplices are " *
                                    "supported"))
            end
            nn = dim + 1
            for k in 1:div(length(nodes), nn)
                push!(elements, nodes[(k-1)*nn+1:k*nn])
            end
        end
    end
    return elements
end

function ConleyDynamics.gmsh_to_simplicial(; dim::Int=3, p::Int=2,
                                           physical::Union{Tuple{Integer,Integer},Nothing}=nothing)
    #
    # Convert the mesh of the current Gmsh model into an EuclideanComplex
    #

    _gmsh_check_initialized("gmsh_to_simplicial")
    if !(1 <= dim <= 3)
        throw(ArgumentError("gmsh_to_simplicial: dim needs to be 1, 2, or 3"))
    end

    # Select the entities

    if isnothing(physical)
        entities = [-1]
    else
        if physical[1] != dim
            throw(ArgumentError("gmsh_to_simplicial: the physical group needs " *
                                "to have dimension dim"))
        end
        entities = gmsh.model.getEntitiesForPhysicalGroup(physical[1], physical[2])
    end

    # Nodes and elements. Node tags need not be contiguous.

    nodetags, nodecoords, _ = gmsh.model.mesh.getNodes()
    tag2index = Dict{UInt64,Int}(nodetags[k] => k for k in 1:length(nodetags))
    points = [nodecoords[3k-2:3k] for k in 1:length(nodetags)]

    elements = _gmsh_elements(dim, entities, "gmsh_to_simplicial")
    if isempty(elements)
        throw(ArgumentError("gmsh_to_simplicial: the model contains no mesh " *
                            "elements of dimension $dim. Has the mesh been generated?"))
    end
    cells = [[tag2index[t] for t in el] for el in elements]

    return tetmesh_to_euclidean(points, cells; p=p, droporphans=true)
end

function ConleyDynamics.gmsh_file_to_simplicial(filename::AbstractString;
                                                dim::Int=3, p::Int=2)
    #
    # Read a file with Gmsh and convert its mesh into an EuclideanComplex
    #

    if !isfile(filename)
        throw(ArgumentError("gmsh_file_to_simplicial: file $filename not found"))
    end

    wasinitialized = gmsh.isInitialized() == 1
    if !wasinitialized
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
    end

    try
        gmsh.open(filename)
        etypes, _, _ = gmsh.model.mesh.getElements(dim)
        if isempty(etypes)
            gmsh.model.mesh.generate(dim)
        end
        return gmsh_to_simplicial(dim=dim, p=p)
    finally
        wasinitialized || gmsh.finalize()
    end
end

function ConleyDynamics.gmsh_marker_cells(ec::EuclideanComplex, dim::Integer, tag::Integer;
                                          physical::Bool=true, closure::Bool=true)
    #
    # Find the cells formed by the mesh elements of a physical group
    # or an entity of the current Gmsh model
    #

    _gmsh_check_initialized("gmsh_marker_cells")
    if !(0 <= dim <= 3)
        throw(ArgumentError("gmsh_marker_cells: dim needs to be 0, 1, 2, or 3"))
    end

    entities = physical ? gmsh.model.getEntitiesForPhysicalGroup(dim, tag) : [tag]
    elements = _gmsh_elements(dim, entities, "gmsh_marker_cells")

    # Identify the vertices of ec by their coordinates

    vlabel = Dict{NTuple{3,Float64},String}()
    for k in 1:ec.ncells
        if ec.dimensions[k] == 0
            c = ec.coords[k][1]
            vlabel[(c[1], c[2], c[3])] = ec.labels[k]
        end
    end

    nodetags, nodecoords, _ = gmsh.model.mesh.getNodes()
    tag2index = Dict{UInt64,Int}(nodetags[k] => k for k in 1:length(nodetags))

    cells = Vector{Int}()
    for el in elements
        labels = String[]
        for t in el
            m = tag2index[t]
            lab = get(vlabel, (nodecoords[3m-2], nodecoords[3m-1], nodecoords[3m]), "")
            isempty(lab) && break
            push!(labels, lab)
        end
        length(labels) == length(el) || continue
        k = get(ec.indices, join(sort(labels)), 0)
        k > 0 && push!(cells, k)
    end

    return closure ? lefschetz_closure(ec, cells) : sort(unique(cells))
end
