module TetMesh

using StaticArrays
using WriteVTK: vtk_grid, MeshCell, VTKCellTypes
using Lazy: @as
using Base: setindex

struct NodeInfo
    coords::MVector{2, Float64}
    adjacent_nodes::Vector{Int}
end

struct ElementInfo{NodesPerFace, InternalNodes}
    control_nodes::MVector{3, Int}
    faces::MVector{3, Int}
    face_nodes::MVector{3, MVector{NodesPerFace, Int}}
    internal_nodes::MVector{InternalNodes, Int}
    adjacent_elements::Vector{Int}
end

struct FaceDetails{NodesPerFace}
    number::Int
    element::Int
    nodes::MVector{NodesPerFace, Int}
end

struct BoundaryRepresentation{NodesPerFace}
    nodes::Vector{Int}
    faces::Vector{FaceDetails{NodesPerFace}}
end

struct Mesh{NodesPerFace, InternalNodes, BoundaryNodesPerFace}
    nodes::Vector{NodeInfo}
    elements::Vector{ElementInfo{NodesPerFace, InternalNodes}}
    boundaries::Vector{BoundaryRepresentation{BoundaryNodesPerFace}}
    boundary_tags::Vector{Vector{String}}
end

function deserialize(stream::IO)
    max_element_adjacencies = read(stream, Int)
    max_node_adjacencies = read(stream, Int)
    NodesPerFace = read(stream, Int)
    InternalNodes = read(stream, Int)
    
    # Read nodes
    count = read(stream, Int)
    nodes = [deserialize(stream, NodeInfo) for i in 1:count]

    # Read elements
    count = read(stream, Int)
    elements = [
        deserialize(stream, ElementInfo{NodesPerFace, InternalNodes})
        for i in 1:count
    ]
    
    # Read boundaries
    count = read(stream, Int)
    boundaries = [
        deserialize(stream, BoundaryRepresentation{NodesPerFace+2})
        for i in 1:count
    ]

    # Read tags
    count = read(stream, Int)
    tags = [deserialize(stream, Vector{String}) for i in 1:count]

    Mesh{NodesPerFace, InternalNodes, NodesPerFace+2}(nodes, elements, boundaries, tags)
end

function deserialize(stream::IO, ::Val{NodesPerFace}, ::Val{InternalNodes}) where {NodesPerFace, InternalNodes}
    # Read nodes
    count = read(stream, Int)
    nodes = [deserialize(stream, NodeInfo) for i in 1:count]

    # Read elements
    count = read(stream, Int)
    elements = [
        deserialize(stream, ElementInfo{NodesPerFace, InternalNodes})
        for i in 1:count
    ]
    
    # Read boundaries
    count = read(stream, Int)
    boundaries = [
        deserialize(stream, BoundaryRepresentation{NodesPerFace+2})
        for i in 1:count
    ]

    # Read tags
    count = read(stream, Int)
    tags = [deserialize(stream, Vector{String}) for i in 1:count]

    Mesh{NodesPerFace, InternalNodes, NodesPerFace+2}(nodes, elements, boundaries, tags)
end

function deserialize(stream::IO, ::Type{NodeInfo})
    coords = deserialize(stream, MVector{2, Float64})
    adjacent = deserialize(stream, Vector{Int})
    NodeInfo(coords, adjacent)
end

function deserialize(stream::IO, ::Type{MVector{N, T}}) where {N, T}
    x = @MVector zeros(T, N)
    read!(stream, x)
end

function deserialize(stream::IO, ::Type{Vector{T}}) where T
    sz = read(stream, Int)
    x = Vector{T}(undef, sz)
    read!(stream, x)
end

function deserialize(stream::IO, ::Type{ElementInfo{NodesPerFace, InternalNodes}}) where {NodesPerFace, InternalNodes}
    control_nodes = deserialize(stream, MVector{3, Int})
    faces = deserialize(stream, MVector{3, Int})
    face_nodes = @MVector [deserialize(stream, MVector{NodesPerFace, Int}) for i in 1:3]
    internal_nodes = deserialize(stream, MVector{InternalNodes, Int})
    adjacent_elements = deserialize(stream, Vector{Int})
    ElementInfo(control_nodes, faces, face_nodes, internal_nodes, adjacent_elements)
end

function deserialize(stream::IO, ::Type{BoundaryRepresentation{N}}) where N
    nodes = deserialize(stream, Vector{Int})
    faces = Vector{FaceDetails{N}}(undef, read(stream, Int))
    for i in eachindex(faces)
        faces[i] = deserialize(stream, FaceDetails{N})
    end
    BoundaryRepresentation{N}(nodes, faces)
end

function deserialize(stream::IO, ::Type{FaceDetails{N}}) where N
    number = read(stream, Int)
    element = read(stream, Int)
    nodes = deserialize(stream, MVector{N, Int})
    FaceDetails{N}(number, element, nodes)
end

function deserialize(stream::IO, ::Type{Vector{String}})
    count = read(stream, Int)
    [deserialize(stream, String) for i in 1:count]
end

function deserialize(stream::IO, ::Type{String})
    count = read(stream, Int)
    String(read!(stream, Vector{UInt8}(undef, count)))
end

"""
Returns three items:
  - points: An array of the coordinates of vertices of triangles in the mesh;
  - elements: An array of WriteVTK.MeshCell with indices of the vertices of
    triangles in `points`;
  - node_map: `node_map[i]` maps a node `i` in the triangle mesh to a node index
    in the full (possibly higher-order) mesh. 
"""
function points_and_cells(msh::Mesh)
    node_indices = zeros(Int, length(msh.nodes))
    coords = Vector{SVector{2, Float64}}(undef, length(msh.nodes))
    vertices = NTuple{3, Int}[]
    node_map = Int[]

    node_index = 1
    for element in msh.elements
        verts = (0, 0, 0)
        i = 1
        for node in (element.control_nodes .+ 1)
            if node_indices[node] == 0
                node_indices[node] = node_index
                coords[node_index] = msh.nodes[node].coords
                node_index += 1
                push!(node_map, node)
            end
            verts = setindex(verts, node_indices[node], i)
            i += 1
        end
        push!(vertices, verts)
    end

    nnodes = node_index - 1
    resize!(coords, nnodes)

    (@as v coords reinterpret(Float64, v) reshape(v, 2, :)),
    [MeshCell(VTKCellTypes.VTK_TRIANGLE, verts) for verts in vertices],
    node_map
end

function initialize_vtk(target, msh)
    points, cells, node_map = points_and_cells(msh)
    vtk_grid(target, points, cells), node_map
end

function add_scalar_data(target, msh::Mesh, node_map::Vector{Int}, data::AbstractVector, name)
    @assert length(data) == length(msh.nodes)
    target[name] = @view data[node_map]
    target
end

end # module TetMesh
