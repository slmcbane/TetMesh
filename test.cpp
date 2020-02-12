#include "TetMesh.hpp"

#include <iostream>

template <class CoordT, size_t MaxElementAdjacencies, size_t MaxNodeAdjacencies,
          size_t NodesPerFace, size_t InternalNodes>
std::ostream &operator<<(std::ostream &stream,
                         const msh::TetMesh<CoordT, MaxElementAdjacencies, MaxNodeAdjacencies,
                                            NodesPerFace, InternalNodes> &mesh)
{
    stream << "Mesh Information\n"
           << "================\n"
           << "Mesh contains " << mesh.num_nodes() << " nodes:\n";

    for (size_t i = 0; i < mesh.num_nodes(); ++i)
    {
        const auto &coord = mesh.coord(i);
        stream << "  (" << coord[0] << ", " << coord[1] << ')';
        stream << (i < mesh.num_nodes() - 1 ? ",\n" : "\n");
    }

    stream << "\nMesh contains " << mesh.num_elements() << " elements:\n";

    for (size_t i = 0; i < mesh.num_elements(); ++i)
    {
        const auto &el = mesh.element(i);
        stream << "  Element " << i << ":\n  -----------------\n";
        stream << "    Control nodes: (" << el.control_nodes[0] << ", "
               << el.control_nodes[1] << ", " << el.control_nodes[2] << ")\n";
        stream << "    Adjacent elements: (";
        for (size_t j = 0; j < el.adjacent_elements.size(); ++j)
        {
            stream << el.adjacent_elements[j] << (j + 1 == el.adjacent_elements.size() ? ")\n" : ", ");
        }
        stream << "    Element's faces: (" << el.faces[0] << ", "
               << el.faces[1] << ", " << el.faces[2] << ")\n";
        if (NodesPerFace != 0)
        {
            stream << "    Nodes on faces: ";
            for (size_t j = 0; j < 3; ++j)
            {
                stream << '(';
                for (size_t k = 0; k < NodesPerFace; ++k)
                {
                    stream << el.face_nodes[j][k] << (k == NodesPerFace - 1 ? ")" : ", ");
                }
                stream << (j < 2 ? ", " : "\n");
            }
        }

        if (InternalNodes != 0)
        {
            stream << "    Internal nodes: (";
            for (size_t j = 0; j < InternalNodes; ++j)
            {
                stream << el.internal_nodes[j] << (j + 1 < InternalNodes ? ", " : ")\n");
            }
        }
        stream << '\n';
    }

    stream << "  Node adjacency information\n"
           << "  --------------------------\n";
    
    for (size_t i = 0; i < mesh.num_nodes(); ++i)
    {
        stream << "    Nodes sharing an element with node " << i << ": (";
        const auto &adj = mesh.adjacent_nodes(i);
        for (size_t j = 0; j < adj.size(); ++j)
        {
            stream << adj[j] << (j == adj.size()-1 ? ")\n" : ", ");
        }
    }

    stream << '\n' << "  Boundary information\n"
                   << "  --------------------\n";
    for (size_t i = 0; i < mesh.num_boundaries(); ++i)
    {
        stream << "    Boundary defined by nodes: (";
        const auto &nodes = mesh.boundary(i).nodes;
        for (size_t j = 0; j < nodes.size(); ++j)
        {
            stream << nodes[j] << (j == nodes.size()-1 ? ")\n" : ", ");
        }
    }

    return stream;
}

const std::array<std::array<double, 2>, 12> nodes = {
    -1, 1,
    -1, -1,
    1, -1,
    1, 1,
    0, 1,
    1, 0,
    0, -1,
    -1, 0,
    -0.5, -0.5,
    0.25, -0.25,
    -0.25, 0.25,
    0.5, 0.5
};

const std::array<std::array<size_t, 3>, 14> tets = {
    0, 4, 10,
    5, 2, 9,
    9, 2, 6,
    0, 10, 7,
    4, 3, 11,
    3, 5, 11,
    1, 8, 6,
    7, 8, 1,
    4, 11, 10,
    8, 9, 6,
    10, 9, 8,
    11, 5, 9,
    10, 11, 9,
    7, 10, 8
};

const std::array<std::array<size_t, 5>, 1> boundaries = {
    0, 4, 3, 5, 2
};

int main()
{
    msh::TetMesh<double, 11, 36> mesh(nodes, tets, boundaries);
    std::cout << mesh;
    std::cout << "Average bandwidth:  " << mesh.average_bandwidth() << "\n";

    std::cout << "\n\nRenumbering-----------------------------------------\n\n\n";
    mesh.renumber_nodes();
    std::cout << mesh;
    std::cout << "Average bandwidth:  " << mesh.average_bandwidth() << "\n";
    return 0;
}
