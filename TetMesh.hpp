/*
 * Copyright 2020 The University of Texas at Austin.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef TETMESH_HPP
#define TETMESH_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

using std::size_t;
using std::get;

#include "SmallVector.hpp"

namespace msh
{

template <size_t MaxElementAdjacencies, size_t NodesPerFace, size_t InternalNodes>
struct ElementInfo
{
    std::array<size_t, 3> control_nodes;
    std::array<size_t, 3> faces;
    std::array<std::array<size_t, NodesPerFace>, 3> face_nodes;
    std::array<size_t, InternalNodes> internal_nodes;
    smv::SmallVector<size_t, MaxElementAdjacencies> adjacent_elements;

    ElementInfo(size_t a, size_t b, size_t c) noexcept : control_nodes{a, b, c},
        faces{}, face_nodes{}, internal_nodes{}, adjacent_elements()
    {
        for (auto &face: face_nodes)
        {
            std::fill(face.begin(), face.end(), static_cast<size_t>(-1));
        }
        std::fill(internal_nodes.begin(), internal_nodes.end(), static_cast<size_t>(-1));
    }
};

template <class CoordT, size_t MaxNodeAdjacencies>
struct NodeInfo
{
    std::array<CoordT, 2> coords;
    smv::SmallVector<size_t, MaxNodeAdjacencies> adjacent_nodes;

    constexpr NodeInfo(CoordT a, CoordT b) noexcept : coords{a, b}, adjacent_nodes{}
    {}

    constexpr NodeInfo() : coords{}, adjacent_nodes{}
    {}
};

/*
 * This class describes a mesh composed of simple triangles
 *
 * The mesh described by this class is defined by a collection of nodes, which
 * are simply points in 2D space, a collection of triangles, specified by the
 * indices of their 3 vertices in the array of nodes, and a collection of
 * bounding curves specified by an ordered list of indices of nodes that, when
 * connected by line segments, form the curve.
 *
 * I am not implementing much in the way of sanity checking for the mesh;
 * e.g. "bounding curves" could be not on the boundary of the mesh and nothing
 * untoward would happen, there could be orphan nodes, the triangles don't form
 * a partition of the domain, etc. It is primarily intended as a data structure
 * to facilitate finite element simulation and the structure of the mesh will
 * come from an external mesh generator.
 *
 * TEMPLATE PARAMETERS
 * -------------------
 *  - CoordT: Data type for real numbers used when specifying nodes.
 *  - MaxElementAdjacencies: The maximum number of elements that share a node
 *    with a given element, not including itself. A vector with a fixed max
 *    size is used to store adjacencies and its max size is fixed by this
 *    argument.
 *  - MaxNodeAdjacencies: The maximum number of adjacent nodes to any single
 *    element; this is the number of nodes in the set of all nodes on either
 *    the element or any of the elements that neighbor it.
 *
 *
 * When setting `MaxElementAdjacencies` and `MaxNodeAdjacencies`, note that
 * an exception *will* be thrown if these values are exceeded. It is set as
 * a template parameter to optimize storage space; just increase it as
 * needed and recompile. The performance penalty for extra space should be
 * reasonable.
 */
template <class CoordT, size_t MaxElementAdjacencies, size_t MaxNodeAdjacencies,
          size_t NodesPerFace = 0, size_t InternalNodes = 0>
class TetMesh
{
public:
    typedef ElementInfo<MaxElementAdjacencies, NodesPerFace, InternalNodes> el_type;
    typedef NodeInfo<CoordT, MaxNodeAdjacencies> node_type;
    /*
     * Construct the TetMesh from a list of nodes, list of elements, and
     * list of curves.
     *
     * All three container arguments must support `size()` and the STL iterator
     * interface. The elements in `nodes` and `tets` should support `get<i>` for
     * element access, and elements of `bounding_curves` should be lists of
     * node indices in a format that allows indexing with `[]`.
     */
    template <class NodeContainer, class ElContainer, class CurveContainer>
    TetMesh(const NodeContainer &nodes, const ElContainer &tets,
            const CurveContainer &bounding_curves)
    {
        m_nodes.reserve(nodes.size());
        m_elems.reserve(tets.size());
        m_bounding_curves.reserve(bounding_curves.size());

        for (const auto &node: nodes)
        {
            m_nodes.emplace_back(get<0>(node), get<1>(node));
        }

        for (const auto &tet: tets)
        {
            auto n1 = get<0>(tet);
            auto n2 = get<1>(tet);
            auto n3 = get<2>(tet);

            auto det = (m_nodes[n3].coords[0] - m_nodes[n1].coords[0]) *
                       (m_nodes[n2].coords[1] - m_nodes[n1].coords[1]) -
                       (m_nodes[n2].coords[0] - m_nodes[n1].coords[0]) *
                       (m_nodes[n3].coords[1] - m_nodes[n1].coords[1]);
            
            if (det > 0)
            {
                m_elems.emplace_back(n1, n2, n3);
            }
            else
            {
                m_elems.emplace_back(n2, n1, n3);
            }
        }

        for (const auto &curve: bounding_curves)
        {
            m_bounding_curves.push_back(std::vector<size_t>());
            auto &m_curve = *(m_bounding_curves.end() - 1);
            m_curve.reserve(curve.size());
            for (size_t i = 0; i < curve.size(); ++i)
            {
                m_curve.push_back(curve[i]);
            }
        }
        
        build_adjacencies();

        if (NodesPerFace != 0 || InternalNodes != 0)
        {
            assign_face_and_internal_nodes();
        }
    } // constructor

    const el_type &element(size_t i) const noexcept
    {
        return m_elems[i];
    }

    const node_type &node(size_t i) const noexcept
    {
        return m_nodes[i];
    }

private:
    std::vector<node_type> m_nodes;
    std::vector<el_type> m_elems;
    std::vector<std::vector<size_t>> m_bounding_curves;
    size_t m_num_faces;

    void build_adjacencies()
    {
        find_adjacent_nodes_and_elements();
        assign_face_numbers();
    }

    void find_adjacent_nodes_and_elements()
    {
        std::vector<smv::SmallVector<size_t, MaxElementAdjacencies+1>> node_neighbors(m_nodes.size());

        for (size_t el = 0; el < m_elems.size(); ++el)
        {
            for (size_t n: m_elems[el].control_nodes)
            {
                node_neighbors[n].push_back(el);
            }

        }

        for (size_t el = 0; el < m_elems.size(); ++el)
        {
            process_adjacencies(el, node_neighbors);
        }
    }

    void process_adjacencies(size_t el,
        const std::vector<smv::SmallVector<size_t, MaxElementAdjacencies+1>> &node_neighbors)
    {
        for (size_t n: m_elems[el].control_nodes)
        {
            // This might be overly clever, but I think by only adding the
            // adjacency once for each pair I will avoid the cost of scanning
            // adjacencies and avoiding duplicates.
            for (size_t m: m_elems[el].control_nodes)
            {
                if (m < n)
                {
                    add_node_adjacency(m, n);
                    add_node_adjacency(n, m);
                }
            }

            for (size_t neighbor: node_neighbors[n])
            {
                if (neighbor < el)
                {
                    add_element_adjacency(neighbor, el);
                    add_element_adjacency(el, neighbor);
                }
            }
        }
    }

    void add_node_adjacency(size_t m, size_t n)
    {
        auto &adj = m_nodes[m].adjacent_nodes;
        if (std::count(adj.begin(), adj.end(), n) == 0)
        {
            adj.push_back(n);
        }
    }

    void add_element_adjacency(size_t m, size_t n)
    {
        auto &adj = m_elems[m].adjacent_elements;
        if (std::count(adj.begin(), adj.end(), n) == 0)
        {
            adj.push_back(n);
        }
    }

    void assign_face_and_internal_nodes()
    {
        size_t node_number = m_nodes.size();

        for (auto &el: m_elems)
        {
            for (size_t facei = 0; facei < 3; ++facei)
            {
                auto &face = el.face_nodes[facei];
                if (face[0] != static_cast<size_t>(-1)) { continue; }

                for (size_t i = 0; i < NodesPerFace; ++i)
                {
                    face[i] = node_number++;
                    m_nodes.emplace_back();
                }
            }

            for (size_t i = 0; i < InternalNodes; ++i)
            {
                el.internal_nodes[i] = node_number++;
                m_nodes.emplace_back();
            }
            add_face_and_internal_nodes_to_adjacent(el);
        }
    }

    void add_face_and_internal_nodes_to_adjacent(const el_type &el)
    {
        for (size_t facei = 0; facei < 3; ++facei)
        {
            const auto &face_nodes = el.face_nodes[facei];
            add_node_adjacency(face_nodes[0], face_nodes[1]);
            add_node_adjacency(face_nodes[1], face_nodes[0]);
            
            for (size_t facej = facei+1; facej < 3; ++facej)
            {
                const auto &other_face_nodes = el.face_nodes[facej];
                add_node_adjacency(face_nodes[0], other_face_nodes[0]);
                add_node_adjacency(face_nodes[0], other_face_nodes[1]);
                add_node_adjacency(face_nodes[1], other_face_nodes[0]);
                add_node_adjacency(face_nodes[1], other_face_nodes[1]);
                add_node_adjacency(other_face_nodes[0], face_nodes[0]);
                add_node_adjacency(other_face_nodes[0], face_nodes[1]);
                add_node_adjacency(other_face_nodes[1], face_nodes[0]);
                add_node_adjacency(other_face_nodes[1], face_nodes[1]);
            }

            for (size_t node: el.control_nodes)
            {
                add_node_adjacency(face_nodes[0], node);
                add_node_adjacency(face_nodes[1], node);
                add_node_adjacency(node, face_nodes[0]);
                add_node_adjacency(node, face_nodes[1]);
                for (size_t inode: el.internal_nodes)
                {
                    add_node_adjacency(node, inode);
                    add_node_adjacency(inode, node);
                }
            }

            for (size_t node: el.internal_nodes)
            {
                add_node_adjacency(face_nodes[0], node);
                add_node_adjacency(face_nodes[1], node);
                add_node_adjacency(node, face_nodes[0]);
                add_node_adjacency(node, face_nodes[1]);
            }

            for (size_t other_el_index: el.adjacent_elements)
            {
                auto &other_el = m_elems[other_el_index];
                auto it = std::find(other_el.faces.begin(), other_el.faces.end(), el.faces[facei]);
                if (it != other_el.faces.end())
                {
                    auto which = it - other_el.faces.begin();
                    other_el.face_nodes[which][0] = face_nodes[1];
                    other_el.face_nodes[which][1] = face_nodes[0];
                }
            }
        }
    }

    void assign_face_numbers()
    {
        size_t face_number = 0;
        std::vector<smv::SmallVector<std::array<size_t, 2>, MaxNodeAdjacencies>>
            faces_adjoining(m_nodes.size());

        for (size_t eli = 0; eli < m_elems.size(); ++eli)
        {
            el_type &el = m_elems[eli];
            bool face_numbered[] = { false, false, false };

            // Discover already assigned faces
            for (auto face: faces_adjoining[el.control_nodes[0]])
            {
                if (face[1] == el.control_nodes[1])
                {
                    el.faces[0] = face[0];
                    face_numbered[0] = true;
                }
                else if (face[1] == el.control_nodes[2])
                {
                    el.faces[2] = face[0];
                    face_numbered[2] = true;
                }
            }

            for (auto face: faces_adjoining[el.control_nodes[1]])
            {
                if (face[1] == el.control_nodes[2])
                {
                    el.faces[1] = face[0];
                    face_numbered[1] = true;
                }
            }

            // For faces not assigned, assign numbers and do the bookkeeping.
            if (!face_numbered[0])
            {
                el.faces[0] = face_number;
                faces_adjoining[el.control_nodes[0]].push_back(
                    std::array<size_t, 2>{face_number, el.control_nodes[1]}
                );
                faces_adjoining[el.control_nodes[1]].push_back(
                    std::array<size_t, 2>{face_number, el.control_nodes[0]}
                );
                face_number += 1;
            }

            if (!face_numbered[1])
            {
                el.faces[1] = face_number;
                faces_adjoining[el.control_nodes[1]].push_back(
                    std::array<size_t, 2>{face_number, el.control_nodes[2]}
                );
                faces_adjoining[el.control_nodes[2]].push_back(
                    std::array<size_t, 2>{face_number, el.control_nodes[1]}
                );
                face_number += 1;
            }

            if (!face_numbered[2])
            {
                el.faces[2] = face_number;
                faces_adjoining[el.control_nodes[0]].push_back(
                    std::array<size_t, 2>{face_number, el.control_nodes[2]}
                );
                faces_adjoining[el.control_nodes[2]].push_back(
                    std::array<size_t, 2>{face_number, el.control_nodes[0]}
                );
                face_number += 1;
            }
        }
        m_num_faces = face_number;
    }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test constructing a first order tet mesh w/ its adjacencies")
{
    const std::vector<std::array<int, 2>> nodes = {
        std::array<int, 2>{-1, -1},
        std::array<int, 2>{-1, 1},
        std::array<int, 2>{1, 1},
        std::array<int, 2>{1, -1}
    };

    const std::vector<std::array<int, 3>> tets = {
        std::array<int, 3>{0, 1, 3},
        std::array<int, 3>{1, 3, 2}
    };

    const TetMesh<int, 1, 3> mesh(nodes, tets, std::vector<std::vector<int>>());
    REQUIRE(mesh.element(0).control_nodes == std::array<size_t, 3>{0, 1, 3});
    REQUIRE(mesh.element(1).control_nodes == std::array<size_t, 3>{3, 1, 2});

    auto eladj = mesh.element(0).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element(1).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    auto nodeadj = mesh.node(0).adjacent_nodes;
    REQUIRE(nodeadj.size() == 2);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 3);

    nodeadj = mesh.node(1).adjacent_nodes;
    REQUIRE(nodeadj.size() == 3);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 2);
    REQUIRE(nodeadj[2] == 3);

    nodeadj = mesh.node(2).adjacent_nodes;
    REQUIRE(nodeadj.size() == 2);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 3);

    nodeadj = mesh.node(3).adjacent_nodes;
    REQUIRE(nodeadj.size() == 3);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 2);

    REQUIRE(mesh.element(0).faces == std::array<size_t, 3>{0, 1, 2});
    REQUIRE(mesh.element(1).faces == std::array<size_t, 3>{1, 3, 4});
} // TEST_CASE

TEST_CASE("Test constructing a third order mesh")
{
    const std::vector<std::array<int, 2>> nodes = {
        std::array<int, 2>{-1, -1},
        std::array<int, 2>{-1, 1},
        std::array<int, 2>{1, 1},
        std::array<int, 2>{1, -1}
    };

    const std::vector<std::array<int, 3>> tets = {
        std::array<int, 3>{0, 1, 3},
        std::array<int, 3>{1, 3, 2}
    };

    const TetMesh<int, 1, 15, 2, 1> mesh(nodes, tets, std::vector<std::vector<int>>());

    auto eladj = mesh.element(0).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element(1).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    auto nodeadj = mesh.node(1).adjacent_nodes;
    REQUIRE(nodeadj.size() == 15);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    for (size_t i = 1; i < 15; ++i)
    {
        CHECK(nodeadj[i] == i+1);
    }

    nodeadj = mesh.node(15).adjacent_nodes;
    REQUIRE(nodeadj.size() == 9);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 2);
    REQUIRE(nodeadj[2] == 3);
    REQUIRE(nodeadj[3] == 6);
    REQUIRE(nodeadj[4] == 7);
    REQUIRE(nodeadj[5] == 11);
    REQUIRE(nodeadj[6] == 12);
    REQUIRE(nodeadj[7] == 13);
    REQUIRE(nodeadj[8] == 14);

    nodeadj = mesh.node(8).adjacent_nodes;
    REQUIRE(nodeadj.size() == 9);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 3);
    REQUIRE(nodeadj[3] == 4);
    REQUIRE(nodeadj[4] == 5);
    REQUIRE(nodeadj[5] == 6);
    REQUIRE(nodeadj[6] == 7);
    REQUIRE(nodeadj[7] == 9);
    REQUIRE(nodeadj[8] == 10);

    auto face_nodes = mesh.element(0).face_nodes;
    REQUIRE(face_nodes[0] == std::array<size_t, 2>{4, 5});
    REQUIRE(face_nodes[1] == std::array<size_t, 2>{6, 7});
    REQUIRE(face_nodes[2] == std::array<size_t, 2>{8, 9});

    face_nodes = mesh.element(1).face_nodes;
    REQUIRE(face_nodes[0] == std::array<size_t, 2>{7, 6});
    REQUIRE(face_nodes[1] == std::array<size_t, 2>{11, 12});
    REQUIRE(face_nodes[2] == std::array<size_t, 2>{13, 14});

    auto internal_nodes = mesh.element(0).internal_nodes;
    REQUIRE(internal_nodes[0] == 10);

    internal_nodes = mesh.element(1).internal_nodes;
    REQUIRE(internal_nodes[0] == 15);
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace msh

#endif // TETMESH_HPP
