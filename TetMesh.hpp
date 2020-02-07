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

template <size_t MaxElementAdjacencies, size_t MaxNodeAdjacencies,
          size_t NodesPerFace, size_t InternalNodes>
struct ElementInfo
{
    std::array<size_t, 3> control_nodes;
    std::array<size_t, 3> faces;
    std::array<size_t, NodesPerFace> face_nodes;
    std::array<size_t, InternalNodes> internal_nodes;
    smv::SmallVector<size_t, MaxElementAdjacencies> adjacent_elements;
    smv::SmallVector<size_t, MaxNodeAdjacencies> adjacent_nodes;

    ElementInfo(size_t a, size_t b, size_t c) noexcept : control_nodes{a, b, c},
        faces{}, face_nodes{}, internal_nodes{}, adjacent_elements(), adjacent_nodes()
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
    typedef ElementInfo<MaxElementAdjacencies, MaxNodeAdjacencies,
                        NodesPerFace, InternalNodes> el_type;
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
            m_nodes.push_back(std::array<CoordT, 2>{get<0>(node), get<1>(node)});
        }

        for (const auto &tet: tets)
        {
            auto n1 = get<0>(tet);
            auto n2 = get<1>(tet);
            auto n3 = get<2>(tet);

            auto det = (m_nodes[n3][0] - m_nodes[n1][0]) * (m_nodes[n2][1] - m_nodes[n1][1]) -
                (m_nodes[n2][0] - m_nodes[n1][0]) * (m_nodes[n3][1] - m_nodes[n1][1]);
            
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

        if (NodesPerFace != 0 || InternalNodes != 0)
        {
            assign_face_and_internal_nodes();
        }
        
        build_adjacencies();
    } // constructor

    const el_type &element(size_t i) const noexcept
    {
        return m_elems[i];
    }

private:
    std::vector<std::array<CoordT, 2>> m_nodes;
    std::vector<el_type> m_elems;
    std::vector<std::vector<size_t>> m_bounding_curves;

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
        auto &eladj = m_elems[el].adjacent_elements;
        auto &nodeadj = m_elems[el].adjacent_nodes;
        for (size_t n: m_elems[el].control_nodes)
        {
            if (std::count(nodeadj.begin(), nodeadj.end(), n) == 0)
            {
                nodeadj.push_back(n);
            }

            for (size_t neighbor: node_neighbors[n])
            {
                if (neighbor == el) { continue; }
                else if (std::count(eladj.begin(), eladj.end(), neighbor) == 0)
                {
                    eladj.push_back(neighbor);
                }

                for (size_t n2: m_elems[neighbor].control_nodes)
                {
                    if (std::count(nodeadj.begin(), nodeadj.end(), n2) == 0)
                    {
                        nodeadj.push_back(n2);
                    }
                }
            }
        }
    }

    void assign_face_and_internal_nodes()
    {}

    void assign_face_numbers()
    {}
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test constructing a tet mesh w/ its adjacencies")
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

    const TetMesh<int, 1, 4> mesh(nodes, tets, std::vector<std::vector<int>>());
    REQUIRE(mesh.element(0).control_nodes == std::array<size_t, 3>{0, 1, 3});
    REQUIRE(mesh.element(1).control_nodes == std::array<size_t, 3>{3, 1, 2});

    auto eladj = mesh.element(0).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element(1).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    auto nodeadj = mesh.element(0).adjacent_nodes;
    REQUIRE(nodeadj.size() == 4);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 2);
    REQUIRE(nodeadj[3] == 3);

    nodeadj = mesh.element(1).adjacent_nodes;
    REQUIRE(nodeadj.size() == 4);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 2);
    REQUIRE(nodeadj[3] == 3);
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace msh

#endif // TETMESH_HPP
