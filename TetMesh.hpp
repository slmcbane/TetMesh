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
template <class CoordT, size_t MaxElementAdjacencies, size_t MaxNodeAdjacencies>
class TetMesh
{
public:
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
        m_tets.reserve(tets.size());
        m_bounding_curves.reserve(bounding_curves.size());

        for (const auto &node: nodes)
        {
            m_nodes.push_back(std::array<CoordT, 2>{get<0>(node), get<1>(node)});
        }

        for (const auto &tet: tets)
        {
            m_tets.push_back(
                std::array<size_t, 3>{
                    static_cast<size_t>(get<0>(tet)),
                    static_cast<size_t>(get<1>(tet)),
                    static_cast<size_t>(get<2>(tet))
                }
            );
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
    } // constructor

    /*
     * For efficiency, this constructor does the same as above but takes rvalue
     * references to vectors for the arguments so that there isn't a bunch of
     * copying involved.
     *
     * Note that the constructor will only disambiguate to this version if you
     * make sure that ALL THREE arguments are rvalues (use std::move) and the
     * types match exactly.
     */
    TetMesh(std::vector<std::array<CoordT, 2>> &&nodes,
            std::vector<std::array<size_t, 3>> &&tets,
            std::vector<std::vector<size_t>> &&m_bounding_curves) :
        m_nodes(std::move(nodes)), m_tets(std::move(tets)),
        m_bounding_curves(std::move(m_bounding_curves))
    {
        build_adjacencies();
    }

    /*
     * Get a list of elements that are adjacent to element `i`.
     *
     * The returned vector contains the indices of all elements that share a
     * node with the element indexed by `i` in the mesh. This list does not
     * include the index `i` itself.
     */
    const smv::SmallVector<size_t, MaxElementAdjacencies> &element_adjacencies(size_t i) const noexcept
    {
        return m_element_adjacencies[i];
    }

    /*
     * Get a list of nodes that are adjacent to element `i`.
     *
     * The returned vector contains the indices of all nodes that either on
     * this element, or any of its neighbors (may be queried with
     * `element_adjacencies`). The entries in this list are unique.
     */
    const smv::SmallVector<size_t, MaxNodeAdjacencies> &node_adjacencies(size_t i) const noexcept
    {
        return m_node_adjacencies[i];
    }

private:
    std::vector<std::array<CoordT, 2>> m_nodes;
    std::vector<std::array<size_t, 3>> m_tets;
    std::vector<smv::SmallVector<size_t, MaxNodeAdjacencies>> m_node_adjacencies;
    std::vector<smv::SmallVector<size_t, MaxElementAdjacencies>> m_element_adjacencies;
    std::vector<std::vector<size_t>> m_bounding_curves;

    void build_adjacencies()
    {
        std::vector<smv::SmallVector<size_t, MaxElementAdjacencies+1>> els_neighboring(m_nodes.size());
        for (size_t el = 0; el < m_tets.size(); ++el)
        {
            for (size_t n: m_tets[el])
            {
                els_neighboring[n].push_back(el);
            }
        }

        m_node_adjacencies.resize(m_tets.size());
        m_element_adjacencies.resize(m_tets.size());

        for (size_t el = 0; el < m_tets.size(); ++el)
        {
            process_adjacencies(el, els_neighboring);
        }
    }

    void process_adjacencies(size_t el,
        const std::vector<smv::SmallVector<size_t, MaxElementAdjacencies+1>> &els_neighboring)
    {
        auto &eladj = m_element_adjacencies[el];
        auto &nodeadj = m_node_adjacencies[el];
        for (size_t n: m_tets[el])
        {
            if (std::count(nodeadj.begin(), nodeadj.end(), n) == 0)
            {
                nodeadj.push_back(n);
            }

            for (size_t neighbor: els_neighboring[n])
            {
                if (neighbor == el) { continue; }
                else if (std::count(eladj.begin(), eladj.end(), neighbor) == 0)
                {
                    eladj.push_back(neighbor);
                }

                for (size_t n2: m_tets[neighbor])
                {
                    if (std::count(nodeadj.begin(), nodeadj.end(), n2) == 0)
                    {
                        nodeadj.push_back(n2);
                    }
                }
            }
        }
    }
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
        std::array<int, 3>{3, 1, 2}
    };

    const TetMesh<int, 1, 4> mesh(nodes, tets, std::vector<std::vector<int>>());

    auto eladj = mesh.element_adjacencies(0);
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element_adjacencies(1);
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    auto nodeadj = mesh.node_adjacencies(0);
    REQUIRE(nodeadj.size() == 4);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 2);
    REQUIRE(nodeadj[3] == 3);

    nodeadj = mesh.node_adjacencies(1);
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
