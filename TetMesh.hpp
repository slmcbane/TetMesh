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

/*
 * This file contains code adapted from the standard library of the Julia
 * language (https://julialang.org) in the implementation of
 * `TetMesh::ArrayHasher`.
 * 
 * This code is
 * Copyright (c) 2009-2019: Jeff Bezanson, Stefan Karpinski, Viral B. Shah, and other contributors
 * under the same MIT license terms as above.
 */

#ifndef TETMESH_HPP
#define TETMESH_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

using std::size_t;
using std::get;

#include "SmallVector.hpp"

namespace msh
{

template <class Iter>
struct ReversedBoundaryIterator : public std::reverse_iterator<Iter>
{
    using std::reverse_iterator<Iter>::reverse_iterator;
};

template <class Iter>
struct is_reversed_iterator : public std::false_type
{};

template <class Iter>
struct is_reversed_iterator<ReversedBoundaryIterator<Iter>> : public std::true_type
{};

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
            face.fill(static_cast<size_t>(-1));
        }
        internal_nodes.fill(static_cast<size_t>(-1));
        std::fill(internal_nodes.begin(), internal_nodes.end(), static_cast<size_t>(-1));
    }
};

template <class CoordT, size_t NodesPerFace>
struct BoundaryRepresentation
{
    struct NodeDetails
    {
        size_t number;
        std::array<CoordT, 2> coords;

        NodeDetails(size_t n, const std::array<CoordT, 2> &c) noexcept :
            number{n}, coords(c)
        {}
    };

    std::vector<NodeDetails> nodes;

    struct FaceDetails
    {
        size_t number;
        size_t element;
        std::array<size_t, 2+NodesPerFace> nodes;

        FaceDetails(size_t n, size_t e, const std::array<size_t, 2+NodesPerFace> &d) noexcept :
            number{n}, element{e}, nodes(d)
        {}
    };

    std::vector<FaceDetails> faces;
};

enum class BoundaryError
{
    LessThanTwoNodes,
    FaceIsNotValid,
    FaceIsInternal
};

class BoundaryException : public std::exception
{
public:
    BoundaryException(BoundaryError c, const char *msg) noexcept : code(c), m_msg(msg)
    {}

    const char *what() const noexcept
    {
        return m_msg;
    }

    BoundaryError code;

private:
    const char   *m_msg;
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
    typedef BoundaryRepresentation<CoordT, NodesPerFace> brep_type;

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
            const CurveContainer &bounding_curves) : m_boundaries()
    {
        m_coords.reserve(nodes.size());
        m_elems.reserve(tets.size());

        for (const auto &node: nodes)
        {
            m_coords.push_back(
                std::array<CoordT, 2>{ static_cast<CoordT>(get<0>(node)),
                                       static_cast<CoordT>(get<1>(node)) });
        }

        for (const auto &tet: tets)
        {
            auto n1 = get<0>(tet);
            auto n2 = get<1>(tet);
            auto n3 = get<2>(tet);

            auto det = (m_coords[n3][0] - m_coords[n1][0]) *
                       (m_coords[n2][1] - m_coords[n1][1]) -
                       (m_coords[n2][0] - m_coords[n1][0]) *
                       (m_coords[n3][1] - m_coords[n1][1]);
            
            if (det > 0)
            {
                m_elems.emplace_back(n1, n2, n3);
            }
            else
            {
                m_elems.emplace_back(n2, n1, n3);
            }
        }
        
        find_adjacent_nodes_and_elements();
        const auto node_faces = assign_face_numbers();

        if (NodesPerFace != 0 || InternalNodes != 0)
        {
            assign_face_and_internal_nodes();
        }

        for (const auto &curve: bounding_curves)
        {
            if (curve.size() < 2)
            {
                throw BoundaryException(BoundaryError::LessThanTwoNodes,
                                        "Boundary curve specified with only 1 node");
            }
            m_boundaries.push_back(
                build_boundary_representation(curve.begin(), curve.end(), node_faces));
        }
    } // constructor

    const el_type &element(size_t i) const
    {
        return m_elems.at(i);
    }

    const std::array<CoordT, 2> &coordinate(size_t i) const
    {
        return m_coords.at(i);
    }

    const smv::SmallVector<size_t, MaxNodeAdjacencies> &adjacent_nodes(size_t m) const
    {
        return m_node_adjacencies.at(m);
    }

    const BoundaryRepresentation<CoordT, NodesPerFace> &boundary(size_t i) const
    {
        return m_boundaries.at(i);
    }

    size_t num_vertices() const noexcept { return m_coords.size(); }
    size_t num_elements() const noexcept { return m_elems.size(); }
    size_t num_nodes() const noexcept { return m_node_adjacencies.size(); }
    size_t num_boundaries() const noexcept { return m_boundaries.size(); }

    double average_bandwidth() const
    {
        std::vector<size_t> max_dist(num_nodes(), 0);

        for (size_t i = 0; i < num_nodes(); ++i)
        {
            for (size_t j: adjacent_nodes(i))
            {
                size_t diff = i > j ? i - j : j - i;
                max_dist[i] = (max_dist[i] >= diff) * max_dist[i] +
                    (diff > max_dist[i]) * diff;
            }
        }

        size_t total = 0;
        for (size_t d: max_dist)
        {
            total += d;
        }

        return total / static_cast<double>(num_nodes());
    }

private:
    std::vector<std::array<CoordT, 2>> m_coords;
    std::vector<smv::SmallVector<size_t, MaxNodeAdjacencies>> m_node_adjacencies;
    std::vector<el_type> m_elems;
    std::vector<BoundaryRepresentation<CoordT, NodesPerFace>> m_boundaries;

    struct NodeFaceInfo
    {
        size_t number;
        size_t other_node;
        smv::SmallVector<size_t, 2> elements;

        NodeFaceInfo(size_t num, size_t other, size_t el) noexcept :
            number{num}, other_node{other}, elements()
        {
            elements.push_back(el);
        }

        NodeFaceInfo() = default;
    };

    void find_adjacent_nodes_and_elements()
    {
        m_node_adjacencies.resize(m_coords.size());
        std::vector<smv::SmallVector<size_t, MaxElementAdjacencies+1>> node_neighbors(m_coords.size());

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
        auto &adj = m_node_adjacencies.at(m);
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
        size_t node_number = m_coords.size();

        for (auto &el: m_elems)
        {
            for (size_t facei = 0; facei < 3; ++facei)
            {
                auto &face = el.face_nodes[facei];
                if (face[0] != static_cast<size_t>(-1)) { continue; }

                for (size_t i = 0; i < NodesPerFace; ++i)
                {
                    face[i] = node_number++;
                    m_node_adjacencies.emplace_back();
                }
            }

            for (size_t i = 0; i < InternalNodes; ++i)
            {
                el.internal_nodes[i] = node_number++;
                m_node_adjacencies.emplace_back();
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

    std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>>
    assign_face_numbers()
    {
        size_t face_number = 0;

        std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>>
            node_faces(m_coords.size());

        for (size_t eli = 0; eli < m_elems.size(); ++eli)
        {
            el_type &el = m_elems[eli];
            bool face_numbered[] = { false, false, false };

            // Discover already assigned faces
            for (auto &face: node_faces[el.control_nodes[1]])
            {
                if (face.other_node == el.control_nodes[0])
                {
                    el.faces[0] = face.number;
                    face_numbered[0] = true;
                    face.elements.push_back(eli);
                }
            }

            for (auto &face: node_faces[el.control_nodes[2]])
            {
                if (face.other_node == el.control_nodes[1])
                {
                    el.faces[1] = face.number;
                    face_numbered[1] = true;
                    face.elements.push_back(eli);
                }
            }

            for (auto &face: node_faces[el.control_nodes[0]])
            {
                if (face.other_node == el.control_nodes[2])
                {
                    el.faces[2] = face.number;
                    face_numbered[2] = true;
                    face.elements.push_back(eli);
                }
            }

            // For faces not assigned, assign numbers and do the bookkeeping.
            if (!face_numbered[0])
            {
                el.faces[0] = face_number;
                node_faces[el.control_nodes[0]].emplace_back(
                    face_number, el.control_nodes[1], eli
                );
                face_number += 1;
            }

            if (!face_numbered[1])
            {
                el.faces[1] = face_number;
                node_faces[el.control_nodes[1]].emplace_back(
                    face_number, el.control_nodes[2], eli
                );
                face_number += 1;
            }

            if (!face_numbered[2])
            {
                el.faces[2] = face_number;
                node_faces[el.control_nodes[2]].emplace_back(
                    face_number, el.control_nodes[0], eli
                );
                face_number += 1;
            }
        }

        return node_faces;
    }

    template <class Iterator>
    typename std::enable_if<!is_reversed_iterator<Iterator>::value,
                            BoundaryRepresentation<CoordT, NodesPerFace>>::type
    build_boundary_representation(const Iterator &begin, const Iterator &end,
        const std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>> &node_faces) const
    {
        auto reversed = check_for_reversal(begin, node_faces);

        // First face is not a valid face.
        if (reversed.first)
        {
            // But it does exist.
            if (reversed.second)
            {
                return build_boundary_representation(ReversedBoundaryIterator<Iterator>(end),
                                                     ReversedBoundaryIterator<Iterator>(begin),
                                                     node_faces);
            }

            // Face does not exist.
            else
            {
                throw BoundaryException(BoundaryError::FaceIsNotValid,
                    "Given nodes do not bound a face in the mesh");
            }
        }

        // At this point, I know that at least the first face is valid.
        brep_type boundary;
        size_t starting_dof = 0;
        for (auto it = begin; (it+1) != end; ++it)
        {
            add_face_to_boundary(*it, *(it+1), boundary, node_faces, starting_dof);
            starting_dof += (1 + NodesPerFace);
        }

        return boundary;
    }

    template <class Iterator>
    typename std::enable_if<is_reversed_iterator<Iterator>::value,
                            BoundaryRepresentation<CoordT, NodesPerFace>>::type
    build_boundary_representation(const Iterator &begin, const Iterator &end,
        const std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>> &node_faces) const
    {
        auto reversed = check_for_reversal(begin, node_faces);

        // Not a valid face
        if (reversed.first)
        {
            // Face exists, there has to be an internal face for this condition to be true.
            if (reversed.second)
            {
                throw BoundaryException(BoundaryError::FaceIsInternal,
                    "Found face internal to the mesh while constructing boundary representation");
            }

            // Face does not exist.
            else
            {
                throw BoundaryException(BoundaryError::FaceIsNotValid,
                    "Given nodes do not bound a face in the mesh");
            }
        }

        brep_type boundary;
        size_t starting_dof = 0;
        for (auto it = begin; (it+1) != end; ++it)
        {
            add_face_to_boundary(*it, *(it+1), boundary, node_faces, starting_dof);
            starting_dof += (1 + NodesPerFace);
        }

        return boundary;
    }

    template <class Iterator>
    std::pair<bool, bool> check_for_reversal(const Iterator &begin,
        const std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>> &node_faces) const
    {
        size_t first_node = *begin;
        size_t second_node = *(begin+1);
        for (const auto &face: node_faces[first_node])
        {
            if (face.other_node == second_node)
            {
                return std::make_pair(false, true);
            }
        }

        for (const auto &face: node_faces[second_node])
        {
            if (face.other_node == first_node)
            {
                return std::make_pair(true, true);
            }
        }

        return std::make_pair(true, false);
    }

    void add_face_to_boundary(size_t n1, size_t n2, brep_type &boundary,
        const std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>> &node_faces,
        size_t starting_dof) const
    {
        // Exception gets thrown if face doesn't exist, or if it is an internal
        // face.
        const NodeFaceInfo &face_info = find_face(n1, n2, node_faces);

        std::array<size_t, NodesPerFace+2> dofs;
        for (auto &dof: dofs)
        {
            dof = starting_dof++;
        }

        boundary.faces.emplace_back(face_info.number, face_info.elements[0], dofs);

        if (boundary.nodes.size() == 0)
        {
            boundary.nodes.emplace_back(n1, m_coords[n1]);
        }

        if (NodesPerFace != 0)
        {
            const auto &face_nodes = get_face_nodes(face_info);

            std::array<CoordT, 2> delta;
            for (size_t i = 0; i < 2; ++i)
            {
                delta[i] = m_coords[n2][i] - m_coords[n1][i];
            }
        
            std::array<CoordT, 2> coord;
            constexpr auto alpha = 1.0 / (NodesPerFace + 1);
            for (size_t i = 0; i < NodesPerFace; ++i)
            {
                coord[0] = m_coords[n1][0] + (i+1) * alpha * m_coords[n2][0]
                    - (i+1) * alpha * m_coords[n1][0];
                coord[1] = m_coords[n1][1] + (i+1) * alpha * m_coords[n2][1]
                    - (i+1) * alpha * m_coords[n1][1];
                boundary.nodes.emplace_back(face_nodes[i], coord);
            }
        }

        boundary.nodes.emplace_back(n2, m_coords[n2]);
    }

    const NodeFaceInfo &find_face(size_t n1, size_t n2,
        const std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>> &node_faces) const
    {
        const auto &adjacent = node_faces[n1];
        for (const auto &face: adjacent)
        {
            if (face.other_node == n2)
            {
                if (face.elements.size() != 1)
                {
                    throw BoundaryException(
                        BoundaryError::FaceIsInternal,
                        "Found face internal to the mesh while constructing boundary representation"
                    );
                }
                return face;
            }
        }

        throw BoundaryException(
            BoundaryError::FaceIsNotValid,
            "Did not find face bounded by two nodes in boundary spec"
        );
    }

    const std::array<size_t, NodesPerFace> &get_face_nodes(const NodeFaceInfo &face_info) const
    {
        size_t eli = face_info.elements[0];
        const auto &el = m_elems.at(eli);

        if (face_info.number == el.faces[0])
        {
            return el.face_nodes[0];
        }
        else if (face_info.number == el.faces[1])
        {
            return el.face_nodes[1];
        }
        else if (face_info.number == el.faces[2])
        {
            return el.face_nodes[2];
        }
        else
        {
            throw BoundaryException(
                BoundaryError::FaceIsNotValid,
                "Did not find face in faces of adjacent element - this should never happen"
            );
        }
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

    const std::array<std::array<size_t, 3>, 1> boundaries = { 3, 2, 1 };

    const TetMesh<int, 1, 3> mesh(nodes, tets, boundaries);

    REQUIRE(mesh.average_bandwidth() == doctest::Approx(2.25));

    REQUIRE(mesh.element(0).control_nodes == std::array<size_t, 3>{0, 1, 3});
    REQUIRE(mesh.element(1).control_nodes == std::array<size_t, 3>{3, 1, 2});

    auto eladj = mesh.element(0).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element(1).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    auto nodeadj = mesh.adjacent_nodes(0);
    REQUIRE(nodeadj.size() == 2);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 3);

    nodeadj = mesh.adjacent_nodes(1);
    REQUIRE(nodeadj.size() == 3);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 2);
    REQUIRE(nodeadj[2] == 3);

    nodeadj = mesh.adjacent_nodes(2);
    REQUIRE(nodeadj.size() == 2);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 3);

    nodeadj = mesh.adjacent_nodes(3);
    REQUIRE(nodeadj.size() == 3);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 2);

    REQUIRE(mesh.element(0).faces == std::array<size_t, 3>{0, 1, 2});
    REQUIRE(mesh.element(1).faces == std::array<size_t, 3>{1, 3, 4});

    const auto &bound = mesh.boundary(0);
    REQUIRE(bound.nodes.size() == 3);
    REQUIRE(bound.nodes[0].number == 1);
    REQUIRE(bound.nodes[1].number == 2);
    REQUIRE(bound.nodes[2].number == 3);
    REQUIRE(bound.nodes[0].coords == std::array<int, 2>{ -1, 1 });
    REQUIRE(bound.nodes[1].coords == std::array<int, 2>{ 1, 1 });
    REQUIRE(bound.nodes[2].coords == std::array<int, 2>{ 1, -1 });

    REQUIRE(bound.faces.size() == 2);
    REQUIRE(bound.faces[0].number == 3);
    REQUIRE(bound.faces[0].element == 1);
    REQUIRE(bound.faces[0].nodes == std::array<size_t, 2>{ 0, 1 });
    REQUIRE(bound.faces[1].number == 4);
    REQUIRE(bound.faces[1].element == 1);
    REQUIRE(bound.faces[1].nodes == std::array<size_t, 2>{ 1, 2 });

    // Check error conditions to make sure they throw the correct error.
    auto make_mesh = [=](std::array<std::array<size_t, 3>, 1> b)
    {
        return TetMesh<int, 1, 4>(nodes, tets, b);
    };

    // This one is a boundary specification containing a face that does not exist.
    try
    {
        auto m = make_mesh(std::array<std::array<size_t, 3>, 1>{ 1, 2, 0 });
        REQUIRE(false);
    }
    catch(const BoundaryException &exc)
    {
        REQUIRE(exc.code == BoundaryError::FaceIsNotValid);
    }

    // The three below are three different possible specifications of a boundary
    // with a face that is internal to the mesh.
    try
    {
        auto m = make_mesh(std::array<std::array<size_t, 3>, 1>{ 0, 1, 3 });
        REQUIRE(false);
    }
    catch(const BoundaryException& exc)
    {
        REQUIRE(exc.code == BoundaryError::FaceIsInternal);
    }

    try
    {
        auto m = make_mesh(std::array<std::array<size_t, 3>, 1>{ 2, 1, 3 });
        REQUIRE(false);
    }
    catch(const BoundaryException& exc)
    {
        REQUIRE(exc.code == BoundaryError::FaceIsInternal);
    }

    try
    {
        auto m = make_mesh(std::array<std::array<size_t, 3>, 1>{ 1, 3, 2 });
        REQUIRE(false);
    }
    catch(const BoundaryException& exc)
    {
        REQUIRE(exc.code == BoundaryError::FaceIsInternal);
    }
} // TEST_CASE

TEST_CASE("Test constructing a third order mesh")
{
    const std::vector<std::array<double, 2>> nodes = {
        std::array<double, 2>{-1, -1},
        std::array<double, 2>{-1, 1},
        std::array<double, 2>{1, 1},
        std::array<double, 2>{1, -1}
    };

    const std::vector<std::array<int, 3>> tets = {
        std::array<int, 3>{0, 1, 3},
        std::array<int, 3>{1, 3, 2}
    };

    const std::array<std::array<size_t, 3>, 1> boundaries = { 1, 2, 3 };

    const TetMesh<double, 1, 15, 2, 1> mesh(nodes, tets, boundaries);

    REQUIRE(mesh.average_bandwidth() == doctest::Approx(164.0 / 16));

    auto eladj = mesh.element(0).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element(1).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    auto nodeadj = mesh.adjacent_nodes(1);
    REQUIRE(nodeadj.size() == 15);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    for (size_t i = 1; i < 15; ++i)
    {
        CHECK(nodeadj[i] == i+1);
    }

    nodeadj = mesh.adjacent_nodes(15);
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

    nodeadj = mesh.adjacent_nodes(8);
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

    const auto &bound = mesh.boundary(0);
    REQUIRE(bound.nodes.size() == 7);
    REQUIRE(bound.nodes[0].number == 1);
    REQUIRE(bound.nodes[1].number == 11);
    REQUIRE(bound.nodes[2].number == 12);
    REQUIRE(bound.nodes[3].number == 2);
    REQUIRE(bound.nodes[4].number == 13);
    REQUIRE(bound.nodes[5].number == 14);
    REQUIRE(bound.nodes[6].number == 3);

    REQUIRE(bound.nodes[0].coords == std::array<double, 2>{-1, 1});
    REQUIRE(bound.nodes[1].coords[0] == doctest::Approx(-1.0 / 3));
    REQUIRE(bound.nodes[1].coords[1] == doctest::Approx(1.0));
    REQUIRE(bound.nodes[2].coords[0] == doctest::Approx(1.0 / 3));
    REQUIRE(bound.nodes[2].coords[1] == doctest::Approx(1.0));
    REQUIRE(bound.nodes[3].coords == std::array<double, 2>{1, 1});
    REQUIRE(bound.nodes[4].coords[0] == doctest::Approx(1.0));
    REQUIRE(bound.nodes[4].coords[1] == doctest::Approx(1.0 / 3));
    REQUIRE(bound.nodes[5].coords[0] == doctest::Approx(1.0));
    REQUIRE(bound.nodes[5].coords[1] == doctest::Approx(-1.0 / 3));
    REQUIRE(bound.nodes[6].coords == std::array<double, 2>{1, -1});

    REQUIRE(bound.faces.size() == 2);
    REQUIRE(bound.faces[0].number == 3);
    REQUIRE(bound.faces[0].element == 1);
    REQUIRE(bound.faces[0].nodes == std::array<size_t, 4>{ 0, 1, 2, 3 });
    REQUIRE(bound.faces[1].number == 4);
    REQUIRE(bound.faces[1].element == 1);
    REQUIRE(bound.faces[1].nodes == std::array<size_t, 4>{ 3, 4, 5, 6 });
} // TEST_CASE

TEST_CASE("A larger third order mesh boundary test")
{
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

    const std::array<std::array<size_t, 5>, 1> boundaries = { 0, 7, 1, 6, 2 };

    const TetMesh<double, 11, 36, 2, 1> mesh(nodes, tets, boundaries);
    
    REQUIRE(mesh.num_boundaries() == 1);
    const auto &boundary = mesh.boundary(0);
    REQUIRE(boundary.nodes.size() == 13);
    REQUIRE(boundary.nodes[0].number == 2);
    REQUIRE(boundary.nodes[1].number == 26);
    REQUIRE(boundary.nodes[2].number == 27);
    REQUIRE(boundary.nodes[3].number == 6);
    REQUIRE(boundary.nodes[4].number == 52);
    REQUIRE(boundary.nodes[5].number == 53);
    REQUIRE(boundary.nodes[6].number == 1);
    REQUIRE(boundary.nodes[7].number == 57);
    REQUIRE(boundary.nodes[8].number == 58);
    REQUIRE(boundary.nodes[9].number == 7);
    REQUIRE(boundary.nodes[10].number == 33);
    REQUIRE(boundary.nodes[11].number == 34);
    REQUIRE(boundary.nodes[12].number == 0);

    // Selective tests on the coordinates of nodes.
    REQUIRE(boundary.nodes[10].coords[0] == doctest::Approx(-1.0));
    REQUIRE(boundary.nodes[10].coords[1] == doctest::Approx(1.0 / 3));
    REQUIRE(boundary.nodes[4].coords[0] == doctest::Approx(-1.0 / 3));
    REQUIRE(boundary.nodes[5].coords[0] == doctest::Approx(-2.0 / 3));

    try
    {
        auto m = TetMesh<double, 11, 36, 2, 1>(nodes, tets,
            std::array<std::array<size_t, 4>, 1>{4, 10, 7, 0});
        REQUIRE(false);
    }
    catch (const BoundaryException &exc)
    {
        REQUIRE(exc.code == BoundaryError::FaceIsInternal);
    }
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace msh

#endif // TETMESH_HPP
