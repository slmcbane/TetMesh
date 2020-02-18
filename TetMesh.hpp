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
#include <cstdio>
#include <cstdint>
#include <exception>
#include <iterator>
#include <queue>
#include <type_traits>
#include <utility>
#include <vector>

using std::size_t;
using std::get;

#include "SmallVector.hpp"
#include "read_gmsh.hpp"

namespace msh
{

// This wrapper around reverse_iterator is used to select which version
// of build_boundary_representation to call using std::enable_if.
template <class Iter>
struct ReversedBoundaryIterator : public std::reverse_iterator<Iter>
{
    using std::reverse_iterator<Iter>::reverse_iterator;
};

// Trait to identify a ReversedBoundaryIterator
template <class Iter>
struct is_reversed_iterator : public std::false_type
{};
template <class Iter>
struct is_reversed_iterator<ReversedBoundaryIterator<Iter>> : public std::true_type
{};

/*
 * The information stored on a per node basis in the mesh.
 * 
 * - coords: The real number coordinates of the node in the X-Y plane
 * - adjacent_nodes: A list of all nodes that share an element with this one.
 * 
 * Construct either with no coordinate info (zero-initialized), or from two
 * coordinates X and Y; either way, the list of adjacent nodes is empty.
 */
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
 * Data structure encapsulating all informatino the mesh stores about each
 * triangle element.
 * 
 * - control_nodes: The three nodes defining the triangle, in clockwise order.
 * - faces: The indices of the faces on the element in the mesh. Note that the
 *          mesh does not actually store any information about faces, only these
 *          numbers; they are used when building the boundary representation and
 *          assigning numbers to nodes on faces.
 * - face_nodes: face_nodes[i] is a size NodesPerFace array with the indices of
 *               the nodes on that face, in clockwise order. Note that the mesh
 *               assumes these are equispaced on the face and assigns coordinates
 *               to them as such.
 * - internal_nodes: A list of the indices of nodes internal to the element, in
 *                   no particular order. It is up to the user to interpret these.
 * - adjacent_elements: A list of all other elements in the mesh that share at
 *                      least one node with this one.
 * 
 * Construct from 3 indices defining the control nodes; all other fields are
 * default initialized.
 */
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

/*
 * The data structure used by the mesh to identify segments of the boundary.
 * 
 * - nodes: A list of the indices of nodes on the boundary, in clockwise order.
 *          Every adjacent pair of nodes should define a face in the mesh, which
 *          guarantees that none are skipped.
 * - faces: Each face stores *its number in the mesh, *the element it belongs to
 *          (there can be only 1 since this is a boundary), and *an array of the
 *          indices of the nodes on this face IN THE BOUNDARY REPRESENTATION.
 * 
 * Thus for face `i in a representation `repr`, the index of the `j-th` node on
 * this face in the boundary representation is `face[i].nodes[j]`, but the index
 * __in the parent mesh__ is `repr.nodes[face[i].nodes[j]]`.
 */
template <size_t NodesPerFace>
struct BoundaryRepresentation
{
    std::vector<size_t> nodes;

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

// Error codes for exceptions thrown when building boundary representation.
enum class BoundaryError
{
    LessThanTwoNodes,
    FaceIsNotValid,
    FaceIsInternal
};

// Thrown during boundary construction in some cases.
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
 * Note that there is not extensive checking on the geometry; there could be
 * disjoint parts of the mesh, etc. that are not caught when initializing. What
 * is checked is that each triangle is defined by a set of control nodes in
 * clockwise order (even if the given nodes were counterclockwise), and that
 * boundary curves are actually on a boundary and are composed of a connected set
 * of faces, again specified in clockwise order.
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
 * When setting `MaxElementAdjacencies` and `MaxNodeAdjacencies`, note that
 * an exception *will* be thrown if these values are exceeded. It is set as
 * a template parameter to optimize storage space; just increase it as
 * needed and recompile. The performance penalty for extra space should be
 * reasonable.
 * 
 * I have put exception handling code to give you a nice error message letting
 * you know exactly where your mesh went wrong when these exceptions are thrown.
 */
template <class CoordT, size_t MaxElementAdjacencies, size_t MaxNodeAdjacencies,
          size_t NodesPerFace = 0, size_t InternalNodes = 0>
class TetMesh
{
public:
    typedef ElementInfo<MaxElementAdjacencies, NodesPerFace, InternalNodes> el_type;
    typedef BoundaryRepresentation<NodesPerFace> brep_type;
    typedef NodeInfo<CoordT, MaxNodeAdjacencies> node_type;

    /*
     * Construct the TetMesh from a list of nodes, list of elements, and
     * list of curves.
     *
     * All three container arguments must support `size()` and the STL iterator
     * interface. The elements in `nodes` and `tets` should support `get<i>` for
     * element access, and elements of `bounding_curves` should be lists of
     * node indices in a format that allows indexing with `[]`.
     * 
     * The initialization does the following:
     *   - Add all node coordinates to the mesh and instantiate triangles,
     *     checking that triangle vertices are given clockwise. They are switched
     *     to clockwise order if not.
     *   - Build the arrays of nodes adjacent to each node and elements adjacent
     *     to each element
     *   - Assign a single, unique number to each face in the mesh
     *     (exposed in the public interface as element(i).faces)
     *   - If either of NodesPerFace or InternalNodes is not 0, assign numbers
     *     to each of these extra nodes. This is simply done in ascending order
     *     by element. For nodes on faces, the code assumes they are equally
     *     spaced and assigns coordinates as such. For internal nodes there is
     *     not really a corresponding assumption; their coordinates are all set
     *     as (0, 0).
     *   - Finally, build the representation of each boundary given. This
     *     operation throws an exception if
     *       * only one node is given on the boundary
     *       * two adjacent nodes on the boundary don't form a face contained in
     *         the mesh
     *       * The specified nodes don't form a collection of faces on the
     *         boundary in clockwise order, and also don't do so when reversed.
     *         I believe if all are faces in the mesh, this is only possible if
     *         a face is internal to the mesh
     *       * A face internal to the mesh is contained in the 'boundary'.
     */
    template <class NodeContainer, class ElContainer, class CurveContainer>
    TetMesh(const NodeContainer &nodes, const ElContainer &tets,
            const CurveContainer &bounding_curves) : m_boundaries()
    {
        m_nodes.reserve(nodes.size());
        m_elems.reserve(tets.size());

        for (const auto &node: nodes)
        {
            m_nodes.emplace_back(
                static_cast<CoordT>(get<0>(node)),
                static_cast<CoordT>(get<1>(node))
            );
        }

        for (const auto &tet: tets)
        {
            auto n1 = get<0>(tet);
            auto n2 = get<1>(tet);
            auto n3 = get<2>(tet);

            auto det = (coord(n3)[0] - coord(n1)[0]) * (coord(n2)[1] - coord(n1)[1]) -
                       (coord(n2)[0] - coord(n1)[0]) * (coord(n3)[1] - coord(n1)[1]);
            
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

        m_boundaries.reserve(bounding_curves.size());
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

        remove_orphaned_nodes();
    } // constructor

    const el_type &element(size_t i) const
    {
        return m_elems.at(i);
    }

    const node_type &node(size_t i) const
    {
        return m_nodes.at(i);
    }

    // Coordinates of the i-th node
    const std::array<CoordT, 2> &coord(size_t i) const
    {
        return node(i).coords;
    }

    // All nodes adjacent to the i-th node
    const smv::SmallVector<size_t, MaxNodeAdjacencies> &adjacent_nodes(size_t i) const
    {
        return node(i).adjacent_nodes;
    }

    const BoundaryRepresentation<NodesPerFace> &boundary(size_t i) const
    {
        return m_boundaries.at(i);
    }

    size_t num_elements() const noexcept { return m_elems.size(); }
    size_t num_nodes() const noexcept { return m_nodes.size(); }
    size_t num_boundaries() const noexcept { return m_boundaries.size(); }

    // The average bandwidth of a finite element matrix assembled from this mesh.
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

    // Heuristically renumber the nodes of the mesh using a greedy algorithm, to
    // attempt to minimize the average bandwidth of a finite element matrix.
    void renumber_nodes()
    {
        constexpr auto sentinel = static_cast<size_t>(-1);
        std::vector<size_t> new_numbers(m_nodes.size(), sentinel);

        // This block computes mapping of node numbers, as well as updating
        // the node numbers stored in elements.
        {
            std::vector<bool> processed(num_elements(), false);
            std::queue<size_t> to_process(std::deque<size_t>(1, 0));
            processed[0] = true;

            size_t node_number = 0;

            while (!(to_process.size() == 0))
            {
                auto &el = m_elems[to_process.front()];
                to_process.pop();
                for (size_t other: el.adjacent_elements)
                {
                    if (!processed[other])
                    {
                        to_process.push(other);
                        processed[other] = true;
                    }
                }

                for (auto &node: el.control_nodes)
                {
                    // Not reassigned yet.
                    if (new_numbers[node] == sentinel)
                    {
                        new_numbers[node] = node_number++;
                    }
                    node = new_numbers[node];
                }

                for (auto &face: el.face_nodes)
                {
                    for (auto &node: face)
                    {
                        if (new_numbers[node] == sentinel)
                        {
                            new_numbers[node] = node_number++;
                        }
                        node = new_numbers[node];
                    }
                }

                for (auto &node: el.internal_nodes)
                {
                    if (new_numbers[node] == sentinel)
                    {
                        new_numbers[node] = node_number++;
                    }
                    node = new_numbers[node];
                }
            }
        }

        // Now new_numbers stores the new number for every node.
        // Swap node array around.
        {
            std::vector<node_type> new_nodes(m_nodes.size());
            for (size_t i = 0; i < num_nodes(); ++i)
            {
                new_nodes[new_numbers[i]] = node(i);
            }
            std::swap(m_nodes, new_nodes);
        }

        // Loop over nodes and update their adjacencies
        for (auto &node: m_nodes)
        {
            for (size_t &i: node.adjacent_nodes)
            {
                i = new_numbers[i];
            }
        }

        // For each boundary, loop over its node numbers and update. Coordinates
        // do not change.
        for (auto &boundary: m_boundaries)
        {
            for (size_t &node: boundary.nodes)
            {
                node = new_numbers[node];
            }
        }
    }

    std::vector<std::tuple<size_t, size_t, int>> sparsity_pattern() const
    {
        std::vector<std::tuple<size_t, size_t, int>> entries;
        for (const auto &element: m_elems)
        {
            for (size_t node: element.control_nodes)
            {
                for (size_t other: element.control_nodes)
                {
                    entries.emplace_back(node, other, 1);
                }
                for (const auto &face: element.face_nodes)
                {
                    for (size_t other: face)
                    {
                        entries.emplace_back(node, other, 1);
                    }
                }
                for (size_t other: element.internal_nodes)
                {
                    entries.emplace_back(node, other, 1);
                }
            }
            for (const auto &face: element.face_nodes)
            {
                for (size_t node: face)
                {
                    for (size_t other : element.control_nodes)
                    {
                        entries.emplace_back(node, other, 1);
                    }
                    for (const auto &face : element.face_nodes)
                    {
                        for (size_t other : face)
                        {
                            entries.emplace_back(node, other, 1);
                        }
                    }
                    for (size_t other : element.internal_nodes)
                    {
                        entries.emplace_back(node, other, 1);
                    }
                }
            }
            for (size_t node: element.internal_nodes)
            {
                for (size_t other : element.control_nodes)
                {
                    entries.emplace_back(node, other, 1);
                }
                for (const auto &face : element.face_nodes)
                {
                    for (size_t other : face)
                    {
                        entries.emplace_back(node, other, 1);
                    }
                }
                for (size_t other : element.internal_nodes)
                {
                    entries.emplace_back(node, other, 1);
                }
            }
        }
        return entries;
    }

    size_t storage_used() const noexcept
    {
        size_t num_bytes = m_nodes.capacity() * sizeof(node_type);
        num_bytes += m_elems.capacity() * sizeof(el_type);
        num_bytes += m_boundaries.capacity() * sizeof(BoundaryRepresentation<NodesPerFace>);
        for (const auto &boundary: m_boundaries)
        {
            num_bytes += boundary.nodes.capacity() * sizeof(size_t);
            num_bytes += boundary.faces.capacity() * sizeof(
                typename BoundaryRepresentation<NodesPerFace>::FaceDetails
            );
        }
        return num_bytes;
    }

private:
    std::vector<node_type> m_nodes;
    std::vector<el_type> m_elems;
    std::vector<BoundaryRepresentation<NodesPerFace>> m_boundaries;

    node_type &node(size_t i)
    {
        return m_nodes.at(i);
    }

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
        std::vector<smv::SmallVector<size_t, MaxElementAdjacencies+1>> node_neighbors(m_nodes.size());

        for (size_t el = 0; el < m_elems.size(); ++el)
        {
            for (size_t n: m_elems[el].control_nodes)
            {
                try
                {
                    node_neighbors[n].push_back(el);
                }
                catch(const smv::MaxSizeExceeded&)
                {
                    fprintf(stderr, "Exceeded max node neighbor elements at node %ld (MaxElementAdjacencies too small)\n", n);
                    throw smv::MaxSizeExceeded{};
                }
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
        try
        {
            auto &adj = m_nodes.at(m).adjacent_nodes;
            if (std::count(adj.begin(), adj.end(), n) == 0)
            {
                adj.push_back(n);
            }
        }
        catch(const smv::MaxSizeExceeded&)
        {
            fprintf(stderr, "Exceed max node adjacencies at node %ld (MaxNodeAdjacencies too small)",
                    m);
            throw smv::MaxSizeExceeded{};
        }
    }

    void add_element_adjacency(size_t m, size_t n)
    {
        try
        {
            auto &adj = m_elems[m].adjacent_elements;
            if (std::count(adj.begin(), adj.end(), n) == 0)
            {
                adj.push_back(n);
            }
        }
        catch (const smv::MaxSizeExceeded&)
        {
            fprintf(stderr, "Exceeded MaxElementAdjacencies at element %ld\n", m);
            throw smv::MaxSizeExceeded{};
        }
    }

    void remove_orphaned_nodes()
    {
        std::vector<size_t> new_index_offsets(num_nodes(), 0);

        // Sort the nodes to be removed to the end while also updating the
        // offsets of nodes that are after it.
        auto new_node_end = std::remove_if(
            m_nodes.begin(), m_nodes.end(),
            [&](const node_type &node)
            {
                if (node.adjacent_nodes.empty())
                {
                    size_t i = &node - m_nodes.data();
                    std::transform(new_index_offsets.begin()+i,
                                   new_index_offsets.end(),
                                   new_index_offsets.begin()+i,
                                   [](size_t x) { return x - 1; });
                }
                return node.adjacent_nodes.empty();
            }
        );

        m_nodes.erase(new_node_end, m_nodes.end());

        // Update references to nodes.
        auto update_node = [&](size_t &n) { n += new_index_offsets.at(n); };
        for (auto &el: m_elems)
        {
            for (size_t &n: el.control_nodes)
            {
                update_node(n);
            }
            for (auto &face: el.face_nodes)
            {
                for (size_t &n: face)
                {
                    update_node(n);
                }
            }
            for (size_t &n: el.internal_nodes)
            {
                update_node(n);
            }
        }

        for (auto &node: m_nodes)
        {
            for (size_t &n: node.adjacent_nodes)
            {
                update_node(n);
            }
        }

        for (auto &boundary: m_boundaries)
        {
            for (size_t &n: boundary.nodes)
            {
                update_node(n);
            }
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

                // Compute coordinates of the new node.
                size_t n1 = el.control_nodes[facei], n2 = el.control_nodes[(facei+1)%3];
                CoordT x1 = coord(n1)[0], y1 = coord(n1)[1], x2 = coord(n2)[0], y2 = coord(n2)[1];

                for (size_t i = 0; i < NodesPerFace; ++i)
                {
                    face[i] = node_number++;
                    double alpha = static_cast<double>(i+1) / (NodesPerFace+1);
                    m_nodes.emplace_back(x1 + alpha * x2 - alpha * x1,
                                         y1 + alpha * y2 - alpha * y1);
                }
            }

            for (size_t &node: el.internal_nodes)
            {
                node = node_number++;
                // In general can't define the coordinate of nodes internal to
                // the element.
                m_nodes.emplace_back();
            }
            add_face_and_internal_nodes_to_adjacent(el);
        }
        m_nodes.shrink_to_fit();
    }

    void add_face_and_internal_nodes_to_adjacent(const el_type &el)
    {
        for (size_t facei = 0; facei < 3; ++facei)
        {
            const auto &face_nodes = el.face_nodes[facei];
            for (size_t i = 0; i < NodesPerFace; ++i)
            {
                for (size_t j = i+1; j < NodesPerFace; ++j)
                {
                    add_node_adjacency(face_nodes[i], face_nodes[j]);
                    add_node_adjacency(face_nodes[j], face_nodes[i]);
                }
            }
            
            for (size_t facej = facei+1; facej < 3; ++facej)
            {
                const auto &other_face_nodes = el.face_nodes[facej];
                for (size_t i = 0; i < NodesPerFace; ++i)
                {
                    for (size_t j = 0; j < NodesPerFace; ++j)
                    {
                        add_node_adjacency(face_nodes[i], other_face_nodes[j]);
                        add_node_adjacency(other_face_nodes[j], face_nodes[i]);
                    }
                }
            }

            for (size_t node: el.control_nodes)
            {
                for (size_t i = 0; i < NodesPerFace; ++i)
                {
                    add_node_adjacency(face_nodes[i], node);
                    add_node_adjacency(node, face_nodes[i]);
                }

                for (size_t inode: el.internal_nodes)
                {
                    add_node_adjacency(node, inode);
                    add_node_adjacency(inode, node);
                }
            }

            for (size_t node: el.internal_nodes)
            {
                for (size_t i = 0; i < NodesPerFace; ++i)
                {
                    add_node_adjacency(face_nodes[i], node);
                    add_node_adjacency(node, face_nodes[i]);
                }
            }

            for (size_t other_el_index: el.adjacent_elements)
            {
                auto &other_el = m_elems[other_el_index];
                auto it = std::find(other_el.faces.begin(), other_el.faces.end(), el.faces[facei]);
                if (it != other_el.faces.end())
                {
                    auto which = it - other_el.faces.begin();
                    for (size_t i = 0; i < NodesPerFace; ++i)
                    {
                        other_el.face_nodes[which][i] = face_nodes[NodesPerFace-i-1];
                    }
                }
            }
        }
    }

    std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>>
    assign_face_numbers()
    {
        size_t face_number = 0;

        std::vector<smv::SmallVector<NodeFaceInfo, MaxElementAdjacencies+2>>
            node_faces(m_nodes.size());

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
                try
                {
                    node_faces[el.control_nodes[0]].emplace_back(
                        face_number, el.control_nodes[1], eli
                    );
                    face_number += 1;
                }
                catch (const smv::MaxSizeExceeded&)
                {
                    fprintf(stderr, "Exceeded number of faces adjoining a single node at node %ld\n", el.control_nodes[0]);
                    throw smv::MaxSizeExceeded{};
                }
            }

            if (!face_numbered[1])
            {
                el.faces[1] = face_number;
                try
                {
                    node_faces[el.control_nodes[1]].emplace_back(
                        face_number, el.control_nodes[2], eli
                    );
                    face_number += 1;
                }
                catch (const smv::MaxSizeExceeded&)
                {
                    fprintf(stderr, "Exceeded number of faces adjoining a single node at node %ld\n", el.control_nodes[1]);
                    throw smv::MaxSizeExceeded{};
                }
            }

            if (!face_numbered[2])
            {
                el.faces[2] = face_number;
                try
                {
                    node_faces[el.control_nodes[2]].emplace_back(
                        face_number, el.control_nodes[0], eli
                    );
                    face_number += 1;
                }
                catch (const smv::MaxSizeExceeded&)
                {
                    fprintf(stderr, "Exceeded number of faces adjoining a single node at node %ld\n", el.control_nodes[2]);
                    throw smv::MaxSizeExceeded{};
                }
            }
        }

        return node_faces;
    }

    template <class Iterator>
    typename std::enable_if<!is_reversed_iterator<Iterator>::value,
                            BoundaryRepresentation<NodesPerFace>>::type
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
                            BoundaryRepresentation<NodesPerFace>>::type
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
            boundary.nodes.push_back(n1);
        }

        if (NodesPerFace != 0)
        {
            const auto &face_nodes = get_face_nodes(face_info);
            for (size_t i = 0; i < NodesPerFace; ++i)
            {
                boundary.nodes.push_back(face_nodes[i]);
            }
        }

        boundary.nodes.push_back(n2);
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

    TetMesh<int, 1, 3> mesh(nodes, tets, boundaries);

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
    REQUIRE(bound.nodes[0] == 1);
    REQUIRE(bound.nodes[1] == 2);
    REQUIRE(bound.nodes[2] == 3);

    REQUIRE(bound.faces.size() == 2);
    REQUIRE(bound.faces[0].number == 3);
    REQUIRE(bound.faces[0].element == 1);
    REQUIRE(bound.faces[0].nodes == std::array<size_t, 2>{ 0, 1 });
    REQUIRE(bound.faces[1].number == 4);
    REQUIRE(bound.faces[1].element == 1);
    REQUIRE(bound.faces[1].nodes == std::array<size_t, 2>{ 1, 2 });

    mesh.renumber_nodes();

    REQUIRE(mesh.element(0).control_nodes == std::array<size_t, 3>{0, 1, 2});
    REQUIRE(mesh.element(1).control_nodes == std::array<size_t, 3>{2, 1, 3});
    REQUIRE(mesh.coord(0) == std::array<int, 2>{-1, -1});
    REQUIRE(mesh.coord(1) == std::array<int, 2>{-1, 1});
    REQUIRE(mesh.coord(2) == std::array<int, 2>{1, -1});
    REQUIRE(mesh.coord(3) == std::array<int, 2>{1, 1});

    nodeadj = mesh.adjacent_nodes(0);
    REQUIRE(nodeadj.size() == 2);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 2);

    nodeadj = mesh.adjacent_nodes(1);
    REQUIRE(nodeadj.size() == 3);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 2);
    REQUIRE(nodeadj[2] == 3);

    nodeadj = mesh.adjacent_nodes(2);
    REQUIRE(nodeadj.size() == 3);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    REQUIRE(nodeadj[1] == 1);
    REQUIRE(nodeadj[2] == 3);

    nodeadj = mesh.adjacent_nodes(3);
    REQUIRE(nodeadj.size() == 2);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 2);

    REQUIRE(bound.nodes.size() == 3);
    REQUIRE(bound.nodes[0] == 1);
    REQUIRE(bound.nodes[1] == 3);
    REQUIRE(bound.nodes[2] == 2);

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

TEST_CASE("Second order mesh")
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

    TetMesh<int, 1, 8, 1> mesh(nodes, tets, boundaries);
    mesh.renumber_nodes();
    
    auto el = mesh.element(0);
    REQUIRE(el.control_nodes == std::array<size_t, 3>{ 0, 1, 2 });
    REQUIRE(el.face_nodes == std::array<std::array<size_t, 1>, 3>{ 3, 4, 5 });
    REQUIRE(el.adjacent_elements[0] == 1);
    
    el = mesh.element(1);
    REQUIRE(el.control_nodes == std::array<size_t, 3>{ 2, 1, 6 });
    REQUIRE(el.face_nodes == std::array<std::array<size_t, 1>, 3>{ 4, 7, 8 });
    REQUIRE(el.adjacent_elements[0] == 0);
    
    const auto &bound = mesh.boundary(0);
    REQUIRE(bound.nodes.size() == 5);
    REQUIRE(bound.nodes[0] == 1);
    REQUIRE(bound.nodes[1] == 7);
    REQUIRE(bound.nodes[2] == 6);
    REQUIRE(bound.nodes[3] == 8);
    REQUIRE(bound.nodes[4] == 2);
    
    auto nodeadj = mesh.adjacent_nodes(4);
    REQUIRE(nodeadj.size() == 8);
    std::sort(nodeadj.begin(), nodeadj.end());
    for (size_t i = 0; i < 8; ++i)
    {
        REQUIRE(nodeadj[i] == (i + (i >= 4)));
    }
    
    nodeadj = mesh.adjacent_nodes(5);
    REQUIRE(nodeadj.size() == 5);
    std::sort(nodeadj.begin(), nodeadj.end());
    for (size_t i = 0; i < 5; ++i)
    {
        REQUIRE(nodeadj[i] == i);
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

    TetMesh<double, 1, 15, 2, 1> mesh(nodes, tets, boundaries);

    REQUIRE(mesh.average_bandwidth() == doctest::Approx(164.0 / 16));

    auto eladj = mesh.element(0).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 1);

    eladj = mesh.element(1).adjacent_elements;
    REQUIRE(eladj.size() == 1);
    REQUIRE(eladj[0] == 0);

    // Inexhaustive checks on coordinates.
    REQUIRE(mesh.coord(3) == std::array<double, 2>{1, -1});
    REQUIRE(mesh.coord(13)[0] == doctest::Approx(1.0));
    REQUIRE(mesh.coord(13)[1] == doctest::Approx(1.0 / 3));
    REQUIRE(mesh.coord(9)[0] == doctest::Approx(-1.0/3));
    REQUIRE(mesh.coord(9)[1] == doctest::Approx(-1.0));
    REQUIRE(mesh.coord(5)[0] == doctest::Approx(-1.0));
    REQUIRE(mesh.coord(5)[1] == doctest::Approx(1.0 / 3));
    REQUIRE(mesh.coord(11)[0] == doctest::Approx(-1.0 / 3));
    REQUIRE(mesh.coord(11)[1] == doctest::Approx(1.0));

    auto nodeadj = mesh.adjacent_nodes(1);
    REQUIRE(nodeadj.size() == 15);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 0);
    for (size_t i = 1; i < 15; ++i)
    {
        CHECK(nodeadj[i] == i+1);
    }

    nodeadj = mesh.adjacent_nodes(6);
    REQUIRE(nodeadj.size() == 15);

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
    REQUIRE(bound.nodes[0] == 1);
    REQUIRE(bound.nodes[1] == 11);
    REQUIRE(bound.nodes[2] == 12);
    REQUIRE(bound.nodes[3] == 2);
    REQUIRE(bound.nodes[4] == 13);
    REQUIRE(bound.nodes[5] == 14);
    REQUIRE(bound.nodes[6] == 3);

    REQUIRE(bound.faces.size() == 2);
    REQUIRE(bound.faces[0].number == 3);
    REQUIRE(bound.faces[0].element == 1);
    REQUIRE(bound.faces[0].nodes == std::array<size_t, 4>{ 0, 1, 2, 3 });
    REQUIRE(bound.faces[1].number == 4);
    REQUIRE(bound.faces[1].element == 1);
    REQUIRE(bound.faces[1].nodes == std::array<size_t, 4>{ 3, 4, 5, 6 });

    mesh.renumber_nodes();

    // Inexhaustive checks on the node renumbering.
    auto el = mesh.element(0);
    REQUIRE(el.control_nodes == std::array<size_t, 3>{0, 1, 2});
    REQUIRE(el.face_nodes[0] == std::array<size_t, 2>{3, 4});
    REQUIRE(el.face_nodes[1] == std::array<size_t, 2>{5, 6});
    REQUIRE(el.face_nodes[2] == std::array<size_t, 2>{7, 8});
    REQUIRE(el.internal_nodes == std::array<size_t, 1>{9});
    
    el = mesh.element(1);
    REQUIRE(el.control_nodes == std::array<size_t, 3>{2, 1, 10});
    REQUIRE(el.face_nodes[0] == std::array<size_t, 2>{6, 5});
    REQUIRE(el.face_nodes[1] == std::array<size_t, 2>{11, 12});
    REQUIRE(el.face_nodes[2] == std::array<size_t, 2>{13, 14});
    REQUIRE(el.internal_nodes == std::array<size_t, 1>{15});

    nodeadj = mesh.adjacent_nodes(0);
    std::sort(nodeadj.begin(), nodeadj.end());
    for (size_t i = 0; i < 9; ++i)
    {
        REQUIRE(nodeadj[i] == i+1);
    }
    
    nodeadj = mesh.adjacent_nodes(5);
    REQUIRE(nodeadj.size() == 15);
    std::sort(nodeadj.begin(), nodeadj.end());
    for (size_t i = 0; i < 5; ++i)
    {
        REQUIRE(nodeadj[i] == i);
    }
    for (size_t i = 5; i < 15; ++i)
    {
        REQUIRE(nodeadj[i] == i+1);
    }

    REQUIRE(bound.nodes.size() == 7);
    REQUIRE(bound.nodes[0] == 1);
    REQUIRE(bound.nodes[1] == 11);
    REQUIRE(bound.nodes[2] == 12);
    REQUIRE(bound.nodes[3] == 10);
    REQUIRE(bound.nodes[4] == 13);
    REQUIRE(bound.nodes[5] == 14);
    REQUIRE(bound.nodes[6] == 2);

    REQUIRE(mesh.coord(13)[0] == doctest::Approx(1.0));
    REQUIRE(mesh.coord(13)[1] == doctest::Approx(1.0 / 3));
    REQUIRE(mesh.coord(6)[0] == doctest::Approx(1.0 / 3));
    REQUIRE(mesh.coord(6)[1] == doctest::Approx(-1.0 / 3));
    REQUIRE(mesh.coord(4)[0] == doctest::Approx(-1.0));
    REQUIRE(mesh.coord(4)[1] == doctest::Approx(1.0 / 3));
} // TEST_CASE

TEST_CASE("Fourth order mesh")
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

    TetMesh<int, 1, 24, 3, 3> mesh(nodes, tets, boundaries);
    mesh.renumber_nodes();
    
    auto el = mesh.element(0);
    REQUIRE(el.control_nodes == std::array<size_t, 3>{ 0, 1, 2 });
    REQUIRE(el.face_nodes == std::array<std::array<size_t, 3>, 3>{
        3, 4, 5, 6, 7, 8, 9, 10, 11});
    REQUIRE(el.internal_nodes == std::array<size_t, 3>{12, 13, 14});
    
    el = mesh.element(1);
    REQUIRE(el.control_nodes == std::array<size_t, 3>{2, 1, 15});
    REQUIRE(el.face_nodes == std::array<std::array<size_t, 3>, 3>{
        8, 7, 6, 16, 17, 18, 19, 20, 21});
    REQUIRE(el.internal_nodes == std::array<size_t, 3>{22, 23, 24});
    
    auto nodeadj = mesh.adjacent_nodes(8);
    REQUIRE(nodeadj.size() == 24);
    std::sort(nodeadj.begin(), nodeadj.end());
    for (size_t i = 0; i < 24; ++i)
    {
        REQUIRE(nodeadj[i] == (i + (i >= 8)));
    }
    
    nodeadj = mesh.adjacent_nodes(9);
    REQUIRE(nodeadj.size() == 14);
    std::sort(nodeadj.begin(), nodeadj.end());
    for (size_t i = 0; i < 14; ++i)
    {
        REQUIRE(nodeadj[i] == i + (i >= 9));
    }
    
    nodeadj = mesh.adjacent_nodes(15);
    REQUIRE(nodeadj.size() == 14);
    std::sort(nodeadj.begin(), nodeadj.end());
    REQUIRE(nodeadj[0] == 1);
    REQUIRE(nodeadj[1] == 2);
    REQUIRE(nodeadj[2] == 6);
    REQUIRE(nodeadj[3] == 7);
    REQUIRE(nodeadj[4] == 8);
    for (size_t i = 5; i < 13; ++i)
    {
        REQUIRE(nodeadj[i] == i + 11);
    }
}

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
    REQUIRE(boundary.nodes[0] == 2);
    REQUIRE(boundary.nodes[1] == 26);
    REQUIRE(boundary.nodes[2] == 27);
    REQUIRE(boundary.nodes[3] == 6);
    REQUIRE(boundary.nodes[4] == 52);
    REQUIRE(boundary.nodes[5] == 53);
    REQUIRE(boundary.nodes[6] == 1);
    REQUIRE(boundary.nodes[7] == 57);
    REQUIRE(boundary.nodes[8] == 58);
    REQUIRE(boundary.nodes[9] == 7);
    REQUIRE(boundary.nodes[10] == 33);
    REQUIRE(boundary.nodes[11] == 34);
    REQUIRE(boundary.nodes[12] == 0);

    // Selective tests on the coordinates of nodes.
    /*
    REQUIRE(boundary.nodes[10].coords[0] == doctest::Approx(-1.0));
    REQUIRE(boundary.nodes[10].coords[1] == doctest::Approx(1.0 / 3));
    REQUIRE(boundary.nodes[4].coords[0] == doctest::Approx(-1.0 / 3));
    REQUIRE(boundary.nodes[5].coords[0] == doctest::Approx(-2.0 / 3)); */

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

struct MeshConversionException : public std::exception
{
    const char* what() const noexcept
    { 
        return "Failed conversion from gmsh to tetmesh";
    }

    MeshConversionException() noexcept = default;
};

template <size_t MaxElementAdjacencies, size_t MaxNodeAdjacencies,
          size_t NodesPerFace = 0, size_t InternalNodes = 0>
TetMesh<double, MaxElementAdjacencies, MaxNodeAdjacencies,
        NodesPerFace, InternalNodes>
parse_gmsh_to_tetmesh(const char *name)
{
    const auto mesh_data = gmsh::parse_gmsh_file(name);

    // For now I assume that only one surface is defined; check that.
    auto num_surfaces = mesh_data.entities.surfs.size();
    if (num_surfaces != 1)
    {
        throw MeshConversionException();
    }

    std::vector<std::array<double, 2>> nodes;
    nodes.reserve(mesh_data.nodes.size());

    for (const auto &node: mesh_data.nodes)
    {
        nodes.emplace_back(std::array<double, 2>{ node.coords[0], node.coords[1] });
    }

    std::vector<std::array<size_t, 3>> tets;
    for (const auto &element: mesh_data.elements)
    {
        if (element.entity_type == gmsh::EntityType::Surface)
        {
            tets.emplace_back(std::array<size_t, 3>{element.node_tags[0],
                                                    element.node_tags[1],
                                                    element.node_tags[2]});
        }
    }

    std::vector<std::vector<size_t>> boundaries;
    // I assume all curves are boundaries of the domain.
    // This limitation should be removed eventually.
    for (size_t i = 0; i < mesh_data.entities.curves.size(); ++i)
    {
        boundaries.emplace_back();
        auto &boundary = boundaries.back();
        // Get iterator to first line element on this curve.
        auto it = std::find_if(mesh_data.elements.begin(),
            mesh_data.elements.end(),
            [=](const gmsh::ElementData &el)
            {
                return el.entity_type == gmsh::EntityType::Curve &&
                       el.entity_tag  == i;
            });
        if (it == mesh_data.elements.end())
        {
            throw MeshConversionException();
        }
        boundary.push_back(it->node_tags[0]);
        boundary.push_back(it->node_tags[1]);
        it += 1;

        while (it != mesh_data.elements.end() && 
               it->entity_type == gmsh::EntityType::Curve &&
               it->entity_tag  == i)
        {
            if (it->node_tags[0] != boundary.back())
            {
                throw MeshConversionException();
            }
            boundary.push_back(it->node_tags[1]);
            it += 1;
        }
    }

    return TetMesh<double, MaxElementAdjacencies, MaxNodeAdjacencies,
                   NodesPerFace, InternalNodes>(
        nodes, tets, boundaries
    );
}

} // namespace msh

#endif // TETMESH_HPP
