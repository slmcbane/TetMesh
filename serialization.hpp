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

#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <cstddef>
#include <cstdio>
#include <cstdint>

/*
 * This file just separates out serialization code from the main TetMesh header.
 * It won't work included on its own; it is included when you include TetMesh.hpp.
 */

struct SerializationException : public std::exception
{};

namespace ser
{

void checked_write(const void *data, size_t size, size_t count, FILE *out)
{
    size_t written = fwrite(data, size, count, out);
    if (written != count)
    {
        throw SerializationException{};
    }
}

template <class T>
T checked_read(FILE *in)
{
    T x;
    if (fread(&x, sizeof(T), 1, in) != 1)
    {
        throw SerializationException{};
    }
    return x;
}

template <class T>
void checked_read(FILE *in, T *dst, const int64_t count)
{
    if (fread(dst, sizeof(T), count, in) != count)
    {
        throw SerializationException{};
    }
}

template <class T, class T2 = T, size_t MaxNodeAdjacencies>
void serialize(const smv::SmallVector<T, MaxNodeAdjacencies> &vec, FILE *out)
{
    if constexpr (std::is_same_v<T, T2>)
    {
        int64_t count = vec.size();
        checked_write(&count, sizeof(int64_t), 1, out);
        checked_write(vec.data(), sizeof(T), count, out);
    }
    else
    {
        smv::SmallVector<T2, MaxNodeAdjacencies> v2;
        v2.resize(vec.size());
        std::copy(vec.begin(), vec.end(), v2.begin());
        serialize(v2, out);
    }
}

template <class T, class T2 = T, size_t N>
void deserialize(FILE *in, smv::SmallVector<T, N> &vec)
{
    if constexpr (std::is_same_v<T, T2>)
    {
        vec.resize(checked_read<int64_t>(in));
        checked_read(in, vec.data(), vec.size());
    }
    else
    {
        smv::SmallVector<T2, N> v2;
        deserialize(in, v2);
        vec.resize(v2.size());
        std::copy(v2.begin(), v2.end(), vec.begin());
    }
}

template <class T, class T2 = T, size_t N>
void serialize(const std::array<T, N> &x, FILE *out)
{
    if constexpr(std::is_same_v<T, T2>)
    {
        checked_write(x.data(), sizeof(T), N, out);
    }
    else
    {
        std::array<T2, N> y;
        std::copy(x.begin(), x.end(), y.begin());
        serialize(y, out);
    }
}

template <class T, class T2 = T, size_t N>
void deserialize(FILE *in, std::array<T, N> &x)
{
    if constexpr (std::is_same_v<T, T2>)
    {
        checked_read(in, x.data(), N);
    }
    else
    {
        std::array<T2, N> y;
        deserialize(in, y);
        std::copy(y.begin(), y.end(), x.begin());
    }
}

template <class CoordT, size_t MaxNodeAdjacencies>
void serialize(const NodeInfo<CoordT, MaxNodeAdjacencies> &node,
               FILE *out)
{
    serialize(node.coords, out);
    serialize<size_t, int64_t>(node.adjacent_nodes, out);
}

template <class CoordT, size_t MaxNodeAdjacencies>
void deserialize(FILE *in, NodeInfo<CoordT, MaxNodeAdjacencies> &node)
{
    deserialize(in, node.coords);
    deserialize<size_t, int64_t>(in, node.adjacent_nodes);
}

template <size_t MaxElementAdjacencies, size_t NodesPerFace, size_t InternalNodes>
void serialize(const ElementInfo<MaxElementAdjacencies, NodesPerFace, InternalNodes> &el,
               FILE *out)
{
    serialize<size_t, int64_t>(el.control_nodes, out);
    serialize<size_t, int64_t>(el.faces, out);
    for (auto face : el.face_nodes)
    {
        serialize<size_t, int64_t>(face, out);
    }
    serialize<size_t, int64_t>(el.internal_nodes, out);
    serialize<size_t, int64_t>(el.adjacent_elements, out);
}

template <size_t... Params>
void deserialize(FILE *in, ElementInfo<Params...> &el)
{
    deserialize<size_t, int64_t>(in, el.control_nodes);
    deserialize<size_t, int64_t>(in, el.faces);
    for (auto &face : el.face_nodes)
    {
        deserialize<size_t, int64_t>(in, face);
    }
    deserialize<size_t, int64_t>(in, el.internal_nodes);
    deserialize<size_t, int64_t>(in, el.adjacent_elements);
}

template <size_t N>
void serialize(const typename BoundaryRepresentation<N>::FaceDetails &face,
               FILE *out)
{
    std::array<int64_t, 4+N> numbers;
    numbers[0] = face.number;
    numbers[1] = face.element;
    std::copy(face.nodes.begin(), face.nodes.end(), numbers.data() + 2);
    serialize(numbers, out);
}

template <size_t N>
void deserialize(FILE *in, typename BoundaryRepresentation<N>::FaceDetails &face)
{
    std::array<int64_t, 4+N> numbers;
    deserialize(in, numbers);
    face.number = numbers[0];
    face.element = numbers[1];
    std::copy(numbers.begin() + 2, numbers.end(), face.nodes.begin());
}

template <size_t NodesPerFace>
void serialize(const BoundaryRepresentation<NodesPerFace> &boundary,
               FILE *out)
{
    int64_t count = boundary.nodes.size();
    checked_write(&count, sizeof(count), 1, out);

    {
        std::vector<int64_t> nodes(boundary.nodes.begin(), boundary.nodes.end());
        checked_write(nodes.data(), sizeof(int64_t), count, out);
    }

    count = boundary.faces.size();
    checked_write(&count, sizeof(count), 1, out);
    for (const auto &face: boundary.faces)
    {
        serialize<NodesPerFace>(face, out);
    }
}

template <size_t NodesPerFace>
void deserialize(FILE *in, BoundaryRepresentation<NodesPerFace> &boundary)
{
    int64_t count = checked_read<int64_t>(in);
    {
        std::vector<int64_t> nodes(static_cast<typename std::vector<int64_t>::size_type>(count));
        checked_read(in, nodes.data(), count);
        boundary.nodes = std::vector<size_t>(nodes.begin(), nodes.end());
    }

    count = checked_read<int64_t>(in);
    boundary.faces.resize(count);
    for (auto &face: boundary.faces)
    {
        deserialize<NodesPerFace>(in, face);
    }
}

void serialize(const std::string &str, FILE *out)
{
    int64_t count = str.size();
    checked_write(&count, sizeof(count), 1, out);
    checked_write(str.data(), sizeof(char), count, out);
}

void deserialize(FILE *in, std::string &str)
{
    int64_t count = checked_read<int64_t>(in);
    str.resize(count);
    checked_read(in, str.data(), count);
}

void serialize(const std::vector<std::string> &tags, FILE *out)
{
    int64_t count = tags.size();
    checked_write(&count, sizeof(count), 1, out);
    for (const auto &str: tags)
    {
        serialize(str, out);
    }
}

void deserialize(FILE *in, std::vector<std::string> &tags)
{
    int64_t count = checked_read<int64_t>(in);
    tags.resize(count);
    for (auto &str: tags)
    {
        deserialize(in, str);
    }
}

template <class T>
void elementwise_serialize(const std::vector<T> &vec, FILE *out)
{
    int64_t count = vec.size();
    checked_write(&count, sizeof(count), 1, out);
    for (const auto &el: vec)
    {
        serialize(el, out);
    }
}

template <class T>
auto elementwise_deserialize(FILE *in)
{
    std::vector<T> vec(
        static_cast<typename std::vector<T>::size_type>(
            checked_read<int64_t>(in)
        )
    );
    for (auto &el: vec)
    {
        deserialize(in, el);
    }
    return vec;
}

} // namespace ser

#endif // SERIALIZATION_HPP
