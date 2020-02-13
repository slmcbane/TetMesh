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

#ifndef READ_GMSH_HPP
#define READ_GMSH_HPP

#include <cstdio>
#include <cstring>
#include <exception>
#include <utility>

namespace msh
{

namespace gmsh
{

struct ParsingException : public std::exception
{
    const char *msg;
    const char *what() const noexcept { return msg; }

    ParsingException(const char* m) noexcept : msg{m}
    {}
};

static const ParsingException unrecognized_header("Unrecognized section header");

enum class SectionType
{
    MeshFormat, PhysicalNames, Entities, PartitionedEntities,
    Nodes, Elements, Periodic, GhostElements, Parametrizations,
    NodeData, ElementData, ElementNodeData, InterpolationScheme
};

class ParserState
{
    size_t m_offset;
    const unsigned char *m_data;
    size_t m_size;
    size_t int_size;

    static_assert(sizeof(double) == 8, "Assumes 64-bit double");

public:
    ParserState(const unsigned char *ptr, size_t size) noexcept :
        m_offset{0}, m_data{ptr}, m_size{size}, int_size{0}
    {}

    ParserState(const char *ptr, size_t size) noexcept :
        ParserState(reinterpret_cast<const unsigned char *>(ptr), size)
    {}

    size_t offset() const noexcept { return m_offset; }

    void add_offset(size_t i)
    {
        if ((m_offset + i) > m_size)
        {
            throw ParsingException("Requested offset increment past end of read data");
        }
        m_offset += i;
    }

    const unsigned char *data() const
    {
        if (m_offset >= m_size)
        {
            throw ParsingException("Requested access to data past the end of what was read");
        }
        return m_data + m_offset;
    }

    unsigned char current() const
    {
        return *data();
    }

    size_t size() const noexcept { return m_size; }

    bool compare_next(const unsigned char *str, size_t len) const noexcept
    {
        if (len > (m_size - m_offset))
        {
            return false;
        }
        else
        {
            return std::memcmp(data(), str, len) == 0;
        }
    }

    bool compare_next(const char *str, size_t len) const noexcept
    {
        return compare_next(reinterpret_cast<const unsigned char*>(str), len);
    }
};

inline SectionType
get_mesh_format_header(ParserState &state)
{
    if (state.compare_next("MeshFormat", 10))
    {
        state.add_offset(10);
        return SectionType::MeshFormat;
    }

    throw unrecognized_header;
}

inline SectionType
get_ghost_elements_header(ParserState &state)
{
    if (state.compare_next("GhostElements", 13))
    {
        state.add_offset(13);
        return SectionType::GhostElements;
    }

    throw unrecognized_header;
}

inline SectionType
get_interpolation_scheme_header(ParserState &state)
{
    if (state.compare_next("InterpolationScheme", 19))
    {
        state.add_offset(19);
        return SectionType::InterpolationScheme;
    }

    throw unrecognized_header;
}

inline SectionType
get_letter_n_header(ParserState &state)
{
    if (!state.compare_next("Node", 4))
    {
        throw unrecognized_header;
    }

    state.add_offset(4);
    switch(state.current())
    {
        case 's':
            return SectionType::Nodes;
        case 'D':
            if (state.compare_next("Data", 4))
            {
                state.add_offset(4);
                return SectionType::NodeData;
            }
        default:
            throw unrecognized_header;
    }
}

inline SectionType
get_letter_p_header(ParserState &state)
{
    state.add_offset(1);
    switch(state.current())
    {
        case 'a':
            if (state.compare_next("artitionedEntities", 18))
            {
                state.add_offset(18);
                return SectionType::PartitionedEntities;
            }
            else if (state.compare_next("arametrizations", 15))
            {
                state.add_offset(15);
                return SectionType::Parametrizations;
            }
        case 'e':
            if (state.compare_next("eriodic", 7))
            {
                state.add_offset(7);
                return SectionType::Periodic;
            }
        case 'h':
            if (state.compare_next("hysicalNames", 12))
            {
                state.add_offset(12);
                return SectionType::PhysicalNames;
            }
        default:
            throw unrecognized_header;
    }
}

inline SectionType
get_letter_e_header(ParserState &state)
{
    state.add_offset(1);
    switch(state.current())
    {
        case 'l':
            if (!state.compare_next("lement", 6))
            {
                throw unrecognized_header;
            }
            state.add_offset(6);
            switch(state.current())
            {
                case 's':
                    return SectionType::Elements;
                case 'D':
                    if (state.compare_next("Data", 4))
                    {
                        state.add_offset(4);
                        return SectionType::ElementData;
                    }
                case 'N':
                    if (state.compare_next("NodeData", 8))
                    {
                        state.add_offset(8);
                        return SectionType::ElementNodeData;
                    }
                default:
                    throw unrecognized_header;
            }
        case 'n':
            if (state.compare_next("ntities", 7))
            {
                state.add_offset(7);
                return SectionType::Entities;
            }
        default:
            throw unrecognized_header;
    }
}

inline SectionType
parse_section_header(ParserState &state)
{
    if (!(state.current() == '$'))
    {
        throw ParsingException("Expected to read a section header");
    }

    state.add_offset(1);

    switch(state.current())
    {
        case 'M':
            return get_mesh_format_header(state);
        case 'G':
            return get_ghost_elements_header(state);
        case 'I':
            return get_interpolation_scheme_header(state);
        case 'P':
            return get_letter_p_header(state);
        case 'E':
            return get_letter_e_header(state);
        case 'N':
            return get_letter_n_header(state);
        default:
            throw unrecognized_header;
    }
}

/********************************************************************************
 * Test parsing of section headers.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[GMSH] Parse section headers")
{
    ParserState state("$Nodes", 6);
    REQUIRE(parse_section_header(state) == SectionType::Nodes);

    state = ParserState("$MeshFormat", 11);
    REQUIRE(parse_section_header(state) == SectionType::MeshFormat);

    state = ParserState("$PhysicalNames", 14);
    REQUIRE(parse_section_header(state) == SectionType::PhysicalNames);

    state = ParserState("$Entities", 9);
    REQUIRE(parse_section_header(state) == SectionType::Entities);

    state = ParserState("$PartitionedEntities", 20);
    REQUIRE(parse_section_header(state) == SectionType::PartitionedEntities);

    state = ParserState("$Elements", 9);
    REQUIRE(parse_section_header(state) == SectionType::Elements);

    state = ParserState("$Periodic", 9);
    REQUIRE(parse_section_header(state) == SectionType::Periodic);

    state = ParserState("$GhostElements", 14);
    REQUIRE(parse_section_header(state) == SectionType::GhostElements);

    state = ParserState("$Parametrizations", 17);
    REQUIRE(parse_section_header(state) == SectionType::Parametrizations);

    state = ParserState("$NodeData", 9);
    REQUIRE(parse_section_header(state) == SectionType::NodeData);

    state = ParserState("$ElementData", 12);
    REQUIRE(parse_section_header(state) == SectionType::ElementData);

    state = ParserState("$ElementNodeData", 16);
    REQUIRE(parse_section_header(state) == SectionType::ElementNodeData);

    state = ParserState("$InterpolationScheme", 20);
    REQUIRE(parse_section_header(state) == SectionType::InterpolationScheme);

    REQUIRE_THROWS(parse_section_header(state));
    state = ParserState("$Elementz", 9);
    REQUIRE_THROWS_WITH(parse_section_header(state), unrecognized_header.what());
    state = ParserState("blah", 4);
    REQUIRE_THROWS_WITH(parse_section_header(state),
                        "Expected to read a section header");
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace gmsh

} // namespace msh

#endif // READ_GMSH_HPP