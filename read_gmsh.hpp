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
static const ParsingException past_the_end("Requested access to data past the end of what was read");

enum class SectionType
{
    MeshFormat, PhysicalNames, Entities, PartitionedEntities,
    Nodes, Elements, Periodic, GhostElements, Parametrizations,
    NodeData, ElementData, ElementNodeData, InterpolationScheme
};

class ParserState
{
    size_t m_offset;
    const char *m_data;
    size_t m_size;
    size_t m_int_size;

    static_assert(sizeof(double) == 8, "Assumes 64-bit double");

public:
    ParserState(const char *ptr, size_t size) noexcept :
        m_offset{0}, m_data{ptr}, m_size{size}, m_int_size{0}
    {}

    size_t offset() const noexcept { return m_offset; }

    void set_data_size(size_t sz) noexcept { m_int_size = sz; }
    size_t int_size() const noexcept { return m_int_size; }

    void add_offset(size_t i)
    {
        if ((m_offset + i) > m_size)
        {
            throw past_the_end;
        }
        m_offset += i;
    }

    const char *data() const
    {
        if (m_offset >= m_size)
        {
            throw past_the_end;
        }
        return m_data + m_offset;
    }

    char current() const
    {
        return *data();
    }

    size_t size() const noexcept { return m_size; }

    bool compare_next(const char *str, size_t len) const noexcept
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

    double extract_ascii_double()
    {
        char *end;
        double val = std::strtod(data(), &end);
        add_offset(end - data());
        return val;
    }

    long extract_ascii_int()
    {
        char *end;
        long val = std::strtol(data(), &end, 10);
        add_offset(end - data());
        m_offset = end - m_data;
        return val;
    }

    size_t extract_int()
    {
        switch(int_size())
        {
            case 4:
                if (m_size - m_offset < 4)
                {
                    throw past_the_end;
                }
                static_assert(sizeof(uint32_t) == 4, "If 4 bytes isn't 32 bits I can't even");
                uint32_t i;
                std::memcpy(&i, data(), 4);
                add_offset(4);
                return static_cast<size_t>(i);
            case 8:
                if (m_size - m_offset < 8)
                {
                    throw past_the_end;
                }
                static_assert(sizeof(uint64_t) == 8, "If 8 bytes isn't 64 bits I can't even");
                uint64_t j;
                std::memcpy(&j, data(), 8);
                add_offset(8);
                return static_cast<size_t>(j);
            default:
                throw ParsingException("Unimplemented int_size: expected 4 or 8 bytes");
        }
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
/********************************************************************************
 *******************************************************************************/

// Doesn't return anything; I'm not planning to figure out exactly which msh
// versions I support, and the only outcome needed is the setting of data_size
// in the ParserState.
inline void
parse_mesh_format(ParserState &state)
{
    auto header = parse_section_header(state);
    if (header != SectionType::MeshFormat)
    {
        throw ParsingException("Did not read mesh format as first section in file");
    }

    state.extract_ascii_double();
    long is_binary = state.extract_ascii_int();
    
    if (is_binary == 0)
    {
        throw ParsingException("Only binary mode msh files accepted");
    }

    long data_size = state.extract_ascii_int();

    // There is a newline after the data size in binary files.
    if (!(state.current() == '\n'))
    {
        throw ParsingException("Expected \\n after data size in format");
    }
    state.add_offset(1);

    // It looks like gmsh uses a 32-bit int for the endianness check, not
    // what data_size says.
    state.set_data_size(4);
    size_t endian_check = state.extract_int();
    if (endian_check != 1)
    {
        throw ParsingException("Conversion of endianness not implemented");
    }
    state.set_data_size(data_size);
}

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parsing the mesh format")
{
    // This is an example from a real binary .msh file
    const char *valid_format = "$MeshFormat\n4.1 1 8\n\x01\0\0\0\n$EndMeshFormat";
    ParserState state(valid_format, 40);
    REQUIRE_NOTHROW(parse_mesh_format(state));
    REQUIRE(state.int_size() == 8);

    const char *test_format = "$MeshFormat\n4.1 0 8\n$EndMeshFormat";
    state = ParserState(test_format, 34);
    REQUIRE_THROWS_WITH(parse_mesh_format(state), "Only binary mode msh files accepted");
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace gmsh

} // namespace msh

#endif // READ_GMSH_HPP