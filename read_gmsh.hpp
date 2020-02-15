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

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <utility>
#include <vector>

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

static const ParsingException unrecognized_section_name("Unrecognized section name");
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
    size_t data_size() const noexcept { return m_int_size; }

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

    int32_t extract_int()
    {
        if (m_size - m_offset < 4)
        {
            throw past_the_end;
        }
        int32_t i;
        static_assert(sizeof(int32_t) == 4, "where is 4 bytes not 32 bits?");
        std::memcpy(&i, data(), 4);
        add_offset(4);
        return i;
    }
    
    size_t extract_size_t()
    {
        switch(data_size())
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
    
    void extract_doubles(double *dst, size_t n)
    {
        if (n * 8 > m_size - m_offset)
        {
            throw past_the_end;
        }
        memcpy(dst, data(), n*8);
        add_offset(n * 8);
    }

    void skip_whitespace()
    {
        while (offset() < size() && isspace(current()))
        {
            add_offset(1);
        }
    }
};

inline SectionType
get_mesh_format_name(ParserState &state)
{
    if (state.compare_next("MeshFormat", 10))
    {
        state.add_offset(10);
        return SectionType::MeshFormat;
    }

    throw unrecognized_section_name;
}

inline SectionType
get_ghost_elements_name(ParserState &state)
{
    if (state.compare_next("GhostElements", 13))
    {
        state.add_offset(13);
        return SectionType::GhostElements;
    }

    throw unrecognized_section_name;
}

inline SectionType
get_interpolation_scheme_name(ParserState &state)
{
    if (state.compare_next("InterpolationScheme", 19))
    {
        state.add_offset(19);
        return SectionType::InterpolationScheme;
    }

    throw unrecognized_section_name;
}

inline SectionType
get_letter_n_name(ParserState &state)
{
    if (!state.compare_next("Node", 4))
    {
        throw unrecognized_section_name;
    }

    state.add_offset(4);
    switch(state.current())
    {
        case 's':
            state.add_offset(1);
            return SectionType::Nodes;
        case 'D':
            if (state.compare_next("Data", 4))
            {
                state.add_offset(4);
                return SectionType::NodeData;
            }
        default:
            throw unrecognized_section_name;
    }
}

inline SectionType
get_letter_p_name(ParserState &state)
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
            throw unrecognized_section_name;
    }
}

inline SectionType
get_letter_e_name(ParserState &state)
{
    state.add_offset(1);
    switch(state.current())
    {
        case 'l':
            if (!state.compare_next("lement", 6))
            {
                throw unrecognized_section_name;
            }
            state.add_offset(6);
            switch(state.current())
            {
                case 's':
                    state.add_offset(1);
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
                    throw unrecognized_section_name;
            }
        case 'n':
            if (state.compare_next("ntities", 7))
            {
                state.add_offset(7);
                return SectionType::Entities;
            }
        default:
            throw unrecognized_section_name;
    }
}

inline SectionType
parse_section_name(ParserState &state)
{
    SectionType type;
    switch(state.current())
    {
        case 'M':
            type = get_mesh_format_name(state);
            break;
        case 'G':
            type = get_ghost_elements_name(state);
            break;
        case 'I':
            type = get_interpolation_scheme_name(state);
            break;
        case 'P':
            type = get_letter_p_name(state);
            break;
        case 'E':
            type = get_letter_e_name(state);
            break;
        case 'N':
            type = get_letter_n_name(state);
            break;
        default:
            throw unrecognized_section_name;
    }
    
    if (state.current() != '\n')
    {
        throw ParsingException("Expected newline at end of section name");
    }
    state.add_offset(1);
    return type;
}

inline SectionType
parse_section_header(ParserState &state)
{
    if (!(state.current() == '$'))
    {
        throw ParsingException("Expected to read a section header");
    }
    state.add_offset(1);
    return parse_section_name(state);
}

/********************************************************************************
 * Test parsing of section headers.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Parse section headers")
{
    ParserState state("$Nodes\n", 7);
    REQUIRE(parse_section_header(state) == SectionType::Nodes);

    state = ParserState("$MeshFormat\n", 12);
    REQUIRE(parse_section_header(state) == SectionType::MeshFormat);

    state = ParserState("$PhysicalNames\n", 15);
    REQUIRE(parse_section_header(state) == SectionType::PhysicalNames);

    state = ParserState("$Entities\n", 10);
    REQUIRE(parse_section_header(state) == SectionType::Entities);

    state = ParserState("$PartitionedEntities\n", 21);
    REQUIRE(parse_section_header(state) == SectionType::PartitionedEntities);

    state = ParserState("$Elements\n", 10);
    REQUIRE(parse_section_header(state) == SectionType::Elements);

    state = ParserState("$Periodic\n", 10);
    REQUIRE(parse_section_header(state) == SectionType::Periodic);

    state = ParserState("$GhostElements\n", 15);
    REQUIRE(parse_section_header(state) == SectionType::GhostElements);

    state = ParserState("$Parametrizations\n", 18);
    REQUIRE(parse_section_header(state) == SectionType::Parametrizations);

    state = ParserState("$NodeData\n", 10);
    REQUIRE(parse_section_header(state) == SectionType::NodeData);

    state = ParserState("$ElementData\n", 13);
    REQUIRE(parse_section_header(state) == SectionType::ElementData);

    state = ParserState("$ElementNodeData\n", 17);
    REQUIRE(parse_section_header(state) == SectionType::ElementNodeData);

    state = ParserState("$InterpolationScheme\n", 21);
    REQUIRE(parse_section_header(state) == SectionType::InterpolationScheme);

    REQUIRE_THROWS(parse_section_header(state));
    state = ParserState("$Elementz", 9);
    REQUIRE_THROWS_WITH(parse_section_header(state), unrecognized_section_name.what());
    state = ParserState("blah", 4);
    REQUIRE_THROWS_WITH(parse_section_header(state),
                        "Expected to read a section header");
    state = ParserState("$Nodes ", 7);
    REQUIRE_THROWS_WITH(parse_section_header(state),
                        "Expected newline at end of section name");
                       
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

inline void
parse_mesh_format(ParserState &state)
{
    auto header = parse_section_header(state);
    if (header != SectionType::MeshFormat)
    {
        throw ParsingException("Did not read mesh format as first section in file");
    }

    double version = state.extract_ascii_double();
    if (version < 4.1)
    {
        throw ParsingException("Parser implemented for msh format 4.1");
    }
    
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
    size_t endian_check = state.extract_int();
    if (endian_check != 1)
    {
        throw ParsingException("Conversion of endianness not implemented");
    }
    state.set_data_size(data_size);
}

/********************************************************************************
 * Test parsing of mesh format
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parsing the mesh format")
{
    // This is an example from a real binary .msh file
    const char *valid_format = "$MeshFormat\n4.1 1 8\n\x01\0\0\0\n$EndMeshFormat";
    ParserState state(valid_format, 40);
    REQUIRE_NOTHROW(parse_mesh_format(state));
    REQUIRE(state.data_size() == 8);

    const char *test_format = "$MeshFormat\n4.1 0 8\n$EndMeshFormat";
    state = ParserState(test_format, 34);
    REQUIRE_THROWS_WITH(parse_mesh_format(state), "Only binary mode msh files accepted");
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

inline void
parse_section_end(ParserState &state, SectionType what_section)
{
    state.skip_whitespace();
    if (!state.compare_next("$End", 4))
    {
        throw ParsingException("Expected end of section");
    }
    state.add_offset(4);
    
    SectionType what_was_read = parse_section_name(state);
    
    if (what_was_read != what_section)
    {
        throw ParsingException("Read section end that doesn't match the current section");
    }
}

/********************************************************************************
 * Test parse of section end.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parsing end of section")
{
    ParserState state("$EndMeshFormat\n", 15);
    REQUIRE_NOTHROW(parse_section_end(state, SectionType::MeshFormat));
    state = ParserState("$EndMeshFormat\n", 15);
    REQUIRE_THROWS_WITH(parse_section_end(state, SectionType::Nodes),
                        "Read section end that doesn't match the current section");
    
    state = ParserState("$EndEntities\n", 13);
    REQUIRE_NOTHROW(parse_section_end(state, SectionType::Entities));
    state = ParserState("$EndEntities\n", 13);
    REQUIRE_THROWS_WITH(parse_section_end(state, SectionType::MeshFormat),
                        "Read section end that doesn't match the current section");
    
    state = ParserState("$EndNodes\n", 10);
    REQUIRE_NOTHROW(parse_section_end(state, SectionType::Nodes));
    
    state = ParserState("EndNodes", 9);
    REQUIRE_THROWS_WITH(parse_section_end(state, SectionType::Nodes), "Expected end of section");
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

inline std::string
parse_ascii_string(ParserState &state)
{
    if (state.current() != '"')
    {
        throw ParsingException("Expected ASCII string");
    }
    state.add_offset(1);

    const char *end = reinterpret_cast<const char*>(
        std::memchr(state.data(), '"', state.size() - state.offset()));

    if (end == nullptr)
    {
        throw ParsingException("Did not find closing quote for ASCII string");
    }

    std::string str(state.data(), end);
    state.add_offset(1 + end - state.data());
    return str;
}

/********************************************************************************
 * Test parse of ascii string
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parse_ascii_string")
{
    ParserState state("\"bottom_points\"", 16);
    REQUIRE(parse_ascii_string(state) == std::string("bottom_points"));

    state = ParserState("\"bottom_points", 14);
    REQUIRE_THROWS_WITH(parse_ascii_string(state), "Did not find closing quote for ASCII string");

    state = ParserState("bottom_points\"", 14);
    REQUIRE_THROWS_WITH(parse_ascii_string(state), "Expected ASCII string");
}

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

struct PhysicalNames
{
    std::vector<size_t> dims;
    std::vector<std::string> names;
};

inline PhysicalNames
parse_physical_names(ParserState &state)
{
    size_t num_physical_names = state.extract_ascii_int();

    PhysicalNames names;
    names.dims.reserve(num_physical_names);
    names.names.reserve(num_physical_names);

    for (size_t i = 0; i < num_physical_names; ++i)
    {
        names.dims.push_back(state.extract_ascii_int());

        size_t tag_value = state.extract_ascii_int();
        if (tag_value != i+1)
        {
            throw ParsingException("Expected sequential tags 1:N in PhysicalNames");
        }
        state.skip_whitespace();

        names.names.emplace_back(parse_ascii_string(state));
        
        if (state.current() != '\n')
        {
            throw ParsingException("Expected newline after name in PhysicalNames");
        }
        state.add_offset(1);
    }

    return names;
}

/********************************************************************************
 * Test parse_physical_names
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parse_physical_names")
{
    const char data[] = "\n4\n0 1 \"top_points\"\n0 2 \"bottom_points\"\n1 3 \""
                        "ports\"\n2 4 \"domain\"\n$EndPhysicalNames\n";
    
    ParserState state(data, sizeof(data));

    auto names = parse_physical_names(state);
    REQUIRE(names.dims.size() == 4);
    REQUIRE(names.names.size() == 4);

    REQUIRE(names.dims[0] == 0);
    REQUIRE(names.dims[1] == 0);
    REQUIRE(names.dims[2] == 1);
    REQUIRE(names.dims[3] == 2);

    REQUIRE(names.names[0] == "top_points");
    REQUIRE(names.names[1] == "bottom_points");
    REQUIRE(names.names[2] == "ports");
    REQUIRE(names.names[3] == "domain");

    REQUIRE_NOTHROW(parse_section_end(state, SectionType::PhysicalNames));
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

struct PointData
{
    std::array<double, 3> coords;
    std::vector<std::string> physical_tags;
};

inline PointData
parse_point_data(ParserState &state, int32_t expected, const PhysicalNames &physical_names)
{
    PointData pt;
    
    int32_t point_tag = state.extract_int();
    if (point_tag != expected+1)
    {
        throw ParsingException("Got non-sequential point tag");
    }
    state.extract_doubles(pt.coords.data(), 3);
    size_t num_tags = state.extract_size_t();
    for (size_t i = 0; i < num_tags; ++i)
    {
        size_t tag = state.extract_int()-1;
        if (tag >= physical_names.dims.size())
        {
            throw ParsingException("Physical tag out of bounds in parse_point_data");
        }
        if (physical_names.dims.at(tag) != 0)
        {
            throw ParsingException("Physical tag in parse_point_data does not refer to point");
        }
        pt.physical_tags.push_back(physical_names.names.at(tag));
    }
    return pt;
}

inline std::vector<PointData>
parse_all_points(ParserState &state, size_t num_points,
                 const PhysicalNames &physical_names)
{
    std::vector<PointData> points;
    for (size_t i = 0; i < num_points; ++i)
    {
        points.emplace_back(parse_point_data(state, i, physical_names));
    }
    return points;
}

/********************************************************************************
 * Test parsing of list of points
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parse_all_points")
{
    PhysicalNames physical_names{
        std::vector<size_t>{0, 0, 1, 2},
        std::vector<std::string>{"top_points", "bottom_points", "ports", "domain"}
    };

    const char data[] = "\x01\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0\xbf\0\0\0"
                        "\0\0\0\0\0\x01\0\0\0\0\0\0\0\x02\0\0\0\x02\0\0\0\0\0\0\0"
                        "\0\0\xf0\xbf\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\x01\0\0\0"
                        "\0\0\0\0\x01\0\0\0\x03\0\0\0\0\0\0\0\0\0\xf0?\0\0\0\0\0"
                        "\0\xf0?\0\0\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\x01\0\0\0\x04"
                        "\0\0\0\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0"
                        "\0\0\x01\0\0\0\0\0\0\0\x02\0\0\0";

    ParserState state(data, sizeof(data));
    state.set_data_size(8);

    auto pts = parse_all_points(state, 4, physical_names);
    REQUIRE(pts.size() == 4);

    REQUIRE(pts[0].coords[0] == doctest::Approx(-1.0));
    REQUIRE(pts[0].coords[1] == doctest::Approx(-1.0));
    REQUIRE(pts[0].coords[2] == 0);
    REQUIRE(pts[0].physical_tags.size() == 1);
    REQUIRE(pts[0].physical_tags[0] == "bottom_points");

    REQUIRE(pts[1].coords[0] == doctest::Approx(-1.0));
    REQUIRE(pts[1].coords[1] == doctest::Approx(1.0));
    REQUIRE(pts[1].coords[2] == 0);
    REQUIRE(pts[1].physical_tags.size() == 1);
    REQUIRE(pts[1].physical_tags[0] == "top_points");

    REQUIRE(pts[2].coords[0] == doctest::Approx(1.0));
    REQUIRE(pts[2].coords[1] == doctest::Approx(1.0));
    REQUIRE(pts[2].coords[2] == 0);
    REQUIRE(pts[2].physical_tags.size() == 1);
    REQUIRE(pts[2].physical_tags[0] == "top_points");

    REQUIRE(pts[3].coords[0] == doctest::Approx(1.0));
    REQUIRE(pts[3].coords[1] == doctest::Approx(-1.0));
    REQUIRE(pts[3].coords[2] == 0);
    REQUIRE(pts[3].physical_tags.size() == 1);
    REQUIRE(pts[3].physical_tags[0] == "bottom_points");
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

struct CurveData
{
    std::array<double, 3> minima;
    std::array<double, 3> maxima;
    std::vector<std::string> physical_tags;
    std::vector<size_t> bounding_points;
};

inline CurveData
parse_curve_data(ParserState &state, int32_t expected, size_t num_points,
                 const PhysicalNames &physical_names)
{
    CurveData curve;

    int32_t curve_tag = state.extract_int();
    if (curve_tag != expected+1)
    {
        throw ParsingException("Got non-sequential curve tag");
    }
    state.extract_doubles(curve.minima.data(), 3);
    state.extract_doubles(curve.maxima.data(), 3);

    size_t num_tags = state.extract_size_t();
    curve.physical_tags.reserve(num_tags);
    for (size_t i = 0; i < num_tags; ++i)
    {
        size_t tag = state.extract_int() - 1;
        if (tag >= physical_names.dims.size())
        {
            throw ParsingException("Physical tag out of bounds in parse_curve_data");
        }
        if (physical_names.dims.at(tag) != 1)
        {
            throw ParsingException("Physical tag in parse_curve_data does not refer to curve");
        }
        curve.physical_tags.push_back(physical_names.names.at(tag));
    }

    size_t num_bpoints = state.extract_size_t();
    curve.bounding_points.reserve(num_bpoints);
    for (size_t i = 0; i < num_bpoints; ++i)
    {
        int32_t tag = state.extract_int();
        size_t pt_tag = tag < 0 ? -(tag+1) : tag-1;
        if (pt_tag >= num_points)
        {
            throw ParsingException("Point tag in parse_curve_data out of bounds");
        }
        curve.bounding_points.push_back(pt_tag);
    }
    return curve;
}

inline std::vector<CurveData>
parse_all_curves(ParserState &state, size_t num_curves, size_t num_points,
                 const PhysicalNames &physical_names)
{
    std::vector<CurveData> curves;
    for (size_t i = 0; i < num_curves; ++i)
    {
        curves.emplace_back(parse_curve_data(state, i, num_points, physical_names));
    }
    return curves;
}

/********************************************************************************
 * Test parse_all_curves
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parse_all_curves")
{
    PhysicalNames physical_names{
        std::vector<size_t>{0, 0, 1, 2},
        std::vector<std::string>{"top_points", "bottom_points", "ports", "domain"}};

    const char data[] = "\x01\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0\xbf\0\0"
                        "\0\0\0\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0?\0\0"
                        "\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\x03\0\0\0\x02\0\0\0\0\0"
                        "\0\0\x01\0\0\0\xfe\xff\xff\xff\x02\0\0\0\0\0\0\0\0\0"
                        "\xf0\xbf\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
                        "\xf0?\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
                        "\x02\0\0\0\0\0\0\0\x02\0\0\0\xfd\xff\xff\xff\x03\0\0\0"
                        "\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\0\0\0"
                        "\0\0\0\0\0\xf0?\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\x01\0\0"
                        "\0\0\0\0\0\x03\0\0\0\x02\0\0\0\0\0\0\0\x03\0\0\0\xfc\xff"
                        "\xff\xff\x04\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0\xbf"
                        "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\xf0\xbf\0\0"
                        "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\x02\0\0\0\0\0\0\0\x04\0\0\0"
                        "\xff\xff\xff\xff";

    ParserState state(data, sizeof(data));
    state.set_data_size(8);
    const auto curves = parse_all_curves(state, 4, 4, physical_names);

    REQUIRE(curves.size() == 4);

    REQUIRE(curves[0].minima == std::array<double, 3>{-1, -1, 0});
    REQUIRE(curves[0].maxima == std::array<double, 3>{-1, 1, 0});
    REQUIRE(curves[0].physical_tags.size() == 1);
    REQUIRE(curves[0].physical_tags[0] == "ports");
    REQUIRE(curves[0].bounding_points.size() == 2);
    REQUIRE(curves[0].bounding_points[0] == 0);
    REQUIRE(curves[0].bounding_points[1] == 1);

    REQUIRE(curves[1].minima == std::array<double, 3>{-1, 1, 0});
    REQUIRE(curves[1].maxima == std::array<double, 3>{1, 1, 0});
    REQUIRE(curves[1].physical_tags.size() == 0);
    REQUIRE(curves[1].bounding_points.size() == 2);
    REQUIRE(curves[1].bounding_points[0] == 1);
    REQUIRE(curves[1].bounding_points[1] == 2);

    REQUIRE(curves[2].minima == std::array<double, 3>{1, -1, 0});
    REQUIRE(curves[2].maxima == std::array<double, 3>{1, 1, 0});
    REQUIRE(curves[2].physical_tags.size() == 1);
    REQUIRE(curves[2].physical_tags[0] == "ports");
    REQUIRE(curves[2].bounding_points.size() == 2);
    REQUIRE(curves[2].bounding_points[0] == 2);
    REQUIRE(curves[2].bounding_points[1] == 3);

    REQUIRE(curves[3].minima == std::array<double, 3>{-1, -1, 0});
    REQUIRE(curves[3].maxima == std::array<double, 3>{1, -1, 0});
    REQUIRE(curves[3].physical_tags.size() == 0);
    REQUIRE(curves[3].bounding_points.size() == 2);
    REQUIRE(curves[3].bounding_points[0] == 3);
    REQUIRE(curves[3].bounding_points[1] == 0);
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

struct SurfaceData
{
    std::array<double, 3> minima;
    std::array<double, 3> maxima;
    std::vector<std::string> physical_tags;
    std::vector<size_t> bounding_curves;
};

inline SurfaceData
parse_surface_data(ParserState &state, int32_t expected, size_t num_curves,
                   const PhysicalNames &physical_names)
{
    SurfaceData surf;

    int32_t surf_tag = state.extract_int();
    if (surf_tag != expected+1)
    {
        throw ParsingException("Got non-sequential surface tag");
    }
    state.extract_doubles(surf.minima.data(), 3);
    state.extract_doubles(surf.maxima.data(), 3);

    size_t num_tags = state.extract_size_t();
    surf.physical_tags.reserve(num_tags);
    for (size_t i = 0; i < num_tags; ++i)
    {
        size_t tag = state.extract_int() - 1;
        if (tag >= physical_names.dims.size())
        {
            throw ParsingException("Physical tag out of bounds in parse_surface_data");
        }
        if (physical_names.dims.at(tag) != 2)
        {
            throw ParsingException("Physical tag in parse_surface_data does not refer to surface");
        }
        surf.physical_tags.push_back(physical_names.names.at(tag));
    }

    size_t num_bcurves = state.extract_size_t();
    surf.bounding_curves.reserve(num_bcurves);
    for (size_t i = 0; i < num_bcurves; ++i)
    {
        size_t tag = state.extract_int() - 1;
        if (tag >= num_curves)
        {
            throw ParsingException("Curve tag in parse_surface_data out of bounds");
        }
        surf.bounding_curves.push_back(tag);
    }
    return surf;
}

inline std::vector<SurfaceData>
parse_all_surfaces(ParserState &state, size_t num_surfs, size_t num_points,
                   const PhysicalNames &physical_names)
{
    std::vector<SurfaceData> surfs;
    for (size_t i = 0; i < num_surfs; ++i)
    {
        surfs.emplace_back(parse_surface_data(state, i, num_points, physical_names));
    }
    return surfs;
}

/********************************************************************************
 * Test parse_all_surfaces
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parse_all_surfaces")
{
    PhysicalNames physical_names{
        std::vector<size_t>{0, 0, 1, 2},
        std::vector<std::string>{"top_points", "bottom_points", "ports", "domain"}};

    const char data[] = "\x01\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0\xbf\0\0"
                        "\0\0\0\0\0\0\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\xf0?\0\0\0\0"
                        "\0\0\0\0\x01\0\0\0\0\0\0\0\x04\0\0\0\x04\0\0\0\0\0\0\0"
                        "\x01\0\0\0\x02\0\0\0\x03\0\0\0\x04\0\0\0\n$EndEntities\n";

    ParserState state(data, sizeof(data));
    state.set_data_size(8);
    const auto surfaces = parse_all_surfaces(state, 1, 4, physical_names);

    REQUIRE(surfaces.size() == 1);

    REQUIRE(surfaces[0].minima == std::array<double, 3>{-1, -1, 0});
    REQUIRE(surfaces[0].maxima == std::array<double, 3>{1, 1, 0});
    REQUIRE(surfaces[0].physical_tags.size() == 1);
    REQUIRE(surfaces[0].physical_tags[0] == "domain");
    REQUIRE(surfaces[0].bounding_curves.size() == 4);
    REQUIRE(surfaces[0].bounding_curves[0] == 0);
    REQUIRE(surfaces[0].bounding_curves[1] == 1);
    REQUIRE(surfaces[0].bounding_curves[2] == 2);
    REQUIRE(surfaces[0].bounding_curves[3] == 3);
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

struct Entities
{
    std::vector<PointData> points;
    std::vector<CurveData> curves;
    std::vector<SurfaceData> surfs;
};

inline Entities
parse_entities(ParserState &state, const PhysicalNames &physical_names)
{
    size_t num_points = state.extract_size_t();
    size_t num_curves = state.extract_size_t();
    size_t num_surfaces = state.extract_size_t();
    size_t num_volumes = state.extract_size_t();
    
    if (num_volumes != 0)
    {
        throw ParsingException("3D meshes not implemented");
    }

    auto points = parse_all_points(state, num_points, physical_names);
    auto curves = parse_all_curves(state, num_curves, num_points, physical_names);
    auto surfaces = parse_all_surfaces(state, num_surfaces, num_curves, physical_names);
    return Entities{std::move(points), std::move(curves), std::move(surfaces)};
}

/********************************************************************************
 * Test parse_entities
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test parse_entities")
{
    PhysicalNames physical_names{
        std::vector<size_t>{0, 0, 1, 2},
        std::vector<std::string>{"top_points", "bottom_points", "ports", "domain"}};

    const char data[] = "\x04\0\0\0\0\0\0\0\x04\0\0\0\0\0\0\0\x01\0\0\0\0\0\0\0"
                        "\0\0\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0"
                        "\0\0\xf0\xbf\0\0\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\x02\0\0"
                        "\0\x02\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0?\0\0\0"
                        "\0\0\0\0\0\x01\0\0\0\0\0\0\0\x01\0\0\0\x03\0\0\0\0\0\0"
                        "\0\0\0\xf0?\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\x01\0\0\0"
                        "\0\0\0\0\x01\0\0\0\x04\0\0\0\0\0\0\0\0\0\xf0?\0\0\0\0\0"
                        "\0\xf0\xbf\0\0\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\x02\0\0\0"
                        "\x01\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0\xbf\0\0"
                        "\0\0\0\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0?\0\0\0"
                        "\0\0\0\0\0\x01\0\0\0\0\0\0\0\x03\0\0\0\x02\0\0\0\0\0\0"
                        "\0\x01\0\0\0\xfe\xff\xff\xff\x02\0\0\0\0\0\0\0\0\0\xf0"
                        "\xbf\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\0\0\0\0\0\0\xf0?"
                        "\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\x02\0"
                        "\0\0\0\0\0\0\x02\0\0\0\xfd\xff\xff\xff\x03\0\0\0\0\0\0"
                        "\0\0\0\xf0?\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\0\0\0\0\0\0"
                        "\0\0\xf0?\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\x01\0\0\0\0"
                        "\0\0\0\x03\0\0\0\x02\0\0\0\0\0\0\0\x03\0\0\0\xfc\xff"
                        "\xff\xff\x04\0\0\0\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\xf0"
                        "\xbf\0\0\0\0\0\0\0\0\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\xf0"
                        "\xbf\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\x02\0\0\0\0\0\0\0"
                        "\x04\0\0\0\xff\xff\xff\xff\x01\0\0\0\0\0\0\0\0\0\xf0"
                        "\xbf\0\0\0\0\0\0\xf0\xbf\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
                        "\xf0?\0\0\0\0\0\0\xf0?\0\0\0\0\0\0\0\0\x01\0\0\0\0\0\0"
                        "\0\x04\0\0\0\x04\0\0\0\0\0\0\0\x01\0\0\0\x02\0\0\0\x03"
                        "\0\0\0\x04\0\0\0\n$EndEntities\n";

    ParserState state(data, sizeof(data));
    state.set_data_size(8);

    const auto entities = parse_entities(state, physical_names);

    const auto &pts = entities.points;
    {
        REQUIRE(pts.size() == 4);

        REQUIRE(pts[0].coords[0] == doctest::Approx(-1.0));
        REQUIRE(pts[0].coords[1] == doctest::Approx(-1.0));
        REQUIRE(pts[0].coords[2] == 0);
        REQUIRE(pts[0].physical_tags.size() == 1);
        REQUIRE(pts[0].physical_tags[0] == "bottom_points");

        REQUIRE(pts[1].coords[0] == doctest::Approx(-1.0));
        REQUIRE(pts[1].coords[1] == doctest::Approx(1.0));
        REQUIRE(pts[1].coords[2] == 0);
        REQUIRE(pts[1].physical_tags.size() == 1);
        REQUIRE(pts[1].physical_tags[0] == "top_points");

        REQUIRE(pts[2].coords[0] == doctest::Approx(1.0));
        REQUIRE(pts[2].coords[1] == doctest::Approx(1.0));
        REQUIRE(pts[2].coords[2] == 0);
        REQUIRE(pts[2].physical_tags.size() == 1);
        REQUIRE(pts[2].physical_tags[0] == "top_points");

        REQUIRE(pts[3].coords[0] == doctest::Approx(1.0));
        REQUIRE(pts[3].coords[1] == doctest::Approx(-1.0));
        REQUIRE(pts[3].coords[2] == 0);
        REQUIRE(pts[3].physical_tags.size() == 1);
        REQUIRE(pts[3].physical_tags[0] == "bottom_points");
    }

    const auto &curves = entities.curves;
    {
        REQUIRE(curves.size() == 4);

        REQUIRE(curves[0].minima == std::array<double, 3>{-1, -1, 0});
        REQUIRE(curves[0].maxima == std::array<double, 3>{-1, 1, 0});
        REQUIRE(curves[0].physical_tags.size() == 1);
        REQUIRE(curves[0].physical_tags[0] == "ports");
        REQUIRE(curves[0].bounding_points.size() == 2);
        REQUIRE(curves[0].bounding_points[0] == 0);
        REQUIRE(curves[0].bounding_points[1] == 1);

        REQUIRE(curves[1].minima == std::array<double, 3>{-1, 1, 0});
        REQUIRE(curves[1].maxima == std::array<double, 3>{1, 1, 0});
        REQUIRE(curves[1].physical_tags.size() == 0);
        REQUIRE(curves[1].bounding_points.size() == 2);
        REQUIRE(curves[1].bounding_points[0] == 1);
        REQUIRE(curves[1].bounding_points[1] == 2);

        REQUIRE(curves[2].minima == std::array<double, 3>{1, -1, 0});
        REQUIRE(curves[2].maxima == std::array<double, 3>{1, 1, 0});
        REQUIRE(curves[2].physical_tags.size() == 1);
        REQUIRE(curves[2].physical_tags[0] == "ports");
        REQUIRE(curves[2].bounding_points.size() == 2);
        REQUIRE(curves[2].bounding_points[0] == 2);
        REQUIRE(curves[2].bounding_points[1] == 3);

        REQUIRE(curves[3].minima == std::array<double, 3>{-1, -1, 0});
        REQUIRE(curves[3].maxima == std::array<double, 3>{1, -1, 0});
        REQUIRE(curves[3].physical_tags.size() == 0);
        REQUIRE(curves[3].bounding_points.size() == 2);
        REQUIRE(curves[3].bounding_points[0] == 3);
        REQUIRE(curves[3].bounding_points[1] == 0);
    }

    const auto &surfaces = entities.surfs;
    {
        REQUIRE(surfaces.size() == 1);

        REQUIRE(surfaces[0].minima == std::array<double, 3>{-1, -1, 0});
        REQUIRE(surfaces[0].maxima == std::array<double, 3>{1, 1, 0});
        REQUIRE(surfaces[0].physical_tags.size() == 1);
        REQUIRE(surfaces[0].physical_tags[0] == "domain");
        REQUIRE(surfaces[0].bounding_curves.size() == 4);
        REQUIRE(surfaces[0].bounding_curves[0] == 0);
        REQUIRE(surfaces[0].bounding_curves[1] == 1);
        REQUIRE(surfaces[0].bounding_curves[2] == 2);
        REQUIRE(surfaces[0].bounding_curves[3] == 3);
    }
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 *******************************************************************************/

/*
inline void
parse_gmsh_file(ParserState &state)
{
    std::array<bool, 13> section_parsed { false, false, false, false, false,
                                          false, false, false, false, false,
                                          false, false, false };
    parse_mesh_format(state);
    parse_section_end(state, SectionType::MeshFormat);
    section_parsed[static_cast<size_t>(SectionType::MeshFormat)] = true;
    
    while (state.offset() < state.size())
    {
        auto what_section = parse_section_header(state);
        if (section_parsed[static_cast<size_t>(what_section)])
        {
            throw ParsingException("Multiple specification of section");
        }
        section_parsed[static_cast<size_t>(what_section)] = true;
        
        PhysicalNames physical_names;

        switch(what_section)
        {
            case SectionType::PhysicalNames:
                physical_names = parse_physical_names(state);
                break;
            case SectionType::Entities:
                auto entities = parse_entities(state, physical_names);
                break;
            case SectionType::PartitionedEntities:
                throw ParsingException("Nothing is implemented regarding partitioned entities");
            case SectionType::Nodes:
                auto nodes = parse_nodes(state);
                break;
            case SectionType::Elements:
                auto elements = parse_elements(state);
                break;
            case SectionType::Periodic:
                throw ParsingException("Nothing is implemented regarding periodic links");
            case SectionType::GhostElements:
                throw ParsingException("Nothing is implemented regarding ghost elements");
            case SectionType::Parametrizations:
                throw ParsingException("Nothing is implemented regarding parametrizations");
            case SectionType::NodeData:
            case SectionType::ElementData:
            case SectionType::ElementNodeData:
                printf("Informational: nothing is currently done with NodeData, "
                       "ElementData, or ElementNodeData sections\n");
                skip_section(state, what_section);
                break;
            case SectionType::InterpolationScheme:
                throw ParsingException("Nothing is implemented regarding interpolation schemes");
            default:
                throw ParsingException("It should be impossible to hit this exception");
        }
        parse_section_end(state, what_section);
    }
}

*/

} // namespace gmsh

} // namespace msh

#endif // READ_GMSH_HPP
