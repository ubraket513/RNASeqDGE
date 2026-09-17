#include "rnaseq/table.hpp"
#include "csv.hpp"
#include <algorithm>
#include <charconv>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace rnaseq {
namespace {
[[noreturn]] void fail(const std::string& source, std::size_t record,
                       const std::string& column, const std::string& message) {
    throw std::runtime_error(source + ": record " + std::to_string(record) +
                             ", column " + column + ": " + message);
}

bool valid_utf8(std::string_view text) {
    for (std::size_t i = 0; i < text.size();) {
        const auto first = static_cast<unsigned char>(text[i++]);
        if (first < 0x80) continue;
        unsigned length = 0;
        std::uint32_t code = 0;
        if (first >= 0xC2 && first <= 0xDF) { length = 1; code = first & 0x1F; }
        else if (first >= 0xE0 && first <= 0xEF) { length = 2; code = first & 0x0F; }
        else if (first >= 0xF0 && first <= 0xF4) { length = 3; code = first & 0x07; }
        else return false;
        if (i + length > text.size()) return false;
        for (unsigned j = 0; j < length; ++j) {
            const auto next = static_cast<unsigned char>(text[i++]);
            if ((next & 0xC0) != 0x80) return false;
            code = (code << 6) | (next & 0x3F);
        }
        if ((length == 2 && code < 0x800) || (length == 3 && code < 0x10000) ||
            code > 0x10FFFF || (code >= 0xD800 && code <= 0xDFFF)) return false;
    }
    return true;
}

// Enforce the restricted manifest dialect before handing decoding to csv-parser.
// The upstream parser intentionally accepts a wider CSV dialect than v1.
std::size_t check_record(const std::string& line, const std::string& source, std::size_t record) {
    if (line.empty()) fail(source, record, "1", "blank records are forbidden");
    if (!valid_utf8(line)) fail(source, record, "1", "invalid UTF-8");
    if (line.find("\xEF\xBB\xBF") != std::string::npos)
        fail(source, record, "1", "BOM is allowed only at file start");
    if (line[0] == '#') fail(source, record, "1", "comments are forbidden");
    enum class State { start, bare, quoted, closed };
    State state = State::start;
    std::size_t column = 1;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '\0' || c == '\r' || c == '\n')
            fail(source, record, std::to_string(column), "NUL/CR/LF inside a cell");
        if (state == State::quoted) {
            if (c == '\t') fail(source, record, std::to_string(column), "tab inside a cell");
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') ++i;
                else state = State::closed;
            }
        } else if (c == '\t') { ++column; state = State::start; }
        else if (state == State::closed) {
            fail(source, record, std::to_string(column), "characters after closing quote");
        } else if (c == '"') {
            if (state != State::start) fail(source, record, std::to_string(column), "quote in unquoted cell");
            state = State::quoted;
        } else state = State::bare;
    }
    if (state == State::quoted) fail(source, record, std::to_string(column), "unclosed quoted cell");
    return column;
}

bool valid_id(const std::string& value) {
    const auto alnum = [](unsigned char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
    };
    return !value.empty() && alnum(value.front()) &&
        std::all_of(value.begin(), value.end(), [&](unsigned char c) {
            return alnum(c) || c == '_' || c == '.' || c == '-';
        });
}

std::size_t column(const Table& table, const std::string& name) {
    const auto it = std::find(table.header.begin(), table.header.end(), name);
    if (it == table.header.end()) fail(table.source, 1, name, "required column missing");
    return static_cast<std::size_t>(it - table.header.begin());
}
}

Table read_table(std::istream& input, const std::string& source) {
    Table table{source, {}, {}};
    std::string line, normalized;
    std::size_t record = 0, width = 0;
    while (std::getline(input, line)) {
        ++record;
        if (record == 1 && line.starts_with("\xEF\xBB\xBF")) line.erase(0, 3);
        // A terminal CR is legal only as part of CRLF, not at EOF.
        if (!line.empty() && line.back() == '\r' && !input.eof()) line.pop_back();
        const auto fields = check_record(line, source, record);
        if (record == 1) width = fields;
        else if (fields != width) fail(source, record, std::to_string(fields), "ragged row");
        normalized += line + '\n';
    }
    if (input.bad() || (input.fail() && !input.eof())) fail(source, record + 1, "1", "read failure");
    if (record == 0) fail(source, 1, "1", "empty file");
    std::istringstream data(normalized);
    csv::CSVFormat format;
    format.delimiter('\t').quote('"').header_row(0)
          .variable_columns(csv::VariableColumnPolicy::THROW).threading(false);
    csv::CSVReader reader(data, format);
    table.header = reader.get_col_names();
    std::set<std::string> headers;
    for (std::size_t i = 0; i < table.header.size(); ++i) {
        const auto& name = table.header[i];
        if (name.empty() || !headers.insert(name).second)
            fail(source, 1, std::to_string(i + 1), "empty or duplicate header");
    }
    for (auto& row : reader) {
        std::vector<std::string> cells;
        for (auto& field : row) cells.push_back(field.get<std::string>());
        table.rows.push_back(std::move(cells));
    }
    if (table.header.size() != width || table.rows.size() + 1 != record)
        fail(source, 1, "1", "parser did not preserve all records");
    return table;
}

Table read_table(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) fail(path, 1, "1", "cannot open input");
    return read_table(input, path);
}

std::uint64_t count_value(const std::string& value, const std::string& context) {
    std::uint64_t number = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), number);
    if (value.empty() || error != std::errc{} || end != value.data() + value.size() ||
        !std::all_of(value.begin(), value.end(), [](char c) { return c >= '0' && c <= '9'; }))
        throw std::runtime_error(context + ": expected uint64 integer (invalid value or overflow)");
    return number;
}

std::uint64_t checked_add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        throw std::overflow_error("count sum exceeds uint64 range");
    return a + b;
}

void validate_samples(const Table& table) {
    const auto id = column(table, "sample_id"), condition = column(table, "condition");
    if (table.rows.empty()) fail(table.source, 2, "sample_id", "no samples");
    std::set<std::string> seen;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        const auto& row = table.rows[i];
        if (!valid_id(row[id]) || !seen.insert(row[id]).second)
            fail(table.source, i + 2, "sample_id", "invalid or duplicate sample ID");
        if (row[condition].empty()) fail(table.source, i + 2, "condition", "required cell empty");
    }
}

void validate_counts(const Table& table) {
    if (table.header.size() < 2 || table.header[0] != "gene_id")
        fail(table.source, 1, "gene_id", "expected gene_id then sample columns");
    for (std::size_t i = 1; i < table.header.size(); ++i)
        if (!valid_id(table.header[i])) fail(table.source, 1, table.header[i], "invalid sample ID");
    if (table.rows.empty()) fail(table.source, 2, "gene_id", "no genes");
    std::set<std::string> genes;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        const auto& row = table.rows[i];
        if (row[0].empty() || !genes.insert(row[0]).second)
            fail(table.source, i + 2, "gene_id", "empty or duplicate gene ID");
        for (std::size_t j = 1; j < row.size(); ++j)
            count_value(row[j], table.source + ": record " + std::to_string(i + 2) +
                        ", column " + table.header[j]);
    }
}
}
