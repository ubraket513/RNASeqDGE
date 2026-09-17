#pragma once
#include <cstdint>
#include <istream>
#include <string>
#include <vector>

namespace rnaseq {
struct Table {
    std::string source;
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};
Table read_table(std::istream& input, const std::string& source);
Table read_table(const std::string& path);
std::uint64_t count_value(const std::string& value, const std::string& context);
std::uint64_t checked_add(std::uint64_t a, std::uint64_t b);
void validate_samples(const Table& table);
void validate_counts(const Table& table);
}
