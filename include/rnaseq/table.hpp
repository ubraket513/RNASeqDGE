#pragma once
#include <cstddef>
#include <cstdint>
#include <istream>
#include <string>
#include <vector>

namespace rnaseq {
class TsvReader {
public:
    TsvReader(std::istream& input, std::string source, bool leading_comments = false);
    bool next(std::vector<std::string>& row);
    std::size_t record() const;

private:
    std::istream& input_;
    std::string source_;
    bool leading_comments_;
    bool saw_data_ = false;
    std::size_t record_ = 0;
    std::size_t width_ = 0;
};

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
