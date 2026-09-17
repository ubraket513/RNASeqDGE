#include "rnaseq/table.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: rnaseq validate-tsv|validate-samples|validate-counts FILE\n";
        return (std::cout << std::flush) ? 0 : 1;
    }
    if (argc != 3) { std::cerr << "Use rnaseq --help for usage\n"; return 2; }
    try {
        const std::string command = argv[1];
        if (command != "validate-tsv" && command != "validate-samples" && command != "validate-counts")
            throw std::runtime_error("unknown command: " + command);
        const auto table = rnaseq::read_table(argv[2]);
        if (command == "validate-samples") rnaseq::validate_samples(table);
        if (command == "validate-counts") rnaseq::validate_counts(table);
        std::cout << "OK: " << table.rows.size() << " records, " << table.header.size() << " columns\n";
        return (std::cout << std::flush) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
