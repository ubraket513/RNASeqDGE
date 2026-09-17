#pragma once

#include <filesystem>

namespace rnaseq {
struct Table;

void validate_run_mapping(const Table &samples, const Table &runs, bool check_fastq_files = true);

void validate_bundle(const std::filesystem::path &samples, const std::filesystem::path &runs,
                     const std::filesystem::path &references, const std::filesystem::path &analysis,
                     const std::filesystem::path &contrasts);

} // namespace rnaseq
