#pragma once

#include <filesystem>

namespace rnaseq {

void merge_counts(const std::filesystem::path& samples,
                  const std::filesystem::path& runs,
                  const std::filesystem::path& genes,
                  const std::filesystem::path& inputs,
                  const std::filesystem::path& output,
                  const std::filesystem::path& annotation = {});

void validate_counts_for_samples(const std::filesystem::path& samples,
                                 const std::filesystem::path& counts);

void validate_annotation(const std::filesystem::path& genes,
                         const std::filesystem::path& annotation);

}
