#pragma once
#include <filesystem>
#include <string>

namespace rnaseq {
// Hash a readable regular file with sha256sum from PATH, without a shell.
std::string sha256_file(const std::filesystem::path& path);
}
