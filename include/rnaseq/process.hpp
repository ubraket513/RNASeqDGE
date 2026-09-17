#pragma once
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
namespace rnaseq {
struct ProcessError : std::runtime_error {
    int status;
    ProcessError(const std::string& message, int code) : std::runtime_error(message), status(code) {}
};
// Synchronous argv execution; merged stdout/stderr log. No shell or PATH lookup.
// Child PATH prepends the executable directory for staged wrapper interpreters;
// the parent environment is unchanged. Direct executable lookup remains absolute.
// SIGINT/SIGTERM terminate the child's process group and are returned as 128+signal.
// Optional stdout capture keeps scheduler protocol output separate from stderr.
void run_process(const std::vector<std::string>& argv, const std::filesystem::path& log,
                 const std::filesystem::path& stdout_log = {});
}
