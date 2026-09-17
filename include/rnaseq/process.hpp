#pragma once
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
namespace rnaseq {
// Linux single-threaded driver: bounded fork workers, fail-fast cancellation,
// SIGINT/SIGTERM forwarding and descendant reaping across process groups.
// Call only when this process owns no unrelated children.
void run_jobs(std::size_t count, std::size_t concurrency, const std::function<void(std::size_t)>& action);
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
