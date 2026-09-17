#include "rnaseq/file_hash.hpp"
#include <algorithm>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <fcntl.h>
#include <sys/prctl.h>
#include <stdexcept>
#include <string_view>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

namespace rnaseq {
namespace {
struct Descriptor {
    int value = -1;
    Descriptor() = default;
    explicit Descriptor(int fd) : value(fd) {}
    Descriptor(const Descriptor&) = delete;
    Descriptor& operator=(const Descriptor&) = delete;
    ~Descriptor() { if (value >= 0) close(value); }
};

}

std::string sha256_file(const std::filesystem::path& path) {
    const auto fail = [&](const std::string& message) -> void {
        throw std::runtime_error(path.string() + ": " + message);
    };
    // Nonblocking open avoids hanging on FIFO input before fstat rejects it.
    Descriptor opened(open(path.c_str(), O_RDONLY | O_NONBLOCK | O_CLOEXEC));
    if (opened.value < 0) fail(std::string("cannot open input: ") + std::strerror(errno));
    struct stat info{};
    if (fstat(opened.value, &info) != 0 || !S_ISREG(info.st_mode))
        fail("expected readable regular file");
    // Keep sources above stderr even if the caller closed standard descriptors.
    Descriptor input(fcntl(opened.value, F_DUPFD_CLOEXEC, 3));
    if (input.value < 0) fail("cannot duplicate input descriptor");
    int ends[2];
    if (pipe(ends) != 0) fail("cannot create checksum pipe");
    Descriptor raw_read(ends[0]), raw_write(ends[1]);
    Descriptor output(fcntl(raw_read.value, F_DUPFD_CLOEXEC, 3));
    Descriptor writer(fcntl(raw_write.value, F_DUPFD_CLOEXEC, 3));
    if (output.value < 0 || writer.value < 0) fail("cannot duplicate checksum pipe");
    close(raw_read.value); raw_read.value = -1;
    close(raw_write.value); raw_write.value = -1;

    char executable[] = "sha256sum", binary[] = "--binary";
    char* arguments[] = {executable, binary, nullptr};
    const pid_t parent = getpid();
    const pid_t child = fork();
    if (child < 0) fail("cannot fork checksum subprocess");
    if (child == 0) {
        // Large input hashing must stop when its driver is cancelled. Check the
        // parent again after installing the death signal to close the fork race.
        if (prctl(PR_SET_PDEATHSIG, SIGKILL) != 0 || getppid() != parent) _exit(125);
        if (dup2(input.value, STDIN_FILENO) < 0 || dup2(writer.value, STDOUT_FILENO) < 0) _exit(126);
        execvp(arguments[0], arguments);
        _exit(127);
    }
    close(writer.value); writer.value = -1;

    // Fixed output bound; unexpected utility output cannot grow native memory.
    char result[69];
    std::size_t size = 0;
    bool read_failed = false;
    while (size < sizeof(result)) {
        const auto count = read(output.value, result + size, sizeof(result) - size);
        if (count == 0) break;
        if (count < 0) {
            if (errno == EINTR) continue;
            read_failed = true;
            break;
        }
        size += static_cast<std::size_t>(count);
    }
    if (read_failed || size == sizeof(result)) kill(child, SIGKILL);
    int status = 0;
    pid_t waited;
    do { waited = waitpid(child, &status, 0); } while (waited < 0 && errno == EINTR);
    if (waited < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
        fail("sha256sum failed");
    if (read_failed || size != 68 || std::string_view(result + 64, 4) != " *-\n" ||
        !std::all_of(result, result + 64, [](char c) {
            return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
        }))
        fail("invalid sha256sum output");
    return std::string(result, 64);
}
}
