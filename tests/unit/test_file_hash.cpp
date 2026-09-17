#include "doctest.h"
#include "rnaseq/file_hash.hpp"
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace {
struct HashFixture {
    std::filesystem::path directory;
    HashFixture() {
        char pattern[] = "/tmp/rnaseq-hash-XXXXXX";
        const auto created = mkdtemp(pattern);
        if (!created) throw std::runtime_error("cannot create hash fixture");
        directory = created;
    }
    ~HashFixture() { std::error_code error; std::filesystem::remove_all(directory, error); }
    std::filesystem::path write(const std::string& name, const std::string& data) {
        auto path = directory / name;
        std::ofstream output(path, std::ios::binary);
        output << data;
        output.close();
        if (!output) throw std::runtime_error("cannot write hash fixture");
        return path;
    }
};
struct HashPath {
    std::string old;
    bool present;
    explicit HashPath(const std::filesystem::path& directory) {
        const auto value = std::getenv("PATH");
        present = value != nullptr;
        if (present) old = value;
        if (setenv("PATH", directory.c_str(), 1)) throw std::runtime_error("setenv failed");
    }
    ~HashPath() { if (present) setenv("PATH", old.c_str(), 1); else unsetenv("PATH"); }
};
}

TEST_CASE("file hashing matches independent SHA-256 vectors") {
    HashFixture fixture;
    CHECK(rnaseq::sha256_file(fixture.write("empty", "")) ==
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    CHECK(rnaseq::sha256_file(fixture.write("- abc ; $x `x` ' \".txt", "abc")) ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    CHECK(rnaseq::sha256_file(fixture.write("large", std::string(1000000, 'a'))) ==
          "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
    std::filesystem::create_symlink(fixture.directory / "empty", fixture.directory / "link");
    CHECK(rnaseq::sha256_file(fixture.directory / "link") ==
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}

TEST_CASE("file hashing works with closed standard descriptors") {
    HashFixture fixture;
    const auto input = fixture.write("input", "abc");
    int saved[3];
    for (int fd = 0; fd < 3; ++fd) {
        saved[fd] = fcntl(fd, F_DUPFD_CLOEXEC, 3);
        REQUIRE(saved[fd] >= 3);
    }
    for (int fd = 0; fd < 3; ++fd) close(fd);
    std::string result, failure;
    try { result = rnaseq::sha256_file(input); }
    catch (const std::exception& error) { failure = error.what(); }
    bool restored = true;
    for (int fd = 0; fd < 3; ++fd) {
        if (dup2(saved[fd], fd) < 0) restored = false;
        close(saved[fd]);
    }
    REQUIRE(restored);
    CHECK(failure.empty());
    CHECK(result == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}

TEST_CASE("file hashing rejects missing and nonregular inputs without hanging") {
    HashFixture fixture;
    CHECK_THROWS_AS(rnaseq::sha256_file(fixture.directory / "missing"), std::runtime_error);
    CHECK_THROWS_AS(rnaseq::sha256_file(fixture.directory), std::runtime_error);
    const auto fifo = fixture.directory / "fifo";
    REQUIRE(mkfifo(fifo.c_str(), 0600) == 0);
    CHECK_THROWS_AS(rnaseq::sha256_file(fifo), std::runtime_error);
}

TEST_CASE("file hashing propagates missing utility, nonzero status and malformed output") {
    HashFixture fixture;
    const auto input = fixture.write("input", "abc");
    HashPath path(fixture.directory);
    CHECK_THROWS_AS(rnaseq::sha256_file(input), std::runtime_error);
    for (const std::string body : {"exit 9\n", "printf 'invalid digest\\n'\n",
            "printf 'za7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad *-\\n'\n",
            "printf 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad  -\\n'\n",
            "printf 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad *-'\n",
            "printf 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad *-\\nx'\n",
            "printf '%0200d' 0\n",
            "printf 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad *-\\n'; exit 7\n"}) {
        const auto utility = fixture.write("sha256sum", "#!/bin/sh\n" + body);
        REQUIRE(chmod(utility.c_str(), 0700) == 0);
        CHECK_THROWS_AS(rnaseq::sha256_file(input), std::runtime_error);
    }
}
