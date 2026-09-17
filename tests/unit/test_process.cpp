#include "doctest.h"
#include "rnaseq/process.hpp"
#include <filesystem>
#include <fstream>
#include <unistd.h>
TEST_CASE("argv subprocess preserves literals and reports failures") {
    auto log = std::filesystem::temp_directory_path() / ("rnaseq-process-" + std::to_string(getpid()));
    rnaseq::run_process({"/bin/printf", "%s", "space ; $(touch NEVER)"}, log);
    std::ifstream input(log); std::string result; std::getline(input, result);
    CHECK(result == "space ; $(touch NEVER)");
    CHECK_THROWS_AS(rnaseq::run_process({"/bin/sh", "-c", "exit 17"}, log), rnaseq::ProcessError);
    CHECK_THROWS(rnaseq::run_process({"/absent-executable"}, log));
    std::filesystem::remove(log);
}
TEST_CASE("staged env interpreter precedes ambient PATH without parent mutation") {
    char pattern[] = "/tmp/rnaseq-interpreter-XXXXXX";
    const char* created = mkdtemp(pattern);
    REQUIRE(created != nullptr);
    const std::filesystem::path root = created;
    const auto bin = root / "literal ; $(unused)";
    const auto ambient = root / "ambient";
    std::filesystem::create_directory(bin);
    std::filesystem::create_directory(ambient);
    const auto wrapper = bin / "wrapper";
    const auto interpreter = bin / "python";
    { std::ofstream out(wrapper); out << "#!/usr/bin/env python\n"; }
    { std::ofstream out(interpreter); out << "#!/bin/sh\nprintf '%s' \"$1\"\n"; }
    std::filesystem::permissions(wrapper, std::filesystem::perms::owner_all);
    std::filesystem::permissions(interpreter, std::filesystem::perms::owner_all);
    const char* previous = std::getenv("PATH");
    const bool had_path = previous != nullptr;
    const std::string original = previous ? previous : "";
    setenv("PATH", ambient.c_str(), 1);
    CHECK_NOTHROW(rnaseq::run_process({wrapper.string()}, root / "log"));
    CHECK(std::string(std::getenv("PATH")) == ambient.string());
    if (had_path) setenv("PATH", original.c_str(), 1); else unsetenv("PATH");
    std::ifstream input(root / "log"); std::string result; std::getline(input, result);
    CHECK(result == wrapper.string());
    std::filesystem::remove_all(root);
}
