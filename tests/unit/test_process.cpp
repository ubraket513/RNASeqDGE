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
TEST_CASE("bounded fork jobs overlap and propagate worker failure") {
    char pattern[] = "/tmp/rnaseq-jobs-XXXXXX";
    const auto root = std::filesystem::path(mkdtemp(pattern));
    CHECK_NOTHROW(rnaseq::run_jobs(2, 2, [&](std::size_t i) {
        std::ofstream(root / std::to_string(i)) << "ready";
        for(int n=0;n<200 && !std::filesystem::exists(root / std::to_string(1-i));++n)usleep(10000);
        if(!std::filesystem::exists(root / std::to_string(1-i)))throw std::runtime_error("no overlap");
    }));
    CHECK(std::filesystem::exists(root / "0"));
    CHECK(std::filesystem::exists(root / "1"));
    CHECK_THROWS_AS(rnaseq::run_jobs(2, 2, [&](std::size_t i) {
        if(i==0) {
            for(int n=0;n<200 && !std::filesystem::exists(root/"descendant");++n)usleep(10000);
            if(!std::filesystem::exists(root/"descendant"))throw std::runtime_error("sibling did not start");
            throw rnaseq::ProcessError("expected failure",17);
        }
        rnaseq::run_process({"/bin/sh","-c","trap '' TERM; sleep 200 & echo $! > '"+(root/"descendant").string()+"'; wait"},root/"log");
    }),rnaseq::ProcessError);
    std::ifstream input(root/"descendant"); int pid=0;
    REQUIRE(static_cast<bool>(input>>pid));
    CHECK(!std::filesystem::exists("/proc/"+std::to_string(pid)));
    std::filesystem::remove_all(root);
}
