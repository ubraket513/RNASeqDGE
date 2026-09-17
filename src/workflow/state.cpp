#include "internal.hpp"
#include "rnaseq/file_hash.hpp"
#include <algorithm>
#include <chrono>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/file.h>
#include <sys/resource.h>
#include <unistd.h>

namespace rnaseq::workflow_detail {
void fail(const std::string &s) {
    throw std::runtime_error("workflow: " + s);
}
std::string read(const fs::path &p) {
    std::ifstream f(p);
    if (!f)
        fail("cannot read " + p.string());
    return {std::istreambuf_iterator<char>(f), {}};
}
void write(const fs::path &p, const std::string &s) {
    std::ofstream f(p);
    f << s;
    f.close();
    if (!f)
        fail("cannot write " + p.string());
}
void atomic_write(const fs::path &p, const std::string &s) {
    const auto tmp = p.string() + ".tmp-" + std::to_string(getpid());
    write(tmp, s);
    fs::rename(tmp, p);
}

Lock::Lock(const fs::path &p, bool shared, bool wait) {
    fd = open(p.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0600);
    if (fd < 0)
        fail("cannot open lock " + p.string());
    if (flock(fd, (shared ? LOCK_SH : LOCK_EX) | (wait ? 0 : LOCK_NB)) != 0) {
        close(fd);
        fd = -1;
        fail("busy lock " + p.string());
    }
}
Lock::~Lock() {
    if (fd >= 0)
        close(fd);
}
void quarantine(const fs::path &p) {
    if (!fs::exists(p))
        return;
    std::string suffix = ".invalid-" + std::to_string(getpid());
    int n = 0;
    auto target = p.string() + suffix;
    while (fs::exists(target))
        target = p.string() + suffix + "-" + std::to_string(++n);
    fs::rename(p, target);
    std::cout << "quarantine " << target << '\n';
}
// Output manifests include every regular file; additions, removals and changes invalidate.
std::string inventory(const fs::path &dir) {
    std::vector<fs::path> files;
    for (const auto &e : fs::recursive_directory_iterator(dir)) {
        if (e.is_symlink())
            fail("symlink in stage output " + e.path().string());
        if (e.is_regular_file() && e.path() != dir / "STAGE.tsv" &&
            !(e.path().parent_path() == dir &&
              e.path().filename().string().find("STAGE.tsv.tmp-") == 0))
            files.push_back(e.path());
    }
    std::sort(files.begin(), files.end());
    std::string s;
    for (const auto &p : files)
        s += sha256_file(p) + "\t" + p.lexically_relative(dir).string() + "\n";
    if (files.empty())
        fail("empty stage " + dir.string());
    return s;
}
bool valid(const fs::path &dir, const std::string &identity) {
    try {
        return fs::is_regular_file(dir / "STAGE.tsv") &&
               read(dir / "STAGE.tsv") == identity + inventory(dir);
    } catch (const std::exception &) {
        return false;
    }
}
void Workflow::verify_sources() const {
    std::string current = "sha256\tpath\n";
    for (const auto &p : sources)
        current += sha256_file(p) + "\t" + p.string() + "\n";
    if (current != source_identity)
        fail("inputs changed during execution; use a new run directory");
}

void Workflow::snapshot(bool create) {
    const auto dir = run / "snapshot";
    if (fs::exists(dir)) {
        if (!valid(dir, "snapshot-v1\n") || read(dir / "sources.tsv") != source_identity)
            fail("snapshot/source/config/tool hash changed; use a new run directory");
        for (const auto &[name, data] : snapshots)
            if (read(dir / name) != data)
                fail("snapshot differs; use a new run directory");
        return;
    }
    if (!create)
        fail("snapshot is missing; prepare with workflow-local or workflow-submit");
    auto tmp = run / ("snapshot.partial-" + std::to_string(getpid()));
    quarantine(tmp);
    fs::create_directory(tmp);
    for (const auto &[name, data] : snapshots)
        write(tmp / name, data);
    write(tmp / "sources.tsv", source_identity);
    verify_sources();
    atomic_write(tmp / "STAGE.tsv", "snapshot-v1\n" + inventory(tmp));
    fs::rename(tmp, dir);
}

std::string Workflow::base_id() const {
    return "workflow-v1\n" + sha256_file(run / "snapshot" / "STAGE.tsv") + "\n";
}

void Workflow::stage(const fs::path &dir, const std::string &id,
                     const std::function<void()> &action) {
    using Clock = std::chrono::steady_clock;
    const auto started = Clock::now();
    auto elapsed = [](auto start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    rusage self_before{}, child_before{};
    getrusage(RUSAGE_SELF, &self_before);
    getrusage(RUSAGE_CHILDREN, &child_before);
    double validation = 0, execution = 0, sources_time = 0, publication = 0;
    const auto profile = [&](const std::string &outcome) {
        rusage self{}, child{};
        getrusage(RUSAGE_SELF, &self);
        getrusage(RUSAGE_CHILDREN, &child);
        auto seconds = [](timeval t) { return t.tv_sec + t.tv_usec / 1e6; };
        std::ostringstream o;
        o << "stage\tstatus\tvalidation_seconds\texecution_seconds\tsource_validation_"
             "seconds\tpublication_hash_seconds\twall_seconds\tuser_cpu_seconds\tsystem_cpu_"
             "seconds\tprocess_peak_rss_kb\tchildren_peak_rss_kb\tlocal_jobs\tlocal_cpus\tlocal_"
             "mem_mb\tlocal_job_mem_mb\n";
        o << dir.string() << '\t' << outcome << '\t' << validation << '\t' << execution << '\t'
          << sources_time << '\t' << publication << '\t' << elapsed(started) << '\t'
          << seconds(self.ru_utime) + seconds(child.ru_utime) - seconds(self_before.ru_utime) -
                 seconds(child_before.ru_utime)
          << '\t'
          << seconds(self.ru_stime) + seconds(child.ru_stime) - seconds(self_before.ru_stime) -
                 seconds(child_before.ru_stime)
          << '\t' << self.ru_maxrss << '\t' << child.ru_maxrss << '\t' << local_jobs << '\t'
          << local_cpus << '\t' << local_mem_mb << '\t' << local_job_mem_mb << '\n';
        fs::create_directories(run / "profiles");
        atomic_write(run / "profiles" /
                         (dir.filename().string() + "-" + std::to_string(getpid()) + "-" +
                          std::to_string(started.time_since_epoch().count()) + ".tsv"),
                     o.str());
    };
    try {
        Lock lock(dir.string() + ".lock");
        auto start = Clock::now();
        const bool reusable = valid(dir, id);
        validation = elapsed(start);
        if (reusable) {
            std::cout << "reuse " << dir << '\n';
            profile("reused");
            return;
        }
        quarantine(dir);
        std::cout << "run " << dir << '\n';
        start = Clock::now();
        try {
            action();
        } catch (...) {
            execution = elapsed(start);
            throw;
        }
        execution = elapsed(start);
        start = Clock::now();
        verify_sources();
        sources_time = elapsed(start);
        write(dir / "generation.txt",
              std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + "\n");
        start = Clock::now();
        atomic_write(dir / "STAGE.tsv", id + inventory(dir));
        publication = elapsed(start);
        profile("executed");
    } catch (...) {
        try {
            profile("failed");
        } catch (...) {
        }
        throw;
    }
}

std::string Workflow::cache_key() const {
    // SHA256 of complete index identity is produced via the ordinary file hasher.
    std::string s = "genome-only-v1\n" + c.at("backend") + "\n" + sha256_file(fasta) + "\n" +
                    sha256_file(gtf) + "\n";
    for (const auto &k : {"threads", "star_sa_bases", "star_chr_bits", "tool_lock"})
        s += c.at(k) + "\n";
    s += sha256_file(c.at("tool_lock")) + "\n" + sha256_file(exe) + "\n";
    std::vector<fs::path> tools;
    for (const auto &e : fs::directory_iterator(c.at("bin_dir")))
        if (e.is_regular_file())
            tools.push_back(e.path());
    std::sort(tools.begin(), tools.end());
    for (const auto &tool : tools)
        s += tool.filename().string() + "\t" + sha256_file(tool) + "\n";
    return s;
}

fs::path Workflow::index_path() {
    if (!c.count("index_cache"))
        return run / "index";
    const auto key = cache_key();
    const auto temp = run / ("cache-key-" + std::to_string(getpid()));
    write(temp, key);
    auto hash = sha256_file(temp);
    fs::remove(temp);
    return fs::path(c.at("index_cache")) / hash;
}

std::string Workflow::index_id() {
    return c.count("index_cache") ? "cache-v1\n" + cache_key() : base_id() + "index\n";
}

} // namespace rnaseq::workflow_detail
