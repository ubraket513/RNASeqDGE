#include "internal.hpp"
#include "rnaseq/counts.hpp"
#include "rnaseq/file_hash.hpp"
#include "rnaseq/manifests.hpp"
#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sched.h>
#include <set>
#include <sstream>
#include <unistd.h>

namespace rnaseq::workflow_detail {
void clean(const std::string &s) {
    if (s.find_first_of("\t\r\n") != std::string::npos)
        fail("control character in path/value");
}
fs::path resolve(const fs::path &base, const std::string &s) {
    clean(s);
    return fs::weakly_canonical(base / fs::path(s));
}
int number(const std::string &s, const std::string &name) {
    int n = 0;
    auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), n);
    if (ec != std::errc() || end != s.data() + s.size() || n < 1)
        fail(name + " must be a positive integer");
    return n;
}
std::size_t col(const Table &t, const std::string &name) {
    auto i = std::find(t.header.begin(), t.header.end(), name);
    if (i == t.header.end())
        fail("missing column " + name);
    return i - t.header.begin();
}
std::string table_text(const Table &t) {
    std::ostringstream o;
    auto row = [&](const auto &r) {
        bool first = true;
        for (const auto &v : r) {
            if (!first)
                o << '\t';
            first = false;
            clean(v);
            o << v;
        }
        o << '\n';
    };
    row(t.header);
    for (const auto &r : t.rows)
        row(r);
    return o.str();
}
// Validation order is intentional: retain the first reported error when several
// settings are invalid, before capturing the immutable source inventory.
Workflow::Workflow(const fs::path &path, const std::string &command) {
    config = fs::canonical(path);
    const auto base = config.parent_path();
    auto t = read_table(config.string());
    if (t.header != std::vector<std::string>{"key", "value"})
        fail("config header must be key/value");
    const std::set<std::string> allowed{
        "samples",    "runs",         "references",    "analysis",        "contrasts",
        "genes",      "annotation",   "bin_dir",       "rscript",         "r_script",
        "tool_lock",  "r_lock",       "run_dir",       "index_cache",     "backend",
        "threads",    "workers",      "star_sa_bases", "star_chr_bits",   "slurm_bin_dir",
        "slurm_cpus", "slurm_mem_mb", "slurm_time",    "slurm_partition", "slurm_concurrency",
        "local_jobs", "local_cpus",   "local_mem_mb",  "local_job_mem_mb"};
    for (const auto &r : t.rows) {
        if (!allowed.count(r[0]) || r[1].empty() || !c.emplace(r[0], r[1]).second)
            fail("unknown, empty or duplicate config key: " + r[0]);
        clean(r[1]);
    }
    for (const auto &k : {"samples", "runs", "references", "analysis", "contrasts", "genes",
                          "bin_dir", "rscript", "r_script", "tool_lock", "r_lock", "run_dir"})
        if (!c.count(k))
            fail(std::string("required config key: ") + k);
    for (const auto &[k, v] : Settings{{"backend", "hisat2"},
                                       {"threads", "1"},
                                       {"workers", "1"},
                                       {"star_sa_bases", "14"},
                                       {"star_chr_bits", "18"},
                                       {"slurm_concurrency", "1"}})
        if (!c.count(k))
            c[k] = v;
    if (c.at("backend") != "hisat2" && c.at("backend") != "star")
        fail("backend must be hisat2 or star");
    threads = number(c.at("threads"), "threads");
    workers = number(c.at("workers"), "workers");
    number(c.at("slurm_concurrency"), "slurm_concurrency");
    if (number(c.at("star_sa_bases"), "star_sa_bases") > 14 ||
        number(c.at("star_chr_bits"), "star_chr_bits") > 18)
        fail("STAR sizing exceeds allowed bounds");
    slurm_resources =
        command == "workflow-submit" || (command == "workflow-plan" && c.count("slurm_cpus"));
    // Local settings still have a strict schema even when submitting, but
    // the login host's affinity/allocation does not constrain remote jobs.
    for (const auto &key : {"local_jobs", "local_cpus", "local_mem_mb", "local_job_mem_mb"})
        if (c.count(key))
            number(c.at(key), key);
    if (slurm_resources) {
        if (!c.count("slurm_cpus"))
            fail("submission requires slurm_cpus");
        local_cpus = number(c.at("slurm_cpus"), "slurm_cpus");
        local_mem_mb = c.count("slurm_mem_mb") ? number(c.at("slurm_mem_mb"), "slurm_mem_mb") : 0;
    } else {
        if (const char *budget = getenv("SLURM_CPUS_PER_TASK"))
            if (std::max(threads, workers) > number(budget, "SLURM_CPUS_PER_TASK"))
                fail("threads/workers exceeds SLURM_CPUS_PER_TASK");
        cpu_set_t affinity;
        CPU_ZERO(&affinity);
        local_cpus =
            sched_getaffinity(0, sizeof(affinity), &affinity) == 0 ? CPU_COUNT(&affinity) : 1;
        if (c.count("local_cpus"))
            local_cpus = std::min(local_cpus, number(c.at("local_cpus"), "local_cpus"));
        if (const char *budget = getenv("SLURM_CPUS_PER_TASK"))
            local_cpus = std::min(local_cpus, number(budget, "SLURM_CPUS_PER_TASK"));
        if (std::max(threads, workers) > local_cpus)
            fail("threads/workers exceeds resolved local CPU budget");
        local_jobs =
            c.count("local_jobs") ? number(c.at("local_jobs"), "local_jobs") : local_cpus / threads;
        local_jobs = std::min(local_jobs, local_cpus / threads);
        if (c.count("local_mem_mb"))
            local_mem_mb = number(c.at("local_mem_mb"), "local_mem_mb");
        if (c.count("local_job_mem_mb"))
            local_job_mem_mb = number(c.at("local_job_mem_mb"), "local_job_mem_mb");
        if (const char *memory = getenv("SLURM_MEM_PER_NODE")) {
            const int allocated = number(memory, "SLURM_MEM_PER_NODE");
            local_mem_mb = local_mem_mb ? std::min(local_mem_mb, allocated) : allocated;
        }
        if (const char *memory = getenv("SLURM_MEM_PER_CPU")) {
            const auto allocated = std::min<long long>(
                std::numeric_limits<int>::max(),
                static_cast<long long>(number(memory, "SLURM_MEM_PER_CPU")) * local_cpus);
            local_mem_mb = local_mem_mb ? std::min(local_mem_mb, static_cast<int>(allocated))
                                        : static_cast<int>(allocated);
        }
        memory_fallback = !local_mem_mb || !local_job_mem_mb;
        if (memory_fallback)
            local_jobs = 1;
        else {
            if (local_job_mem_mb > local_mem_mb)
                fail("local_job_mem_mb exceeds resolved memory budget");
            local_jobs = std::min(local_jobs, local_mem_mb / local_job_mem_mb);
        }
    }
    if (c.count("slurm_cpus") &&
        number(c.at("slurm_cpus"), "slurm_cpus") < std::max(threads, workers))
        fail("slurm_cpus smaller than threads/workers");
    if (c.count("slurm_mem_mb"))
        number(c.at("slurm_mem_mb"), "slurm_mem_mb");
    if (c.count("slurm_time") &&
        c.at("slurm_time").find_first_not_of("0123456789:-") != std::string::npos)
        fail("slurm_time must use numeric Slurm time syntax");
    for (const auto &k : {"samples", "runs", "references", "analysis", "contrasts", "genes",
                          "annotation", "bin_dir", "rscript", "r_script", "tool_lock", "r_lock",
                          "run_dir", "index_cache", "slurm_bin_dir"})
        if (c.count(k))
            c[k] = resolve(base, c.at(k)).string();
    run = c.at("run_dir");
    exe = fs::canonical("/proc/self/exe");
    if (!fs::is_directory(run.parent_path()))
        fail("run_dir parent must exist");
    if (c.count("index_cache") && !fs::is_directory(fs::path(c.at("index_cache")).parent_path()))
        fail("index_cache parent must exist");
    if (c.at("bin_dir").find(':') != std::string::npos)
        fail("bin_dir contains PATH separator");
    validate_bundle(c.at("samples"), c.at("runs"), c.at("references"), c.at("analysis"),
                    c.at("contrasts"));
    runs = read_table(c.at("runs"));
    refs = read_table(c.at("references"));
    for (auto &r : runs.rows)
        for (const auto &key : {"fastq_1", "fastq_2"}) {
            auto &v = r[col(runs, key)];
            if (v.empty())
                continue;
            v = resolve(fs::path(c.at("runs")).parent_path(), v).string();
            if (fs::path(v).extension() == ".gz")
                fail("P4 requires uncompressed FASTQ");
            sources.emplace_back(v);
        }
    for (auto &r : refs.rows) {
        auto &v = r[col(refs, "path")];
        v = resolve(fs::path(c.at("references")).parent_path(), v).string();
        sources.emplace_back(v);
        if (r[col(refs, "role")] == "genome")
            fasta = v;
        else
            gtf = v;
    }
    auto genes = read_table(c.at("genes"));
    if (genes.header != std::vector<std::string>{"gene_id"} || genes.rows.empty())
        fail("genes must be nonempty gene_id table");
    std::set<std::string> ids;
    for (const auto &r : genes.rows)
        if (r[0].empty() || !ids.insert(r[0]).second)
            fail("empty/duplicate gene_id");
    if (c.count("annotation"))
        validate_annotation(c.at("genes"), c.at("annotation"));
    sources.push_back(config);
    sources.push_back(exe);
    for (const auto &k : {"samples", "runs", "references", "analysis", "contrasts", "genes",
                          "annotation", "rscript", "r_script", "tool_lock", "r_lock"})
        if (c.count(k))
            sources.emplace_back(c.at(k));
    for (const auto &k :
         {"samples", "analysis", "contrasts", "genes", "annotation", "tool_lock", "r_lock"})
        if (c.count(k))
            snapshots[std::string(k) + ".tsv"] = read(c.at(k));
    snapshots["r-preflight.R"] = "stopifnot(getRversion() == '4.5.3')\npins <- c(DESeq2='1.50.2', "
                                 "apeglm='1.32.0', BiocParallel='1.44.0', pheatmap='1.0.13')\nfor "
                                 "(p in names(pins)) stopifnot(as.character(packageVersion(p)) == "
                                 "pins[[p]])\ncat('Pinned R package versions verified\\n')\n";
    snapshots["runs.tsv"] = table_text(runs);
    snapshots["references.tsv"] = table_text(refs);
    std::string resolved = "key\tvalue\n";
    for (const auto &[k, v] : c)
        resolved += k + "\t" + v + "\n";
    snapshots["config.tsv"] = resolved;
    const auto bin = fs::path(c.at("bin_dir"));
    for (const auto &name : std::vector<std::string>{
             c.at("backend") == "star" ? "STAR" : "hisat2-align-s",
             c.at("backend") == "star" ? "STAR" : "hisat2-build-s", "samtools", "featureCounts"})
        if (access((bin / name).c_str(), X_OK) != 0)
            fail("missing executable " + (bin / name).string());
    if (access(c.at("rscript").c_str(), X_OK) != 0)
        fail("rscript not executable");
    // Include wrapper companions/interpreters in the selected tool directory.
    for (const auto &e : fs::directory_iterator(bin))
        if (e.is_regular_file())
            sources.push_back(e.path());
    std::sort(sources.begin(), sources.end());
    sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
    source_identity = "sha256\tpath\n";
    for (const auto &p : sources) {
        if (!fs::is_regular_file(p) || access(p.c_str(), R_OK) != 0)
            fail("unreadable input " + p.string());
        source_identity += sha256_file(p) + "\t" + p.string() + "\n";
    }
}

void Workflow::plan() const {
    std::cout << "DAG: index -> alignment/counting[" << runs.rows.size()
              << "] -> sample merge -> offline R\n";
    for (const auto &[k, v] : c)
        std::cout << k << '\t' << v << '\n';
    std::cout << "reference_fasta\t" << fasta << "\nreference_gtf\t" << gtf << "\nBLAS/OpenMP: 1\n";
    if (slurm_resources) {
        std::cout << "resource_scope\tslurm_requested\nresolved_task_cpus\t" << local_cpus
                  << "\nresolved_task_mem_mb\t" << local_mem_mb << '\n';
    } else {
        std::cout << "resource_scope\tlocal_allocation\n";
        std::cout << "local_resolved_cpus\t" << local_cpus << "\nlocal_resolved_mem_mb\t"
                  << local_mem_mb << "\nlocal_job_mem_mb\t" << local_job_mem_mb
                  << "\nlocal_resolved_jobs\t" << local_jobs << "\nlocal_memory_policy\t"
                  << (memory_fallback ? "serial fallback: total/per-job memory estimate missing"
                                      : "reservation estimates; not a hard RSS limit")
                  << '\n';
    }
    std::cout << table_text(runs);
}

} // namespace rnaseq::workflow_detail
