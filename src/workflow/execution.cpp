#include "internal.hpp"
#include "rnaseq/alignment.hpp"
#include "rnaseq/counts.hpp"
#include "rnaseq/file_hash.hpp"
#include "rnaseq/process.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <unistd.h>

namespace rnaseq::workflow_detail {
namespace {
std::string shell_quote(const std::string &value) {
    std::string s = "'";
    for (char c : value)
        s += c == '\'' ? "'\\''" : std::string(1, c);
    return s + "'";
}
} // namespace

void Workflow::align(const std::string &cmd, const fs::path &out,
                     const std::vector<std::string> &extra) {
    std::vector<std::string> args{"rnaseq",    cmd,
                                  "--backend", c.at("backend"),
                                  "--bin-dir", c.at("bin_dir"),
                                  "--fasta",   fasta.string(),
                                  "--gtf",     gtf.string(),
                                  "--threads", c.at("threads"),
                                  "--output",  out.string()};
    args.insert(args.end(), extra.begin(), extra.end());
    std::vector<char *> ptrs;
    for (auto &a : args)
        ptrs.push_back(a.data());
    alignment_command(ptrs.size(), ptrs.data());
}

void Workflow::index() {
    if (c.count("index_cache"))
        fs::create_directories(c.at("index_cache"));
    const auto dir = index_path();
    const auto id = c.count("index_cache") ? "cache-v1\n" + cache_key() : base_id() + "index\n";
    stage(dir, id, [&] {
        std::vector<std::string> extra;
        if (c.at("backend") == "star")
            extra = {"--star-sa-bases", c.at("star_sa_bases"), "--star-chr-bits",
                     c.at("star_chr_bits")};
        align("index", dir, extra);
    });
}

void Workflow::alignment(std::size_t i) {
    if (i >= runs.rows.size())
        fail("array task outside immutable run mapping");
    const auto idx = index_path();
    // Shared cache lock holds the checked index stable throughout alignment.
    Lock index_lock(idx.string() + ".lock", true);
    if (!valid(idx, index_id()))
        fail("index incomplete or corrupt; rerun index stage");
    const auto &row = runs.rows[i];
    const auto name = row[col(runs, "run_id")];
    const auto dir = run / "align" / name;
    const auto id = base_id() + "align\t" + name + "\n" + sha256_file(idx / "STAGE.tsv") + "\n";
    stage(dir, id, [&] {
        std::vector<std::string> extra{"--index",        idx.string(),
                                       "--reads1",       row[col(runs, "fastq_1")],
                                       "--layout",       row[col(runs, "layout")],
                                       "--strandedness", row[col(runs, "strandedness")]};
        if (!row[col(runs, "fastq_2")].empty())
            extra.insert(extra.end(), {"--reads2", row[col(runs, "fastq_2")]});
        align("align-count", dir, extra);
    });
}

void Workflow::finish() {
    const auto idx = index_path();
    Lock index_lock(idx.string() + ".lock", true);
    if (!valid(idx, index_id()))
        fail("invalid index before merge");
    std::vector<std::unique_ptr<Lock>> alignment_locks;
    std::string merge_id = base_id() + "merge\n", inputs = "run_id\tformat\tpath\tcolumn\n";
    for (const auto &row : runs.rows) {
        const auto name = row[col(runs, "run_id")];
        const auto dir = run / "align" / name;
        alignment_locks.push_back(std::make_unique<Lock>(dir.string() + ".lock", true));
        const auto id = base_id() + "align\t" + name + "\n" + sha256_file(idx / "STAGE.tsv") + "\n";
        if (!valid(dir, id))
            fail("incomplete or corrupt alignment " + name);
        merge_id += sha256_file(dir / "STAGE.tsv") + "\n";
        std::ifstream f(dir / "counts.txt");
        TsvReader reader(f, (dir / "counts.txt").string(), true);
        std::vector<std::string> header;
        if (!reader.next(header) || header.size() != 7 || header[0] != "Geneid")
            fail("featureCounts header must contain exactly one BAM column");
        inputs +=
            name + "\tfeaturecounts\t" + (dir / "counts.txt").string() + "\t" + header[6] + "\n";
    }
    const auto snap = run / "snapshot", merged = run / "merge";
    stage(merged, merge_id, [&] {
        fs::create_directory(merged);
        write(merged / "inputs.tsv", inputs);
        merge_counts(snap / "samples.tsv", snap / "runs.tsv", snap / "genes.tsv",
                     merged / "inputs.tsv", merged / "counts.tsv",
                     c.count("annotation") ? snap / "annotation.tsv" : fs::path{});
    });
    const auto out = run / "analysis";
    stage(out, base_id() + "R\n" + sha256_file(merged / "STAGE.tsv") + "\n", [&] {
        // R owns its own atomic publication; wrapper logs live outside its destination.
        std::vector<std::string> args{c.at("rscript"),
                                      "--vanilla",
                                      c.at("r_script"),
                                      "--counts",
                                      (merged / "counts.tsv").string(),
                                      "--samples",
                                      (snap / "samples.tsv").string(),
                                      "--analysis",
                                      (snap / "analysis.tsv").string(),
                                      "--contrasts",
                                      (snap / "contrasts.tsv").string(),
                                      "--out",
                                      out.string(),
                                      "--workers",
                                      c.at("workers")};
        if (c.count("annotation"))
            args.insert(args.end(), {"--annotation", (snap / "annotation.tsv").string()});
        for (const auto &name : {"R_HOME", "R_LIBS", "LD_LIBRARY_PATH", "LD_PRELOAD"})
            unsetenv(name);
        for (const auto &name : {"R_LIBS_USER", "R_LIBS_SITE"})
            setenv(name, "/dev/null", 1);
        const auto preflight_log = run / ("R-preflight-" + std::to_string(getpid()) + ".log");
        quarantine(preflight_log);
        run_process({c.at("rscript"), "--vanilla", (snap / "r-preflight.R").string()},
                    preflight_log);
        const auto log = run / ("R-" + std::to_string(getpid()) + ".log");
        quarantine(log);
        run_process(args, log);
        if (!fs::is_directory(out))
            fail("R returned success without output directory");
        fs::copy_file(log, out / "workflow-R.log");
        fs::copy_file(preflight_log, out / "workflow-R-preflight.log");
    });
}

void Workflow::submit() {
    if (getenv("SLURM_JOB_ID"))
        fail("nested submission is forbidden");
    for (const auto &k : {"slurm_bin_dir", "slurm_cpus", "slurm_mem_mb", "slurm_time"})
        if (!c.count(k))
            fail(std::string("submission requires ") + k);
    const auto bin = fs::path(c.at("slurm_bin_dir"));
    for (const auto &name : {"sbatch", "scancel"})
        if (access((bin / name).c_str(), X_OK) != 0)
            fail(std::string("missing scheduler executable ") + name);
    if (fs::exists(run / "submission"))
        fail(
            "submission already attempted; inspect recorded jobs before using a new run directory");
    fs::create_directory(run / "submission");
    const auto dir = run / "submission";
    for (const auto &task : {"index", "array", "finish"})
        write(dir / (std::string(task) + ".sh"),
              "#!/bin/sh\nset -eu\nexec " + shell_quote(exe.string()) + " workflow-task " +
                  shell_quote(config.string()) + " " + task + "\n");
    std::vector<std::string> accepted;
    try {
        for (const auto &task : {"index", "array", "finish"}) {
            std::vector<std::string> args{(bin / "sbatch").string(),
                                          "--parsable",
                                          "--cpus-per-task=" + c.at("slurm_cpus"),
                                          "--mem=" + c.at("slurm_mem_mb"),
                                          "--time=" + c.at("slurm_time"),
                                          "--output=" +
                                              (dir / (std::string(task) + "-%A_%a.log")).string()};
            if (c.count("slurm_partition"))
                args.push_back("--partition=" + c.at("slurm_partition"));
            if (!accepted.empty())
                args.push_back("--dependency=afterok:" + accepted.back());
            if (std::string(task) == "array")
                args.push_back("--array=0-" + std::to_string(runs.rows.size() - 1) + "%" +
                               c.at("slurm_concurrency"));
            args.push_back((dir / (std::string(task) + ".sh")).string());
            std::string record;
            for (const auto &a : args)
                record += a + "\n";
            write(dir / (std::string(task) + ".argv"), record);
            const auto log = dir / (std::string(task) + ".submit.stdout");
            run_process(args, dir / (std::string(task) + ".submit.stderr"), log);
            std::string id = read(log);
            while (!id.empty() && (id.back() == '\n' || id.back() == '\r'))
                id.pop_back();
            const auto semicolon = id.find(';');
            if (semicolon != std::string::npos)
                id = id.substr(0, semicolon);
            if (id.empty() || id.find_first_not_of("0123456789") != std::string::npos)
                fail("sbatch did not return a parsable numeric job ID; inspect " + log.string());
            accepted.push_back(id);
            std::string ids;
            for (const auto &job : accepted)
                ids += job + "\n";
            atomic_write(dir / "accepted.txt", ids);
            std::cout << "accepted " << task << " job " << id << '\n';
        }
        atomic_write(dir / "SUBMITTED", "jobs accepted; live execution is not yet verified\n");
    } catch (...) {
        for (const auto &id : accepted)
            try {
                run_process({(bin / "scancel").string(), id}, dir / ("cancel-" + id + ".log"));
            } catch (const std::exception &e) {
                std::cerr << "rollback cancellation failed: " << e.what() << '\n';
            }
        write(dir / "FAILED", "submission failed; inspect accepted IDs and cancellation logs\n");
        throw;
    }
}

} // namespace rnaseq::workflow_detail
