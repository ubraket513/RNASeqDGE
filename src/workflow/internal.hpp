#pragma once

#include "rnaseq/table.hpp"
#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <vector>

// Private workflow implementation. Only workflow.hpp is a public interface.
namespace rnaseq::workflow_detail {
namespace fs = std::filesystem;
using Settings = std::map<std::string, std::string>;

void fail(const std::string &message);
std::string read(const fs::path &path);
void write(const fs::path &path, const std::string &contents);
void atomic_write(const fs::path &path, const std::string &contents);
void clean(const std::string &value);
std::size_t col(const Table &table, const std::string &name);
std::string table_text(const Table &table);
void quarantine(const fs::path &path);
bool valid(const fs::path &directory, const std::string &identity);

// Keep each lock alive throughout verification and consumption/publication.
class Lock {
  public:
    explicit Lock(const fs::path &path, bool shared = false, bool wait = false);
    ~Lock();
    Lock(const Lock &) = delete;

  private:
    int fd = -1;
};

struct Workflow {
    Settings c;
    fs::path config, run, exe, fasta, gtf;
    Table runs, refs;
    std::string source_identity;
    std::map<std::string, std::string> snapshots;
    std::vector<fs::path> sources;
    int threads = 1, workers = 1, local_jobs = 1, local_cpus = 1, local_mem_mb = 0,
        local_job_mem_mb = 0;
    bool memory_fallback = true, slurm_resources = false;

    // config.cpp: validate and resolve inputs and resource budgets.
    explicit Workflow(const fs::path &path, const std::string &command);
    void plan() const;

    // state.cpp: publish or verify the immutable run snapshot.
    void snapshot(bool create);

    // execution.cpp: coordinate the existing native and R stage implementations.
    void index();
    void alignment(std::size_t i);
    void finish();
    void submit();

  private:
    // State helpers retain locks through validation, execution and publication.
    void verify_sources() const;
    std::string base_id() const;
    void stage(const fs::path &dir, const std::string &id, const std::function<void()> &action);
    std::string cache_key() const;
    fs::path index_path();
    std::string index_id();
    void align(const std::string &cmd, const fs::path &out, const std::vector<std::string> &extra);
};
} // namespace rnaseq::workflow_detail
