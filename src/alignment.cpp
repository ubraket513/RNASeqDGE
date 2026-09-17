#include "rnaseq/alignment.hpp"
#include "rnaseq/file_hash.hpp"
#include "rnaseq/process.hpp"
#include <cerrno>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <linux/fs.h>
#include <map>
#include <set>
#include <sstream>
#include <sys/syscall.h>
#include <unistd.h>
namespace rnaseq {
namespace {
namespace fs = std::filesystem;
using Options = std::map<std::string, std::string>;
int positive(const std::string &value, const std::string &name) {
    int result = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc() || end != value.data() + value.size() || result < 1)
        throw std::runtime_error(name + " must be a positive integer");
    return result;
}
void readable(const fs::path &file) {
    if (!fs::is_regular_file(file) || access(file.c_str(), R_OK) != 0 || fs::file_size(file) == 0)
        throw std::runtime_error("expected nonempty readable regular file: " + file.string());
}
std::string contents(const fs::path &path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot read " + path.string());
    return {std::istreambuf_iterator<char>(input), {}};
}
void write(const fs::path &path, const std::string &data) {
    std::ofstream stream(path);
    stream << data;
    stream.close();
    if (!stream)
        throw std::runtime_error("cannot write " + path.string());
}
bool pinned_version(const std::string &executable, const std::string &output,
                    const std::string &expected) {
    std::istringstream lines(output);
    std::string line;
    while (std::getline(lines, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.find_first_not_of(" \t") == std::string::npos)
            continue;
        if (executable == "STAR")
            return line == expected;
        if (executable == "samtools")
            return line == "samtools " + expected;
        if (executable == "featureCounts")
            return line == "featureCounts v" + expected;
        const auto delimiter = line.rfind(" version ");
        if (delimiter == std::string::npos || line.substr(delimiter + 9) != expected)
            return false;
        // HISAT2 prints either its wrapper name or a full native executable path.
        const auto name = fs::path(line.substr(0, delimiter)).filename().string();
        if (executable == "hisat2-build-s")
            return name == "hisat2-build" || name == "hisat2-build-s" || name == "hisat2-build-l";
        return name == "hisat2" || name == "hisat2-align-s" || name == "hisat2-align-l";
    }
    return false;
}
std::vector<std::string> index_files(const std::string &backend) {
    if (backend == "star")
        return {"Genome",       "SA",
                "SAindex",      "genomeParameters.txt",
                "chrName.txt",  "chrLength.txt",
                "chrStart.txt", "chrNameLength.txt"};
    std::vector<std::string> result;
    for (int i = 1; i <= 8; ++i)
        result.push_back("genome." + std::to_string(i) + ".ht2");
    return result;
}
} // namespace
void alignment_command(int argc, char **argv) {
    const std::string command = argv[1];
    const bool indexing = command == "index";
    std::set<std::string> allowed{"--backend", "--bin-dir", "--fasta",
                                  "--gtf",     "--threads", "--output"};
    if (indexing)
        allowed.insert({"--star-sa-bases", "--star-chr-bits"});
    else
        allowed.insert({"--index", "--reads1", "--reads2", "--layout", "--strandedness"});
    Options options;
    for (int i = 2; i < argc; i += 2) {
        if (i + 1 == argc || !allowed.count(argv[i]) ||
            !options.emplace(argv[i], argv[i + 1]).second)
            throw std::runtime_error("unknown, duplicate, or incomplete option: " +
                                     std::string(argv[i]));
        if (std::string(argv[i + 1]).find_first_of("\t\r\n") != std::string::npos)
            throw std::runtime_error("control character in option");
    }
    auto required = [&](const std::string &key) -> std::string {
        if (!options.count(key) || options[key].empty())
            throw std::runtime_error("required option: " + key);
        return options[key];
    };
    const auto backend = required("--backend");
    if (backend != "star" && backend != "hisat2")
        throw std::runtime_error("backend must be star or hisat2");
    const int threads =
        positive(options.count("--threads") ? options.at("--threads") : "1", "threads");
    if (const char *allocation = std::getenv("SLURM_CPUS_PER_TASK"))
        if (threads > positive(allocation, "SLURM_CPUS_PER_TASK"))
            throw std::runtime_error("threads exceeds SLURM_CPUS_PER_TASK");
    const fs::path bin = fs::absolute(required("--bin-dir"));
    const fs::path fasta = fs::absolute(required("--fasta")), gtf = fs::absolute(required("--gtf"));
    const fs::path output = fs::absolute(required("--output")).lexically_normal();
    readable(fasta);
    readable(gtf);
    if (fs::exists(fs::symlink_status(output)) || !fs::is_directory(output.parent_path()))
        throw std::runtime_error("output must be absent with an existing parent");
    std::vector<std::string> executables{backend == "star" ? "STAR"
                                         : indexing        ? "hisat2-build-s"
                                                           : "hisat2-align-s"};
    if (!indexing) {
        executables.push_back("samtools");
        executables.push_back("featureCounts");
    }
    for (const auto &executable : executables)
        if (!fs::is_regular_file(bin / executable) || access((bin / executable).c_str(), X_OK))
            throw std::runtime_error("missing executable: " + (bin / executable).string());
    int sa = 14, chr = 18;
    fs::path index, reads1, reads2;
    std::string layout, strand;
    if (indexing) {
        if (backend != "star" &&
            (options.count("--star-sa-bases") || options.count("--star-chr-bits")))
            throw std::runtime_error("STAR parameters require star backend");
        if (options.count("--star-sa-bases"))
            sa = positive(options.at("--star-sa-bases"), "star-sa-bases");
        if (options.count("--star-chr-bits"))
            chr = positive(options.at("--star-chr-bits"), "star-chr-bits");
        if (sa > 14 || chr > 18)
            throw std::runtime_error("STAR sa-bases must be <=14 and chr-bits <=18");
    } else {
        index = fs::absolute(required("--index"));
        reads1 = fs::absolute(required("--reads1"));
        layout = required("--layout");
        strand = required("--strandedness");
        if (layout != "single" && layout != "paired")
            throw std::runtime_error("layout must be single or paired");
        if (strand != "unstranded" && strand != "forward" && strand != "reverse")
            throw std::runtime_error("invalid strandedness");
        readable(reads1);
        if (reads1.extension() == ".gz")
            throw std::runtime_error("P4 requires uncompressed FASTQ");
        if (layout == "paired") {
            reads2 = fs::absolute(required("--reads2"));
            readable(reads2);
            if (reads2.extension() == ".gz" || fs::equivalent(reads1, reads2))
                throw std::runtime_error("paired reads must be distinct uncompressed FASTQ files");
        } else if (options.count("--reads2"))
            throw std::runtime_error("single layout cannot have reads2");
        for (const auto &file : index_files(backend))
            readable(index / file);
    }
    const std::string identity = "backend\t" + backend + "\npolicy\tgenome-only\nfasta_sha256\t" +
                                 sha256_file(fasta) + "\ngtf_sha256\t" + sha256_file(gtf) + "\n";
    if (!indexing && contents(index / "index.tsv").find(identity) != 0)
        throw std::runtime_error("index reference/backend identity mismatch");
    auto pattern =
        (output.parent_path() / (output.filename().string() + ".staging-XXXXXX")).string();
    std::vector<char> buffer(pattern.begin(), pattern.end());
    buffer.push_back(0);
    const char *created = mkdtemp(buffer.data());
    if (!created)
        throw std::runtime_error("cannot create staging directory");
    const fs::path stage = created;
    auto run = [&](const std::string &label, std::vector<std::string> args) {
        std::ofstream record(stage / "commands.tsv", std::ios::app);
        record << label;
        for (const auto &arg : args)
            record << '\t' << arg;
        record << '\n';
        record.close();
        if (!record)
            throw std::runtime_error("cannot record command");
        run_process(args, stage / (label + ".log"));
    };
    try {
        for (const auto &executable : executables) {
            run("version-" + executable,
                {(bin / executable).string(), executable == "featureCounts" ? "-v" : "--version"});
            const auto version = contents(stage / ("version-" + executable + ".log"));
            const std::string expected = executable == "STAR"            ? "2.7.10b"
                                         : executable == "samtools"      ? "1.17"
                                         : executable == "featureCounts" ? "2.0.6"
                                                                         : "2.2.1";
            if (!pinned_version(executable, version, expected))
                throw std::runtime_error(executable + " expected pinned version " + expected);
        }
        const auto t = std::to_string(threads);
        if (indexing) {
            if (backend == "star")
                run("index",
                    {(bin / "STAR").string(), "--runMode", "genomeGenerate", "--runThreadN", t,
                     "--genomeDir", stage.string(), "--genomeFastaFiles", fasta.string(),
                     "--genomeSAindexNbases", std::to_string(sa), "--genomeChrBinNbits",
                     std::to_string(chr), "--outFileNamePrefix", (stage / "star-").string()});
            else
                run("index", {(bin / "hisat2-build-s").string(), "--wrapper", "basic-0", "-p", t,
                              fasta.string(), (stage / "genome").string()});
            for (const auto &file : index_files(backend))
                readable(stage / file);
            write(stage / "index.tsv", identity + "threads\t" + t + "\nstar_sa_bases\t" +
                                           std::to_string(sa) + "\nstar_chr_bits\t" +
                                           std::to_string(chr) + "\n");
        } else {
            std::vector<std::string> args;
            fs::path unsorted;
            if (backend == "star") {
                unsorted = stage / "star-Aligned.out.bam";
                args = {(bin / "STAR").string(), "--runThreadN", t, "--genomeDir", index.string(),
                        "--readFilesIn",         reads1.string()};
                if (layout == "paired")
                    args.push_back(reads2.string());
                args.insert(args.end(), {"--outFileNamePrefix", (stage / "star-").string(),
                                         "--outSAMtype", "BAM", "Unsorted"});
            } else {
                unsorted = stage / "aligned.sam";
                args = {(bin / "hisat2-align-s").string(), "--wrapper", "basic-0", "-p", t, "-x",
                        (index / "genome").string()};
                if (layout == "paired")
                    args.insert(args.end(), {"-1", reads1.string(), "-2", reads2.string()});
                else
                    args.insert(args.end(), {"-U", reads1.string()});
                if (strand != "unstranded")
                    args.insert(args.end(),
                                {"--rna-strandness", layout == "paired"
                                                         ? (strand == "forward" ? "FR" : "RF")
                                                         : (strand == "forward" ? "F" : "R")});
                args.insert(args.end(), {"-S", unsorted.string()});
            }
            run("align", args);
            readable(unsorted);
            const auto bam = (stage / "aligned.bam").string();
            run("sort",
                {(bin / "samtools").string(), "sort", "-@", std::to_string(threads - 1), "-m",
                 "256M", "-T", (stage / "sort-tmp").string(), "-o", bam, unsorted.string()});
            run("quickcheck", {(bin / "samtools").string(), "quickcheck", bam});
            args = {(bin / "featureCounts").string(),
                    "-T",
                    t,
                    "-s",
                    strand == "unstranded" ? "0"
                    : strand == "forward"  ? "1"
                                           : "2",
                    "-Q",
                    "0",
                    "-t",
                    "exon",
                    "-g",
                    "gene_id"};
            if (layout == "paired")
                args.insert(args.end(), {"-p", "--countReadPairs"});
            args.insert(args.end(),
                        {"-a", gtf.string(), "-o", (stage / "counts.txt").string(), bam});
            run("count", args);
            readable(stage / "counts.txt");
            readable(stage / "counts.txt.summary");
            write(stage / "run.tsv", identity + "threads\t" + t + "\nlayout\t" + layout +
                                         "\nstrandedness\t" + strand + "\nreads1_sha256\t" +
                                         sha256_file(reads1) + "\nreads2_sha256\t" +
                                         (layout == "paired" ? sha256_file(reads2) : "") +
                                         "\nsort_memory_per_thread\t256M\n");
        }
        write(stage / "COMPLETE", "success\n");
        if (fs::exists(fs::symlink_status(output)))
            throw std::runtime_error("output appeared during execution");
        // Linux no-replace rename closes the absent-output race without a P5 lock.
        if (syscall(SYS_renameat2, AT_FDCWD, stage.c_str(), AT_FDCWD, output.c_str(),
                    RENAME_NOREPLACE) != 0)
            throw std::runtime_error(std::string("cannot publish output: ") + std::strerror(errno));
    } catch (const ProcessError &error) {
        std::error_code ignored;
        fs::remove(stage / "COMPLETE", ignored);
        throw ProcessError(std::string(error.what()) + "; preserved stage: " + stage.string(),
                           error.status);
    } catch (const std::exception &error) {
        std::error_code ignored;
        fs::remove(stage / "COMPLETE", ignored);
        throw std::runtime_error(std::string(error.what()) +
                                 "; preserved stage: " + stage.string());
    }
}
} // namespace rnaseq
