#include "rnaseq/counts.hpp"

#include "rnaseq/manifests.hpp"
#include "rnaseq/table.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <unistd.h>

namespace rnaseq {
namespace {
namespace fs = std::filesystem;

[[noreturn]] void fail(const std::string& source, std::size_t record,
                       const std::string& column, const std::string& message) {
    throw std::runtime_error(source + ": record " + std::to_string(record) +
                             ", column " + column + ": " + message);
}

[[noreturn]] void fail(const Table& table, std::size_t record,
                       const std::string& column, const std::string& message) {
    fail(table.source, record, column, message);
}

std::size_t column(const Table& table, const std::string& name) {
    const auto found = std::find(table.header.begin(), table.header.end(), name);
    if (found == table.header.end()) fail(table, 1, name, "required column missing");
    return static_cast<std::size_t>(found - table.header.begin());
}

void fixed_schema(const Table& table, const std::vector<std::string>& expected) {
    for (const auto& name : expected) column(table, name);
    const std::set<std::string> allowed(expected.begin(), expected.end());
    for (const auto& name : table.header)
        if (!allowed.contains(name)) fail(table, 1, name, "unexpected column");
}

void readable_regular(const fs::path& path, const std::string& source,
                      std::size_t record, const std::string& field) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error)
        fail(source, record, field, "path does not name a regular readable file");
    std::ifstream input(path, std::ios::binary);
    if (!input) fail(source, record, field, "path does not name a regular readable file");
}

Table read_regular_table(const fs::path& path) {
    readable_regular(path, path.string(), 1, "1");
    return read_table(path.string());
}

fs::path data_path(const Table& table, std::size_t record,
                   const std::string& field, const std::string& value) {
    if (value.empty()) fail(table, record, field, "path is empty");
    return fs::path(table.source).parent_path() / fs::path(value);
}

bool bytewise_less(const std::string& left, const std::string& right) {
    return std::lexicographical_compare(
        left.begin(), left.end(), right.begin(), right.end(),
        [](char a, char b) {
            return static_cast<unsigned char>(a) < static_cast<unsigned char>(b);
        });
}

struct GeneUniverse {
    std::vector<std::string> ids;
    std::map<std::string, std::size_t, decltype(&bytewise_less)> index{&bytewise_less};
};

GeneUniverse read_genes(const fs::path& path) {
    const auto genes = read_regular_table(path);
    fixed_schema(genes, {"gene_id"});
    if (genes.rows.empty()) fail(genes, 2, "gene_id", "no genes");

    GeneUniverse result;
    result.ids.reserve(genes.rows.size());
    std::set<std::string, decltype(&bytewise_less)> seen(&bytewise_less);
    for (std::size_t row = 0; row < genes.rows.size(); ++row) {
        const auto& id = genes.rows[row][0];
        if (id.empty()) fail(genes, row + 2, "gene_id", "gene ID is empty");
        if (!seen.insert(id).second)
            fail(genes, row + 2, "gene_id", "duplicate gene ID");
        result.ids.push_back(id);
    }
    std::sort(result.ids.begin(), result.ids.end(), bytewise_less);
    for (std::size_t i = 0; i < result.ids.size(); ++i)
        result.index.emplace(result.ids[i], i);
    return result;
}

void reject_output_alias(const fs::path& output, const fs::path& input) {
    std::error_code input_error;
    if (!fs::exists(input, input_error) || input_error) return;
    std::error_code output_error;
    if (!fs::exists(output, output_error) || output_error) return;
    std::error_code equivalent_error;
    const bool same = fs::equivalent(output, input, equivalent_error);
    if (equivalent_error)
        throw std::runtime_error("cannot compare output and input paths: " +
                                 equivalent_error.message());
    if (same)
        throw std::runtime_error("output aliases input file: " + input.string());
}

void reject_output_aliases(const fs::path& output,
                           const std::vector<fs::path>& inputs) {
    for (const auto& input : inputs) reject_output_alias(output, input);
}

std::ifstream open_regular(const fs::path& path, const std::string& source,
                           std::size_t record, const std::string& field) {
    readable_regular(path, source, record, field);
    std::ifstream input(path, std::ios::binary);
    if (!input) fail(source, record, field, "cannot open input");
    return input;
}

void require_header(const std::vector<std::string>& actual,
                    const std::vector<std::string>& expected,
                    const std::string& source) {
    if (actual != expected)
        fail(source, 1, "1", "unexpected header");
}

std::size_t unique_header_column(const std::vector<std::string>& header,
                                 const std::string& selected,
                                 const std::string& source) {
    std::set<std::string> seen;
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i].empty() || !seen.insert(header[i]).second)
            fail(source, 1, std::to_string(i + 1), "empty or duplicate header");
    }
    const auto found = std::find(header.begin(), header.end(), selected);
    if (found == header.end()) fail(source, 1, selected, "selected column missing");
    return static_cast<std::size_t>(found - header.begin());
}

void add_count(const std::vector<std::string>& row, std::size_t gene_column,
               std::size_t count_column, const std::string& count_name,
               const std::string& source, std::size_t record,
               const GeneUniverse& genes, std::vector<bool>& seen,
               std::vector<std::uint64_t>& totals, std::size_t sample) {
    const auto found = genes.index.find(row[gene_column]);
    if (row[gene_column].empty() || found == genes.index.end())
        fail(source, record, "gene_id", "gene is absent from reference universe");
    const auto gene = found->second;
    if (seen[gene]) fail(source, record, "gene_id", "duplicate gene ID");
    seen[gene] = true;
    const auto value = count_value(
        row[count_column], source + ": record " + std::to_string(record) +
                               ", column " + count_name);
    const auto offset = gene * (totals.size() / genes.ids.size()) + sample;
    try {
        totals[offset] = checked_add(totals[offset], value);
    } catch (const std::overflow_error&) {
        throw std::overflow_error(source + ": record " + std::to_string(record) +
                                  ", column " + count_name +
                                  ": count sum exceeds uint64 range");
    }
}

void require_complete_universe(const std::vector<bool>& seen,
                               const GeneUniverse& genes,
                               const std::string& source) {
    const auto missing = std::find(seen.begin(), seen.end(), false);
    if (missing != seen.end()) {
        const auto index = static_cast<std::size_t>(missing - seen.begin());
        fail(source, 1, "gene_id", "reference gene missing: " + genes.ids[index]);
    }
}

void read_native_counts(const fs::path& path, const GeneUniverse& genes,
                        std::vector<std::uint64_t>& totals, std::size_t sample,
                        const std::string& manifest_source, std::size_t manifest_record) {
    auto input = open_regular(path, manifest_source, manifest_record, "path");
    TsvReader reader(input, path.string());
    std::vector<std::string> row;
    if (!reader.next(row)) fail(path.string(), 1, "1", "empty file");
    require_header(row, {"gene_id", "count"}, path.string());
    std::vector<bool> seen(genes.ids.size());
    while (reader.next(row))
        add_count(row, 0, 1, "count", path.string(), reader.record(), genes,
                  seen, totals, sample);
    require_complete_universe(seen, genes, path.string());
}

void read_legacy_counts(const fs::path& path, const GeneUniverse& genes,
                        std::vector<std::uint64_t>& totals, std::size_t sample,
                        const std::string& manifest_source, std::size_t manifest_record) {
    auto input = open_regular(path, manifest_source, manifest_record, "path");
    TsvReader reader(input, path.string());
    std::vector<bool> seen(genes.ids.size());
    std::vector<std::string> row;
    bool any = false;
    while (reader.next(row)) {
        any = true;
        if (row.size() != 2)
            fail(path.string(), reader.record(), std::to_string(row.size()),
                 "legacy count row must have exactly two columns");
        add_count(row, 0, 1, "count", path.string(), reader.record(), genes,
                  seen, totals, sample);
    }
    if (!any) fail(path.string(), 1, "1", "empty file");
    require_complete_universe(seen, genes, path.string());
}

void validate_featurecounts_summary(const fs::path& path,
                                    const std::string& selected,
                                    const std::string& manifest_source,
                                    std::size_t manifest_record) {
    auto input = open_regular(path, manifest_source, manifest_record, "path");
    TsvReader reader(input, path.string());
    std::vector<std::string> row;
    if (!reader.next(row)) fail(path.string(), 1, "1", "empty summary");
    if (row.size() < 2 || row[0] != "Status")
        fail(path.string(), 1, "Status", "expected Status then count columns");
    if (unique_header_column(row, selected, path.string()) == 0)
        fail(path.string(), 1, selected, "selected column is not a count column");

    std::set<std::string> statuses;
    bool assigned = false;
    while (reader.next(row)) {
        if (row[0].empty() || !statuses.insert(row[0]).second)
            fail(path.string(), reader.record(), "Status", "empty or duplicate status");
        assigned = assigned || row[0] == "Assigned";
        for (std::size_t i = 1; i < row.size(); ++i)
            count_value(row[i], path.string() + ": record " +
                                      std::to_string(reader.record()) + ", column " +
                                      std::to_string(i + 1));
    }
    if (!assigned) fail(path.string(), 1, "Status", "Assigned status missing");
}

void read_featurecounts(const fs::path& path, const std::string& selected,
                        const GeneUniverse& genes,
                        std::vector<std::uint64_t>& totals, std::size_t sample,
                        const std::string& manifest_source, std::size_t manifest_record) {
    auto input = open_regular(path, manifest_source, manifest_record, "path");
    TsvReader reader(input, path.string(), true);
    std::vector<std::string> row;
    if (!reader.next(row)) fail(path.string(), 1, "1", "empty file");
    const std::vector<std::string> prefix{"Geneid", "Chr", "Start", "End", "Strand", "Length"};
    if (row.size() < prefix.size() + 1 ||
        !std::equal(prefix.begin(), prefix.end(), row.begin()))
        fail(path.string(), reader.record(), "1",
             "expected six featureCounts annotation columns then count columns");
    const auto selected_column = unique_header_column(row, selected, path.string());
    if (selected_column < prefix.size())
        fail(path.string(), reader.record(), selected, "selected column is not a count column");

    std::vector<bool> seen(genes.ids.size());
    while (reader.next(row)) {
        for (std::size_t i = prefix.size(); i < row.size(); ++i)
            count_value(row[i], path.string() + ": record " +
                                      std::to_string(reader.record()) + ", column " +
                                      std::to_string(i + 1));
        add_count(row, 0, selected_column, selected, path.string(), reader.record(),
                  genes, seen, totals, sample);
    }
    require_complete_universe(seen, genes, path.string());
    validate_featurecounts_summary(fs::path(path.string() + ".summary"), selected,
                                   manifest_source, manifest_record);
}

std::string output_field(std::string_view value) {
    if ((value.empty() || value.front() != '#') &&
        value.find('"') == std::string_view::npos)
        return std::string(value);
    std::string result{"\""};
    for (const char c : value) {
        result.push_back(c);
        if (c == '"') result.push_back('"');
    }
    result.push_back('"');
    return result;
}

class TemporaryOutput {
public:
    explicit TemporaryOutput(const fs::path& output) {
        auto parent = output.parent_path();
        if (parent.empty()) parent = ".";
        const auto base = output.filename().string();
        std::string pattern = (parent / ("." + base + ".tmp.XXXXXX")).string();
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        const int descriptor = ::mkstemp(buffer.data());
        if (descriptor == -1)
            throw std::runtime_error("cannot create temporary output beside " +
                                     output.string() + ": " + std::strerror(errno));
        path_ = buffer.data();
        stream_ = ::fdopen(descriptor, "wb");
        if (!stream_) {
            const auto error = errno;
            ::close(descriptor);
            std::error_code ignored;
            fs::remove(path_, ignored);
            throw std::runtime_error("cannot open temporary output: " +
                                     std::string(std::strerror(error)));
        }
    }

    TemporaryOutput(const TemporaryOutput&) = delete;
    TemporaryOutput& operator=(const TemporaryOutput&) = delete;

    ~TemporaryOutput() {
        if (stream_) std::fclose(stream_);
        if (!published_) {
            std::error_code ignored;
            fs::remove(path_, ignored);
        }
    }

    void write(std::string_view text) {
        if (std::fwrite(text.data(), 1, text.size(), stream_) != text.size())
            throw std::runtime_error("failed writing temporary output");
    }

    void publish(const fs::path& output) {
        if (std::fflush(stream_) != 0)
            throw std::runtime_error("failed flushing temporary output");
        if (std::fclose(stream_) != 0) {
            stream_ = nullptr;
            throw std::runtime_error("failed closing temporary output");
        }
        stream_ = nullptr;
        std::error_code error;
        fs::rename(path_, output, error);
        if (error) throw std::runtime_error("cannot publish output: " + error.message());
        published_ = true;
    }

private:
    fs::path path_;
    std::FILE* stream_ = nullptr;
    bool published_ = false;
};

void write_counts(const fs::path& output, const GeneUniverse& genes,
                  const std::vector<std::string>& samples,
                  const std::vector<std::uint64_t>& totals) {
    TemporaryOutput temporary(output);
    temporary.write("gene_id");
    for (const auto& sample : samples) temporary.write("\t" + output_field(sample));
    temporary.write("\n");
    for (std::size_t gene = 0; gene < genes.ids.size(); ++gene) {
        temporary.write(output_field(genes.ids[gene]));
        for (std::size_t sample = 0; sample < samples.size(); ++sample)
            temporary.write("\t" + std::to_string(totals[gene * samples.size() + sample]));
        temporary.write("\n");
    }
    temporary.publish(output);
}

void validate_annotation_table(const GeneUniverse& genes, const Table& annotation) {
    fixed_schema(annotation, {"gene_id", "gene_symbol", "biotype", "chromosome"});
    const auto id = column(annotation, "gene_id");
    std::set<std::string, decltype(&bytewise_less)> seen(&bytewise_less);
    for (std::size_t row = 0; row < annotation.rows.size(); ++row) {
        const auto& gene = annotation.rows[row][id];
        if (gene.empty() || !genes.index.contains(gene))
            fail(annotation, row + 2, "gene_id", "gene is absent from reference universe");
        if (!seen.insert(gene).second)
            fail(annotation, row + 2, "gene_id", "duplicate gene ID");
    }
}
}

void validate_counts_for_samples(const fs::path& samples_path,
                                 const fs::path& counts_path) {
    const auto samples = read_regular_table(samples_path);
    const auto counts = read_regular_table(counts_path);
    validate_samples(samples);
    validate_counts(counts);
    const auto sample_id = column(samples, "sample_id");
    std::vector<std::string> expected{"gene_id"};
    for (const auto& row : samples.rows) expected.push_back(row[sample_id]);
    if (counts.header != expected)
        fail(counts, 1, "gene_id", "count columns must exactly match sample manifest order");
}

void validate_annotation(const fs::path& genes_path,
                         const fs::path& annotation_path) {
    const auto genes = read_genes(genes_path);
    const auto annotation = read_regular_table(annotation_path);
    validate_annotation_table(genes, annotation);
}

void merge_counts(const fs::path& samples_path, const fs::path& runs_path,
                  const fs::path& genes_path, const fs::path& inputs_path,
                  const fs::path& output, const fs::path& annotation_path) {
    if (output.empty() || output.filename().empty())
        throw std::runtime_error("output filename is empty");
    reject_output_aliases(output, {samples_path, runs_path, genes_path, inputs_path});
    const auto samples = read_regular_table(samples_path);
    const auto runs = read_regular_table(runs_path);
    const auto genes = read_genes(genes_path);
    const auto inputs = read_regular_table(inputs_path);
    validate_samples(samples);
    validate_run_mapping(samples, runs, false);

    if (!annotation_path.empty()) {
        reject_output_alias(output, annotation_path);
        const auto annotation = read_regular_table(annotation_path);
        validate_annotation_table(genes, annotation);
    }

    fixed_schema(inputs, {"run_id", "format", "path", "column"});
    const auto input_run = column(inputs, "run_id");
    const auto format = column(inputs, "format");
    const auto path_column = column(inputs, "path");
    const auto selected_column = column(inputs, "column");
    const auto run_id = column(runs, "run_id");
    const auto run_sample = column(runs, "sample_id");
    const auto sample_id = column(samples, "sample_id");

    std::vector<std::string> sample_ids;
    std::map<std::string, std::size_t> sample_indices;
    for (std::size_t row = 0; row < samples.rows.size(); ++row) {
        const auto& id = samples.rows[row][sample_id];
        if (id == "gene_id")
            fail(samples, row + 2, "sample_id",
                 "gene_id is reserved for the canonical count header");
        sample_indices.emplace(id, sample_ids.size());
        sample_ids.push_back(id);
    }
    std::map<std::string, std::size_t> run_samples;
    for (const auto& row : runs.rows)
        run_samples.emplace(row[run_id], sample_indices.at(row[run_sample]));

    struct CountInput {
        std::string format;
        fs::path path;
        std::string selected;
        std::size_t sample;
        std::size_t record;
    };
    std::vector<CountInput> count_inputs;
    std::set<std::string> covered_runs;
    std::vector<fs::path> protected_paths;
    for (std::size_t row = 0; row < inputs.rows.size(); ++row) {
        const auto record = row + 2;
        const auto& values = inputs.rows[row];
        const auto mapped = run_samples.find(values[input_run]);
        if (mapped == run_samples.end())
            fail(inputs, record, "run_id", "unknown run ID");
        if (!covered_runs.insert(values[input_run]).second)
            fail(inputs, record, "run_id", "duplicate run ID");
        if (values[format] != "tsv" && values[format] != "legacy" &&
            values[format] != "featurecounts")
            fail(inputs, record, "format", "expected tsv, legacy, or featurecounts");
        if (values[format] == "featurecounts") {
            if (values[selected_column].empty())
                fail(inputs, record, "column", "featureCounts column is empty");
        } else if (!values[selected_column].empty()) {
            fail(inputs, record, "column", "column must be empty for tsv and legacy");
        }
        const auto resolved = data_path(inputs, record, "path", values[path_column]);
        readable_regular(resolved, inputs.source, record, "path");
        protected_paths.push_back(resolved);
        if (values[format] == "featurecounts") {
            const fs::path summary(resolved.string() + ".summary");
            readable_regular(summary, inputs.source, record, "path");
            protected_paths.push_back(summary);
        }
        count_inputs.push_back({values[format], resolved, values[selected_column],
                                mapped->second, record});
    }
    for (const auto& [id, unused] : run_samples) {
        static_cast<void>(unused);
        if (!covered_runs.contains(id))
            fail(inputs, 1, "run_id", "run input missing: " + id);
    }
    reject_output_aliases(output, protected_paths);

    // The aggregate matrix is the only gene-by-sample structure retained;
    // individual run tables are consumed one row at a time.
    if (sample_ids.size() > std::numeric_limits<std::size_t>::max() / genes.ids.size())
        throw std::overflow_error("gene-by-sample matrix dimensions exceed size_t range");
    std::vector<std::uint64_t> totals(genes.ids.size() * sample_ids.size());
    for (const auto& source : count_inputs) {
        if (source.format == "tsv")
            read_native_counts(source.path, genes, totals, source.sample,
                               inputs.source, source.record);
        else if (source.format == "legacy")
            read_legacy_counts(source.path, genes, totals, source.sample,
                               inputs.source, source.record);
        else
            read_featurecounts(source.path, source.selected, genes, totals,
                               source.sample, inputs.source, source.record);
    }
    write_counts(output, genes, sample_ids, totals);
}

}
