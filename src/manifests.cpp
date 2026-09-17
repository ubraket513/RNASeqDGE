#include "rnaseq/manifests.hpp"

#include "rnaseq/file_hash.hpp"
#include "rnaseq/table.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rnaseq {
namespace {
namespace fs = std::filesystem;

[[noreturn]] void fail(const Table &table, std::size_t record, const std::string &column,
                       const std::string &message) {
    throw std::runtime_error(table.source + ": record " + std::to_string(record) + ", column " +
                             column + ": " + message);
}

std::size_t column(const Table &table, const std::string &name) {
    const auto found = std::find(table.header.begin(), table.header.end(), name);
    if (found == table.header.end())
        fail(table, 1, name, "required column missing");
    return static_cast<std::size_t>(found - table.header.begin());
}

void fixed_schema(const Table &table, const std::vector<std::string> &expected) {
    for (const auto &name : expected)
        column(table, name);
    const std::set<std::string> allowed(expected.begin(), expected.end());
    for (const auto &name : table.header)
        if (!allowed.contains(name))
            fail(table, 1, name, "unexpected column");
}

bool valid_id(std::string_view value) {
    const auto alnum = [](unsigned char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
    };
    return !value.empty() && alnum(static_cast<unsigned char>(value.front())) &&
           std::all_of(value.begin(), value.end(), [&](unsigned char c) {
               return alnum(c) || c == '_' || c == '.' || c == '-';
           });
}

bool valid_identifier(std::string_view value) {
    const auto alpha = [](unsigned char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    };
    const auto digit = [](unsigned char c) { return c >= '0' && c <= '9'; };
    return !value.empty() && alpha(static_cast<unsigned char>(value.front())) &&
           std::all_of(value.begin() + 1, value.end(),
                       [&](unsigned char c) { return alpha(c) || digit(c) || c == '_'; });
}

bool decimal(std::string_view value, double &result) {
    if (value.empty())
        return false;
    std::size_t position = 0;
    if (value[position] == '+' || value[position] == '-') {
        if (++position == value.size())
            return false;
    }
    bool digits = false;
    while (position < value.size() && value[position] >= '0' && value[position] <= '9') {
        digits = true;
        ++position;
    }
    if (position < value.size() && value[position] == '.') {
        ++position;
        while (position < value.size() && value[position] >= '0' && value[position] <= '9') {
            digits = true;
            ++position;
        }
    }
    if (!digits)
        return false;
    if (position < value.size() && (value[position] == 'e' || value[position] == 'E')) {
        ++position;
        if (position < value.size() && (value[position] == '+' || value[position] == '-'))
            ++position;
        const auto exponent = position;
        while (position < value.size() && value[position] >= '0' && value[position] <= '9')
            ++position;
        if (position == exponent)
            return false;
    }
    if (position != value.size())
        return false;
    std::string normalized(value);
    if (normalized.front() == '+')
        normalized.erase(0, 1);
    const auto parsed = std::from_chars(normalized.data(), normalized.data() + normalized.size(),
                                        result, std::chars_format::general);
    return parsed.ec == std::errc{} && parsed.ptr == normalized.data() + normalized.size() &&
           std::isfinite(result);
}

fs::path data_path(const Table &table, std::size_t record, const std::string &column_name,
                   const std::string &value) {
    if (value.empty())
        fail(table, record, column_name, "path is empty");
    const fs::path supplied(value);
    return fs::path(table.source).parent_path() / supplied;
}

void readable_regular(const Table &table, std::size_t record, const std::string &column_name,
                      const fs::path &path) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error)
        fail(table, record, column_name, "path does not name a regular readable file");
    std::ifstream input(path, std::ios::binary);
    if (!input)
        fail(table, record, column_name, "path does not name a regular readable file");
}

std::vector<std::string> split_design(const Table &analysis, std::size_t record,
                                      const std::string &value) {
    std::vector<std::string> result;
    std::set<std::string> seen;
    std::size_t start = 0;
    while (start <= value.size()) {
        const auto comma = value.find(',', start);
        const auto term = value.substr(start, comma == std::string::npos ? comma : comma - start);
        if (!valid_identifier(term) || !seen.insert(term).second)
            fail(analysis, record, "value", "design terms must be unique identifiers");
        result.push_back(term);
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }
    return result;
}

std::string lower_hex(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        if (c >= 'A' && c <= 'F')
            return static_cast<char>(c - 'A' + 'a');
        return static_cast<char>(c);
    });
    return value;
}

bool valid_digest(std::string_view value) {
    return value.size() == 64 && std::all_of(value.begin(), value.end(), [](unsigned char c) {
               return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
           });
}

struct Setting {
    std::string value;
    std::size_t record;
};
} // namespace

void validate_run_mapping(const Table &samples, const Table &runs, bool check_fastq_files) {
    validate_samples(samples);
    const auto sample_id = column(samples, "sample_id");
    std::map<std::string, std::size_t> sample_records;
    for (std::size_t row = 0; row < samples.rows.size(); ++row)
        sample_records.emplace(samples.rows[row][sample_id], row + 2);

    fixed_schema(runs, {"run_id", "sample_id", "fastq_1", "fastq_2", "layout", "strandedness"});
    if (runs.rows.empty())
        fail(runs, 2, "run_id", "no runs");
    const auto run_id = column(runs, "run_id"), run_sample = column(runs, "sample_id");
    const auto fastq_1 = column(runs, "fastq_1"), fastq_2 = column(runs, "fastq_2");
    const auto layout = column(runs, "layout"), strand = column(runs, "strandedness");
    std::set<std::string> run_ids, covered;
    std::map<std::string, std::pair<std::string, std::string>> sample_library;
    for (std::size_t row = 0; row < runs.rows.size(); ++row) {
        const auto record = row + 2;
        const auto &values = runs.rows[row];
        if (!valid_id(values[run_id]) || !run_ids.insert(values[run_id]).second)
            fail(runs, record, "run_id", "invalid or duplicate run ID");
        if (!sample_records.contains(values[run_sample]))
            fail(runs, record, "sample_id", "unknown sample ID");
        covered.insert(values[run_sample]);
        if (values[layout] != "single" && values[layout] != "paired")
            fail(runs, record, "layout", "expected single or paired");
        if (values[strand] != "unstranded" && values[strand] != "forward" &&
            values[strand] != "reverse")
            fail(runs, record, "strandedness", "expected unstranded, forward, or reverse");
        const auto [existing, inserted] =
            sample_library.emplace(values[run_sample], std::pair{values[layout], values[strand]});
        if (!inserted && existing->second.first != values[layout])
            fail(runs, record, "layout", "runs for a sample must use one layout");
        if (!inserted && existing->second.second != values[strand])
            fail(runs, record, "strandedness", "runs for a sample must use one strandedness");
        const auto first = data_path(runs, record, "fastq_1", values[fastq_1]);
        if (check_fastq_files)
            readable_regular(runs, record, "fastq_1", first);
        if (values[layout] == "single") {
            if (!values[fastq_2].empty())
                fail(runs, record, "fastq_2", "single-end run must not have a second FASTQ");
        } else {
            const auto second = data_path(runs, record, "fastq_2", values[fastq_2]);
            if (!check_fastq_files)
                continue;
            readable_regular(runs, record, "fastq_2", second);
            std::error_code error;
            if (fs::equivalent(first, second, error) && !error)
                fail(runs, record, "fastq_2", "paired FASTQ files must be distinct");
            if (error)
                fail(runs, record, "fastq_2", "cannot compare paired FASTQ files");
        }
    }
    for (const auto &[id, record] : sample_records)
        if (!covered.contains(id))
            fail(samples, record, "sample_id", "sample has no run");
}

void validate_bundle(const fs::path &samples_path, const fs::path &runs_path,
                     const fs::path &references_path, const fs::path &analysis_path,
                     const fs::path &contrasts_path) {
    const auto samples = read_table(samples_path.string());
    const auto runs = read_table(runs_path.string());
    const auto references = read_table(references_path.string());
    const auto analysis = read_table(analysis_path.string());
    const auto contrasts = read_table(contrasts_path.string());

    validate_run_mapping(samples, runs);

    fixed_schema(references, {"role", "path", "sha256", "source", "release"});
    const auto role = column(references, "role"), reference_path = column(references, "path");
    const auto digest = column(references, "sha256"), source = column(references, "source");
    const auto release = column(references, "release");
    std::set<std::string> roles;
    for (std::size_t row = 0; row < references.rows.size(); ++row) {
        const auto record = row + 2;
        const auto &values = references.rows[row];
        if ((values[role] != "genome" && values[role] != "annotation") ||
            !roles.insert(values[role]).second)
            fail(references, record, "role", "expected genome and annotation exactly once");
        if (values[source].empty())
            fail(references, record, "source", "required cell empty");
        if (values[release].empty())
            fail(references, record, "release", "required cell empty");
        if (!valid_digest(values[digest]))
            fail(references, record, "sha256", "expected 64 hexadecimal digits");
        const auto path = data_path(references, record, "path", values[reference_path]);
        readable_regular(references, record, "path", path);
        std::string actual;
        try {
            actual = sha256_file(path);
        } catch (const std::exception &error) {
            fail(references, record, "path", error.what());
        }
        if (actual != lower_hex(values[digest]))
            fail(references, record, "sha256", "checksum does not match file contents");
    }
    if (!roles.contains("genome"))
        fail(references, 1, "role", "genome role missing");
    if (!roles.contains("annotation"))
        fail(references, 1, "role", "annotation role missing");

    fixed_schema(analysis, {"key", "value"});
    const auto key = column(analysis, "key"), setting_value = column(analysis, "value");
    std::map<std::string, Setting> settings;
    for (std::size_t row = 0; row < analysis.rows.size(); ++row) {
        const auto record = row + 2;
        const auto &name = analysis.rows[row][key];
        if (name.empty() ||
            !settings.emplace(name, Setting{analysis.rows[row][setting_value], record}).second)
            fail(analysis, record, "key", "empty or duplicate setting key");
    }
    const auto require = [&](const std::string &name) -> const Setting & {
        const auto found = settings.find(name);
        if (found == settings.end())
            fail(analysis, 1, name, "required setting missing");
        return found->second;
    };
    if (require("version").value != "1")
        fail(analysis, require("version").record, "value", "version must be 1");
    if (require("filter").value != "zero_total")
        fail(analysis, require("filter").record, "value", "filter must be zero_total");
    if (require("shrinkage").value != "apeglm")
        fail(analysis, require("shrinkage").record, "value", "shrinkage must be apeglm");
    double alpha = 0;
    if (!decimal(require("alpha").value, alpha) || !(alpha > 0 && alpha < 1))
        fail(analysis, require("alpha").record, "value",
             "alpha must be a finite decimal between 0 and 1");
    const auto design =
        split_design(analysis, require("design_terms").record, require("design_terms").value);

    std::map<std::string, std::size_t> sample_columns;
    for (std::size_t i = 0; i < samples.header.size(); ++i)
        sample_columns.emplace(samples.header[i], i);
    std::map<std::string, std::string> types{{"condition", "categorical"}};
    std::map<std::string, Setting> references_by_column;
    for (const auto &[name, setting] : settings) {
        if (name == "version" || name == "design_terms" || name == "alpha" || name == "filter" ||
            name == "shrinkage")
            continue;
        if (name.starts_with("type.")) {
            const auto field = name.substr(5);
            if (field == "sample_id" || !sample_columns.contains(field))
                fail(analysis, setting.record, "key", "type declaration names absent covariate");
            if (setting.value != "categorical" && setting.value != "numeric")
                fail(analysis, setting.record, "value", "type must be categorical or numeric");
            if (field == "condition" && setting.value != "categorical")
                fail(analysis, setting.record, "value", "condition must be categorical");
            types[field] = setting.value;
        } else if (name.starts_with("reference.")) {
            const auto field = name.substr(10);
            if (field == "sample_id" || !sample_columns.contains(field))
                fail(analysis, setting.record, "key", "reference names absent covariate");
            references_by_column.emplace(field, setting);
        } else {
            fail(analysis, setting.record, "key", "unknown analysis setting");
        }
    }
    for (std::size_t i = 0; i < samples.header.size(); ++i) {
        const auto &name = samples.header[i];
        if (name == "sample_id")
            continue;
        if (name != "condition" && !types.contains(name))
            fail(samples, 1, name, "covariate requires type declaration");
        for (std::size_t row = 0; row < samples.rows.size(); ++row) {
            const auto &value = samples.rows[row][i];
            if (types.at(name) == "categorical") {
                if (value.empty())
                    fail(samples, row + 2, name, "categorical value is empty");
            } else {
                double parsed = 0;
                if (!decimal(value, parsed))
                    fail(samples, row + 2, name, "numeric value must be a finite decimal");
            }
        }
    }
    std::set<std::string> design_set;
    for (const auto &term : design) {
        if (term == "sample_id" || !sample_columns.contains(term))
            fail(analysis, require("design_terms").record, "value",
                 "design term is not a sample covariate");
        design_set.insert(term);
    }
    for (const auto &[field, setting] : references_by_column) {
        if (!types.contains(field) || types.at(field) != "categorical")
            fail(analysis, setting.record, "key", "references require categorical covariates");
        const auto index = sample_columns.at(field);
        const bool observed =
            std::any_of(samples.rows.begin(), samples.rows.end(),
                        [&](const auto &row) { return row[index] == setting.value; });
        if (!observed)
            fail(analysis, setting.record, "value", "reference level is not observed");
    }

    fixed_schema(contrasts, {"contrast_id", "factor", "numerator", "denominator"});
    if (contrasts.rows.empty())
        fail(contrasts, 2, "contrast_id", "no contrasts");
    const auto contrast_id = column(contrasts, "contrast_id"), factor = column(contrasts, "factor");
    const auto numerator = column(contrasts, "numerator"),
               denominator = column(contrasts, "denominator");
    std::set<std::string> contrast_ids, contrasted;
    std::map<std::string, std::set<std::string>> denominators;
    for (std::size_t row = 0; row < contrasts.rows.size(); ++row) {
        const auto record = row + 2;
        const auto &values = contrasts.rows[row];
        if (!valid_id(values[contrast_id]) || !contrast_ids.insert(values[contrast_id]).second)
            fail(contrasts, record, "contrast_id", "invalid or duplicate contrast ID");
        if (!valid_identifier(values[factor]) || !design_set.contains(values[factor]))
            fail(contrasts, record, "factor", "contrast factor must be a design term");
        if (!types.contains(values[factor]) || types.at(values[factor]) != "categorical")
            fail(contrasts, record, "factor", "contrast factor must be categorical");
        if (values[numerator].empty() || values[numerator] == values[denominator])
            fail(contrasts, record, "numerator", "contrast levels must be nonempty and distinct");
        if (values[denominator].empty())
            fail(contrasts, record, "denominator", "contrast level must be nonempty");
        const auto index = sample_columns.at(values[factor]);
        const auto observed = [&](const std::string &level) {
            return std::any_of(samples.rows.begin(), samples.rows.end(),
                               [&](const auto &sample) { return sample[index] == level; });
        };
        if (!observed(values[numerator]))
            fail(contrasts, record, "numerator", "categorical level is not observed");
        if (!observed(values[denominator]))
            fail(contrasts, record, "denominator", "categorical level is not observed");
        contrasted.insert(values[factor]);
        denominators[values[factor]].insert(values[denominator]);
    }
    for (const auto &[field, requested] : denominators) {
        const auto explicit_reference = references_by_column.find(field);
        if (explicit_reference != references_by_column.end() &&
            (requested.size() != 1 || !requested.contains(explicit_reference->second.value)))
            fail(analysis, explicit_reference->second.record, "value",
                 "explicit reference must agree with all contrast denominators");
    }
    for (const auto &term : design)
        if (types.at(term) == "categorical" && !contrasted.contains(term) &&
            !references_by_column.contains(term))
            fail(analysis, require("design_terms").record, "value",
                 "categorical design term requires a reference or contrast");
}

} // namespace rnaseq
