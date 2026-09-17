#include "rnaseq/alignment.hpp"
#include "rnaseq/counts.hpp"
#include "rnaseq/manifests.hpp"
#include "rnaseq/process.hpp"
#include "rnaseq/table.hpp"
#include "rnaseq/workflow.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: rnaseq workflow-plan|workflow-local|workflow-submit CONFIG\n"
                     "       rnaseq workflow-task CONFIG index|array|finish\n"
                     "       rnaseq validate-tsv|validate-samples|validate-counts FILE\n"
                     "       rnaseq validate-bundle SAMPLES RUNS REFERENCES ANALYSIS CONTRASTS\n"
                     "       rnaseq merge-counts SAMPLES RUNS GENES INPUTS OUTPUT [ANNOTATION]\n"
                     "       rnaseq validate-counts-for-samples SAMPLES COUNTS\n"
                     "       rnaseq validate-annotation GENES ANNOTATION\n"
                     "       rnaseq index|align-count --backend star|hisat2 --bin-dir DIR --fasta "
                     "FILE --gtf FILE --output DIR [--threads N]\n"
                     "       index: [--star-sa-bases N --star-chr-bits N] (genome-only policy)\n"
                     "       align-count: --index DIR --reads1 FILE [--reads2 FILE] --layout "
                     "single|paired --strandedness unstranded|forward|reverse\n";
        return (std::cout << std::flush) ? 0 : 1;
    }
    if (argc < 2) {
        std::cerr << "Use rnaseq --help for usage\n";
        return 2;
    }
    try {
        const std::string command = argv[1];
        if (command.starts_with("workflow-")) {
            rnaseq::workflow_command(argc, argv);
            return 0;
        }
        if (command == "index" || command == "align-count") {
            rnaseq::alignment_command(argc, argv);
            std::cout << "OK: " << command << " completed\n";
            return (std::cout << std::flush) ? 0 : 1;
        }
        if (command == "merge-counts") {
            if (argc != 7 && argc != 8) {
                std::cerr << "Use rnaseq --help for usage\n";
                return 2;
            }
            rnaseq::merge_counts(argv[2], argv[3], argv[4], argv[5], argv[6],
                                 argc == 8 ? argv[7] : "");
            std::cout << "OK: counts merged\n";
            return (std::cout << std::flush) ? 0 : 1;
        }
        if (command == "validate-counts-for-samples" || command == "validate-annotation") {
            if (argc != 4) {
                std::cerr << "Use rnaseq --help for usage\n";
                return 2;
            }
            if (command == "validate-counts-for-samples")
                rnaseq::validate_counts_for_samples(argv[2], argv[3]);
            else
                rnaseq::validate_annotation(argv[2], argv[3]);
            std::cout << "OK: inputs valid\n";
            return (std::cout << std::flush) ? 0 : 1;
        }
        if (command == "validate-bundle") {
            if (argc != 7) {
                std::cerr << "Use rnaseq --help for usage\n";
                return 2;
            }
            rnaseq::validate_bundle(argv[2], argv[3], argv[4], argv[5], argv[6]);
            std::cout << "OK: bundle valid\n";
            return (std::cout << std::flush) ? 0 : 1;
        }
        if (argc != 3) {
            std::cerr << "Use rnaseq --help for usage\n";
            return 2;
        }
        if (command != "validate-tsv" && command != "validate-samples" &&
            command != "validate-counts")
            throw std::runtime_error("unknown command: " + command);
        const auto table = rnaseq::read_table(argv[2]);
        if (command == "validate-samples")
            rnaseq::validate_samples(table);
        if (command == "validate-counts")
            rnaseq::validate_counts(table);
        std::cout << "OK: " << table.rows.size() << " records, " << table.header.size()
                  << " columns\n";
        return (std::cout << std::flush) ? 0 : 1;
    } catch (const rnaseq::ProcessError &error) {
        std::cerr << error.what() << '\n';
        return error.status;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
