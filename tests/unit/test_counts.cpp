#include "doctest.h"
#include "rnaseq/counts.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include <sys/stat.h>

namespace {
namespace fs = std::filesystem;

class CountsFixture {
public:
    CountsFixture() {
        char pattern[] = "/tmp/rnaseq-counts-XXXXXX";
        const auto created = mkdtemp(pattern);
        if (!created) throw std::runtime_error("cannot create counts fixture");
        root = created;
        fs::create_directories(root / "data");
        write("samples.tsv",
              "sample_id\tcondition\n"
              "sample_b\ttreated\n"
              "sample_a\tcontrol\n");
        write("runs.tsv",
              "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
              "run_tsv\tsample_b\tarchived/one.fastq\t\tsingle\tunstranded\n"
              "run_legacy\tsample_a\tarchived/two.fastq\t\tsingle\tunstranded\n"
              "run_fc\tsample_a\tarchived/three.fastq\t\tsingle\tunstranded\n");
        write("genes.tsv", "gene_id\ngene_z\ngene_a\ngene_b\n");
        write("inputs.tsv",
              "run_id\tformat\tpath\tcolumn\n"
              "run_fc\tfeaturecounts\tdata/fc.txt\tselected.bam\n"
              "run_tsv\ttsv\tdata/native.tsv\t\n"
              "run_legacy\tlegacy\tdata/legacy.tsv\t\n");
        write("data/native.tsv",
              "gene_id\tcount\n"
              "gene_b\t3\n"
              "gene_z\t4\n"
              "gene_a\t2\n");
        write("data/legacy.tsv", "gene_z\t7\ngene_a\t5\ngene_b\t6\n");
        write("data/fc.txt",
              "# Program:featureCounts v2.0.6\n"
              "Geneid\tChr\tStart\tEnd\tStrand\tLength\tother.bam\tselected.bam\n"
              "gene_a\tchr1\t1\t10\t+\t10\t900\t11\n"
              "gene_z\tchr1\t21\t30\t+\t10\t901\t13\n"
              "gene_b\tchr1\t11\t20\t-\t10\t902\t12\n");
        write("data/fc.txt.summary",
              "Status\tother.bam\tselected.bam\n"
              "Assigned\t999\t36\n"
              "Unassigned_Ambiguity\t0\t1\n");
        write("annotation.tsv",
              "gene_id\tgene_symbol\tbiotype\tchromosome\n"
              "gene_b\tSHARED\tprotein_coding\tchr1\n"
              "gene_a\tSHARED\t\t\n");
    }

    ~CountsFixture() {
        std::error_code error;
        fs::remove_all(root, error);
    }

    void write(const fs::path& relative, const std::string& contents) const {
        const auto path = root / relative;
        if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
        std::ofstream output(path, std::ios::binary);
        output << contents;
        if (!output) throw std::runtime_error("cannot write counts fixture");
    }

    std::string read(const fs::path& relative) const {
        std::ifstream input(root / relative, std::ios::binary);
        return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
    }

    void merge(const fs::path& output = "counts.tsv", const fs::path& annotation = {}) const {
        rnaseq::merge_counts(root / "samples.tsv", root / "runs.tsv", root / "genes.tsv",
                             root / "inputs.tsv", root / output,
                             annotation.empty() ? fs::path{} : root / annotation);
    }

    fs::path root;
};
}

TEST_CASE("all count formats aggregate technical runs in canonical order") {
    CountsFixture fixture;
    CHECK_NOTHROW(fixture.merge("counts.tsv", "annotation.tsv"));
    CHECK(fixture.read("counts.tsv") ==
          "gene_id\tsample_b\tsample_a\n"
          "gene_a\t2\t16\n"
          "gene_b\t3\t18\n"
          "gene_z\t4\t20\n");
}

TEST_CASE("canonical output quotes a gene ID beginning with a comment marker") {
    CountsFixture fixture;
    fixture.write("genes.tsv", "gene_id\ngene_z\n\"#gene\"\ngene_b\n");
    fixture.write("data/native.tsv",
                  "gene_id\tcount\n"
                  "gene_b\t3\n"
                  "gene_z\t4\n"
                  "\"#gene\"\t2\n");
    fixture.write("data/legacy.tsv", "gene_z\t7\n\"#gene\"\t5\ngene_b\t6\n");
    fixture.write("data/fc.txt",
                  "Geneid\tChr\tStart\tEnd\tStrand\tLength\tother.bam\tselected.bam\n"
                  "\"#gene\"\tchr1\t1\t10\t+\t10\t900\t11\n"
                  "gene_z\tchr1\t21\t30\t+\t10\t901\t13\n"
                  "gene_b\tchr1\t11\t20\t-\t10\t902\t12\n");

    REQUIRE_NOTHROW(fixture.merge());
    CHECK(fixture.read("counts.tsv") ==
          "gene_id\tsample_b\tsample_a\n"
          "\"#gene\"\t2\t16\n"
          "gene_b\t3\t18\n"
          "gene_z\t4\t20\n");
    CHECK_NOTHROW(rnaseq::validate_counts_for_samples(fixture.root / "samples.tsv",
                                                       fixture.root / "counts.tsv"));
}

TEST_CASE("sample IDs cannot collide with the canonical gene column") {
    CountsFixture fixture;
    fixture.write("samples.tsv",
                  "sample_id\tcondition\n"
                  "gene_id\ttreated\n"
                  "sample_a\tcontrol\n");
    fixture.write("runs.tsv",
                  "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
                  "run_tsv\tgene_id\tarchived/one.fastq\t\tsingle\tunstranded\n"
                  "run_legacy\tsample_a\tarchived/two.fastq\t\tsingle\tunstranded\n"
                  "run_fc\tsample_a\tarchived/three.fastq\t\tsingle\tunstranded\n");
    fixture.write("counts.tsv", "previous output\n");

    CHECK_THROWS_WITH(fixture.merge(),
                      doctest::Contains("samples.tsv: record 2, column sample_id"));
    CHECK(fixture.read("counts.tsv") == "previous output\n");
}

TEST_CASE("canonical counts require the exact manifest sample columns") {
    CountsFixture fixture;
    fixture.write("counts.tsv",
                  "gene_id\tsample_b\tsample_a\n"
                  "gene_a\t2\t16\n");
    CHECK_NOTHROW(rnaseq::validate_counts_for_samples(fixture.root / "samples.tsv",
                                                       fixture.root / "counts.tsv"));

    for (const std::string header : {
             "gene_id\tsample_a\tsample_b\n", "gene_id\tsample_b\n",
             "gene_id\tsample_b\tsample_a\textra\n"}) {
        CAPTURE(header);
        fixture.write("counts.tsv", header + "gene_a\t2\t16\n");
        CHECK_THROWS_AS(rnaseq::validate_counts_for_samples(fixture.root / "samples.tsv",
                                                            fixture.root / "counts.tsv"),
                        std::runtime_error);
    }
}

TEST_CASE("annotation is an independently validated gene subset") {
    CountsFixture fixture;
    CHECK_NOTHROW(rnaseq::validate_annotation(fixture.root / "genes.tsv",
                                               fixture.root / "annotation.tsv"));

    SUBCASE("unknown gene") {
        fixture.write("annotation.tsv",
                      "gene_id\tgene_symbol\tbiotype\tchromosome\n"
                      "missing\tX\t\t\n");
        CHECK_THROWS_WITH(rnaseq::validate_annotation(fixture.root / "genes.tsv",
                                                       fixture.root / "annotation.tsv"),
                          doctest::Contains("annotation.tsv: record 2, column gene_id"));
    }
    SUBCASE("duplicate gene") {
        fixture.write("annotation.tsv",
                      "gene_id\tgene_symbol\tbiotype\tchromosome\n"
                      "gene_a\tX\t\t\n"
                      "gene_a\tX\t\t\n");
        CHECK_THROWS_AS(rnaseq::validate_annotation(fixture.root / "genes.tsv",
                                                     fixture.root / "annotation.tsv"),
                        std::runtime_error);
    }
    SUBCASE("wrong schema") {
        fixture.write("annotation.tsv", "gene_id\tgene_symbol\n gene_a\tX\n");
        CHECK_THROWS_AS(rnaseq::validate_annotation(fixture.root / "genes.tsv",
                                                     fixture.root / "annotation.tsv"),
                        std::runtime_error);
    }
}

TEST_CASE("mapping and gene-universe failures preserve an existing output") {
    CountsFixture fixture;
    fixture.write("counts.tsv", "previous output\n");

    SUBCASE("missing run input") {
        fixture.write("inputs.tsv",
                      "run_id\tformat\tpath\tcolumn\n"
                      "run_tsv\ttsv\tdata/native.tsv\t\n"
                      "run_legacy\tlegacy\tdata/legacy.tsv\t\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("duplicate run input") {
        fixture.write("inputs.tsv",
                      "run_id\tformat\tpath\tcolumn\n"
                      "run_tsv\ttsv\tdata/native.tsv\t\n"
                      "run_tsv\tlegacy\tdata/legacy.tsv\t\n"
                      "run_fc\tfeaturecounts\tdata/fc.txt\tselected.bam\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("duplicate reference gene") {
        fixture.write("genes.tsv", "gene_id\ngene_z\ngene_z\ngene_a\n");
        CHECK_THROWS_WITH(fixture.merge(),
                          doctest::Contains("genes.tsv: record 3, column gene_id"));
    }
    SUBCASE("run gene does not match universe") {
        fixture.write("data/native.tsv",
                      "gene_id\tcount\n"
                      "gene_a\t2\n"
                      "gene_b\t3\n"
                      "gene_other\t4\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    CHECK(fixture.read("counts.tsv") == "previous output\n");
}

TEST_CASE("malformed counts and checked aggregation preserve an existing output") {
    CountsFixture fixture;
    fixture.write("counts.tsv", "previous output\n");

    SUBCASE("fractional count") {
        fixture.write("data/native.tsv",
                      "gene_id\tcount\n"
                      "gene_a\t1.5\n"
                      "gene_b\t3\n"
                      "gene_z\t4\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("duplicate run gene") {
        fixture.write("data/legacy.tsv", "gene_z\t7\ngene_a\t5\ngene_a\t6\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("technical-run overflow") {
        fixture.write("data/legacy.tsv",
                      "gene_z\t7\n"
                      "gene_a\t18446744073709551615\n"
                      "gene_b\t6\n");
        CHECK_THROWS_AS(fixture.merge(), std::overflow_error);
    }
    SUBCASE("late malformed row") {
        fixture.write("data/native.tsv",
                      "gene_id\tcount\n"
                      "gene_a\t2\n"
                      "gene_b\t3\n"
                      "gene_z\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    CHECK(fixture.read("counts.tsv") == "previous output\n");
}

TEST_CASE("featureCounts selected column and summary are strict") {
    CountsFixture fixture;
    fixture.write("counts.tsv", "previous output\n");

    SUBCASE("selected count column is absent") {
        fixture.write("inputs.tsv",
                      "run_id\tformat\tpath\tcolumn\n"
                      "run_fc\tfeaturecounts\tdata/fc.txt\tmissing.bam\n"
                      "run_tsv\ttsv\tdata/native.tsv\t\n"
                      "run_legacy\tlegacy\tdata/legacy.tsv\t\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("summary is absent") {
        fs::remove(fixture.root / "data/fc.txt.summary");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("Assigned is absent") {
        fixture.write("data/fc.txt.summary",
                      "Status\tother.bam\tselected.bam\n"
                      "Unassigned_Ambiguity\t0\t1\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("summary status repeats") {
        fixture.write("data/fc.txt.summary",
                      "Status\tother.bam\tselected.bam\n"
                      "Assigned\t999\t36\n"
                      "Assigned\t999\t36\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("summary value is fractional") {
        fixture.write("data/fc.txt.summary",
                      "Status\tother.bam\tselected.bam\n"
                      "Assigned\t999\t35.5\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("selected summary column must follow Status") {
        fixture.write("inputs.tsv",
                      "run_id\tformat\tpath\tcolumn\n"
                      "run_fc\tfeaturecounts\tdata/fc.txt\tStatus\n"
                      "run_tsv\ttsv\tdata/native.tsv\t\n"
                      "run_legacy\tlegacy\tdata/legacy.tsv\t\n");
        fixture.write("data/fc.txt",
                      "Geneid\tChr\tStart\tEnd\tStrand\tLength\tStatus\n"
                      "gene_a\tchr1\t1\t10\t+\t10\t11\n"
                      "gene_z\tchr1\t21\t30\t+\t10\t13\n"
                      "gene_b\tchr1\t11\t20\t-\t10\t12\n");
        fixture.write("data/fc.txt.summary",
                      "Status\tother.bam\n"
                      "Assigned\t36\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    CHECK(fixture.read("counts.tsv") == "previous output\n");
}

TEST_CASE("manifest-relative paths retain symlink parent traversal semantics") {
    CountsFixture fixture;
    fs::create_directories(fixture.root / "actual/subdir");
    fixture.write("actual/native.tsv", fixture.read("data/native.tsv"));
    fs::create_directory_symlink(fixture.root / "actual/subdir", fixture.root / "link");
    fixture.write("inputs.tsv",
                  "run_id\tformat\tpath\tcolumn\n"
                  "run_fc\tfeaturecounts\tdata/fc.txt\tselected.bam\n"
                  "run_tsv\ttsv\tlink/../native.tsv\t\n"
                  "run_legacy\tlegacy\tdata/legacy.tsv\t\n");
    CHECK_NOTHROW(fixture.merge());
}

TEST_CASE("non-regular inputs and output aliases are rejected") {
    CountsFixture fixture;

    SUBCASE("FIFO count input") {
        const auto fifo = fixture.root / "data/fifo.tsv";
        REQUIRE(::mkfifo(fifo.c_str(), 0600) == 0);
        fixture.write("inputs.tsv",
                      "run_id\tformat\tpath\tcolumn\n"
                      "run_fc\tfeaturecounts\tdata/fc.txt\tselected.bam\n"
                      "run_tsv\ttsv\tdata/fifo.tsv\t\n"
                      "run_legacy\tlegacy\tdata/legacy.tsv\t\n");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("hard-linked output") {
        fs::create_hard_link(fixture.root / "data/native.tsv", fixture.root / "counts.tsv");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
        CHECK(fixture.read("data/native.tsv") ==
              "gene_id\tcount\n"
              "gene_b\t3\n"
              "gene_z\t4\n"
              "gene_a\t2\n");
    }
    SUBCASE("symlinked output") {
        fs::create_symlink(fixture.root / "data/native.tsv", fixture.root / "counts.tsv");
        CHECK_THROWS_AS(fixture.merge(), std::runtime_error);
    }
    SUBCASE("output path has no filename") {
        fs::create_directories(fixture.root / "output-dir");
        const fs::path output(fixture.root.string() + "/output-dir/");
        CHECK_THROWS_WITH(fixture.merge(output), doctest::Contains("output filename is empty"));
    }
}
