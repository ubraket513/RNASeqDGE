#include "doctest.h"
#include "rnaseq/manifests.hpp"
#include "rnaseq/table.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>

namespace {
namespace fs = std::filesystem;

class BundleFixture {
public:
    BundleFixture() {
        char pattern[] = "/tmp/rnaseq-manifests-XXXXXX";
        const auto created = mkdtemp(pattern);
        if (!created) throw std::runtime_error("cannot create manifest fixture");
        root = created;
        fs::create_directories(root / "reads");
        write("reads/s1.fastq", "");
        write("reads/s2.fastq", "");
        write("genome.fa", "");
        write("genes.gtf", "");
        write("samples.tsv",
              "sample_id\tcondition\tage\tbatch\n"
              "s1\ttreated\t41\tb1\n"
              "s2\tcontrol\t52.5\tb2\n");
        write("runs.tsv",
              "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
              "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n"
              "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        write("references.tsv",
              "role\tpath\tsha256\tsource\trelease\n"
              "genome\tgenome.fa\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n"
              "annotation\tgenes.gtf\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n");
        write("analysis.tsv",
              "key\tvalue\n"
              "version\t1\n"
              "design_terms\tcondition,age,batch\n"
              "alpha\t0.05\n"
              "filter\tzero_total\n"
              "shrinkage\tapeglm\n"
              "type.age\tnumeric\n"
              "type.batch\tcategorical\n"
              "reference.batch\tb1\n");
        write("contrasts.tsv",
              "contrast_id\tfactor\tnumerator\tdenominator\n"
              "treated_vs_control\tcondition\ttreated\tcontrol\n");
    }

    ~BundleFixture() { std::error_code error; fs::remove_all(root, error); }

    void write(const fs::path& relative, const std::string& contents) const {
        std::ofstream output(root / relative, std::ios::binary);
        output << contents;
    }

    void validate() const {
        rnaseq::validate_bundle(root / "samples.tsv", root / "runs.tsv",
                                root / "references.tsv", root / "analysis.tsv",
                                root / "contrasts.tsv");
    }

    fs::path root;
};
}

TEST_CASE("merge run mapping permits archived reads but preserves metadata validation") {
    const auto parse = [](const std::string& value, const std::string& source) {
        std::istringstream input(value);
        return rnaseq::read_table(input, source);
    };
    const auto samples = parse("sample_id\tcondition\ns1\tcontrol\n", "samples.tsv");
    auto runs = parse("run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
                      "r1\ts1\t/nonexistent-rnaseq-archived-read.fastq\t\tsingle\tunstranded\n", "runs.tsv");
    CHECK_NOTHROW(rnaseq::validate_run_mapping(samples, runs, false));
    CHECK_THROWS_WITH(rnaseq::validate_run_mapping(samples, runs), doctest::Contains("fastq_1"));
    runs.rows[0][2].clear();
    CHECK_THROWS_WITH(rnaseq::validate_run_mapping(samples, runs, false), doctest::Contains("fastq_1"));
    runs.rows[0][2] = "archived.fastq";
    runs.rows[0][1] = "unknown";
    CHECK_THROWS_WITH(rnaseq::validate_run_mapping(samples, runs, false), doctest::Contains("sample_id"));
    runs.rows[0][1] = "s1";
    runs.rows[0][4] = "paired";
    CHECK_THROWS_WITH(rnaseq::validate_run_mapping(samples, runs, false), doctest::Contains("fastq_2"));
}

TEST_CASE("a complete manifest bundle is accepted") {
    BundleFixture bundle;
    CHECK_NOTHROW(bundle.validate());
}

TEST_CASE("run mappings are complete, unique, and internally consistent") {
    BundleFixture bundle;
    SUBCASE("duplicate run ID") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n"
            "r1\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 3, column run_id"));
    }
    SUBCASE("unknown sample") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n"
            "r2\tmissing\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 3, column sample_id"));
    }
    SUBCASE("sample has no run") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("samples.tsv: record 3, column sample_id"));
    }
    SUBCASE("same sample has mixed layout") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n"
            "r1b\ts1\treads/s1.fastq\treads/s2.fastq\tpaired\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 3, column layout"));
    }
    SUBCASE("same sample has mixed strandedness") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n"
            "r1b\ts1\treads/s1.fastq\t\tsingle\tforward\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 3, column strandedness"));
    }
}

TEST_CASE("run layout controls FASTQ paths and all paths name readable regular files") {
    BundleFixture bundle;
    SUBCASE("single has second mate") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\treads/s2.fastq\tsingle\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 2, column fastq_2"));
    }
    SUBCASE("paired lacks second mate") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tpaired\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 2, column fastq_2"));
    }
    SUBCASE("missing FASTQ") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/missing.fastq\t\tsingle\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("runs.tsv: record 2, column fastq_1"));
    }
    SUBCASE("directory is not a FASTQ file") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads\t\tsingle\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("regular readable file"));
    }
    SUBCASE("absolute paths are accepted unchanged") {
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\t" + (bundle.root / "reads/s1.fastq").string() + "\t\tsingle\tunstranded\n"
            "r2\ts2\t" + (bundle.root / "reads/s2.fastq").string() + "\t\tsingle\tunstranded\n");
        CHECK_NOTHROW(bundle.validate());
    }
    SUBCASE("paired aliases to the same file are rejected") {
        fs::create_symlink(bundle.root / "reads/s1.fastq", bundle.root / "reads/s1-link.fastq");
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\treads/s1-link.fastq\tpaired\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("paired FASTQ files must be distinct"));
    }
    SUBCASE("parent traversal after a symlink follows filesystem resolution") {
        fs::create_directories(bundle.root / "actual/subdir");
        bundle.write("actual/input.fastq", "");
        fs::create_directory_symlink(bundle.root / "actual/subdir", bundle.root / "link");
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\tlink/../input.fastq\t\tsingle\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n");
        CHECK_NOTHROW(bundle.validate());
    }
}

TEST_CASE("reference roles, metadata, digests, and file contents are validated") {
    BundleFixture bundle;
    SUBCASE("checksum mismatch") {
        bundle.write("genome.fa", "changed\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("references.tsv: record 2, column sha256"));
    }
    SUBCASE("malformed checksum") {
        bundle.write("references.tsv",
            "role\tpath\tsha256\tsource\trelease\n"
            "genome\tgenome.fa\tnot-a-digest\tsynthetic\tv1\n"
            "annotation\tgenes.gtf\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("references.tsv: record 2, column sha256"));
    }
    SUBCASE("duplicate role") {
        bundle.write("references.tsv",
            "role\tpath\tsha256\tsource\trelease\n"
            "genome\tgenome.fa\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n"
            "genome\tgenes.gtf\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("references.tsv: record 3, column role"));
    }
    SUBCASE("empty source") {
        bundle.write("references.tsv",
            "role\tpath\tsha256\tsource\trelease\n"
            "genome\tgenome.fa\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\t\tv1\n"
            "annotation\tgenes.gtf\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("references.tsv: record 2, column source"));
    }
    SUBCASE("reference path traversal preserves symlink semantics") {
        fs::create_directories(bundle.root / "actual/subdir");
        bundle.write("actual/genome.fa", "abc");
        fs::create_directory_symlink(bundle.root / "actual/subdir", bundle.root / "link");
        bundle.write("references.tsv",
            "role\tpath\tsha256\tsource\trelease\n"
            "genome\tlink/../genome.fa\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tsynthetic\tv1\n"
            "annotation\tgenes.gtf\te3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\tsynthetic\tv1\n");
        CHECK_NOTHROW(bundle.validate());
    }
}

TEST_CASE("fixed manifest schemas reject extra columns") {
    BundleFixture bundle;
    bundle.write("analysis.tsv",
        "key\tvalue\textra\n"
        "version\t1\tx\n");
    CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 1, column extra"));
}

TEST_CASE("analysis settings are strict and designs are identifier lists") {
    BundleFixture bundle;
    SUBCASE("unknown setting") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\nsurprise\tyes\n"
            "type.age\tnumeric\ntype.batch\tcategorical\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 7, column key"));
    }
    SUBCASE("formula syntax") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition+age\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 3, column value"));
    }
    SUBCASE("duplicate design term") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition,condition\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 3, column value"));
    }
    SUBCASE("undeclared sample covariate type") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\ntype.age\tnumeric\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("samples.tsv: record 1, column batch"));
    }
    SUBCASE("typed metadata outside the design need not be an identifier") {
        bundle.write("samples.tsv",
            "sample_id\tcondition\tcollection-site\n"
            "s1\ttreated\tSeoul\n"
            "s2\tcontrol\tBusan\n");
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.collection-site\tcategorical\n");
        CHECK_NOTHROW(bundle.validate());
    }
}

TEST_CASE("numeric settings and covariates must be finite decimals") {
    SUBCASE("signed exponent forms are accepted") {
        BundleFixture bundle;
        bundle.write("samples.tsv",
            "sample_id\tcondition\tage\tbatch\n"
            "s1\ttreated\t+4.1e1\tb1\n"
            "s2\tcontrol\t-.525E+2\tb2\n");
        CHECK_NOTHROW(bundle.validate());
    }
    for (const std::string value : {"NaN", "Inf", "-Inf", "", "1x"}) {
        CAPTURE(value);
        BundleFixture bundle;
        bundle.write("samples.tsv",
            "sample_id\tcondition\tage\tbatch\n"
            "s1\ttreated\t" + value + "\tb1\n"
            "s2\tcontrol\t52.5\tb2\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("samples.tsv: record 2, column age"));
    }
    for (const std::string alpha : {"0", "1", "NaN", "0.2x"}) {
        CAPTURE(alpha);
        BundleFixture bundle;
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition,age,batch\nalpha\t" + alpha + "\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\nreference.batch\tb1\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 4, column value"));
    }
}

TEST_CASE("contrasts use observed categorical design levels and explicit references") {
    BundleFixture bundle;
    SUBCASE("numeric contrast factor") {
        bundle.write("contrasts.tsv",
            "contrast_id\tfactor\tnumerator\tdenominator\n"
            "old_vs_young\tage\t52.5\t41\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("contrasts.tsv: record 2, column factor"));
    }
    SUBCASE("factor absent from design") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition,age\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\n");
        bundle.write("contrasts.tsv",
            "contrast_id\tfactor\tnumerator\tdenominator\n"
            "b2_vs_b1\tbatch\tb2\tb1\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("contrasts.tsv: record 2, column factor"));
    }
    SUBCASE("missing categorical level") {
        bundle.write("contrasts.tsv",
            "contrast_id\tfactor\tnumerator\tdenominator\n"
            "treated_vs_missing\tcondition\ttreated\tmissing\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("contrasts.tsv: record 2, column denominator"));
    }
    SUBCASE("same numerator and denominator") {
        bundle.write("contrasts.tsv",
            "contrast_id\tfactor\tnumerator\tdenominator\n"
            "same\tcondition\tcontrol\tcontrol\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("contrasts.tsv: record 2, column numerator"));
    }
    SUBCASE("missing other categorical reference") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition,age,batch\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 3, column value"));
    }
    SUBCASE("reference is unobserved") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition,age,batch\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\nreference.batch\tmissing\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 9, column value"));
    }
    SUBCASE("explicit condition reference disagrees with denominator") {
        bundle.write("analysis.tsv",
            "key\tvalue\nversion\t1\ndesign_terms\tcondition,age,batch\nalpha\t0.05\n"
            "filter\tzero_total\nshrinkage\tapeglm\n"
            "type.age\tnumeric\ntype.batch\tcategorical\nreference.batch\tb1\n"
            "reference.condition\ttreated\n");
        CHECK_THROWS_WITH(bundle.validate(), doctest::Contains("analysis.tsv: record 10, column value"));
    }
    SUBCASE("multiple contrasts may request different denominators without a fixed reference") {
        bundle.write("reads/s3.fastq", "");
        bundle.write("samples.tsv",
            "sample_id\tcondition\tage\tbatch\n"
            "s1\ttreated\t41\tb1\n"
            "s2\tcontrol\t52.5\tb2\n"
            "s3\tvehicle\t38\tb1\n");
        bundle.write("runs.tsv",
            "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n"
            "r1\ts1\treads/s1.fastq\t\tsingle\tunstranded\n"
            "r2\ts2\treads/s2.fastq\t\tsingle\tunstranded\n"
            "r3\ts3\treads/s3.fastq\t\tsingle\tunstranded\n");
        bundle.write("contrasts.tsv",
            "contrast_id\tfactor\tnumerator\tdenominator\n"
            "treated_vs_control\tcondition\ttreated\tcontrol\n"
            "treated_vs_vehicle\tcondition\ttreated\tvehicle\n");
        CHECK_NOTHROW(bundle.validate());
    }
}
