#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "rnaseq/table.hpp"
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
rnaseq::Table parse(const std::string& text) {
    std::istringstream input(text);
    return rnaseq::read_table(input, "fixture.tsv");
}
}

TEST_CASE("quoted cells, CRLF, initial BOM, empty optional cells and UTF-8") {
    auto table = parse("\xEF\xBB\xBFsample_id\tcondition\tnote\r\n"
                       "s1\t\"treated\"\t\"say \"\"hello\"\"\"\r\n"
                       "s2\t대조군\t\r\n");
    REQUIRE(table.rows.size() == 2);
    CHECK(table.rows[0][2] == "say \"hello\"");
    CHECK(table.rows[1][1] == "대조군");
    CHECK(table.rows[1][2].empty());
    CHECK_NOTHROW(rnaseq::validate_samples(table));
    CHECK(parse("x\ny").rows[0][0] == "y");
}

TEST_CASE("streaming TSV records are consumed lazily") {
    std::istringstream input("gene_id\tcount\nvalid\t1\n\"unterminated\t2\n");
    rnaseq::TsvReader reader(input, "lazy.tsv");
    std::vector<std::string> row;

    CHECK(reader.record() == 0);
    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"gene_id", "count"});
    CHECK(reader.record() == 1);
    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"valid", "1"});
    CHECK(reader.record() == 2);
    CHECK_THROWS_WITH(reader.next(row), doctest::Contains("lazy.tsv: record 3, column 1"));
    CHECK(reader.record() == 3);
}

TEST_CASE("streaming TSV records preserve strict decoded cells") {
    std::istringstream input("\xEF\xBB\xBFgene_id\tnote\tcount\r\n"
                             "유전자\t\"say \"\"hello\"\"\"\t0\r\n");
    rnaseq::TsvReader reader(input, "stream.tsv");
    std::vector<std::string> row;

    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"gene_id", "note", "count"});
    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"유전자", "say \"hello\"", "0"});
    CHECK_FALSE(reader.next(row));
    CHECK(reader.record() == 2);
}

TEST_CASE("streaming TSV reader supports headerless records") {
    std::istringstream input("gene-a\t4\ngene-b\t0\n");
    rnaseq::TsvReader reader(input, "legacy.tsv");
    std::vector<std::string> row;

    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"gene-a", "4"});
    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"gene-b", "0"});
    CHECK_FALSE(reader.next(row));
}

TEST_CASE("featureCounts comments are allowed only before data") {
    std::istringstream input("\xEF\xBB\xBF# Program:featureCounts\n"
                             "# Command:featureCounts -a genes.gtf\n"
                             "Geneid\tcount\n"
                             "#late\t1\n");
    rnaseq::TsvReader reader(input, "featurecounts.txt", true);
    std::vector<std::string> row;

    REQUIRE(reader.next(row));
    CHECK(row == std::vector<std::string>{"Geneid", "count"});
    CHECK(reader.record() == 3);
    CHECK_THROWS_WITH(reader.next(row),
                      doctest::Contains("featurecounts.txt: record 4, column 1: comments are forbidden"));
}

TEST_CASE("streaming TSV reader rejects late ragged records") {
    std::istringstream input("gene_id\tcount\ngene-a\t1\ngene-b\n");
    rnaseq::TsvReader reader(input, "ragged.tsv");
    std::vector<std::string> row;

    REQUIRE(reader.next(row));
    REQUIRE(reader.next(row));
    CHECK_THROWS_WITH(reader.next(row),
                      doctest::Contains("ragged.tsv: record 3, column 1: ragged row"));
}

TEST_CASE("malformed records never disappear") {
    for (const std::string text : {"", "x\tx\n1\t2\n", "x\t\n1\t2\n",
        "x\ty\n1\n", "x\ty\n1\t2\t3\n", "x\ny\n\n", "x\n#comment\n",
        "x\n\"unclosed\n", "x\n\"closed\"suffix\n", "x\nba\"re\n",
        "x\n\"a\tb\"\n", "x\n\"a\nb\"\n", "x\ny\r", "x\ny\rz\n",
        "x\n\xEF\xBB\xBFy\n"}) {
        CAPTURE(text);
        CHECK_THROWS_AS(parse(text), std::runtime_error);
    }
    CHECK_THROWS_AS(parse(std::string("x\ny\0z\n", 6)), std::runtime_error);
    CHECK_THROWS_AS(parse("x\n\xC0\xAF\n"), std::runtime_error);
    CHECK_THROWS_AS(parse("x\n\xED\xA0\x80\n"), std::runtime_error);
    CHECK_THROWS_AS(parse("x\n\xF4\x90\x80\x80\n"), std::runtime_error);
    CHECK_THROWS_AS(parse("x\n\xE2\x82\n"), std::runtime_error);
}

TEST_CASE("sample IDs and required cells") {
    for (const std::string text : {"sample_id\tcondition\n", "sample_id\tother\ns\tc\n",
        "sample_id\tcondition\ns\tc\ns\tt\n", "sample_id\tcondition\n..\tc\n",
        "sample_id\tcondition\n-s\tc\n", "sample_id\tcondition\ns x\tc\n",
        "sample_id\tcondition\ns\t\n"}) {
        CAPTURE(text);
        CHECK_THROWS_AS(rnaseq::validate_samples(parse(text)), std::runtime_error);
    }
    CHECK_NOTHROW(rnaseq::validate_samples(parse("sample_id\tcondition\ns.1-a_b\tc\n")));
}

TEST_CASE("strict unsigned counts and checked aggregation") {
    CHECK(rnaseq::count_value("18446744073709551615", "test") == std::numeric_limits<std::uint64_t>::max());
    CHECK(rnaseq::count_value("0", "test") == 0);
    for (const std::string value : {"", "-1", "+1", "1.0", "1e3", " 1", "1 ",
                                     "18446744073709551616", "NaN"}) {
        CAPTURE(value);
        CHECK_THROWS_AS(rnaseq::count_value(value, "test"), std::runtime_error);
    }
    CHECK(rnaseq::checked_add(2, 4) == 6);
    CHECK(rnaseq::checked_add(std::numeric_limits<std::uint64_t>::max(), 0) ==
          std::numeric_limits<std::uint64_t>::max());
    CHECK_THROWS_AS(rnaseq::checked_add(std::numeric_limits<std::uint64_t>::max(), 1), std::overflow_error);
}

TEST_CASE("count schemas preserve identifiers and reject duplicates") {
    CHECK_NOTHROW(rnaseq::validate_counts(parse("gene_id\ts2\ts1\ngene.2\t0\t4\n")));
    for (const std::string text : {"Geneid\ts\ng\t1\n", "gene_id\ng\n", "gene_id\ts\n",
        "gene_id\ts\ng\t1\ng\t2\n", "gene_id\ts\n\t1\n", "gene_id\ts\ng\t-1\n"}) {
        CAPTURE(text);
        CHECK_THROWS_AS(rnaseq::validate_counts(parse(text)), std::runtime_error);
    }
}

TEST_CASE("actionable errors include source, record and column") {
    CHECK_THROWS_WITH(rnaseq::validate_counts(parse("gene_id\ts\ng\t1.5\n")),
                     doctest::Contains("fixture.tsv: record 2, column s"));
    CHECK_THROWS_WITH(parse("x\ty\nz\n"), doctest::Contains("fixture.tsv: record 2, column 1"));
    std::istringstream broken("x\ny\n");
    broken.setstate(std::ios::badbit);
    CHECK_THROWS_WITH(rnaseq::read_table(broken, "broken.tsv"), doctest::Contains("read failure"));
}
