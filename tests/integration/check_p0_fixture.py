#!/usr/bin/env python3
"""Independent, stdlib-only arithmetic check of the hand-authored P0 oracle.

This is fixture verification, not the production TSV validator or count adapter.
"""
import csv
import hashlib
from pathlib import Path
import unittest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "p0"


def rows(name):
    with (FIXTURE / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class FixtureChecks(unittest.TestCase):
    def test_sequence_manifests(self):
        references = rows("references.tsv")
        self.assertEqual([row["role"] for row in references], ["genome", "annotation"])
        for row in references:
            self.assertEqual(hashlib.sha256((FIXTURE / row["path"]).read_bytes()).hexdigest(),
                             row["sha256"])
        runs = rows("runs.tsv")
        self.assertEqual([(r["run_id"], r["sample_id"]) for r in runs],
                         [(r["run_id"], r["sample_id"]) for r in rows("count_inputs.tsv")])
        by_sample = {}
        for run in runs:
            props = (run["layout"], run["strandedness"])
            self.assertEqual(by_sample.setdefault(run["sample_id"], props), props)
            self.assertEqual(bool(run["fastq_2"]), run["layout"] == "paired")
            for key in ("fastq_1", "fastq_2"):
                if not run[key]:
                    continue
                lines = (FIXTURE / run[key]).read_text().splitlines()
                self.assertEqual(len(lines) % 4, 0)
                for i in range(0, len(lines), 4):
                    self.assertTrue(lines[i].startswith("@"))
                    self.assertEqual(lines[i+2], "+")
                    self.assertEqual(len(lines[i+1]), 75)
                    self.assertEqual(len(lines[i+1]), len(lines[i+3]))
        sequence = "".join((FIXTURE / "reference.fa").read_text().splitlines()[1:])
        reads = {}
        for path in sorted((FIXTURE / "reads").glob("*.fastq")):
            lines = path.read_text().splitlines()
            reads.update((lines[i][1:], lines[i+1]) for i in range(0, len(lines), 4))
        origins = rows("read_origins.tsv")
        self.assertEqual(set(reads), {r["read_id"] for r in origins})
        for origin in origins:
            expected = ""
            for interval in origin["intervals_1based_inclusive"].split(","):
                start, end = map(int, interval.split("-"))
                expected += sequence[start-1:end]
            if origin["orientation"] == "-":
                expected = expected.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            self.assertEqual(reads[origin["read_id"]], expected)
        self.assertEqual(sequence[300:302], "GT")
        self.assertEqual(sequence[498:500], "AG")

    def test_hand_calculated_counts(self):
        samples = [row["sample_id"] for row in rows("samples.tsv")]
        self.assertEqual(len(samples), len(set(samples)))
        genes = [row["gene_id"] for row in rows("annotation.tsv")]
        self.assertEqual(len(genes), len(set(genes)))
        totals = {gene: dict.fromkeys(samples, 0) for gene in genes}
        inputs = rows("count_inputs.tsv")
        self.assertEqual(len(inputs), len({row["run_id"] for row in inputs}))
        self.assertEqual(set(samples), {row["sample_id"] for row in inputs})
        for run in inputs:
            counts = rows(run["path"])
            self.assertEqual(len(counts), len(genes))
            self.assertEqual({row["gene_id"] for row in counts}, set(genes))
            for row in counts:
                self.assertRegex(row["count"], r"^[0-9]+$")
                totals[row["gene_id"]][run["sample_id"]] += int(row["count"])
        expected = rows("expected/counts.tsv")
        self.assertEqual(list(expected[0]), ["gene_id", *samples])
        self.assertEqual([row["gene_id"] for row in expected], sorted(genes))
        for row in expected:
            self.assertEqual({sample: int(row[sample]) for sample in samples},
                             totals[row["gene_id"]])

    def test_explicit_comparison_and_annotation(self):
        samples = rows("samples.tsv")
        self.assertEqual([row["condition"] for row in samples].count("control"), 2)
        self.assertEqual([row["condition"] for row in samples].count("treated"), 2)
        self.assertEqual(rows("contrasts.tsv"), [{
            "contrast_id": "treated_vs_control", "factor": "condition",
            "numerator": "treated", "denominator": "control"}])
        annotation = rows("annotation.tsv")
        self.assertEqual(annotation[0]["gene_symbol"], annotation[1]["gene_symbol"])
        self.assertEqual(annotation[2]["gene_symbol"], "")
        self.assertEqual(dict((r["key"], r["value"]) for r in rows("analysis.tsv")),
                         {"version": "1", "design_terms": "condition", "alpha": "0.05",
                          "filter": "zero_total", "shrinkage": "apeglm"})


if __name__ == "__main__":
    unittest.main()
