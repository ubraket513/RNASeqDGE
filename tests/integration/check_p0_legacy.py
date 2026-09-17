#!/usr/bin/env python3
"""Exercise the unchanged merger with local data; prohibit network access."""
from pathlib import Path
import socket
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pandas as pd
import merge_transcripts as legacy


class LegacyEvidence(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.old_connect = socket.socket.connect
        socket.socket.connect = lambda *args: self.fail("Unexpected network access")
        self.addCleanup(setattr, socket.socket, "connect", self.old_connect)

    def annotation(self, ids, symbols):
        return pd.DataFrame({"Ensembl_ID": ids, "GeneSymbol": symbols,
                             "Biotype": "", "Chromosome": "chrSynthetic"})

    def record_original(self, name, result):
        evidence = ROOT / "tests/output/p0/original"
        evidence.mkdir(parents=True, exist_ok=True)
        result.to_csv(evidence / f"{name}.csv", index=False)

    def test_gene_id_symbol_mismatch_loses_all_genes(self):
        (self.path / "run.csv").write_text("gene_A.1\t2\ngene_B.2\t3\n")
        result = legacy.build_count_matrix(str(self.path),
            self.annotation(["gene_A.1", "gene_B.2"], ["SHARED", "SHARED"]),
            {"run": "sample"})
        self.assertEqual(len(result), 0)  # observed defect, not desired behavior
        self.record_original("gene_id_mismatch", result)

    def test_missing_annotation_and_sample_are_silently_lost(self):
        (self.path / "run.csv").write_text("gene_A.1\t2\ngene_C\t7\n")
        (self.path / "unmapped.csv").write_text("gene_A.1\t9\ngene_C\t8\n")
        result = legacy.build_count_matrix(str(self.path),
            self.annotation(["gene_A.1"], ["gene_A.1"]), {"run": "sample"})
        self.assertEqual(result["Ensembl_ID"].tolist(), ["gene_A.1"])
        self.assertEqual(result.columns.tolist()[-1:], ["sample"])
        self.assertNotIn("unmapped", result.columns)
        self.record_original("missing_gene_and_sample", result)

    def test_technical_runs_are_not_aggregated(self):
        (self.path / "run_a.csv").write_text("gene_A.1\t2\n")
        (self.path / "run_b.csv").write_text("gene_A.1\t4\n")
        result = legacy.build_count_matrix(str(self.path),
            self.annotation(["gene_A.1"], ["gene_A.1"]),
            {"run_a": "sample", "run_b": "sample"})
        self.assertNotIn("sample", result.columns)
        self.assertEqual(result[["sample_x", "sample_y"]].values.tolist(), [[2, 4]])
        self.record_original("technical_runs", result)

    def test_reviewed_gene_id_merge_matches_independent_oracle(self):
        fixture = ROOT / "tests/fixtures/p0"
        samples = pd.read_csv(fixture / "samples.tsv", sep="\t").sample_id.tolist()
        mapping = pd.read_csv(fixture / "count_inputs.tsv", sep="\t")
        genes = sorted(pd.read_csv(fixture / "annotation.tsv", sep="\t").gene_id)
        # Reviewed reference implementation only; not a replacement production CLI.
        output = pd.DataFrame(0, index=pd.Index(genes, name="gene_id"), columns=samples)
        for row in mapping.itertuples(index=False):
            counts = pd.read_csv(fixture / row.path, sep="\t", index_col="gene_id")
            self.assertEqual(set(counts.index), set(genes))
            self.assertTrue(counts.index.is_unique)
            output[row.sample_id] += counts["count"].reindex(genes)
        expected = pd.read_csv(fixture / "expected/counts.tsv", sep="\t", index_col="gene_id")
        pd.testing.assert_frame_equal(output, expected)
        evidence = ROOT / "tests/output/p0/reviewed"
        evidence.mkdir(parents=True, exist_ok=True)
        output.to_csv(evidence / "counts.tsv", sep="\t")


if __name__ == "__main__":
    unittest.main()
