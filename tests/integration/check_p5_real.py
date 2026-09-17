#!/usr/bin/env python3
"""Real P0 reads through P5 alignment/merge, expecting the genuine tiny-data R failure."""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile

repo = Path(__file__).resolve().parents[2]
work = Path(tempfile.mkdtemp(prefix="p5-real-", dir=repo / "tests/output"))
config = work / "workflow.tsv"
settings = {
    "samples": repo / "tests/fixtures/p0/samples.tsv",
    "runs": repo / "tests/fixtures/p0/runs.tsv",
    "references": repo / "tests/fixtures/p0/references.tsv",
    "analysis": repo / "tests/fixtures/p0/analysis.tsv",
    "contrasts": repo / "tests/fixtures/p0/contrasts.tsv",
    "genes": repo / "tests/fixtures/p2/genes.tsv",
    "bin_dir": repo / ".deps/alignment/bin",
    "rscript": repo / ".deps/p0/bin/Rscript",
    "r_script": repo / "run_deg_analysis_offline.R",
    "tool_lock": repo / "config/alignment-linux-64.explicit.txt",
    "r_lock": repo / "config/p0-linux-64.explicit.txt",
    "run_dir": work / "run",
}
config.write_text("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in settings.items()))
with (work / "workflow.log").open("w") as log:
    result = subprocess.run([str(repo / "build/rnaseq"), "workflow-local", str(config)],
                            stdout=log, stderr=subprocess.STDOUT, env={**os.environ, "PATH": "/usr/bin:/bin"})
assert result.returncode != 0, "tiny fixture must not be represented as valid DESeq2 analysis"
run = work / "run"
assert (run / "merge/STAGE.tsv").is_file()
assert not (run / "analysis/STAGE.tsv").exists()
with (run / "merge/counts.tsv").open() as handle:
    counts = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
for sample in ("treated_2", "control_1", "treated_1", "control_2"):
    assert int(counts["gene_A.1"][sample]) == 2
    assert int(counts["gene_B.2"][sample]) == (0 if sample == "control_1" else 1)
    assert int(counts["gene_C"][sample]) == 0
    assert int(counts["gene_zero"][sample]) == 0
logs = [p for p in run.glob("R-*.log") if not p.name.startswith("R-preflight-")]
assert len(logs) == 1
message = logs[0].read_text()
assert "dispersion" in message or "too few" in message or "fewer" in message, message
(work / "verified.json").write_text(json.dumps({
    "alignment": "HISAT2 2.2.1, samtools 1.17, featureCounts 2.0.6",
    "run_count": 5, "sample_count": 4,
    "counts": "independent P0 SE and PE expectations passed",
    "R": "real tiny-data failure retained; no analysis success manifest",
    "exit_status": result.returncode,
}, indent=2) + "\n")
print(f"PASS real P5 alignment/count merge and honest R failure: {work}")
