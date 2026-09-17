#!/usr/bin/env python3
"""Capture bounded legacy evidence without downloads, writes to inputs, or jobs.

Extracted expressions are evaluated by their actual interpreters. This does not
claim a full Snakemake or DESeq2 run. Outputs belong outside golden fixtures.
"""
import ast
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(argv, stdin=None):
    try:
        result = subprocess.run(argv, input=stdin, cwd=ROOT, text=True,
                                capture_output=True, timeout=30)
        return {"argv": argv, "exit_code": result.returncode,
                "stdout": result.stdout, "stderr": result.stderr}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"argv": argv, "error": str(exc)}


def main():
    report = {"scope": "Legacy startup and isolated expressions; no full pipeline run"}
    report["merge_startup"] = run([sys.executable, "merge_transcripts.py", "--help"])
    # Reconstruct the Python string used for the counting shell block. Formatting
    # only resolves escaped braces; no workflow or external tool is executed.
    source = (ROOT / "pipeline.smk").read_text()
    block = source.split("rule count_reads:", 1)[1].split("shell:", 1)[1]
    literal = block[block.index('"""'):]
    literal = literal[:literal.index('"""', 3) + 3]
    shell = ast.literal_eval(literal)
    awk_line = next(line for line in shell.splitlines() if line.strip().startswith("awk "))
    program = awk_line.split("'", 2)[1].replace("{{", "{").replace("}}", "}")
    report["count_delimiter"] = run(["awk", program], "# fixture\nGeneid Count\ngene_A.1 7\n")
    report["count_delimiter"]["matches_tab_oracle"] = (
        report["count_delimiter"].get("stdout") == "gene_A.1\t7\n")
    if shutil.which("Rscript"):
        report["condition_expression"] = run(["Rscript", "--vanilla", "-e", '''
x <- as.list(parse("run_deg_analysis.R"))
e <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
            identical(e[[2]], quote(coldata$condition)), x)[[1]]
condition_raw <- c("control-2", "bipolar-4")
coldata <- data.frame(title=condition_raw)
eval(e)
stopifnot(identical(as.character(coldata$condition), c("control", "bipolar")))
cat("labels:", as.character(coldata$condition), "\\n")
cat("alphabetical reference:", levels(coldata$condition)[1], "\\n")
'''])
    else:
        report["condition_expression"] = {"skipped": "Rscript unavailable"}
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
