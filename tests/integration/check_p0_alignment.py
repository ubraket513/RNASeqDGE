#!/usr/bin/env python3
"""Real HISAT2/featureCounts check with an independent read-origin oracle."""
import csv
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/fixtures/p0"
OUTPUT = ROOT / "tests/output/p0/alignment"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    commands = []

    def run(argv, stdout=None):
        commands.append(list(map(str, argv)))
        with (OUTPUT / "commands.json").open("w") as stream:
            json.dump(commands, stream, indent=2)
        with (OUTPUT / "tools.log").open("a") as log:
            if stdout:
                with Path(stdout).open("w") as out:
                    subprocess.run(argv, stdout=out, stderr=log, check=True, timeout=180)
            else:
                subprocess.run(argv, stdout=log, stderr=log, check=True, timeout=180)

    run(["hisat2-build", "-p", "1", str(FIXTURE / "reference.fa"), str(OUTPUT / "index")])
    # Use the legacy genome-only index and default discovery, not a changed
    # transcript-aware method. Junction anchors are 37/38 bases and GT/AG.
    expected = {"single": {"gene_A.1": 2, "gene_B.2": 1, "gene_C": 0, "gene_zero": 0},
                "paired": {"gene_A.1": 1, "gene_B.2": 0, "gene_C": 0, "gene_zero": 0}}
    for layout in ("single", "paired"):
        sam = OUTPUT / f"{layout}.sam"
        bam = OUTPUT / f"{layout}.bam"
        counts = OUTPUT / f"{layout}.counts.txt"
        reads = (["-U", str(FIXTURE / "reads/single.fastq")] if layout == "single" else
                 ["-1", str(FIXTURE / "reads/paired_1.fastq"),
                  "-2", str(FIXTURE / "reads/paired_2.fastq")])
        run(["hisat2", "-p", "1", "-x", str(OUTPUT / "index"), *reads, "-S", str(sam)])
        run(["samtools", "sort", "-@", "0", "-o", str(bam), str(sam)])
        run(["samtools", "quickcheck", str(bam)])
        run(["featureCounts", "-T", "1", "-s", "0", "-Q", "0", "-t", "exon", "-g", "gene_id",
             *(["-p", "--countReadPairs"] if layout == "paired" else []),
             "-a", str(FIXTURE / "reference.gtf"), "-o", str(counts), str(bam)])
        with counts.open() as stream:
            table = list(csv.reader((line for line in stream if not line.startswith("#")), delimiter="\t"))
        actual = {row[0]: int(row[-1]) for row in table[1:]}
        assert actual == expected[layout], (layout, actual, expected[layout])
        with Path(str(counts) + ".summary").open() as stream:
            summary = {row[0]: int(row[1]) for row in list(csv.reader(stream, delimiter="\t"))[1:]}
        assert summary["Assigned"] == (3 if layout == "single" else 1), summary
        assert summary["Unassigned_Ambiguity"] == (1 if layout == "single" else 0), summary
        if layout == "single":
            records = [line.split("\t") for line in sam.read_text().splitlines() if not line.startswith("@")]
            junction = [row for row in records if row[0] == "junction_plus"]
            assert len(junction) == 1 and junction[0][5] == "37M200N38M", junction
        print(f"PASS {layout}: counts, assignment summary, BAM integrity")
    (OUTPUT / "verified.json").write_text(json.dumps(expected, indent=2) + "\n")


if __name__ == "__main__":
    main()
