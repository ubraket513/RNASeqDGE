#!/usr/bin/env python3
"""Stage checked local GEO/ENA downloads into explicit study-analysis inputs."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import tempfile

SUBSET = ("GSM2124739", "GSM2124742", "GSM2124744", "GSM2124760", "GSM2124763", "GSM2124766")
ROOT = Path(__file__).resolve().parents[1]


def metadata(text):
    samples = {}
    current = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            current = line.split(" = ", 1)[1]
            if current in samples:
                raise ValueError("duplicate GEO sample")
            samples[current] = {"sample_id": current}
        elif current and " = " in line:
            key, value = line.split(" = ", 1)
            if key == "!Sample_title":
                match = re.fullmatch(r"(control|bipolar)-(\d+)", value)
                if not match:
                    raise ValueError(f"unexpected study title: {value}")
                samples[current].update(condition=match[1], title=value,
                    count_column=("C_" if match[1] == "control" else "BD_") + match[2])
            elif key == "!Sample_relation" and "SRX" in value:
                samples[current]["experiment"] = re.search(r"SRX\d+", value)[0]
            elif key == "!Sample_characteristics_ch1":
                label, item = value.split(": ", 1)
                column = {"age (years)": "age", "Sex": "sex", "postmortem interval (hours)": "pmi", "rin": "rin"}[label]
                samples[current][column] = item
    if len(samples) != 36 or "PRJNA318642" not in text or "SRP073382" not in text:
        raise ValueError("wrong study or incomplete GEO sample set")
    for sample in samples.values():
        if set(sample) != {"sample_id", "condition", "title", "count_column", "experiment", "age", "sex", "pmi", "rin"}:
            raise ValueError("incomplete GEO sample metadata")
    return samples


def write_tsv(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def stage(source, destination):
    if destination.exists() or destination.is_symlink() or not destination.parent.is_dir():
        raise ValueError("output must be absent with an existing parent")
    source_info = json.loads((ROOT / "config/gse80336-sources.json").read_text())
    for name, item in source_info["files"].items():
        actual = hashlib.sha256((source / name).read_bytes()).hexdigest()
        if actual != item["sha256"]:
            raise ValueError(f"source checksum mismatch: {name}")
    text = gzip.decompress((source / "GSE80336_family.soft.gz").read_bytes()).decode("utf-8")
    samples = metadata(text)
    with (source / "ena-runs.tsv").open() as stream:
        ena = list(csv.DictReader(stream, delimiter="\t"))
    experiments = {r["experiment_accession"]: r for r in ena}
    if len(experiments) != 36 or len(ena) != 36 or set(experiments) != {s["experiment"] for s in samples.values()}:
        raise ValueError("GEO/ENA mapping is not one-to-one for all 36 samples")
    with gzip.open(source / "GSE80336_Counts.txt.gz", "rt") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        expected = ["Ensembl_ID", "GeneSymbol", "Biotype", "Chromosome"]
        if reader.fieldnames[:4] != expected or set(reader.fieldnames[4:]) != {s["count_column"] for s in samples.values()}:
            raise ValueError("supplementary count columns do not match GEO cohort")
        rows = sorted(reader, key=lambda row: row["Ensembl_ID"].encode("utf-8"))
    if len(rows) != 47886 or len({r["Ensembl_ID"] for r in rows}) != len(rows):
        raise ValueError("unexpected or duplicate gene identities")
    columns = [s["count_column"] for s in samples.values()]
    if any(not re.fullmatch(r"[0-9]+", r[c]) or int(r[c]) > 2147483647 for r in rows for c in columns):
        raise ValueError("counts must be DESeq2-compatible raw integers")
    tmp = Path(tempfile.mkdtemp(prefix=".gse80336-", dir=destination.parent))
    # Keep failed stage for diagnosis; no successful output is published on errors.
    mapping = []
    for sid, s in sorted(samples.items()):
        run = experiments[s["experiment"]]
        if run["study_accession"] != "PRJNA318642" or run["library_layout"] != "SINGLE" or run["library_strategy"] != "RNA-Seq":
            raise ValueError("unexpected ENA study/layout/strategy")
        mapping.append([sid, s["title"], s["count_column"], s["experiment"], run["run_accession"],
                        "single", "reverse", run["fastq_ftp"], run["fastq_md5"], run["fastq_bytes"], run["read_count"]])
    write_tsv(tmp / "sample_run_mapping.tsv",
        ["sample_id", "title", "count_column", "experiment_id", "run_id", "layout", "strandedness", "fastq_ftp", "fastq_md5", "fastq_bytes", "read_count"], mapping)
    cohorts = {"full36": sorted(samples), "without_C28": sorted(set(samples) - {"GSM2124750"}), "subset6": list(SUBSET)}
    for name, ids in cohorts.items():
        out = tmp / name
        out.mkdir()
        if {samples[s]["condition"] for s in ids} != {"control", "bipolar"}:
            raise ValueError("cohort requires both conditions")
        header = ["sample_id", "condition", "age", "sex", "pmi", "rin"]
        write_tsv(out / "samples.tsv", header, [[samples[s][h] for h in header] for s in ids])
        write_tsv(out / "counts.tsv", ["gene_id", *ids],
                  [[r["Ensembl_ID"], *[r[samples[s]["count_column"]] for s in ids]] for r in rows])
        write_tsv(out / "annotation.tsv", ["gene_id", "gene_symbol", "biotype", "chromosome"],
                  [[r[c] for c in expected] for r in rows])
        write_tsv(out / "analysis.tsv", ["key", "value"],
                  [("version", "1"), ("design_terms", "condition"), ("alpha", "0.05"), ("filter", "zero_total"),
                   ("shrinkage", "apeglm"), ("type.age", "numeric"), ("type.sex", "categorical"),
                   ("type.pmi", "numeric"), ("type.rin", "numeric")])
        write_tsv(out / "contrasts.tsv", ["contrast_id", "factor", "numerator", "denominator"],
                  [("bipolar_vs_control", "condition", "bipolar", "control")])
    report = dict(source_info, genes=len(rows), cohort_sizes={k: len(v) for k, v in cohorts.items()},
                  compressed_fastq_bytes=sum(int(r["fastq_bytes"]) for r in ena),
                  design="~ condition; covariates retained as metadata, not included in report-model replication",
                  source_count_scope="Published TopHat/HTSeq processed counts, including upstream symbol collapse; not new reference-gene featureCounts output",
                  subset_selection="Three per condition, selected explicitly with similar ages and matched sex distribution, not first N runs",
                  direction="bipolar/control; reverse of displayed report results contrast")
    (tmp / "provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    tmp.rename(destination)
    print(f"Staged {len(rows)} genes and {len(samples)} mapped samples: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    stage(args.source, args.out)
