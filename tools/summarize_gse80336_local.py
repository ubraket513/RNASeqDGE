#!/usr/bin/env python3
"""Summarize completed reduced workflows without claiming whole-genome parity."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def correlation(x, y):
    if len(x) < 2 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    return statistics.correlation(x, y)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--expected-repeats", type=int, default=3)
    args = parser.parse_args()
    root = args.directory
    records = json.loads((root / "runs.json").read_text())
    wanted = records.get("requested_replicates", args.expected_repeats)
    expected = {(b, n) for b in ("hisat2", "star") for n in range(1, wanted + 1)}
    observed = [(r["backend"], r["repeat"]) for r in records["runs"]]
    if wanted < 1 or set(observed) != expected or len(observed) != len(expected):
        raise SystemExit("incomplete or duplicate requested backend/repetition pairs")
    records["requested_replicates"] = wanted
    records["complete"] = True
    for run in records["runs"]:
        if run["status"] != 0:
            raise SystemExit("refusing completed comparison with failed workflows")
        name = f"{run['backend']}-{run['repeat']}"
        path = root / name
        run["run_directory_bytes_excluding_shared_index"] = sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
        run["resource_usage_by_attempt"] = {}
        for filename in run.get("timing_files", [f"{name}.time.txt"]):
            usage = {}
            for line in (root / filename).read_text().splitlines():
                if ": " in line:
                    key, value = line.strip().split(": ", 1)
                    usage[key] = value
            run["resource_usage_by_attempt"][filename] = usage
        run["assignment"] = {}
        for summary in sorted((path / "align").glob("*/counts.txt.summary")):
            values = table(summary)
            totals = {r["Status"]: int(list(r.values())[1]) for r in values}
            run["assignment"][summary.parent.name] = totals
        run["analysis"] = {r["key"]: r["value"] for r in table(
            path / "analysis/contrasts/bipolar_vs_control/summary.tsv")}
    data = {backend: table(root / f"{backend}-1/merge/counts.tsv") for backend in ("hisat2", "star")}
    h, s = data["hisat2"], data["star"]
    if [r["gene_id"] for r in h] != [r["gene_id"] for r in s] or list(h[0]) != list(s[0]):
        raise SystemExit("incomparable gene/sample identity or order")
    ids = list(h[0])[1:]
    hv = [int(row[col]) for row in h for col in ids]
    sv = [int(row[col]) for row in s for col in ids]
    records["counts_comparison"] = dict(genes=len(h), samples=len(ids),
        differing_cells=sum(a != b for a, b in zip(hv, sv)),
        total_absolute_difference=sum(abs(a - b) for a, b in zip(hv, sv)),
        hisat2_total=sum(hv), star_total=sum(sv), pearson=correlation(hv, sv))
    for backend, total in (("hisat2", sum(hv)), ("star", sum(sv))):
        first = next(r for r in records["runs"] if r["backend"] == backend and r["repeat"] == 1)
        if sum(item["Assigned"] for item in first["assignment"].values()) != total:
            raise SystemExit("merged counts disagree with featureCounts assigned totals")
    differences = [(sum(abs(int(a[c]) - int(b[c])) for c in ids), a["gene_id"])
                   for a, b in zip(h, s)]
    records["counts_comparison"]["largest_gene_differences"] = [
        dict(gene_id=gene, summed_absolute_difference=delta)
        for delta, gene in sorted(differences, reverse=True)[:10]]
    # Mapping/junction counts use primary records; assignment summaries above
    # may count multiple alignments and therefore use a different denominator.
    config = {r["key"]: r["value"] for r in table(root / "hisat2-1.tsv")}
    source = json.loads((Path(config["runs"]).parent / "provenance.json").read_text())
    read_totals = {r["run_id"]: r["prefix_records"] for r in source["reads"]}
    records["mapping_primary_records"] = {}
    samtools = Path(config["bin_dir"]) / "samtools"
    for backend in ("hisat2", "star"):
        records["mapping_primary_records"][backend] = {}
        for run_id, total in read_totals.items():
            metrics = dict(input_reads=total, mapped=0, unique=0, multimapped=0, spliced=0)
            process = subprocess.Popen([str(samtools), "view", str(root / f"{backend}-1/align" / run_id / "aligned.bam")],
                                       stdout=subprocess.PIPE, text=True)
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                flag = int(fields[1])
                if flag & (4 | 256 | 2048):
                    continue
                metrics["mapped"] += 1
                nh = [int(tag[5:]) for tag in fields[11:] if tag.startswith("NH:i:")]
                if len(nh) != 1:
                    process.kill(); process.wait()
                    raise SystemExit("missing/duplicate NH tag in mapped primary record")
                metrics["unique" if nh[0] == 1 else "multimapped"] += 1
                metrics["spliced"] += "N" in fields[5]
            if process.wait() != 0 or metrics["mapped"] > total:
                raise SystemExit("invalid BAM mapping summary")
            metrics["mapping_rate"] = metrics["mapped"] / total
            records["mapping_primary_records"][backend][run_id] = metrics
    result = {backend: {r["gene_id"]: r for r in table(root / f"{backend}-1/analysis/contrasts/bipolar_vs_control/results.tsv")}
              for backend in ("hisat2", "star")}
    genes = sorted(set(result["hisat2"]) & set(result["star"]))
    finite = lambda x: x != "NA" and math.isfinite(float(x))
    genes = [g for g in genes if all(finite(result[b][g]["log2FoldChange"]) for b in result)]
    lfc = {b: [float(result[b][g]["log2FoldChange"]) for g in genes] for b in result}
    deg = {b: {g for g, r in rows.items() if finite(r["padj"]) and float(r["padj"]) < .05}
           for b, rows in result.items()}
    union = deg["hisat2"] | deg["star"]
    records["statistics_comparison"] = dict(common_finite_lfc_genes=len(genes),
        lfc_pearson=correlation(lfc["hisat2"], lfc["star"]),
        median_absolute_lfc_difference=statistics.median(abs(a - b) for a, b in zip(lfc["hisat2"], lfc["star"])) if genes else None,
        lfc_sign_disagreements=sum((a > 0) != (b > 0) for a, b in zip(lfc["hisat2"], lfc["star"]) if a and b),
        deg={b: sorted(v) for b, v in deg.items()},
        deg_jaccard=len(deg["hisat2"] & deg["star"]) / len(union) if union else None)
    records["repeat_equality"] = {}
    for backend in ("hisat2", "star"):
        repeats = [r["repeat"] for r in records["runs"] if r["backend"] == backend]
        records["repeat_equality"][backend] = {}
        for artifact in ("merge/counts.tsv", "analysis/contrasts/bipolar_vs_control/results.tsv",
                         "analysis/contrasts/bipolar_vs_control/shrunken.tsv"):
            hashes = [hashlib.sha256((root / f"{backend}-{n}" / artifact).read_bytes()).hexdigest() for n in repeats]
            records["repeat_equality"][backend][artifact] = dict(repeats=len(repeats),
                identical=len(set(hashes)) == 1 if len(repeats) > 1 else None,
                assessment="checked" if len(repeats) > 1 else "not assessed: single completed run")
    records["limitations"] = ["Restricted chromosome reference can misassign reads from omitted chromosomes.",
        "Prefix sampling is not random and preserves sequencing-order bias.",
        "Low-depth DEG/LFC differences are diagnostic, not biological findings.",
        "OS cache uncontrolled; first run includes indexing, later runs reuse index.",
        "Timing includes provenance hashes and R; no independent per-stage resource profiler."]
    (root / "comparison.json").write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: records[k] for k in ("counts_comparison", "statistics_comparison", "repeat_equality")}, indent=2))


if __name__ == "__main__":
    main()
