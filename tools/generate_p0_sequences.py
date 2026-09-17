#!/usr/bin/env python3
"""Recreate synthetic P0 sequence inputs; never generates expected count tables."""
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "p0"


def main():
    # SHA-256 counter bytes give a deterministic sequence independent of PRNG version.
    sequence = "".join("ACGT"[byte % 4] for i in range(100)
                       for byte in hashlib.sha256(f"RNASeqDGE-P0-{i}".encode()).digest())
    sequence = list(sequence[:3000])
    sequence[300:302] = "GT"
    sequence[498:500] = "AG"
    sequence = "".join(sequence)
    (ROOT / "reference.fa").write_text(
        ">chrSynthetic\n" + "\n".join(sequence[i:i+60] for i in range(0, len(sequence), 60)) + "\n")
    exons = [("gene_A.1", "+", 101, 300), ("gene_A.1", "+", 501, 700),
             ("gene_B.2", "-", 1001, 1300), ("gene_C", "+", 1101, 1400),
             ("gene_zero", "+", 2001, 2300)]
    (ROOT / "reference.gtf").write_text("".join(
        f'chrSynthetic\tP0\texon\t{start}\t{end}\t.\t{strand}\t.\t'
        f'gene_id "{gene}"; transcript_id "{gene}.tx";\n'
        for gene, strand, start, end in exons))
    reads_dir = ROOT / "reads"
    reads_dir.mkdir(exist_ok=True)

    def reverse_complement(seq):
        return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def fastq(name, records):
        (reads_dir / name).write_text("".join(
            f"@{identifier}\n{seq}\n+\n{'I' * len(seq)}\n" for identifier, seq in records))

    # 1-based inclusive origins are recorded separately for manual review.
    fastq("single.fastq", [("exon_plus", sequence[150:225]),
                           ("junction_plus", sequence[263:300] + sequence[500:538]),
                           ("exon_minus", reverse_complement(sequence[1000:1075])),
                           ("ambiguous_overlap", sequence[1150:1225])])
    fastq("paired_1.fastq", [("fragment_A/1", sequence[150:225])])
    fastq("paired_2.fastq", [("fragment_A/2", reverse_complement(sequence[575:650]))])
    (ROOT / "read_origins.tsv").write_text(
        "read_id\tintervals_1based_inclusive\torientation\tcase\n"
        "exon_plus\t151-225\t+\texonic_gene_A.1\n"
        "junction_plus\t264-300,501-538\t+\tspliced_gene_A.1\n"
        "exon_minus\t1001-1075\t-\texonic_gene_B.2\n"
        "ambiguous_overlap\t1151-1225\t+\tgene_B.2_and_gene_C_unstranded\n"
        "fragment_A/1\t151-225\t+\tpaired_gene_A.1\n"
        "fragment_A/2\t576-650\t-\tpaired_gene_A.1\n")
    (ROOT / "runs.tsv").write_text(
        "run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\n" + "".join(
            f"{run}\t{sample}\t" + (
                "reads/paired_1.fastq\treads/paired_2.fastq\tpaired\tunstranded\n"
                if sample == "control_1" else "reads/single.fastq\t\tsingle\tunstranded\n")
            for run, sample in [("run_t2", "treated_2"), ("run_c1a", "control_1"),
                                ("run_c1b", "control_1"), ("run_t1", "treated_1"),
                                ("run_c2", "control_2")]))
    (ROOT / "references.tsv").write_text(
        "role\tpath\tsha256\tsource\trelease\n" + "".join(
            f"{role}\t{name}\t{hashlib.sha256((ROOT / name).read_bytes()).hexdigest()}"
            "\tsynthetic:tools/generate_p0_sequences.py\tp0-v1\n"
            for role, name in [("genome", "reference.fa"), ("annotation", "reference.gtf")]))


if __name__ == "__main__":
    main()
