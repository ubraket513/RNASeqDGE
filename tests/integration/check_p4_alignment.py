#!/usr/bin/env python3
"""Real-tool gate. Expectations follow P0 read origins, never captured outputs."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'tests/fixtures/p0'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, default=ROOT / 'build/rnaseq')
    parser.add_argument('--bin-dir', type=Path, default=ROOT / '.deps/alignment/bin')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    out = args.output or Path(tempfile.mkdtemp(prefix='p4-real-', dir=ROOT / 'tests/output'))
    out.mkdir(parents=True, exist_ok=True)
    exe, tools = args.exe.resolve(), args.bin_dir.resolve()
    common = ['--bin-dir', tools, '--fasta', FIXTURE / 'reference.fa', '--gtf', FIXTURE / 'reference.gtf']
    evidence = []

    def run(argv):
        subprocess.run(list(map(str, argv)), check=True, timeout=180)

    # Two plus-strand gene_A reads (one junction), one minus-strand gene_B read,
    # one overlap read excluded as ambiguous. One paired gene_A fragment.
    expected = {'single': {'gene_A.1': 2, 'gene_B.2': 1, 'gene_C': 0, 'gene_zero': 0},
                'paired': {'gene_A.1': 1, 'gene_B.2': 0, 'gene_C': 0, 'gene_zero': 0}}
    for backend in ('hisat2', 'star'):
        for threads in (1, 2):
            index = out / f'{backend}-{threads}-index'
            run([exe, 'index', '--backend', backend, *common, '--threads', threads, '--output', index,
                 *(['--star-sa-bases', 3, '--star-chr-bits', 10] if backend == 'star' else [])])
            for layout in ('single', 'paired'):
                dest = out / f'{backend}-{threads}-{layout}'
                reads = ['--reads1', FIXTURE / ('reads/single.fastq' if layout == 'single' else 'reads/paired_1.fastq')]
                if layout == 'paired':
                    reads += ['--reads2', FIXTURE / 'reads/paired_2.fastq']
                run([exe, 'align-count', '--backend', backend, *common, '--threads', threads,
                     '--index', index, '--layout', layout, '--strandedness', 'unstranded',
                     *reads, '--output', dest])
                with (dest / 'counts.txt').open() as stream:
                    rows = list(csv.reader((line for line in stream if not line.startswith('#')), delimiter='\t'))
                counts = {row[0]: int(row[-1]) for row in rows[1:]}
                assert counts == expected[layout], (backend, threads, layout, counts, expected[layout])
                with (dest / 'counts.txt.summary').open() as stream:
                    summary = {row[0]: int(row[1]) for row in list(csv.reader(stream, delimiter='\t'))[1:]}
                assert summary['Assigned'] == (3 if layout == 'single' else 1), summary
                assert summary['Unassigned_Ambiguity'] == (1 if layout == 'single' else 0), summary
                assert sum(summary.values()) == (4 if layout == 'single' else 1), summary
                records = subprocess.check_output([tools / 'samtools', 'view', dest / 'aligned.bam'], text=True)
                if layout == 'single':
                    junction = [line.split('\t') for line in records.splitlines() if line.split('\t')[0] == 'junction_plus']
                    assert len(junction) == 1 and junction[0][5] == '37M200N38M', junction
                assert (dest / 'COMPLETE').is_file()
                evidence.append(dict(backend=backend, threads=threads, layout=layout, counts=counts, summary=summary))
                print(f'PASS {backend} threads={threads} {layout}: counts, summary, junction, BAM')
    (out / 'verified.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print(f'Evidence: {out}')


if __name__ == '__main__':
    main()
