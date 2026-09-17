#!/usr/bin/env python3
"""Stage an explicitly biased, restricted-reference GSE80336 smoke-test dataset.

Downloads are preparation only. Never presents prefix samples as random samples,
never claims full ENA MD5 verification for truncated gzip streams, and verifies
pinned full NCBI reference MD5s before extracting exact NC_000022.11.
"""
import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = 'NC_000022.11'
ASSEMBLY = 'GCF_000001405.40_GRCh38.p14'


def copy_fastq_prefix(source, destination, records, check_deadline):
    """Copy exactly N complete, structurally validated four-line FASTQ records."""
    for record in range(records):
        check_deadline()
        lines = [source.readline() for _ in range(4)]
        if any(not line or not line.endswith(b'\n') for line in lines):
            raise ValueError(f'incomplete FASTQ record {record + 1}')
        header, sequence, plus, qualities = [line.rstrip(b'\r\n') for line in lines]
        if not header.startswith(b'@') or not plus.startswith(b'+'):
            raise ValueError(f'invalid FASTQ delimiters at record {record + 1}')
        if not sequence or len(sequence) != len(qualities):
            raise ValueError(f'FASTQ sequence/quality length mismatch at record {record + 1}')
        if any(base not in b'ACGTNacgtn' for base in sequence):
            raise ValueError(f'invalid FASTQ sequence at record {record + 1}')
        if any(value < 33 or value > 126 for value in qualities):
            raise ValueError(f'invalid FASTQ quality at record {record + 1}')
        destination.write(b'\n'.join((header, sequence, plus, qualities)) + b'\n')


def digest(path, algorithm='sha256'):
    result = hashlib.new(algorithm)
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def write_tsv(path, fields, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'tests/output/p6-local-data')
    parser.add_argument('--samples', type=Path, default=ROOT / 'tests/output/p6-inputs/subset6/samples.tsv')
    parser.add_argument('--mapping', type=Path, default=ROOT / 'tests/output/p6-inputs/sample_run_mapping.tsv')
    parser.add_argument('--sources', type=Path, default=ROOT / 'config/gse80336-sources.json')
    parser.add_argument('--records', type=int, default=50000)
    parser.add_argument('--ena-scheme', choices=['https', 'http', 'ftp'], default='https',
                        help='Public ENA transport; explicit http/ftp fallback is recorded')
    parser.add_argument('--deadline', required=True, help='Absolute UTC deadline, e.g. 2026-09-17T13:25:00Z')
    parser.add_argument('--timeout', type=int, default=30, help='Network inactivity timeout seconds')
    parser.add_argument('--reuse-reference-downloads', action='store_true', help='Use only complete existing reference files, still verifying MD5')
    args = parser.parse_args()
    if args.records < 1 or args.timeout < 1:
        parser.error('records and timeout must be positive')
    deadline = datetime.fromisoformat(args.deadline.replace('Z', '+00:00')).timestamp()
    start = time.time()

    def check_deadline():
        if time.time() >= deadline:
            raise TimeoutError('staging deadline exceeded')

    def fetch(url):
        check_deadline()
        request = urllib.request.Request(url, headers={'User-Agent': 'RNASeqDGE-local-reduced-P6/1'})
        return urllib.request.urlopen(request, timeout=max(1, min(args.timeout, deadline - time.time())))

    out = args.output.resolve()
    downloads = out / 'downloads'
    reads_dir = out / 'reads'
    downloads.mkdir(parents=True, exist_ok=True)
    reads_dir.mkdir(exist_ok=True)
    if (out / 'STAGED').exists():
        raise ValueError('output already staged; choose an absent staging result')
    sources = json.loads(args.sources.read_text())
    if sources['reference']['assembly'] != ASSEMBLY:
        raise ValueError('reference assembly mismatch')
    with args.samples.open() as stream:
        samples = list(csv.DictReader(stream, delimiter='\t'))
    with args.mapping.open() as stream:
        mapping = {row['sample_id']: row for row in csv.DictReader(stream, delimiter='\t')}
    if len(samples) != 6 or len({row['sample_id'] for row in samples}) != 6:
        raise ValueError('exactly six distinct explicitly selected samples required')
    selected = [mapping[sample['sample_id']] for sample in samples]
    if any(row['layout'] != 'single' or row['strandedness'] != 'reverse' for row in selected):
        raise ValueError('selected mapping must be reverse single-end')
    provenance = {'study': 'GSE80336', 'assembly': ASSEMBLY, 'accession': ACCESSION,
                  'started_utc': datetime.fromtimestamp(start, timezone.utc).isoformat(),
                  'deadline_utc': args.deadline, 'records_per_run': args.records,
                  'scope': 'restricted-reference smoke check, not a full method benchmark',
                  'sampling': 'first N complete FASTQ records in archive order; nonrandom prefix sampling may be biased',
                  'source_reads_full_md5_verified': False,
                  'references': [], 'reads': []}

    def reference_file(name, info):
        path = downloads / name
        if not args.reuse_reference_downloads:
            temporary = path.with_suffix(path.suffix + '.partial')
            with fetch(info['url']) as response, temporary.open('wb') as stream:
                while block := response.read(1024 * 1024):
                    check_deadline()
                    stream.write(block)
            temporary.replace(path)
        actual = digest(path, 'md5')
        if actual != info['md5']:
            raise ValueError(f'NCBI full reference MD5 mismatch: {name}: {actual}')
        return {'file': name, 'url': info['url'], 'md5': actual, 'md5_verified': True,
                'compressed_bytes': path.stat().st_size}

    def read_prefix(row):
        url = args.ena_scheme + '://' + row['fastq_ftp']
        if ';' in url:
            raise ValueError('expected exactly one single-end FASTQ source')
        path = reads_dir / (row['run_id'] + '.fastq')
        temporary = path.with_suffix('.fastq.partial')
        begin = time.time()
        with fetch(url) as response, gzip.GzipFile(fileobj=response) as stream, temporary.open('wb') as target:
            copy_fastq_prefix(stream, target, args.records, check_deadline)
        temporary.replace(path)
        result = {'sample_id': row['sample_id'], 'run_id': row['run_id'], 'url': url,
                  'source_full_md5': row['fastq_md5'], 'source_full_md5_verified': False,
                  'source_full_gzip_crc_verified': False, 'source_total_records': int(row['read_count']),
                  'prefix_records': args.records, 'path': str(path.relative_to(out)),
                  'sha256': digest(path), 'bytes': path.stat().st_size,
                  'seconds': round(time.time() - begin, 3)}
        (reads_dir / (row['run_id'] + '.provenance.json')).write_text(json.dumps(result, indent=2) + '\n')
        print(f"STAGED {row['run_id']}: {args.records} records in {result['seconds']}s", flush=True)
        return result

    # Three bounded network workers keep reference preparation and read prefixes
    # moving without ever materializing full read files or whole references in RAM.
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        read_jobs = [pool.submit(read_prefix, row) for row in selected]
        reference_jobs = [pool.submit(reference_file, name, info)
                          for name, info in sources['reference']['files'].items()]
        provenance['reads'] = [job.result() for job in read_jobs]
        provenance['references'] = [job.result() for job in reference_jobs]
    fasta_source = downloads / (ASSEMBLY + '_genomic.fna.gz')
    gtf_source = downloads / (ASSEMBLY + '_genomic.gtf.gz')
    bases = 0
    found = 0
    active = False
    with gzip.open(fasta_source, 'rb') as stream, (out / 'reference.fa').open('wb') as target:
        for line in stream:
            check_deadline()
            if line.startswith(b'>'):
                active = line[1:].split()[0].decode() == ACCESSION
                if active:
                    found += 1
                    target.write(line)
            elif active:
                bases += len(line.strip())
                target.write(line)
    if found != 1 or bases != 50818468:
        raise ValueError(f'NC_000022.11 FASTA mismatch: records={found}, bases={bases}')
    genes = {}
    exon_rows = 0
    with gzip.open(gtf_source, 'rt') as stream, (out / 'reference.gtf').open('w') as target:
        for line in stream:
            check_deadline()
            if line.startswith('#'):
                target.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('invalid reference GTF row')
            if fields[0] != ACCESSION:
                continue
            if int(fields[3]) < 1 or int(fields[4]) > bases:
                raise ValueError('GTF coordinates outside selected accession')
            target.write(line)
            if fields[2] == 'exon':
                attributes = dict(re.findall(r'(\w+) "([^"]*)";', fields[8]))
                gene = attributes['gene_id']
                genes[gene] = {'gene_id': gene, 'gene_symbol': attributes.get('gene', ''),
                               'biotype': attributes.get('gene_biotype', ''), 'chromosome': ACCESSION}
                exon_rows += 1
    if not exon_rows or not genes:
        raise ValueError('selected accession has no annotated exons')
    shutil.copyfile(args.samples, out / 'samples.tsv')
    write_tsv(out / 'runs.tsv', ['run_id', 'sample_id', 'fastq_1', 'fastq_2', 'layout', 'strandedness'],
              [dict(run_id=row['run_id'], sample_id=row['sample_id'], fastq_1=f"reads/{row['run_id']}.fastq",
                    fastq_2='', layout='single', strandedness='reverse') for row in selected])
    write_tsv(out / 'references.tsv', ['role', 'path', 'sha256', 'source', 'release'],
              [dict(role=role, path=name, sha256=digest(out / name),
                    source=sources['reference']['files'][ASSEMBLY + suffix]['url'] + '#' + ACCESSION,
                    release=ASSEMBLY) for role, name, suffix in
               [('genome', 'reference.fa', '_genomic.fna.gz'), ('annotation', 'reference.gtf', '_genomic.gtf.gz')]])
    analysis = [('version', '1'), ('design_terms', 'condition'), ('alpha', '0.05'),
                ('filter', 'zero_total'), ('shrinkage', 'apeglm'), ('reference.condition', 'control'),
                ('type.age', 'numeric'), ('type.sex', 'categorical'), ('type.pmi', 'numeric'), ('type.rin', 'numeric')]
    write_tsv(out / 'analysis.tsv', ['key', 'value'], [dict(key=key, value=value) for key, value in analysis])
    write_tsv(out / 'contrasts.tsv', ['contrast_id', 'factor', 'numerator', 'denominator'],
              [dict(contrast_id='bipolar_vs_control', factor='condition', numerator='bipolar', denominator='control')])
    write_tsv(out / 'genes.tsv', ['gene_id'], [dict(gene_id=gene) for gene in sorted(genes)])
    write_tsv(out / 'annotation.tsv', ['gene_id', 'gene_symbol', 'biotype', 'chromosome'],
              [genes[gene] for gene in sorted(genes)])
    provenance.update(reference_bases=bases, exon_rows=exon_rows, genes=len(genes),
                      elapsed_seconds=round(time.time() - start, 3),
                      completed_utc=datetime.now(timezone.utc).isoformat())
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    (out / 'STAGED').write_text('restricted-reference smoke data prepared; not full-method validation\n')
    print(f'STAGED {out}: {len(genes)} genes, {bases} reference bases', flush=True)


if __name__ == '__main__':
    main()
