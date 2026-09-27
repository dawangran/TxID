#!/usr/bin/env python3
"""Locked, bounded real-read case: whole-genome alignment, two chr22 callers.

This is an interoperability case, not a discovery-accuracy benchmark. All outputs
must go to a new directory. The separately downloaded StringTie binary is supplied
through PATH.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def sha256(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reads', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--annotation', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    for name in ('inputs', 'alignment', 'callers'):
        (out / name).mkdir(exist_ok=True)
    manifest = out / 'inputs' / 'selection.json'
    if manifest.exists() or any((out / 'alignment').iterdir()) or any((out / 'callers').iterdir()):
        raise SystemExit('Refusing to overwrite an existing case; use a new output directory.')
    expected = 'a335dc75d3356e6e16dc1097f856a050b1f3d978732935ab5b6b0f14848a3427'
    digest = sha256(args.reads)
    if digest != expected:
        raise SystemExit('ENCFF105WIJ source checksum differs from locked public manifest')
    selected = out / 'inputs' / 'ENCFF105WIJ.first100000.fastq.gz'
    total_bases = 0
    with gzip.open(args.reads, 'rb') as source, selected.open('xb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as target:
            for index in range(100000):
                record = [source.readline() for _ in range(4)]
                if not record[0].startswith(b'@') or not record[2].startswith(b'+'):
                    raise ValueError(f'Invalid or incomplete FASTQ record {index + 1}')
                if not record[3] or len(record[1].rstrip()) != len(record[3].rstrip()):
                    raise ValueError(f'Sequence/quality mismatch at record {index + 1}')
                target.writelines(record)
                total_bases += len(record[1].rstrip())
    annotation = out / 'inputs' / 'gencode-v29.identity.chr22.gtf'
    records = 0
    with gzip.open(args.annotation, 'rt') as source, annotation.open('x') as target:
        for line in source:
            if line.startswith('chr22\t'):
                target.write(line)
                records += 1
    if not records:
        raise ValueError('No chr22 annotation records')
    selection = {
        'accession': 'ENCFF105WIJ',
        'url': 'https://www.encodeproject.org/files/ENCFF105WIJ/@@download/ENCFF105WIJ.fastq.gz',
        'source_sha256': digest, 'selection': 'first 100000 complete FASTQ records in file order',
        'reads': 100000, 'bases': total_bases, 'selected_sha256': sha256(selected),
        'annotation_source_sha256': sha256(args.annotation),
        'annotation_sha256': sha256(annotation), 'annotation_records': records,
        'annotation_selection': 'chr22 records from retained GENCODE v29 identity view',
        'reference': str(args.reference.resolve()), 'reference_sha256': sha256(args.reference),
        'region': 'chr22', 'primary_alignment_exclusion_mask': 2308,
        'scripts': {str(q.relative_to(root)): sha256(q) for q in [
            Path(__file__), root / 'workflows/publication_benchmark/scripts/align_reads.py',
            root / 'workflows/publication_benchmark/scripts/run_caller.py',
            root / 'workflows/publication_benchmark/scripts/isoquant_entry.py']},
        'commands': [],
    }

    def run(argv):
        selection['commands'].append(argv)
        manifest.write_text(json.dumps(selection, indent=2) + '\n')
        subprocess.run(argv, check=True)

    bam = out / 'alignment' / 'whole-genome.bam'
    run([sys.executable, str(root / 'workflows/publication_benchmark/scripts/align_reads.py'),
         '--reference', str(args.reference.resolve()), '--reads', str(selected),
         '--read-type', 'pacbio_ccs', '--bam', str(bam), '--bai', str(bam) + '.bai',
         '--flagstat', str(out / 'alignment/flagstat.txt'),
         '--run-json', str(out / 'alignment/run.json'), '--threads', str(args.threads)])
    chr_bam = out / 'alignment' / 'chr22.primary.bam'
    run(['samtools', 'view', '-bh', '-F', '2308', '-o', str(chr_bam), str(bam), 'chr22'])
    run(['samtools', 'index', str(chr_bam)])
    selection['chr22_primary_alignments'] = int(subprocess.check_output(
        ['samtools', 'view', '-c', str(chr_bam)], text=True).strip())
    selection['chr22_bam_sha256'] = sha256(chr_bam)
    caller_commands = []
    for caller in ('isoquant', 'stringtie'):
        caller_commands.append([
            sys.executable, str(root / 'workflows/publication_benchmark/scripts/run_caller.py'),
            '--caller', caller, '--dataset', 'ENCFF105WIJ_chr22',
            '--reference', str(args.reference.resolve()), '--annotation', str(annotation),
            '--bam', str(chr_bam), '--read-type', 'pacbio_ccs',
            '--output-gtf', str(out / 'callers' / f'{caller}.gtf'),
            '--work-dir', str(out / 'callers' / caller),
            '--run-json', str(out / 'callers' / f'{caller}.run.json'),
            '--threads', str(max(1, args.threads // 2)),
            *(['--infer-annotation-meta-features'] if caller == 'isoquant' else [])])
    selection['commands'].extend(caller_commands)
    manifest.write_text(json.dumps(selection, indent=2) + '\n')
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda argv: subprocess.run(argv, check=False).returncode, caller_commands))
    selection['caller_exit_codes'] = dict(zip(('isoquant', 'stringtie'), results))
    manifest.write_text(json.dumps(selection, indent=2) + '\n')
    if any(results):
        raise SystemExit(f'Caller failure: {results}; inspect retained logs')
    print(json.dumps({'reads': 100000, 'chr22_primary_alignments': selection['chr22_primary_alignments'],
                      'caller_exit_codes': selection['caller_exit_codes']}))


if __name__ == '__main__':
    main()
