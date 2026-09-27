#!/usr/bin/env python3
"""Package the local scientific review materials and their compact source data.

No publication or upload is performed. Large raw reads, alignments, reference
FASTA files, databases and historical draft ZIPs are deliberately excluded.
"""
import hashlib
import io
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ROOT = PACKAGE.parents[1]
STEM = 'TxID-submission-materials-2026-09-23'


def main():
    contents = {}
    source = HERE / 'TxID-0.1.3-submission-source.zip'
    with zipfile.ZipFile(source) as z:
        assert z.testzip() is None
        for name in z.namelist():
            relative = name.split('/', 1)[1]
            if relative in {'RELEASE-MANIFEST.json', 'SHA256SUMS'}:
                continue
            contents[relative] = z.read(name)
    extensions = {'.md', '.json', '.tsv', '.bib', '.py', '.txt', '.yaml', '.gtf',
                  '.log', '.docx', '.png', '.svg', '.pdf', '.tif'}
    # These post-build reports contain this archive's digest and must remain
    # sidecars; including them would create a circular checksum dependency.
    sidecar_reports = {
        PACKAGE / 'reproducibility/standalone-materials-validation.json',
        PACKAGE / 'reproducibility/final-integrity-check.json',
    }
    candidates = [p for p in PACKAGE.rglob('*') if p.is_file()
                  and p.suffix in extensions and '__pycache__' not in p.parts
                  and p not in sidecar_reports
                  and not p.name.startswith(STEM)]
    manifest = json.loads((PACKAGE / 'reproducibility/source-manifest.json').read_text())
    candidates += [ROOT / row['path'] for row in manifest]
    candidates.append(source)
    for p in candidates:
        if p.is_symlink() or not p.resolve().is_relative_to(ROOT):
            raise ValueError(f'Unapproved archive source: {p}')
        if p.stat().st_size > 20 * 1024 * 1024:
            raise ValueError(f'Unexpected large source: {p}')
        relative = p.relative_to(ROOT).as_posix()
        data = p.read_bytes()
        if relative in contents and contents[relative] != data:
            raise ValueError(f'Source archive differs from current file: {relative}')
        contents[relative] = data
    instructions = '''# TxID local scientific review and reproduction package

This is a local, versioned working snapshot, not a public DOI deposit or a
submission to Bioinformatics. Use the English manuscript for submission
preparation; manuscript.md/docx pair English and Chinese for author review.
Author metadata, declarations and public software/data access still require
completion. Initial submission permits Format-Free presentation; journal-style
pagination is checked when preparing the formatted version. Read the manuscript
package's README and editorial-notes.zh-CN.md. Archive-digest integrity and
standalone-rebuild reports are supplied as external sidecars to avoid circular
checksum dependencies.

From the extracted root, install the software and plotting dependencies:

    python -m pip install .
    python -m pip install -r deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/reproducibility/requirements-artifacts.txt

Rebuild and validate the paper from included compact reports:

    python deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/reproducibility/build_artifacts.py
    python deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/reproducibility/verify_manuscript.py

The package includes the small real-case caller GTFs, reference annotation,
coordinate evidence, matrices, hashes and execution logs. It excludes raw FASTQ,
BAM, the full reference FASTA, and SQLite registries. To rerun registry evaluation,
obtain the full reference matching the recorded checksum and run:

    python benchmarks/evaluate_real_interop.py --case-dir deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/reproducibility/real-interop --reference /path/to/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --output /new/evaluation/directory

To repeat alignment and caller reconstruction, use the source accession and
pinned tools in the selection/tool manifests, then benchmarks/run_real_interop_callers.py
with a new output directory. The executed caller environment used Python 3.12.13;
the benchmark selection script requires Python >=3.11. Software runtime remains
Python >=3.10. No historical 0.1.0 benchmark is relabelled as a 0.1.3 run.
'''
    contents['REPRODUCE.md'] = instructions.encode()
    records = [{'path': n, 'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
               for n, b in sorted(contents.items())]
    contents['MATERIALS-MANIFEST.json'] = (json.dumps({
        'format': 'txid-local-scientific-materials-v1', 'public_deposit_verified': False,
        'source_archive_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'excluded': ['raw FASTQ', 'BAM', 'reference FASTA', 'SQLite databases', 'old draft ZIPs',
                     'post-build archive-digest sidecar reports'],
        'files': records}, indent=2) + '\n').encode()
    contents['SHA256SUMS'] = ''.join(f'{hashlib.sha256(b).hexdigest()}  {n}\n'
                                     for n, b in sorted(contents.items())).encode()
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for name, data in sorted(contents.items()):
            info = zipfile.ZipInfo(f'{STEM}/{name}', (2026, 9, 23, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            z.writestr(info, data)
    payload = buffer.getvalue()
    output = HERE / f'{STEM}.zip'
    output.write_bytes(payload)
    with zipfile.ZipFile(output) as z:
        assert z.testzip() is None
        for row in records:
            assert hashlib.sha256(z.read(f'{STEM}/' + row['path'])).hexdigest() == row['sha256']
    report = {'archive': output.name, 'sha256': hashlib.sha256(payload).hexdigest(),
              'bytes': len(payload), 'files': len(records), 'crc_and_all_file_hashes_valid': True,
              'public_deposit_verified': False}
    (HERE / f'{STEM}.manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    (HERE / f'{STEM}.zip.sha256').write_text(report['sha256'] + '  ' + output.name + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
