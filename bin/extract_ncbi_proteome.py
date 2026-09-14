#!/usr/bin/env python3
"""Extract a matched NCBI protein/GFF3 pair using the package assembly report."""
import argparse
import json
from pathlib import Path
import re
import shutil

from filter_isoforms_gff3 import checksum


def extract(directory, accession, protein_output, annotation_output, metadata_output):
    if not re.fullmatch(r'GC[AF]_\d+(\.\d+)?', accession):
        raise ValueError(f'Invalid NCBI assembly accession: {accession}')
    root = directory / 'ncbi_dataset' / 'data'
    with (root / 'assembly_data_report.jsonl').open() as handle:
        reports = [json.loads(line) for line in handle if line.strip()]
    matches = [report for report in reports if report.get('accession') == accession or
               ('.' not in accession and report.get('accession', '').split('.')[0] == accession)]
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one assembly report for {accession}; found {len(matches)}')
    resolved = matches[0]['accession']
    if not re.fullmatch(r'GC[AF]_\d+\.\d+', resolved):
        raise ValueError(f'Invalid resolved assembly accession: {resolved}')
    protein = root / resolved / 'protein.faa'
    annotation = root / resolved / 'genomic.gff'
    for path in (protein, annotation):
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f'Missing/empty annotation package file: {path}')
    metadata = dict(requested_accession=accession, resolved_accession=resolved,
                    assembly_report=matches[0], protein_sha256=checksum(protein),
                    annotation_sha256=checksum(annotation))
    shutil.copyfile(protein, protein_output)
    shutil.copyfile(annotation, annotation_output)
    metadata_output.write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--accession', required=True)
    for flag in ('protein-output', 'annotation-output', 'metadata-output'):
        parser.add_argument('--' + flag, type=Path, required=True)
    try:
        extract(**vars(parser.parse_args()))
    except (ValueError, OSError) as error:
        parser.exit(1, f'ERROR: {error}\n')
