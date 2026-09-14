#!/usr/bin/env python3
"""Select longest proteins using exact CDS protein_id -> gene ancestry in GFF3.

Never infer gene identity from symbols, sequence identity, or accession prefixes.
Unmapped, duplicate and ambiguously assigned FASTA IDs are fatal and audited.
"""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import unquote

from normalize_proteome_fasta import open_text


def checksum(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path):
    header, residues = None, []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(residues)
                header, residues = line[1:].strip(), []
                if not header:
                    raise ValueError(f'Empty FASTA identifier at line {number}')
            elif line.strip():
                if header is None:
                    raise ValueError(f'Sequence before FASTA header at line {number}')
                residues.append(''.join(line.split()))
    if header is not None:
        yield header, ''.join(residues)


def parse_attributes(text):
    attributes = {}
    if text == '.':
        return attributes
    for field in text.rstrip(';').split(';'):
        key, separator, value = field.partition('=')
        if not separator or not key or key in attributes:
            raise ValueError(f'Malformed/duplicate GFF3 attribute: {field!r}')
        # Split BEFORE decoding: %2C is part of an ID, not a list separator.
        attributes[key] = tuple(unquote(item) for item in value.split(','))
    return attributes


def protein_gene_mapping(path):
    features, cds, metadata = {}, [], []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#!'):
                metadata.append(line.strip())
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'GFF3 line {number}: expected nine columns')
            attributes = parse_attributes(fields[8])
            identifiers = attributes.get('ID', ())
            if len(identifiers) > 1:
                raise ValueError(f'GFF3 line {number}: multiple feature IDs')
            identifier = identifiers[0] if identifiers else None
            gene_ids = frozenset(x for x in attributes.get('Dbxref', ()) if x.startswith('GeneID:') and x != 'GeneID:')
            proteins = attributes.get('protein_id', ())
            feature = (fields[2], attributes.get('Parent', ()), gene_ids, fields[0], fields[6])
            if identifier:
                if identifier in features and features[identifier] != feature:
                    raise ValueError(f'Conflicting definitions for GFF3 feature {identifier}')
                features[identifier] = feature
            if fields[2] == 'CDS':
                cds.append((feature, proteins))

    memo = {}

    def ancestry(identifier, visiting):
        if identifier in visiting:
            raise ValueError(f'Cycle in GFF3 ancestry at {identifier}')
        if identifier not in features:
            raise ValueError(f'Missing GFF3 parent {identifier}')
        if identifier in memo:
            return memo[identifier]
        kind, parents, own_ids, _, _ = features[identifier]
        gene_ids = set(own_ids)
        genes = {identifier} if kind in ('gene', 'pseudogene') else set()
        for parent in parents:
            parent_ids, parent_genes = ancestry(parent, visiting | {identifier})
            gene_ids.update(parent_ids)
            genes.update(parent_genes)
        memo[identifier] = gene_ids, genes
        return gene_ids, genes

    mapping, problems = defaultdict(set), defaultdict(set)
    for feature, proteins in cds:
        if not proteins:
            continue
        _, parents, own_ids, _, _ = feature
        gene_ids, genes = set(own_ids), set()
        try:
            for parent in parents:
                parent_ids, parent_genes = ancestry(parent, set())
                gene_ids.update(parent_ids)
                genes.update(parent_genes)
            if len(gene_ids) == 1:
                gene = next(iter(gene_ids))
            elif not gene_ids and len(genes) == 1:
                gene = 'gene_feature:' + next(iter(genes))
            else:
                raise ValueError(f'Missing or conflicting gene ancestry: {sorted(gene_ids or genes)}')
            for protein in proteins:
                mapping[protein].add(gene)
        except ValueError as error:
            for protein in proteins:
                problems[protein].add(str(error))
    return mapping, problems, metadata


def filter_proteome(fasta, gff3, output, ledger, qc, species):
    paths = [fasta, gff3, output, ledger, qc]
    if len({path.resolve() for path in paths}) != len(paths):
        raise ValueError('Input and output paths must be distinct')
    mapping, problems, metadata = protein_gene_mapping(gff3)
    annotation_sha256 = checksum(gff3)
    rows, best, seen, errors = [], {}, set(), []
    for header, sequence in fasta_records(fasta):
        protein = header.split()[0]
        length = len(sequence[:-1] if sequence.endswith('*') else sequence)
        genes = mapping.get(protein, set())
        reasons = sorted(problems.get(protein, set()))
        if protein in seen:
            reasons.append('duplicate_fasta_id')
        seen.add(protein)
        if not genes:
            reasons.append('unmapped_protein_id')
        elif len(genes) != 1:
            reasons.append('ambiguous_gene_assignment: ' + ', '.join(sorted(genes)))
        if not length:
            reasons.append('empty_protein')
        gene = next(iter(genes)) if len(genes) == 1 else ''
        row = dict(species=species, protein_id=protein, gene_id=gene, length=length,
                   mapping_source=('NCBI_GeneID' if gene.startswith('GeneID:') else 'GFF3_gene_feature') if gene else 'unresolved',
                   annotation_sha256=annotation_sha256,
                   representative='', decision='', reason='; '.join(reasons))
        rows.append(row)
        if reasons:
            errors.append(f'{protein}: {row["reason"]}')
            row['decision'] = 'error'
        elif gene not in best or (-length, protein) < best[gene]:
            best[gene] = (-length, protein)
    if not rows:
        errors.append('No protein records in FASTA')
    for row in rows:
        if row['decision'] == 'error':
            continue
        row['representative'] = best[row['gene_id']][1]
        row['decision'] = 'retained' if row['protein_id'] == row['representative'] else 'removed_isoform'
        row['reason'] = 'longest_protein_then_lexical_id'
    with ledger.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['species', 'protein_id', 'gene_id', 'length',
                               'representative', 'decision', 'reason', 'mapping_source', 'annotation_sha256'],
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    report = dict(species=species, status='failed' if errors else 'passed',
                  protein_sha256=checksum(fasta), annotation_sha256=annotation_sha256,
                  annotation_metadata=metadata, input_proteins=len(rows), mapped_genes=len(best),
                  retained=sum(row['decision'] == 'retained' for row in rows),
                  removed_isoforms=sum(row['decision'] == 'removed_isoform' for row in rows),
                  errors=errors)
    qc.write_text(json.dumps(report, indent=2) + '\n')
    if errors:
        raise ValueError(f'{len(errors)} mapping errors; see {ledger}: ' + '; '.join(errors[:10]))
    assert report['input_proteins'] == report['retained'] + report['removed_isoforms']
    assert report['retained'] == report['mapped_genes']
    selected = {protein for _, protein in best.values()}
    with output.open('w') as handle:
        for header, sequence in fasta_records(fasta):
            if header.split()[0] in selected:
                handle.write(f'>{header}\n{sequence}\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('fasta', 'gff3', 'output', 'ledger', 'qc'):
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--species', required=True)
    args = parser.parse_args()
    try:
        filter_proteome(**vars(args))
    except (ValueError, OSError) as error:
        print(f'ERROR: GFF3 isoform filtering: {error}', file=sys.stderr)
        # Parse failures occur before the per-protein ledger can be constructed.
        if args.qc.resolve() not in {args.fasta.resolve(), args.gff3.resolve()} and not args.qc.exists():
            args.qc.write_text(json.dumps(dict(status='failed', errors=[str(error)]), indent=2) + '\n')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
