"""Small structural checks and per-case QC summaries for existing stage boundaries."""
import json
import os
from pathlib import Path
import re

import pandas as pd

KEY = ['Chr_Pos', 'Variant_REF_ALT']
ALIGNED = ['Gene_ID', 'Transcript_biotype', 'Consequences']


def is_mitochondrial_variant(row):
    chrom = str(row.get('Chr_Pos', '')).split(':', 1)[0].casefold().removeprefix('chr')
    flag = str(row.get('is_mitochondrial', '')).casefold() in {'true', '1', 'yes'}
    routes = re.split(r'\s*(?:\||,)\s*', str(row.get('Evidence_routes', '')).upper())
    return chrom in {'m', 'mt'} or flag or 'MITOCHONDRIAL' in routes


def require_columns(frame, columns, label):
    missing = set(columns) - set(frame.columns)
    if missing or frame.columns.duplicated().any():
        raise ValueError(f'{label}: missing or duplicate columns: {sorted(missing)}')


def variant_keys(frame, label):
    require_columns(frame, KEY, label)
    locations = frame['Chr_Pos'].fillna('').astype(str)
    alleles = frame['Variant_REF_ALT'].fillna('').astype(str)
    if not locations.str.match(r'^\S+:[1-9][0-9]*$').all():
        raise ValueError(f'{label}: malformed chromosome/position')
    if not alleles.str.match(r'^[^\s>,]+>[^\s>,]+$').all():
        raise ValueError(f'{label}: malformed REF/ALT')
    return locations.str.replace(r'^chr', '', regex=True).str.replace(r'^M:', 'MT:', regex=True) + '|' + alleles.str.upper()


def validate_table(frame, label, required=(), unique=True, aligned=False):
    require_columns(frame, required, label)
    keys = variant_keys(frame, label)
    if unique and keys.duplicated().any():
        raise ValueError(f'{label}: duplicate canonical variant keys')
    if aligned:
        fields = ['Transcript_ID', *ALIGNED]
        gene = 'Gene_symbol_transcript_aligned' if 'Gene_symbol_transcript_aligned' in frame else 'Gene_symbol'
        fields.append(gene)
        require_columns(frame, fields, label)
        lengths = frame[fields].fillna('').astype(str).apply(lambda col: col.str.split(' | ', regex=False).str.len())
        if lengths.nunique(axis=1).gt(1).any():
            raise ValueError(f'{label}: transcript/gene fields lost alignment')
    return keys


def record_qc(label, input_rows, output_rows, input_variants, output_variants,
              operation, warnings=(), **details):
    report = dict(stage=label, status='PASS', operation=operation,
                  input_rows=int(input_rows), output_rows=int(output_rows),
                  input_variants=int(input_variants), output_variants=int(output_variants),
                  warnings=list(warnings), **details)
    directory = Path('results/qc')
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f'{label}.json'
    temporary = target.with_suffix(f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    os.replace(temporary, target)
    print(f'QC {label}: {input_rows} rows/{input_variants} variants -> {output_rows} rows/{output_variants} variants ({operation})')
    return report


def check_transition(before, after, label, operation='enrichment', unique=True,
                     required=(), aligned=False, changed_fields=(), warnings=()):
    old = validate_table(before, label + ' input', unique=unique)
    new = validate_table(after, label + ' output', required, unique=True if operation == 'collapse' else unique, aligned=aligned)
    if operation == 'enrichment':
        require_columns(after, before.columns, label)
        if len(before) != len(after) or old.tolist() != new.tolist():
            raise ValueError(f'{label}: evidence enrichment changed variant identities/counts/order')
        for column in before.columns:
            if column in changed_fields:
                continue
            a = before[column].astype(str).where(before[column].notna(), '').reset_index(drop=True)
            b = after[column].astype(str).where(after[column].notna(), '').reset_index(drop=True)
            if not a.equals(b):
                raise ValueError(f'{label}: original evidence changed: {column}')
    elif operation == 'collapse':
        if set(old) != set(new):
            raise ValueError(f'{label}: collapse lost or added a variant')
    elif operation == 'filter':
        if not set(new).issubset(set(old)):
            raise ValueError(f'{label}: filtering added variants')
    else:
        raise ValueError(f'Unknown QC operation: {operation}')
    return record_qc(label, len(before), len(after), old.nunique(), new.nunique(),
                     operation, warnings, removed_variants=int(old.nunique() - new.nunique()))


class VCFCounts:
    """Count sorted, normalized VCF records without retaining a genome-wide key set."""
    def __init__(self, contigs):
        self.ranks = {chrom: i for i, chrom in enumerate(contigs)}
        self.previous = None
        self.alleles = set()
        self.count = 0

    def observe(self, chrom, pos, ref, alt):
        if chrom not in self.ranks or int(pos) < 1 or not ref or not alt or ',' in alt:
            raise ValueError('Malformed normalized VCF identity')
        locus = (self.ranks[chrom], int(pos))
        if self.previous is not None and locus < self.previous:
            raise ValueError('Normalized VCF records are out of order')
        if locus != self.previous:
            self.alleles.clear()
        allele = (str(ref).upper(), str(alt).upper())
        if allele in self.alleles:
            raise ValueError('Duplicate normalized VCF variant')
        self.alleles.add(allele)
        self.previous = locus
        self.count += 1
