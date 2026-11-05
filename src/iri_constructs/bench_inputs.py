"""Reviewed local T05 inputs linked to the unchanged T02 CSV manifest."""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlparse

from .runtime import InputError, digest_file, read_json, require_text
from .validate_manifest import validate_manifest

AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
DNA = frozenset('ACGTRYSWKMBDHVN')


def validate_sequence(sequence: str, molecule: str) -> None:
    alphabet = AA if molecule == 'protein' else DNA if molecule == 'dna' else None
    if alphabet is None or not sequence or set(sequence) - alphabet:
        raise InputError('Empty sequence or unsupported alphabet; protein inference requires canonical amino acids')


@dataclass(frozen=True)
class ReviewedSource:
    input_id: str
    path: Path
    format: str
    molecule: str
    sequence_id: str
    review: dict = field(repr=False)


@dataclass(frozen=True)
class ProteinInput:
    input_id: str
    construct_id: str
    sequence: str = field(repr=False)
    source_sha256: str
    feature_id: str | None = None

    @property
    def sequence_sha256(self):
        return hashlib.sha256(self.sequence.encode('ascii')).hexdigest()

    def provenance(self):
        return {'input_id': self.input_id, 'construct_id': self.construct_id,
                'source_sha256': self.source_sha256, 'sequence_sha256': self.sequence_sha256,
                'length_aa': len(self.sequence), 'feature_id': self.feature_id}


def reviewed_sources(config_path: Path) -> tuple[dict, list[ReviewedSource], dict]:
    config_path = Path(config_path).resolve()
    config = read_json(config_path)
    if config.get('schema_version') != 't05-inputs-1':
        raise InputError('Expected schema_version t05-inputs-1')
    manifest = (config_path.parent / require_text(config, 'review_manifest')).resolve()
    try:
        errors = validate_manifest(manifest)
        if errors:
            raise InputError('T02 manifest validation failed: ' + '; '.join(errors))
        with manifest.open(encoding='utf-8-sig', newline='') as handle:
            rows = {row['record_id']: row for row in csv.DictReader(handle)}
    except OSError as exc:
        raise InputError('Cannot read T02 review manifest') from exc
    entries = config.get('records')
    if not isinstance(entries, list) or not entries:
        raise InputError('records must be a nonempty list of reviewed local inputs')
    sources, seen = [], set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise InputError('Every records entry must be an object')
        iid = require_text(entry, 'input_id')
        if iid in seen:
            raise InputError('Duplicate input_id')
        seen.add(iid)
        row = rows.get(require_text(entry, 'record_id'))
        if row is None:
            raise InputError('Input has no matching T02 review record')
        if any(row[k] != 'approved' for k in ('review_status', 'licensing_status', 'safety_review')):
            raise InputError('Input is not approved under all T02 record, license and safety gates')
        for k in ('notes', 'checksum_sha256'):
            require_text(row, k)
        source_url = urlparse(row['source_url'])
        if source_url.scheme not in ('https', 'http') or not source_url.netloc:
            raise InputError('Reviewed source_url must be a public HTTP(S) source reference')
        path = (config_path.parent / require_text(entry, 'path')).resolve()
        if entry.get('format') not in ('fasta', 'genbank') or entry.get('molecule_type') not in ('dna', 'protein'):
            raise InputError('Input format must be fasta/genbank and molecule_type dna/protein')
        try:
            checksum = digest_file(path)
        except OSError as exc:
            raise InputError('Reviewed local source file is missing or unreadable') from exc
        if checksum.lower() != row['checksum_sha256'].lower():
            raise InputError('Source checksum does not match the reviewed T02 manifest')
        sources.append(ReviewedSource(iid, path, entry['format'], entry['molecule_type'],
                                      require_text(entry, 'sequence_id'), row))
    return config, sources, {'input_config_sha256': digest_file(config_path),
                             'review_manifest_sha256': digest_file(manifest),
                             'sources': [{k: s.review[k] for k in
                                 ('record_id', 'source', 'source_url', 'accession', 'accession_version',
                                  'retrieved_on', 'checksum_sha256', 'archetype', 'annotation_status',
                                  'review_status', 'licensing_status', 'safety_review')} for s in sources]}


def read_source(source: ReviewedSource):
    from Bio import SeqIO
    import warnings
    try:
        # Parser warnings can contain raw annotation/sequence text. Treat them as controlled failures.
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with source.path.open(encoding='utf-8-sig') as handle:
                records = list(SeqIO.parse(handle, source.format))
    except (ValueError, Warning) as exc:
        raise InputError('Malformed sequence record or annotation; inspect the local source file') from exc
    ids = [r.id for r in records]
    if len(ids) != len(set(ids)):
        raise InputError('Duplicate sequence IDs in source file')
    matches = [r for r in records if r.id == source.sequence_id]
    if len(matches) != 1:
        raise InputError('sequence_id must match exactly one source record')
    record = matches[0]
    try:
        sequence = str(record.seq).upper()
    except Exception as exc:
        raise InputError('Source contains undefined sequence content') from exc
    validate_sequence(sequence, source.molecule)
    if source.molecule == 'dna' and len(sequence) != int(source.review['sequence_length_nt']):
        raise InputError('DNA length does not match T02 source metadata')
    if source.format == 'genbank':
        accession = source.review['accession']
        expected = accession if accession.endswith('.' + source.review['accession_version']) else accession + '.' + source.review['accession_version']
        if record.id != expected:
            raise InputError('GenBank accession/version does not match reviewed source metadata')
        if len(record.features) != int(source.review['feature_count']):
            raise InputError('GenBank feature count does not match reviewed source metadata')
    return record


def load_proteins(config_path: Path) -> tuple[dict, dict[str, ProteinInput], dict, list[dict]]:
    config, sources, provenance = reviewed_sources(config_path)
    proteins, exclusions = {}, []
    for source in sources:
        record = read_source(source)
        checksum = digest_file(source.path)
        if source.review['annotation_status'] != 'approved':
            exclusions.append({'input_id': source.input_id, 'reason': 'T02 source annotation review is not approved'})
            continue
        if source.molecule == 'protein':
            if source.format != 'fasta':
                raise InputError('Protein input must be a FASTA record')
            proteins[source.input_id] = ProteinInput(source.input_id, source.review['record_id'],
                                                     str(record.seq).upper(), checksum)
        else:
            from .annotations import normalize_record
            document, audit, extracted = normalize_record(record, source.review, checksum)
            for unit_id, sequence in extracted.items():
                iid = source.input_id + '/' + unit_id
                unit = next(u for u in document['units'] if u['unit_id'] == unit_id)
                proteins[iid] = ProteinInput(iid, source.review['record_id'], sequence, checksum,
                                             unit['feature_ids'][0])
            exclusions.extend({'input_id': source.input_id, 'unit_id': u['unit_id'],
                               'reason': u['representation_note']} for u in document['units']
                              if u['unit_id'] not in extracted)
    return config, proteins, provenance, exclusions
