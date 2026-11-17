"""Public preprocessing API: annotated GenBank files to a traceable output bundle.

    from iri_constructs import preprocess_files
    result = preprocess_files(['plasmid.gb'], 'out/')

Nothing here needs a review manifest, family labels, a reference panel, model
packages, the Git checkout or network access. Biopython is imported lazily so
that importing this module (and ``--help``) stays dependency-light.

Coordinates are zero-based and half-open, with strand and *ordered* source
segments. A joined or origin-spanning feature keeps its segment list; it is
never flattened into a bounding interval.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import shutil
import tempfile
import time
from typing import Mapping, Sequence

from . import __version__
from .annotations import (CDS_OUTCOMES, LOCATION_ISSUES, construct_annotation_status,
                          json_safe, normalized_features)
from .runtime import InputError
from .sequences import DNA

BUNDLE_SCHEMA = 'iri-constructs-bundle-1'
BUNDLE_FILES = ('normalized.json', 'annotations.json', 'proteins.fasta', 'features.tsv',
                'units.tsv', 'report.json', 'run.json')
INPUT_SUFFIXES = ('.gb', '.gbk', '.gbff', '.genbank', '.fasta', '.fa', '.fna')
COORDINATES = ('Zero-based, half-open [start0, end0) on the forward strand of the parent record; '
               'strand is 1 or -1; segments are listed in the source-annotated order.')
PROVENANCE_FIELDS = ('source', 'source_url', 'retrieved_on', 'terms_note')

# Every reason a unit can carry, with a one-line explanation for reports.
REASONS = {
    'translation_validated': 'Source /translation equals the translation of the documented coordinates, frame and table.',
    'noncoding_metadata': 'Non-CDS feature kept as metadata; no protein unit is expected.',
    'missing_translation': 'CDS has no /translation; nothing is inferred or synthesized.',
    'exceptional_cds_annotation': 'CDS carries /pseudo, /pseudogene, /transl_except or /exception.',
    'conflicting_cds_qualifiers': 'CDS repeats /translation, /product, /gene, /codon_start or /transl_table.',
    'translation_mismatch': 'Source /translation disagrees with the documented coordinates, frame and table.',
    'empty_translation': 'CDS /translation is present but empty.',
    'noncanonical_translation_residues': 'Source /translation contains residues outside the 20 canonical amino acids.',
    'invalid_codon_start': '/codon_start is not 1, 2 or 3.',
    'invalid_transl_table': '/transl_table is not a known translation table.',
    'incomplete_cds': 'Documented coordinates are not a complete CDS (start codon, whole codons, single terminal stop).',
    **{code: text for code, text in LOCATION_ISSUES.items()},
}
NOTES = ('Extraction validates agreement between a source translation and the source coordinates. '
         'It does not establish that the annotation is biologically correct, that the protein is functional, '
         'or anything about safety.',
         'Withheld units are expected, documented outcomes; record failures are parsing problems.',
         'Model compatibility is reported separately and never changes what is extracted.')

_PLACEHOLDER_IDS = {'', '.', 'unknown', 'exported', 'unnamed', 'untitled', 'record'}


class PreprocessError(Exception):
    """Invalid invocation or a fatal input/output problem; nothing was written."""


class OutputError(PreprocessError):
    """The output location is unsafe: it exists, is not a bundle, or contains an input."""


@dataclass(frozen=True)
class PreprocessOptions:
    """``overwrite`` replaces an existing bundle only. ``model_max_residues`` adds a
    compatibility column for one chosen model limit; it never filters or truncates.
    ``provenance`` maps an input file name or SHA-256 to caller-supplied facts."""
    overwrite: bool = False
    model_max_residues: int | None = None
    provenance: Mapping[str, Mapping] | None = None


@dataclass(frozen=True)
class Issue:
    code: str
    message: str

    def as_dict(self) -> dict:
        return {'code': self.code, 'message': self.message}


@dataclass(frozen=True)
class FeatureResult:
    feature_id: str
    construct_id: str
    unit_id: str
    source_feature_index: int
    source_feature_type: str
    feature_label: str
    display_name: str | None
    strand: str
    segments: tuple
    location_status: str
    source_location: str
    annotation_status: str
    names: Mapping[str, str | None] = field(default_factory=dict)


@dataclass(frozen=True)
class UnitResult:
    unit_id: str
    construct_id: str
    feature_id: str
    molecule_type: str           # 'protein' | 'metadata_only'
    extraction_state: str        # 'exported' | 'withheld' | 'metadata_only'
    reason_code: str
    detail: Mapping = field(default_factory=dict)
    length_aa: int | None = None
    sequence_sha256: str | None = None
    model_compatibility: str = 'not_evaluated'
    sequence: str | None = field(default=None, repr=False)


@dataclass(frozen=True)
class RecordResult:
    input_index: int
    record_index: int            # 1-based position inside its input file; 0 for a file-level failure
    input_name: str
    status: str                  # 'ok' | 'failed' | 'skipped_duplicate'
    source_record_id: str | None = None
    construct_id: str | None = None
    failure: Issue | None = None
    warnings: tuple = ()
    duplicate_of: str | None = None
    construct: Mapping | None = None
    features: tuple = ()
    units: tuple = ()
    annotations: Mapping | None = field(default=None, repr=False)


@dataclass(frozen=True)
class InputFile:
    index: int
    path: Path
    name: str
    sha256: str | None
    size_bytes: int | None
    format: str | None


@dataclass(frozen=True)
class BatchResult:
    options: PreprocessOptions
    inputs: tuple
    records: tuple
    output_dir: Path | None = None

    @property
    def ok_records(self) -> tuple:
        return tuple(r for r in self.records if r.status == 'ok')

    @property
    def failed_records(self) -> tuple:
        return tuple(r for r in self.records if r.status == 'failed')

    @property
    def units(self) -> tuple:
        return tuple(u for r in self.ok_records for u in r.units)

    @property
    def proteins(self) -> dict:
        """Validated source proteins keyed by unit ID."""
        return {u.unit_id: u.sequence for u in self.units if u.extraction_state == 'exported'}

    @property
    def exit_code(self) -> int:
        """0 all records processed; 2 partial batch; 1 nothing usable."""
        if not self.ok_records:
            return 1
        return 2 if self.failed_records else 0

    def summary(self) -> dict:
        units = self.units
        withheld = Counter(u.reason_code for u in units if u.extraction_state == 'withheld')
        exported = [u for u in units if u.extraction_state == 'exported']
        summary = {
            'inputs': len(self.inputs), 'records_total': len(self.records),
            'records_ok': len(self.ok_records), 'records_failed': len(self.failed_records),
            'records_skipped_duplicate': sum(r.status == 'skipped_duplicate' for r in self.records),
            'features': len(units),
            'cds_features': sum(u.molecule_type == 'protein' for u in units),
            'noncoding_metadata_features': sum(u.molecule_type == 'metadata_only' for u in units),
            'proteins_exported': len(exported), 'cds_withheld': sum(withheld.values()),
            'cds_withheld_by_reason': dict(sorted(withheld.items())),
            'failures_by_code': dict(sorted(Counter(r.failure.code for r in self.failed_records).items())),
        }
        limit = self.options.model_max_residues
        summary['model_compatibility'] = {
            'evaluated': limit is not None, 'max_residues': limit,
            'compatible': sum(u.model_compatibility == 'compatible' for u in exported) if limit is not None else None,
            'exceeds_max_residues': sum(u.model_compatibility == 'exceeds_max_residues' for u in exported) if limit is not None else None}
        return summary


# ---------------------------------------------------------------- reading inputs

def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _first_line(text: object, limit: int = 160) -> str:
    """Parser messages can quote annotation payload after the first line; keep the headline only."""
    lines = str(text).strip().splitlines()
    return (lines[0] if lines else '')[:limit]


def _slug(text: str | None) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', text or '').strip('._-')[:64]


def expand_inputs(paths: Sequence) -> list[Path]:
    """Files are taken as given; a directory contributes its sequence files (not recursive), sorted."""
    files = []
    for raw in paths:
        path = Path(raw)
        if path.is_dir():
            files.extend(sorted(p for p in path.iterdir() if p.is_file() and p.suffix.lower() in INPUT_SUFFIXES))
        elif path.is_file():
            files.append(path)
        else:
            raise PreprocessError(f'Input does not exist or is not a regular file: {path}')
    if not files:
        raise PreprocessError('No input files found')
    return files


def _split_records(text: str, fmt: str) -> list[str]:
    """Split on record boundaries so one malformed record cannot hide its neighbours."""
    chunks, current = [], []
    for line in text.splitlines():
        starts = line.startswith('LOCUS') if fmt == 'genbank' else line.startswith('>')
        if starts and any(x.strip() for x in current):
            chunks.append(current)
            current = []
        current.append(line)
        if fmt == 'genbank' and line.rstrip() == '//':
            chunks.append(current)
            current = []
    if any(x.strip() for x in current):
        chunks.append(current)
    return ['\n'.join(x.rstrip() for x in chunk).strip() + '\n' for chunk in chunks]


def _parse_chunk(chunk: str, fmt: str):
    """Return (SeqRecord | None, sequence | None, warnings, failure)."""
    from Bio import SeqIO
    import warnings
    notes = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            records = list(SeqIO.parse(io.StringIO(chunk), fmt))
        except Exception as exc:
            return None, None, notes, Issue('malformed_record', f'Parser rejected the record ({type(exc).__name__}): {_first_line(exc)}')
    notes = [Issue('parser_warning', _first_line(w.message)) for w in caught]
    if len(records) != 1:
        return None, None, notes, Issue('malformed_record', 'No complete sequence record could be read from this text')
    record = records[0]
    try:
        sequence = str(record.seq).upper()
    except Exception:
        return None, None, notes, Issue('sequence_undefined', 'Record has no sequence data (for example a CONTIG-only or annotation-only record)')
    if not sequence:
        return None, None, notes, Issue('sequence_undefined', 'Record has an empty sequence')
    if 'protein' in str(record.annotations.get('molecule_type', '')).lower():
        return None, None, notes, Issue('not_nucleotide_record', 'Protein record; this toolkit reads annotated nucleotide constructs')
    if set(sequence) - DNA:
        return None, None, notes, Issue('unsupported_sequence_alphabet', 'Sequence is not IUPAC DNA (protein, RNA or other characters present)')
    if fmt == 'genbank':
        declared = re.match(r'LOCUS\s+\S+\s+(\d+)\s+bp\b', chunk)
        if declared and int(declared.group(1)) != len(sequence):
            return None, None, notes, Issue('declared_length_mismatch',
                f'LOCUS declares {declared.group(1)} bp but {len(sequence)} were read; the record may be truncated')
    return record, sequence, notes, None


def _accession_version(record_id: str) -> str | None:
    return record_id if re.fullmatch(r'[A-Za-z]{1,6}_?\d{3,}\.\d+', record_id or '') else None


# ---------------------------------------------------------------- normalization

def _strand_text(segments: list, code: str | None) -> str:
    strands = {s['strand'] for s in segments}
    if code == 'location_mixed_strand' or len(strands) > 1:
        return 'mixed'
    return {1: '+', -1: '-'}.get(next(iter(strands), None), '')


def _build_record(parsed: dict, construct_id: str, options: PreprocessOptions, source_input: dict) -> dict:
    record, sequence = parsed['record'], parsed['sequence']
    entries = normalized_features(record, len(sequence), strict=False)
    features, units, originals = [], [], {}
    limit = options.model_max_residues
    for entry in entries:
        suffix = entry['signature'][:20] + ('' if entry['duplicate_index'] == 1 else f"-dup{entry['duplicate_index']}")
        fid, uid = f'{construct_id}-f-{suffix}', f'{construct_id}-u-{suffix}'
        q = entry['original']['qualifiers']
        features.append(FeatureResult(
            feature_id=fid, construct_id=construct_id, unit_id=uid,
            source_feature_index=entry['source_index'], source_feature_type=entry['type'],
            feature_label=entry['label'], display_name=entry['display_name'],
            strand=_strand_text(entry['segments'], entry['location_code']),
            segments=tuple(dict(s) for s in entry['segments']),
            location_status=entry['location_code'] or 'exact',
            source_location=entry['original']['location'], annotation_status=entry['status'],
            names={k: (q.get(k) or [None])[0] for k in ('label', 'gene', 'product', 'protein_id')}))
        protein = entry['protein']
        state = 'metadata_only' if not entry['is_cds'] else 'exported' if protein else 'withheld'
        compatibility = 'not_evaluated'
        if protein and limit is not None:
            compatibility = 'compatible' if len(protein) <= limit else 'exceeds_max_residues'
        units.append(UnitResult(
            unit_id=uid, construct_id=construct_id, feature_id=fid,
            molecule_type='protein' if entry['is_cds'] else 'metadata_only',
            extraction_state=state, reason_code=entry['code'], detail=dict(entry['detail']),
            length_aa=len(protein) if protein else None,
            sequence_sha256=_sha256(protein.encode('ascii')) if protein else None,
            model_compatibility=compatibility, sequence=protein))
        originals[fid] = {'type': entry['type'], 'location': entry['original']['location'],
                          'operator': entry['original']['operator'],
                          'source_feature_id': entry['original']['source_feature_id'],
                          'parts': entry['parts'], 'qualifier_order': [str(k) for k in q],
                          'qualifiers': json_safe(q), 'duplicate_index': entry['duplicate_index']}
    warnings = list(parsed['warnings'])
    duplicates = sum(e['duplicate_index'] > 1 for e in entries)
    if duplicates:
        warnings.append(Issue('duplicate_feature', f'{duplicates} identical duplicate source feature(s) kept with a -dupN suffix'))
    topology = record.annotations.get('topology')
    accessions = [a for a in record.annotations.get('accessions', []) if _slug(a).lower() not in _PLACEHOLDER_IDS]
    supplied = None
    for key in (source_input['sha256'], source_input['name']):
        if options.provenance and key in options.provenance:
            supplied = json_safe(dict(options.provenance[key]))
            break
    construct = {
        'construct_id': construct_id, 'source_record_id': record.id, 'source_record_name': record.name,
        'description': record.description or None,
        'accession': accessions[0] if accessions else None,
        'accession_version': _accession_version(record.id),
        'organism': record.annotations.get('organism'), 'source_date': record.annotations.get('date'),
        'molecule_type': record.annotations.get('molecule_type'),
        'topology': topology if topology in ('linear', 'circular') else 'unknown',
        'sequence_length_nt': len(sequence), 'sequence_sha256': _sha256(sequence.encode('ascii')),
        'record_text_sha256': parsed['digest'], 'record_index': parsed['record_index'],
        'source_input': source_input,
        'annotation_status': construct_annotation_status(f.annotation_status for f in features),
        'provenance': {'supplied': supplied,
                       'missing_fields': [k for k in PROVENANCE_FIELDS if not (supplied or {}).get(k)]},
    }
    annotations = {'source_record_id': record.id, 'record_text_sha256': parsed['digest'],
                   'record_annotations': json_safe(dict(record.annotations)), 'features': originals}
    return {'construct': construct, 'features': tuple(features), 'units': tuple(units),
            'annotations': annotations, 'warnings': tuple(warnings)}


def process_inputs(paths: Sequence, options: PreprocessOptions | None = None) -> BatchResult:
    """Parse and normalize without writing anything. Raises PreprocessError for unusable invocations."""
    options = options or PreprocessOptions()
    if options.model_max_residues is not None and (type(options.model_max_residues) is not int or options.model_max_residues < 1):
        raise PreprocessError('model_max_residues must be a positive integer')
    files = expand_inputs(paths)
    inputs, parsed_records, results = [], [], {}
    for index, path in enumerate(files, 1):
        try:
            data = path.read_bytes()
        except OSError as exc:
            inputs.append(InputFile(index, path, path.name, None, None, None))
            results[(index, 0)] = RecordResult(index, 0, path.name, 'failed',
                failure=Issue('unreadable_input', f'Cannot read the file ({type(exc).__name__})'))
            continue
        file_notes = []
        try:
            text = data.decode('utf-8-sig')
        except UnicodeDecodeError:
            text = data.decode('latin-1')
            file_notes.append(Issue('non_utf8_input', 'File is not UTF-8; decoded as Latin-1 so no byte is dropped'))
        first = next((line.strip() for line in text.splitlines() if line.strip()), '')
        fmt = 'genbank' if first.startswith('LOCUS') else 'fasta' if first.startswith('>') else None
        inputs.append(InputFile(index, path, path.name, _sha256(data), len(data), fmt))
        if fmt is None:
            code, message = ('empty_input', 'File is empty') if not first else (
                'unrecognized_format', 'First line is neither a GenBank LOCUS line nor a FASTA header')
            results[(index, 0)] = RecordResult(index, 0, path.name, 'failed', failure=Issue(code, message))
            continue
        for record_index, chunk in enumerate(_split_records(text, fmt), 1):
            record, sequence, notes, failure = _parse_chunk(chunk, fmt)
            if failure:
                results[(index, record_index)] = RecordResult(
                    index, record_index, path.name, 'failed', failure=failure, warnings=tuple(file_notes + notes),
                    source_record_id=(re.match(r'(?:LOCUS\s+|>)(\S+)', chunk) or [None, None])[1])
                continue
            parsed_records.append({'input_index': index, 'record_index': record_index, 'path': path,
                                   'record': record, 'sequence': sequence, 'warnings': file_notes + notes,
                                   'digest': _sha256(chunk.encode('utf-8')), 'fmt': fmt})
    # Identity: exact duplicate record text is processed once; different records that
    # share an identifier are all disambiguated by content digest, never overwritten.
    first_by_digest, bases = {}, {}
    for item in parsed_records:
        item['duplicate_of'] = first_by_digest.get(item['digest'])
        first_by_digest.setdefault(item['digest'], item)
        base, placeholder = _slug(item['record'].id), False
        if base.lower() in _PLACEHOLDER_IDS or item['record'].id == '<unknown id>':
            base = _slug(item['record'].name)
            if base.lower() in _PLACEHOLDER_IDS or item['record'].name == '<unknown name>':
                base = _slug(item['path'].stem) or 'record'
            placeholder = True
        item['base'], item['placeholder'] = base, placeholder
        if item['duplicate_of'] is None:
            bases.setdefault(base, []).append(item)
    for base, group in bases.items():
        for item in group:
            item['shared_id'] = len(group) > 1
            item['construct_id'] = base if len(group) == 1 and not item['placeholder'] else f"{base}~{item['digest'][:10]}"
    for item in parsed_records:
        key = (item['input_index'], item['record_index'])
        name = item['path'].name
        if item['duplicate_of'] is not None:
            results[key] = RecordResult(item['input_index'], item['record_index'], name, 'skipped_duplicate',
                source_record_id=item['record'].id, construct_id=item['duplicate_of']['construct_id'],
                duplicate_of=item['duplicate_of']['construct_id'],
                warnings=(Issue('duplicate_record', 'Identical record text was already processed in this batch'),))
            continue
        source = inputs[item['input_index'] - 1]
        try:
            built = _build_record(item, item['construct_id'], options,
                                  {'input_index': source.index, 'name': source.name, 'sha256': source.sha256, 'format': item['fmt']})
        except InputError as exc:
            results[key] = RecordResult(item['input_index'], item['record_index'], name, 'failed',
                                        source_record_id=item['record'].id, failure=Issue('invalid_record', str(exc)))
            continue
        except Exception as exc:  # defensive: one unexpected record must not sink the batch
            results[key] = RecordResult(item['input_index'], item['record_index'], name, 'failed',
                source_record_id=item['record'].id,
                failure=Issue('normalization_error', f'Unexpected {type(exc).__name__} while normalizing this record'))
            continue
        warnings = list(built['warnings'])
        if item['shared_id']:
            warnings.append(Issue('shared_record_identifier', 'Another, different record in this batch has the same identifier; IDs carry a content digest'))
        if item['placeholder']:
            warnings.append(Issue('placeholder_record_identifier', 'Record has no usable accession/name; ID uses the name or file stem plus a content digest'))
        if item['fmt'] == 'fasta':
            warnings.append(Issue('sequence_only_input', 'FASTA carries no annotations; no features or proteins are invented'))
        results[key] = RecordResult(item['input_index'], item['record_index'], name, 'ok',
            source_record_id=item['record'].id, construct_id=item['construct_id'], warnings=tuple(warnings),
            construct=built['construct'], features=built['features'], units=built['units'], annotations=built['annotations'])
    return BatchResult(options, tuple(inputs), tuple(results[k] for k in sorted(results)))


# ---------------------------------------------------------------- output bundle

def _tsv(value) -> str:
    if value is None:
        return ''
    return str(value).replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n').replace('\r', '\\r')


def _segments_text(segments) -> str:
    """Biopython-style zero-based half-open spans in source order, e.g. ``[6:12](+);[0:6](+)``."""
    return ';'.join(f"[{s['start0']}:{s['end0']}]({ {1: '+', -1: '-'}.get(s['strand'], '?') })" for s in segments)


def _table(header: Sequence[str], rows) -> str:
    return ''.join('\t'.join(_tsv(v) for v in row) + '\n' for row in [header, *rows])


def _dump(value) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'


def bundle_documents(result: BatchResult) -> dict[str, str]:
    """Deterministic text of every bundle file except the variable ``run.json``."""
    tool = {'name': 'iri-constructs', 'version': __version__}
    constructs, annotations, fasta, feature_rows, unit_rows = [], {}, [], [], []
    for record in result.ok_records:
        cid = record.construct_id
        units = {u.feature_id: u for u in record.units}
        constructs.append({
            **record.construct,
            'features': [{
                'feature_id': f.feature_id, 'construct_id': cid, 'unit_id': f.unit_id,
                'source_feature_index': f.source_feature_index, 'source_feature_type': f.source_feature_type,
                'feature_label': f.feature_label, 'display_name': f.display_name, 'strand': f.strand,
                'segments': list(f.segments), 'location_status': f.location_status,
                'source_location': f.source_location, 'annotation_status': f.annotation_status,
                'annotation_ref': f'annotations.json#/constructs/{cid}/features/{f.feature_id}'} for f in record.features],
            'units': [{
                'unit_id': u.unit_id, 'construct_id': cid, 'feature_ids': [u.feature_id],
                'molecule_type': u.molecule_type, 'extraction_state': u.extraction_state,
                'reason_code': u.reason_code, 'detail': dict(u.detail),
                'protein': None if u.sequence is None else {'length_aa': u.length_aa, 'sequence_sha256': u.sequence_sha256,
                                                             'fasta_id': u.unit_id},
                'model_compatibility': u.model_compatibility} for u in record.units]})
        annotations[cid] = dict(record.annotations)
        for f in record.features:
            u = units[f.feature_id]
            feature_rows.append([cid, f.feature_id, f.unit_id, f.source_feature_index, f.source_feature_type,
                                 f.feature_label, f.display_name, f.names.get('gene'), f.names.get('product'),
                                 f.names.get('protein_id'), f.strand, len(f.segments), _segments_text(f.segments),
                                 f.location_status, f.source_location, f.annotation_status])
            unit_rows.append([cid, u.unit_id, u.feature_id, u.molecule_type, u.extraction_state, u.reason_code,
                              f.strand, _segments_text(f.segments), u.length_aa, u.sequence_sha256,
                              u.model_compatibility, record.source_record_id, record.input_name])
            if u.sequence:
                fasta.append(f'>{u.unit_id} construct={cid} feature={u.feature_id} len={u.length_aa} sha256={u.sequence_sha256}\n'
                             + '\n'.join(u.sequence[i:i+60] for i in range(0, len(u.sequence), 60)) + '\n')
    options = {'overwrite': result.options.overwrite, 'model_max_residues': result.options.model_max_residues,
               'provenance_supplied_for': sorted(result.options.provenance or {})}
    report = {
        'schema_version': BUNDLE_SCHEMA, 'tool': tool, 'configuration': options,
        'coordinate_convention': COORDINATES, 'exit_code': result.exit_code,
        'status': {0: 'complete', 2: 'partial', 1: 'no_usable_records'}[result.exit_code],
        'summary': result.summary(), 'reason_codes': REASONS, 'notes': list(NOTES),
        'inputs': [{'input_index': i.index, 'name': i.name, 'sha256': i.sha256, 'size_bytes': i.size_bytes,
                    'format': i.format} for i in result.inputs],
        'records': [{
            'input_index': r.input_index, 'record_index': r.record_index, 'input_name': r.input_name,
            'status': r.status, 'source_record_id': r.source_record_id, 'construct_id': r.construct_id,
            'duplicate_of': r.duplicate_of, 'failure': r.failure.as_dict() if r.failure else None,
            'warnings': [w.as_dict() for w in r.warnings],
            'counts': None if r.status != 'ok' else dict(sorted(Counter(u.extraction_state for u in r.units).items())),
            'units': [{'unit_id': u.unit_id, 'extraction_state': u.extraction_state, 'reason_code': u.reason_code,
                       'model_compatibility': u.model_compatibility} for u in r.units if u.molecule_type == 'protein'],
        } for r in result.records],
    }
    return {
        'normalized.json': _dump({'schema_version': BUNDLE_SCHEMA, 'tool': tool,
                                  'coordinate_convention': COORDINATES, 'constructs': constructs}),
        'annotations.json': _dump({'schema_version': BUNDLE_SCHEMA, 'constructs': annotations}),
        'proteins.fasta': ''.join(fasta),
        'features.tsv': _table(['construct_id', 'feature_id', 'unit_id', 'source_feature_index', 'source_feature_type',
                                'feature_label', 'display_name', 'gene', 'product', 'protein_id', 'strand',
                                'segment_count', 'segments', 'location_status', 'source_location', 'annotation_status'], feature_rows),
        'units.tsv': _table(['construct_id', 'unit_id', 'feature_id', 'molecule_type', 'extraction_state', 'reason_code',
                             'strand', 'segments', 'length_aa', 'protein_sha256', 'model_compatibility',
                             'source_record_id', 'input_name'], unit_rows),
        'report.json': _dump(report),
    }


def _check_output(output_dir: Path, inputs: Sequence[Path], overwrite: bool) -> None:
    resolved = output_dir.resolve()
    for path in inputs:
        source = Path(path).resolve()
        if source == resolved or resolved in source.parents:
            raise OutputError(f'Output directory would contain an input file ({path}); choose a different --out')
    if output_dir.exists():
        if not output_dir.is_dir():
            raise OutputError(f'Output path exists and is not a directory: {output_dir}')
        existing = sorted(p.name for p in output_dir.iterdir())
        if existing:
            is_bundle = (output_dir / 'report.json').is_file() and set(existing) <= set(BUNDLE_FILES)
            if not overwrite:
                raise OutputError(f'Output directory is not empty: {output_dir}. Choose a new directory'
                                  + (' or pass overwrite to replace this bundle.' if is_bundle else '.'))
            if not is_bundle:
                raise OutputError(f'Refusing to overwrite {output_dir}: it contains files that are not a previous bundle')


def write_bundle(result: BatchResult, output_dir, *, run_metadata: Mapping | None = None) -> BatchResult:
    """Write the bundle atomically (temporary sibling directory, then rename)."""
    output_dir = Path(output_dir)
    _check_output(output_dir, [i.path for i in result.inputs], result.options.overwrite)
    documents = bundle_documents(result)
    documents['run.json'] = _dump({'schema_version': BUNDLE_SCHEMA, **dict(run_metadata or {})})
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f'.{output_dir.name}.tmp-', dir=output_dir.parent))
    try:
        for name, text in documents.items():
            (staging / name).write_text(text, encoding='utf-8', newline='\n')
        if output_dir.exists():
            for name in BUNDLE_FILES:
                (output_dir / name).unlink(missing_ok=True)
            for name in documents:
                os.replace(staging / name, output_dir / name)
            staging.rmdir()
        else:
            os.replace(staging, output_dir)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return BatchResult(result.options, result.inputs, result.records, output_dir)


def _run_metadata(result: BatchResult, started: float, started_at: str) -> dict:
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    except Exception:
        peak = None
    try:
        import Bio
        biopython = Bio.__version__
    except Exception:
        biopython = None
    return {'tool': {'name': 'iri-constructs', 'version': __version__},
            'started_at_utc': started_at, 'finished_at_utc': datetime.now(timezone.utc).isoformat(),
            'wall_seconds': round(time.perf_counter() - started, 6),
            'peak_rss_platform_units': peak, 'python': platform.python_version(),
            'platform': platform.platform(), 'biopython': biopython,
            'input_paths': [str(i.path.resolve()) for i in result.inputs],
            'note': 'Variable run facts live here so that every other bundle file is reproducible byte for byte.'}


def preprocess_files(input_paths: Sequence, output_dir, *, options: PreprocessOptions | None = None) -> BatchResult:
    """Process files and/or directories into one bundle at ``output_dir``.

    Returns a BatchResult whose ``exit_code`` follows the CLI contract; library
    callers get the structured result instead of a process exit. Raises
    PreprocessError (or OutputError) when nothing can safely be written.
    """
    options = options or PreprocessOptions()
    started, started_at = time.perf_counter(), datetime.now(timezone.utc).isoformat()
    files = expand_inputs(input_paths)
    _check_output(Path(output_dir), files, options.overwrite)  # fail before doing any work
    result = process_inputs(files, options)
    return write_bundle(result, output_dir, run_metadata=_run_metadata(result, started, started_at))


def preprocess_file(input_path, output_dir, *, options: PreprocessOptions | None = None) -> BatchResult:
    """Process one file (which may hold several records) or one directory."""
    return preprocess_files([input_path], output_dir, options=options)
