"""Biopython records to approved T03 v0.1 objects plus a separate local audit sidecar."""
from __future__ import annotations

import copy
import hashlib
import json
from urllib.parse import quote

from .bench_inputs import validate_sequence
from .runtime import InputError
from .schema_check import FEATURE_TYPES

UNKNOWN_LOCATION = {'start0': None, 'end0': None, 'strand': None}


def feature_location(feature, length: int) -> tuple[dict, list[dict], str | None]:
    from Bio.SeqFeature import ExactPosition
    loc = feature.location
    if loc is None:
        return dict(UNKNOWN_LOCATION), [], 'Missing location'
    segments, issue = [], None
    for part in loc.parts:
        if part.ref is not None or part.ref_db is not None:
            issue = 'Remote location withheld; no sequence reconstruction'
            continue
        if type(part.start) is not ExactPosition or type(part.end) is not ExactPosition:
            issue = 'Fuzzy or unknown boundary; original location retained in audit'
            continue
        start, end = int(part.start), int(part.end)
        if not 0 <= start < end <= length:
            raise InputError('Feature coordinates are outside the parent or have zero length')
        strand = part.strand if part.strand in (-1, 1) else None
        segments.append({'start0': start, 'end0': end, 'strand': strand})
        if strand is None:
            issue = 'Unknown feature strand'
    if len({s['strand'] for s in segments}) > 1:
        issue = 'Mixed strands; original segment order preserved without extraction'
    if getattr(loc, 'operator', None) not in (None, 'join'):
        issue = 'Uncertain compound-location operator; no extraction'
    # A joined annotation is never represented by its misleading bounding interval.
    location = segments[0].copy() if len(loc.parts) == 1 and len(segments) == 1 else dict(UNKNOWN_LOCATION)
    return location, segments, issue


def annotated_translation(record, feature, location_issue: str | None) -> tuple[str | None, str, str]:
    q = feature.qualifiers
    if location_issue:
        return None, 'ambiguous', location_issue
    if any(k in q for k in ('pseudo', 'pseudogene', 'transl_except', 'exception')):
        return None, 'ambiguous', 'Exceptional/pseudo CDS requires manual review'
    translations = q.get('translation', [])
    if not translations:
        return None, 'unknown', 'Missing annotated translation; no ORF inference or translation synthesis'
    if len(translations) != 1 or any(len(q.get(k, [])) > 1 for k in ('product', 'gene', 'codon_start', 'transl_table')):
        return None, 'ambiguous', 'Multiple conflicting CDS qualifiers'
    sequence = ''.join(translations[0].split()).upper()
    if sequence.endswith('*'):
        sequence = sequence[:-1]
    try:
        validate_sequence(sequence, 'protein')
        offset = int(q.get('codon_start', ['1'])[0])
        table = int(q.get('transl_table', ['1'])[0])
        if offset not in (1, 2, 3):
            raise ValueError
        # This validates the provided translation; it never creates a new protein unit.
        translated = str(feature.extract(record.seq)[offset-1:].translate(table=table, cds=True))
        if translated != sequence:
            return None, 'ambiguous', 'Annotated translation disagrees with the documented CDS coordinates/table'
    except Exception:
        return None, 'ambiguous', 'Invalid/incomplete CDS, translation alphabet, frame or translation table'
    return sequence, 'known', 'Source annotated translation verified against exact CDS coordinates/table'


def normalize_record(record, review: dict, source_sha256: str) -> tuple[dict, dict, dict]:
    """Pure normalization; callers enforce the T02 gate before reading source records.

    Tests use explicitly fictional in-memory records, never approve a source record.
    Protein outputs are only existing /translation qualifiers validated against the CDS.
    """
    sequence = str(record.seq).upper()
    validate_sequence(sequence, 'dna')
    cid = review['record_id']
    parent = {
        'construct_id': cid, 'source': review['source'], 'source_url': review['source_url'],
        'accession': review['accession'], 'accession_version': review['accession_version'],
        'sequence_length_nt': len(sequence),
        'topology': record.annotations.get('topology', 'unknown'),
        'primary_archetype': review['archetype'], 'secondary_archetypes': [],
        'annotation_status': 'unknown', 'license_note': review['notes'], 'review_note': review['notes'],
    }
    if parent['topology'] not in ('linear', 'circular'):
        parent['topology'] = 'unknown'
    document = {'schema_version': '0.1', 'construct': parent, 'features': [], 'units': []}
    audit = {'audit_schema_version': 't05-annotation-audit-1', 'construct_id': cid,
             'source_sha256': source_sha256, 'source_record_id': record.id,
             'record_annotations': copy.deepcopy(record.annotations), 'features': {}}
    extracted = {}
    signatures = set()
    for feature in record.features:
        original = {'type': feature.type, 'location': str(feature.location),
                    'operator': getattr(feature.location, 'operator', None),
                    'source_feature_id': feature.id, 'qualifiers': copy.deepcopy(feature.qualifiers)}
        signature = hashlib.sha256(json.dumps(original, sort_keys=True).encode()).hexdigest()
        if signature in signatures:
            raise InputError('Identical duplicate source features require manual disambiguation')
        signatures.add(signature)
        fid, uid = cid + '-f-' + signature[:20], cid + '-u-' + signature[:20]
        location, segments, issue = feature_location(feature, len(sequence))
        label = feature.type if feature.type in FEATURE_TYPES else 'other'
        display = feature.qualifiers.get('label') or feature.qualifiers.get('product') or feature.qualifiers.get('gene')
        status, note, protein = 'known', 'Noncoding/source feature retained as unsupported metadata', None
        if feature.type == 'CDS':
            protein, status, note = annotated_translation(record, feature, issue)
        elif issue:
            status, note = 'ambiguous', issue
        normalized = {'feature_id': fid, 'construct_id': cid, 'feature_label': label,
                      'display_name': display[0] if display else None, 'location': location,
                      'segments': segments, 'annotation_source': 'annotation-audit.json#/' + quote(cid, safe='') + '/features/' + quote(fid, safe=''),
                      'annotation_status': status}
        if status in ('unknown', 'ambiguous'):
            normalized['ambiguity_note'] = note
        document['features'].append(normalized)
        document['units'].append({'unit_id': uid, 'construct_id': cid, 'feature_ids': [fid],
                                 'source_segments': segments, 'molecule_type': 'protein' if feature.type == 'CDS' else 'metadata_only',
                                 'model_input_status': ('pending_model_selection' if protein else 'needs_review') if feature.type == 'CDS' else 'unsupported',
                                 'representation_note': note})
        audit['features'][fid] = original
        if protein:
            extracted[uid] = protein
    document['features'].sort(key=lambda f: f['feature_id'])
    document['units'].sort(key=lambda u: u['unit_id'])
    statuses = {f['annotation_status'] for f in document['features']}
    parent['annotation_status'] = 'ambiguous' if 'ambiguous' in statuses else 'unknown' if not statuses or 'unknown' in statuses else 'known'
    return document, audit, extracted


def json_safe(value):
    """Preserve audit strings/containers, including Biopython reference objects locally."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)
