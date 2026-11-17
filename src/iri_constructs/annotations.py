"""Biopython records to normalized construct/feature/unit objects.

One normalizer serves two callers:

* ``normalize_record`` keeps the historical T03 v0.1 document and audit sidecar
  byte-for-byte for the frozen experiment runners. It still needs their review row.
* ``normalized_features`` is the review-free core. The installed preprocessing
  toolkit (``iri_constructs.preprocess``) calls it directly with ``strict=False``.

Protein outputs are only existing ``/translation`` qualifiers validated against
the documented CDS coordinates, frame and translation table. Nothing is inferred.
"""
from __future__ import annotations

import copy
import hashlib
import json
from urllib.parse import quote

from .runtime import InputError
from .schema_check import FEATURE_TYPES
from .sequences import validate_sequence

UNKNOWN_LOCATION = {'start0': None, 'end0': None, 'strand': None}

# Stable location reason codes and the historical free-text they replace.
LOCATION_ISSUES = {
    'location_missing': 'Missing location',
    'location_remote': 'Remote location withheld; no sequence reconstruction',
    'location_fuzzy': 'Fuzzy or unknown boundary; original location retained in audit',
    'location_unknown_strand': 'Unknown feature strand',
    'location_mixed_strand': 'Mixed strands; original segment order preserved without extraction',
    'location_unsupported_operator': 'Uncertain compound-location operator; no extraction',
    # Reported only by the non-strict toolkit path; the historical path raises instead.
    'location_out_of_bounds': 'Feature coordinates lie outside the parent sequence; no extraction',
    'location_zero_length': 'Zero-length (between-bases) location; retained as a site, no extraction',
}

# Stable CDS outcome codes -> (annotation status, historical note).
_INVALID = 'Invalid/incomplete CDS, translation alphabet, frame or translation table'
CDS_OUTCOMES = {
    'translation_validated': ('known', 'Source annotated translation verified against exact CDS coordinates/table'),
    'missing_translation': ('unknown', 'Missing annotated translation; no ORF inference or translation synthesis'),
    'exceptional_cds_annotation': ('ambiguous', 'Exceptional/pseudo CDS requires manual review'),
    'conflicting_cds_qualifiers': ('ambiguous', 'Multiple conflicting CDS qualifiers'),
    'translation_mismatch': ('ambiguous', 'Annotated translation disagrees with the documented CDS coordinates/table'),
    'empty_translation': ('ambiguous', _INVALID),
    'noncanonical_translation_residues': ('ambiguous', _INVALID),
    'invalid_codon_start': ('ambiguous', _INVALID),
    'invalid_transl_table': ('ambiguous', _INVALID),
    'incomplete_cds': ('ambiguous', _INVALID),
}
NONCODING_NOTE = 'Noncoding/source feature retained as unsupported metadata'

_POSITION_KINDS = {'ExactPosition': 'exact', 'BeforePosition': 'before', 'AfterPosition': 'after',
                   'WithinPosition': 'within', 'BetweenPosition': 'between', 'OneOfPosition': 'one_of',
                   'UnknownPosition': 'unknown', 'UncertainPosition': 'uncertain'}


def _position(position) -> tuple[int | None, str, str]:
    kind = _POSITION_KINDS.get(type(position).__name__, 'other')
    try:
        value = int(position)
    except (TypeError, ValueError):
        value = None
    return value, kind, str(position)


def location_details(feature, length: int, *, strict: bool = True) -> tuple[dict, list[dict], str | None, list[dict]]:
    """Return (single-span location, ordered exact segments, issue code, original parts).

    ``parts`` describes every source location part, including fuzzy/remote ones
    that are deliberately absent from ``segments``. With ``strict`` an
    out-of-bounds or zero-length part raises (historical behavior); otherwise it
    becomes a reason code and the remaining record is still processed.
    """
    from Bio.SeqFeature import ExactPosition
    loc = feature.location
    if loc is None:
        return dict(UNKNOWN_LOCATION), [], 'location_missing', []
    segments, parts, issue = [], [], None
    for part in loc.parts:
        start, start_kind, start_text = _position(part.start)
        end, end_kind, end_text = _position(part.end)
        strand = part.strand if part.strand in (-1, 1) else None
        ref = None if part.ref is None and part.ref_db is None else ':'.join(
            str(x) for x in (part.ref_db, part.ref) if x is not None)
        parts.append({'start0': start, 'end0': end, 'start_kind': start_kind, 'end_kind': end_kind,
                      'start_text': start_text, 'end_text': end_text, 'strand': strand, 'ref': ref})
        if ref is not None:
            issue = 'location_remote'
            continue
        if type(part.start) is not ExactPosition or type(part.end) is not ExactPosition:
            issue = 'location_fuzzy'
            continue
        if not 0 <= start < end <= length:
            if strict:
                raise InputError('Feature coordinates are outside the parent or have zero length')
            issue = 'location_zero_length' if 0 <= start == end <= length else 'location_out_of_bounds'
            continue
        segments.append({'start0': start, 'end0': end, 'strand': strand})
        if strand is None:
            issue = 'location_unknown_strand'
    if len({s['strand'] for s in segments}) > 1:
        issue = 'location_mixed_strand'
    if getattr(loc, 'operator', None) not in (None, 'join'):
        issue = 'location_unsupported_operator'
    # A joined annotation is never represented by its misleading bounding interval.
    location = segments[0].copy() if len(loc.parts) == 1 and len(segments) == 1 else dict(UNKNOWN_LOCATION)
    return location, segments, issue, parts


def feature_location(feature, length: int) -> tuple[dict, list[dict], str | None]:
    """Historical interface: free-text issue, raising on out-of-bounds coordinates."""
    location, segments, code, _ = location_details(feature, length, strict=True)
    return location, segments, LOCATION_ISSUES.get(code)


def _cds_defects(coding: str, table_id: int) -> list[str]:
    """Explain why a strict complete-CDS translation was refused (diagnostic only)."""
    from Bio.Data import CodonTable
    table = CodonTable.ambiguous_dna_by_id[table_id]
    defects = []
    if len(coding) % 3:
        defects.append('length_not_multiple_of_three')
    if coding[:3] not in table.start_codons:
        defects.append('no_start_codon')
    in_frame = coding[:len(coding) - len(coding) % 3]
    if len(in_frame) < 6 or in_frame[-3:] not in table.stop_codons:
        defects.append('no_terminal_stop_codon')
    return defects


def check_cds_translation(record, feature, location_code: str | None) -> dict:
    """Validate an existing ``/translation``; never create a protein.

    Returns ``protein`` (or None), a stable ``code`` and diagnostic ``detail``.
    ``detail['in_frame_translation_agrees']`` reports, for an incomplete CDS,
    whether the source translation still matches the documented frame. That is a
    diagnostic for the user; such a unit is withheld all the same.
    """
    q = feature.qualifiers
    if location_code:
        return {'protein': None, 'code': location_code, 'detail': {}}
    exceptional = sorted(k for k in ('pseudo', 'pseudogene', 'transl_except', 'exception') if k in q)
    if exceptional:
        return {'protein': None, 'code': 'exceptional_cds_annotation', 'detail': {'qualifiers': exceptional}}
    translations = q.get('translation', [])
    if not translations:
        return {'protein': None, 'code': 'missing_translation', 'detail': {}}
    repeated = sorted(k for k in ('product', 'gene', 'codon_start', 'transl_table') if len(q.get(k, [])) > 1)
    if len(translations) != 1 or repeated:
        return {'protein': None, 'code': 'conflicting_cds_qualifiers',
                'detail': {'qualifiers': (['translation'] if len(translations) != 1 else []) + repeated}}
    sequence = ''.join(translations[0].split()).upper()
    if sequence.endswith('*'):
        sequence = sequence[:-1]
    if not sequence:
        return {'protein': None, 'code': 'empty_translation', 'detail': {}}
    try:
        validate_sequence(sequence, 'protein')
    except InputError:
        return {'protein': None, 'code': 'noncanonical_translation_residues', 'detail': {}}
    try:
        offset = int(q.get('codon_start', ['1'])[0])
        if offset not in (1, 2, 3):
            raise ValueError
    except (TypeError, ValueError):
        return {'protein': None, 'code': 'invalid_codon_start', 'detail': {}}
    try:
        table = int(q.get('transl_table', ['1'])[0])
    except (TypeError, ValueError):
        return {'protein': None, 'code': 'invalid_transl_table', 'detail': {}}
    try:
        coding = feature.extract(record.seq)[offset-1:]
        # This validates the provided translation; it never creates a new protein unit.
        translated = str(coding.translate(table=table, cds=True))
    except Exception:
        detail = {}
        try:
            text = str(coding).upper()
            detail['cds_defects'] = _cds_defects(text, table)
            loose = str(coding[:len(coding) - len(coding) % 3].translate(table=table))
            if loose.endswith('*'):
                loose = loose[:-1]
            if '*' in loose:
                detail['cds_defects'].append('internal_stop_codon')
            detail['in_frame_translation_agrees'] = loose == sequence or 'M' + loose[1:] == sequence
        except KeyError:
            return {'protein': None, 'code': 'invalid_transl_table', 'detail': {}}
        except Exception:
            detail = {}
        return {'protein': None, 'code': 'incomplete_cds', 'detail': detail}
    if translated != sequence:
        return {'protein': None, 'code': 'translation_mismatch', 'detail': {}}
    return {'protein': sequence, 'code': 'translation_validated', 'detail': {}}


def annotated_translation(record, feature, location_issue: str | None) -> tuple[str | None, str, str]:
    """Historical interface: (protein, status, free-text note) for a free-text location issue."""
    if location_issue:
        return None, 'ambiguous', location_issue
    result = check_cds_translation(record, feature, None)
    status, note = CDS_OUTCOMES[result['code']]
    return result['protein'], status, note


def feature_signature(feature) -> tuple[dict, str]:
    """Original source payload and its content digest (the basis of stable IDs)."""
    original = {'type': feature.type, 'location': str(feature.location),
                'operator': getattr(feature.location, 'operator', None),
                'source_feature_id': feature.id, 'qualifiers': copy.deepcopy(feature.qualifiers)}
    return original, hashlib.sha256(json.dumps(original, sort_keys=True).encode()).hexdigest()


def normalized_features(record, sequence_length: int, *, strict: bool = True) -> list[dict]:
    """Review-free core: one entry per source feature, in source order.

    Each entry carries the content ``signature``, the ``original`` payload, exact
    ordered ``segments``, all source location ``parts``, a stable ``code`` with
    its ``status``/``note`` and, for a validated CDS, the source ``protein``.
    ``duplicate_index`` is 1 for the first occurrence of identical content. With
    ``strict`` an identical duplicate raises (historical behavior); otherwise
    later copies are kept and numbered so nothing is silently merged.
    """
    entries, seen = [], {}
    for index, feature in enumerate(record.features):
        original, signature = feature_signature(feature)
        seen[signature] = seen.get(signature, 0) + 1
        if strict and seen[signature] > 1:
            raise InputError('Identical duplicate source features require manual disambiguation')
        location, segments, location_code, parts = location_details(feature, sequence_length, strict=strict)
        display = feature.qualifiers.get('label') or feature.qualifiers.get('product') or feature.qualifiers.get('gene')
        entry = {'source_index': index, 'signature': signature, 'duplicate_index': seen[signature],
                 'original': original, 'type': feature.type,
                 'label': feature.type if feature.type in FEATURE_TYPES else 'other',
                 'display_name': display[0] if display else None,
                 'location': location, 'segments': segments, 'parts': parts,
                 'location_code': location_code, 'is_cds': feature.type == 'CDS',
                 'protein': None, 'detail': {}}
        if entry['is_cds']:
            result = check_cds_translation(record, feature, location_code)
            entry.update(protein=result['protein'], code=result['code'], detail=result['detail'])
            if location_code:
                entry.update(status='ambiguous', note=LOCATION_ISSUES[location_code])
            else:
                entry['status'], entry['note'] = CDS_OUTCOMES[result['code']]
        elif location_code:
            entry.update(code=location_code, status='ambiguous', note=LOCATION_ISSUES[location_code])
        else:
            entry.update(code='noncoding_metadata', status='known', note=NONCODING_NOTE)
        entries.append(entry)
    return entries


def construct_annotation_status(statuses) -> str:
    statuses = set(statuses)
    return 'ambiguous' if 'ambiguous' in statuses else 'unknown' if not statuses or 'unknown' in statuses else 'known'


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
    for entry in normalized_features(record, len(sequence), strict=True):
        fid, uid = cid + '-f-' + entry['signature'][:20], cid + '-u-' + entry['signature'][:20]
        status, note, protein, segments = entry['status'], entry['note'], entry['protein'], entry['segments']
        normalized = {'feature_id': fid, 'construct_id': cid, 'feature_label': entry['label'],
                      'display_name': entry['display_name'], 'location': entry['location'],
                      'segments': segments, 'annotation_source': 'annotation-audit.json#/' + quote(cid, safe='') + '/features/' + quote(fid, safe=''),
                      'annotation_status': status}
        if status in ('unknown', 'ambiguous'):
            normalized['ambiguity_note'] = note
        document['features'].append(normalized)
        document['units'].append({'unit_id': uid, 'construct_id': cid, 'feature_ids': [fid],
                                 'source_segments': segments, 'molecule_type': 'protein' if entry['is_cds'] else 'metadata_only',
                                 'model_input_status': ('pending_model_selection' if protein else 'needs_review') if entry['is_cds'] else 'unsupported',
                                 'representation_note': note})
        audit['features'][fid] = entry['original']
        if protein:
            extracted[uid] = protein
    document['features'].sort(key=lambda f: f['feature_id'])
    document['units'].sort(key=lambda u: u['unit_id'])
    parent['annotation_status'] = construct_annotation_status(f['annotation_status'] for f in document['features'])
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
