"""Sequence-free T06 candidate discovery and manual review queue, not admission."""
from __future__ import annotations

from collections import Counter
from datetime import date
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse

from .runtime import InputError, read_json, require_text
from .schema_check import ARCHETYPES

GATES = ('annotation_status', 'review_status', 'licensing_status', 'safety_review')
UNREVIEWED = {'pending', 'uncertain', 'rejected'}
TEXT_FIELDS = ('record_id', 'construct_name', 'source', 'source_url', 'accession',
               'accession_version', 'record_version_note', 'archetype_evidence',
               'source_metadata_status', 'source_reuse_note', 'terms_url',
               'annotation_completeness', 'inclusion_reason', 'exclusion_reason',
               'reviewer_notes', 'independence_group')
LIST_FIELDS = ('secondary_archetypes', 'documented_cds', 'other_annotations', 'evidence_urls')
CHECKLIST = (
    'Confirm exact public material/accession identity and immutable source-file version.',
    'Researcher records ordinary benign rationale and source-specific safety review.',
    'Researcher checks sequence-data access/reuse/redistribution terms; material MTA alone is insufficient.',
    'Acquire only a permitted reviewed local file; record retrieval date and SHA256.',
    'Review complete feature table, parent nucleotide length, original CDS/noncoding qualifiers and ambiguity.',
    'Resolve variant/tag/boundary discrepancies; document unknown states rather than guessing.',
    'Assign independent reference groups and evidence-supported family labels separately.',
    'Create a distinct complete T02 admission row and t05-inputs-1 source mapping; do not relabel this inventory as approved.',
)


def _public_url(value: str) -> bool:
    parsed = urlparse(value)
    return parsed.scheme in ('http', 'https') and bool(parsed.hostname) and not parsed.username and not parsed.password


def validate_candidate_records(records: list) -> None:
    """Validate discovery metadata, rejecting approvals or sequence payloads.

    Review outcomes belong in the separate T02 admission manifest. This schema
    deliberately cannot authorize a record for parsing or scored experiments.
    """
    if not isinstance(records, list) or not records:
        raise InputError('Candidate inventory records must be a nonempty list')
    seen_ids, seen_sources, priorities = set(), set(), set()
    for row in records:
        if not isinstance(row, dict):
            raise InputError('Candidate entries must be objects')
        required = set(TEXT_FIELDS) | set(LIST_FIELDS) | set(GATES) | {
            'primary_archetype', 'candidate_status', 'retrieved_on', 'sequence_length_nt',
            'feature_count', 'checksum_sha256', 'pilot_priority',
        }
        if required - row.keys():
            raise InputError('Candidate metadata fields are missing; unknown values must be explicit')
        for field in TEXT_FIELDS:
            require_text(row, field)
        rid = row['record_id']
        if any(c.isspace() for c in rid):
            raise InputError('Candidate record_id must not contain whitespace')
        if rid in seen_ids:
            raise InputError('Duplicate candidate record_id')
        seen_ids.add(rid)
        source_key = (row['source'], row['accession'], row['accession_version'])
        if source_key in seen_sources:
            raise InputError('Duplicate candidate source identity/version')
        seen_sources.add(source_key)
        if row.get('primary_archetype') not in ARCHETYPES:
            raise InputError('Candidate primary_archetype is invalid')
        for field in LIST_FIELDS:
            if not isinstance(row.get(field), list) or any(not isinstance(x, str) or not x.strip() for x in row[field]):
                raise InputError('Candidate list metadata is missing or invalid: ' + field)
        if any(x not in ARCHETYPES for x in row['secondary_archetypes']):
            raise InputError('Candidate secondary_archetypes are invalid')
        if len(set(row['secondary_archetypes'])) != len(row['secondary_archetypes']):
            raise InputError('Repeated candidate secondary archetype')
        if not row['evidence_urls'] or any(not _public_url(url) for url in row['evidence_urls']):
            raise InputError('Candidate evidence requires public HTTP(S) URLs without credentials')
        if not _public_url(row['source_url']) or not _public_url(row['terms_url']):
            raise InputError('Candidate source/terms URLs must be public HTTP(S) references without credentials')
        for field in GATES:
            if row.get(field) not in UNREVIEWED:
                raise InputError('Discovery inventory cannot contain approved/invalid review gates: ' + field)
        if row.get('candidate_status') not in UNREVIEWED:
            raise InputError('Candidate status must remain pending, uncertain or rejected')
        try:
            date.fromisoformat(row.get('retrieved_on', ''))
        except (TypeError, ValueError) as exc:
            raise InputError('Candidate retrieved_on must be YYYY-MM-DD') from exc
        for field, minimum in (('sequence_length_nt', 1), ('feature_count', 0)):
            value = row.get(field)
            if value is not None and (type(value) is not int or value < minimum):
                raise InputError('Candidate count must be a valid integer or explicit null: ' + field)
        if row.get('checksum_sha256') is not None:
            digest = row['checksum_sha256']
            if not isinstance(digest, str) or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
                raise InputError('Candidate checksum must be SHA256 or explicit null')
        priority = row.get('pilot_priority')
        if priority is not None:
            if type(priority) is not int or priority < 1 or priority in priorities:
                raise InputError('Pilot review priorities must be unique positive integers or null')
            priorities.add(priority)
        forbidden = {'sequence', 'protein', 'dna', 'translation', 'embeddings', 'qualifiers', 'features'} & row.keys()
        if forbidden:
            raise InputError('Candidate inventory is metadata-only; raw sequence/feature payload fields are forbidden')


def load_candidate_inventory(path: Path) -> list[dict]:
    """Read the separately versioned candidate schema without contacting sources."""
    doc = read_json(path)
    if doc.get('schema_version') != 't06-candidates-1':
        raise InputError('Expected candidate schema_version t06-candidates-1')
    for key in ('inventory_version', 'metadata_retrieved_on', 'scope', 'provenance_method'):
        require_text(doc, key)
    try:
        date.fromisoformat(doc['metadata_retrieved_on'])
    except ValueError as exc:
        raise InputError('Inventory metadata_retrieved_on must be YYYY-MM-DD') from exc
    rows = doc.get('records')
    validate_candidate_records(rows)
    return rows


def review_queue(records: list[dict]) -> dict:
    """Return explicit blockers and priorities, never an approved source manifest."""
    validate_candidate_records(records)
    entries = []
    for row in sorted(records, key=lambda r: (r['pilot_priority'] is None, r['pilot_priority'] or 0, r['record_id'])):
        blockers = [field + ':' + row[field] for field in GATES if row[field] != 'approved']
        for field in ('sequence_length_nt', 'feature_count', 'checksum_sha256'):
            if row.get(field) is None:
                blockers.append('missing_reviewed_' + field)
        if row['accession_version'] == 'unknown':
            blockers.append('immutable_source_version_not_selected')
        entries.append({'record_id': row['record_id'], 'construct_name': row['construct_name'],
                        'source_url': row['source_url'], 'primary_archetype': row['primary_archetype'],
                        'pilot_priority': row['pilot_priority'], 'independence_group': row['independence_group'],
                        'status': 'excluded_from_scoring', 'blockers': blockers,
                        'reason': row['exclusion_reason'], 'reviewer_notes': row['reviewer_notes']})
    return {'review_schema_version': 't06-review-queue-1', 'status': 'pending_manual_review',
            'candidate_count': len(records), 'approved_for_scoring': 0,
            'proposed_pilot_count': sum(r['pilot_priority'] is not None for r in records),
            'proposed_pilot_ids': [e['record_id'] for e in entries if e['pilot_priority'] is not None],
            'primary_archetype_counts': dict(sorted(Counter(r['primary_archetype'] for r in records).items())),
            'manual_checklist': list(CHECKLIST), 'entries': entries,
            'limitations': ['Discovery metadata does not admit any biological record.',
                            'Candidate relatedness groups are conservative review aids, not validated homology clusters.',
                            'Unknown counts/versions are explicit and cannot be used to fabricate a T02 admission row.']}


def build_pilot_review_queue(inventory_path: Path) -> dict:
    """Derive a shareable pending-priority queue, never a scored biological dataset."""
    queue = review_queue(load_candidate_inventory(inventory_path))
    proposed = [entry for entry in queue['entries'] if entry['pilot_priority'] is not None]
    canonical = json.dumps(read_json(inventory_path),sort_keys=True,separators=(',', ':'),allow_nan=False).encode('utf-8')
    return {'schema_version': 't08-pilot-review-1', 'inventory_sha256': hashlib.sha256(canonical).hexdigest(),
            'inventory_hash_policy': 'canonical-json-v1: sorted keys, compact separators, UTF-8; independent of Git line endings',
            'status': 'pending_manual_review', 'record_count': len(proposed),
            'approved_record_count': 0, 'experimental_input_ready': False,
            'primary_archetype_counts': dict(sorted(Counter(e['primary_archetype'] for e in proposed).items())),
            'records': proposed, 'manual_checklist': queue['manual_checklist'],
            'release_scope': 'Public discovery metadata only; no source sequence/feature tables or approved dataset',
            'limitations': queue['limitations']}


def validate_pilot_review_queue(path: Path, inventory_path: Path) -> None:
    if read_json(path) != build_pilot_review_queue(inventory_path):
        raise InputError('Pilot queue is stale or altered; regenerate metadata from the candidate inventory, not approvals')
