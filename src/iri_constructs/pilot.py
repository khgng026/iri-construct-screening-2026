"""Reviewed T08 pilot preparation and Experiment A source-annotation agreement.

This audits preservation of supplied annotations, not biological correctness.
All sequence/qualifier artifacts remain inside the existing ignored output roots.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

from .annotations import json_safe, normalize_record
from .bench_inputs import read_source, reviewed_sources
from .runtime import (InputError, Measure, digest_file, failure, local_output,
                      read_json, require_text, run_metadata, write_json)
from .schema_check import FEATURE_TYPES, check_schema_example

COUNTS = (
    'source_features', 'source_cds', 'source_noncoding_features',
    'retained_features', 'retained_noncoding_features', 'original_qualifiers_preserved',
    'original_locations_preserved', 'coordinate_agreement_features', 'strand_agreement_features',
    'feature_parent_links_valid', 'unit_parent_links_valid', 'annotated_proteins_emitted',
    'annotated_translations_verified', 'missing_translation_cds', 'ambiguous_cds',
    'withheld_cds', 'unknown_features', 'ambiguous_features',
)


def _ratio(numerator: int, denominator: int):
    return numerator / denominator if denominator else None


def _rates(counts: dict) -> dict:
    return {
        'feature_retention': _ratio(counts['retained_features'], counts['source_features']),
        'qualifier_preservation': _ratio(counts['original_qualifiers_preserved'], counts['source_features']),
        'coordinate_agreement': _ratio(counts['coordinate_agreement_features'], counts['source_features']),
        'strand_agreement': _ratio(counts['strand_agreement_features'], counts['source_features']),
        'source_cds_usable_translation_coverage': _ratio(counts['annotated_proteins_emitted'], counts['source_cds']),
        'exported_translation_agreement': _ratio(counts['annotated_translations_verified'], counts['annotated_proteins_emitted']),
        'missing_translation_cds': _ratio(counts['missing_translation_cds'], counts['source_cds']),
        'withheld_cds': _ratio(counts['withheld_cds'], counts['source_cds']),
        'ambiguous_features': _ratio(counts['ambiguous_features'], counts['source_features']),
        'unknown_features': _ratio(counts['unknown_features'], counts['source_features']),
    }


def _source_signature(feature) -> dict:
    return {'type': feature.type, 'location': str(feature.location),
            'operator': getattr(feature.location, 'operator', None),
            'source_feature_id': feature.id, 'qualifiers': feature.qualifiers}


def _expected_spans(feature, parent_length: int) -> tuple[list[dict], dict]:
    """Independently inspect source positions without rounding uncertain bounds."""
    from Bio.SeqFeature import ExactPosition
    parts = [] if feature.location is None else feature.location.parts
    spans = []
    for part in parts:
        if (part.ref is not None or part.ref_db is not None or
                type(part.start) is not ExactPosition or type(part.end) is not ExactPosition):
            continue
        start, end = int(part.start), int(part.end)
        if not 0 <= start < end <= parent_length:
            raise InputError('Source contains invalid feature coordinates')
        spans.append({'start0': start, 'end0': end,
                      'strand': part.strand if part.strand in (-1, 1) else None})
    bounds = spans[0] if len(parts) == len(spans) == 1 else dict(start0=None, end0=None, strand=None)
    return spans, bounds


def audit_annotation_preservation(record, document: dict, audit: dict, extracted: dict) -> dict:
    """Check a normalized record against the original supplied feature semantics.

    Returned values contain counts/IDs only. Unknown coordinates can agree through
    explicit nulls plus preserved original location text; they are not exact bounds.
    """
    parent = document['construct']
    cid = parent['construct_id']
    normalized = {f['feature_id']: f for f in document['features']}
    units = document['units']
    by_feature = defaultdict(list)
    for unit in units:
        for fid in unit.get('feature_ids', []):
            by_feature[fid].append(unit)
    counts = {name: 0 for name in COUNTS}
    counts['source_features'] = len(record.features)
    counts['source_cds'] = sum(f.type == 'CDS' for f in record.features)
    counts['source_noncoding_features'] = len(record.features) - counts['source_cds']
    failed = []
    expected_ids = set()
    for original in record.features:
        signature = _source_signature(original)
        digest = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
        fid = cid + '-f-' + digest[:20]
        expected_ids.add(fid)
        feature = normalized.get(fid)
        saved = audit.get('features', {}).get(fid)
        if feature is None or saved is None:
            failed.append('source_feature_missing')
            continue
        counts['retained_features'] += 1
        if original.type != 'CDS':
            counts['retained_noncoding_features'] += 1
        if saved.get('qualifiers') == original.qualifiers:
            counts['original_qualifiers_preserved'] += 1
        else:
            failed.append('qualifier_preservation_failed')
        if all(saved.get(k) == signature[k] for k in ('type', 'location', 'operator', 'source_feature_id')):
            counts['original_locations_preserved'] += 1
        else:
            failed.append('original_feature_semantics_changed')
        spans, location = _expected_spans(original, len(record))
        if feature.get('location') == location and feature.get('segments', []) == spans:
            counts['coordinate_agreement_features'] += 1
        else:
            failed.append('coordinate_agreement_failed')
        if ([s.get('strand') for s in feature.get('segments', [])] == [s['strand'] for s in spans]
                and feature.get('location', {}).get('strand') == location['strand']):
            counts['strand_agreement_features'] += 1
        else:
            failed.append('strand_agreement_failed')
        if feature.get('construct_id') == cid:
            counts['feature_parent_links_valid'] += 1
        else:
            failed.append('feature_parent_link_failed')
        if feature.get('feature_label') != (original.type if original.type in FEATURE_TYPES else 'other'):
            failed.append('unsupported_feature_label_introduced')
        linked = by_feature[fid]
        if (len(linked) != 1 or linked[0].get('feature_ids') != [fid] or
                linked[0].get('construct_id') != cid or linked[0].get('source_segments') != spans):
            failed.append('unit_traceability_failed')
            continue
        unit = linked[0]
        counts['unit_parent_links_valid'] += 1
        state = feature.get('annotation_status')
        counts['unknown_features'] += state == 'unknown'
        counts['ambiguous_features'] += state == 'ambiguous'
        if original.type != 'CDS':
            if unit.get('molecule_type') != 'metadata_only' or unit.get('model_input_status') != 'unsupported':
                failed.append('noncoding_feature_misrepresented')
            if unit['unit_id'] in extracted:
                failed.append('noncoding_protein_invented')
            continue
        if unit.get('molecule_type') != 'protein':
            failed.append('cds_model_molecule_type_changed')
        counts['missing_translation_cds'] += not original.qualifiers.get('translation')
        counts['ambiguous_cds'] += state == 'ambiguous'
        sequence = extracted.get(unit['unit_id'])
        if sequence is None:
            counts['withheld_cds'] += 1
            if unit.get('model_input_status') != 'needs_review':
                failed.append('withheld_cds_marked_ready')
            continue
        counts['annotated_proteins_emitted'] += 1
        from Bio.SeqFeature import ExactPosition
        parts = [] if original.location is None else original.location.parts
        unresolved = (not parts or any(type(p.start) is not ExactPosition or type(p.end) is not ExactPosition
            or p.ref is not None or p.ref_db is not None or p.strand not in (-1, 1) for p in parts)
            or len({p.strand for p in parts}) != 1
            or getattr(original.location, 'operator', None) not in (None, 'join')
            or any(k in original.qualifiers for k in ('pseudo', 'pseudogene', 'transl_except', 'exception')))
        if unresolved:
            failed.append('unresolved_cds_exported')
        annotations = original.qualifiers.get('translation', [])
        expected = ''.join(annotations[0].split()).upper() if len(annotations) == 1 else None
        if expected and expected.endswith('*'):
            expected = expected[:-1]
        try:
            offset = int(original.qualifiers.get('codon_start', ['1'])[0]) - 1
            table = int(original.qualifiers.get('transl_table', ['1'])[0])
            validation = str(original.extract(record.seq)[offset:].translate(table=table, cds=True))
        except Exception:
            validation = None
        if not unresolved and expected is not None and sequence == expected == validation and state == 'known':
            counts['annotated_translations_verified'] += 1
        else:
            failed.append('annotated_translation_agreement_failed')
    if (len(normalized) != len(document['features']) or set(normalized) != expected_ids
            or set(audit.get('features', {})) != expected_ids):
        failed.append('extra_or_missing_feature_ids')
    unit_ids = [u['unit_id'] for u in units]
    if len(unit_ids) != len(set(unit_ids)) or len(units) != len(record.features):
        failed.append('extra_or_duplicate_units')
    if not set(extracted).issubset(unit_ids):
        failed.append('untraceable_protein_export')
    if (parent.get('sequence_length_nt') != len(record) or audit.get('construct_id') != cid
            or audit.get('source_record_id') != record.id):
        failed.append('parent_source_binding_failed')
    if json_safe(audit.get('record_annotations')) != json_safe(record.annotations):
        failed.append('record_annotation_metadata_changed')
    return {'construct_id': cid, 'primary_archetype': parent['primary_archetype'],
            'status': 'agrees_with_source' if not failed else 'failed',
            'counts': counts, 'rates': _rates(counts), 'failed_checks': sorted(set(failed)),
            'source_annotation_observable': bool(record.features)}


def _summarize(records: list[dict]) -> dict:
    counts = {name: sum(r['counts'][name] for r in records) for name in COUNTS}
    return {'records': len(records), 'counts': counts, 'rates': _rates(counts)}


def _validate_entries(config: dict, config_path: Path) -> list[dict]:
    if config.get('schema_version') != 't05-inputs-1':
        raise InputError('Expected schema_version t05-inputs-1')
    entries = config.get('records')
    if not isinstance(entries, list):
        raise InputError('records must be a list')
    parents, physical = set(), set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise InputError('Every records entry must be an object')
        for key in ('input_id', 'record_id', 'sequence_id'):
            value = require_text(entry, key)
            if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._:/|+\-]*', value) is None:
                raise InputError('Pilot identifiers must be nonempty metadata tokens without whitespace')
        parent = entry['record_id']
        identity = ((config_path.parent / require_text(entry, 'path')).resolve(), entry['sequence_id'])
        if parent in parents or identity in physical:
            raise InputError('Duplicate parent or source-record alias in pilot inputs')
        parents.add(parent)
        physical.add(identity)
        if entry.get('molecule_type') != 'dna':
            raise InputError('T08 annotation pilot requires parent DNA records, not protein-only inputs')
    return entries


def prepare_pilot(inputs_path: Path, output_path: Path) -> dict:
    """Prepare a reviewed pilot, returning and writing a truthful sequence-free report.

    An empty records list yields ``no_inputs`` without reading a manifest or source.
    All review gates are checked for the full panel before any parsing. Any pending
    source blocks admission; restrict the input config to approved records to run.
    """
    inputs_path = Path(inputs_path).resolve()
    output = local_output(output_path)
    if output.exists() and any(output.iterdir()):
        raise InputError('Pilot output directory must be new or empty')
    report = run_metadata('t08_pilot_annotation_preservation', {'input_contract': 't05-inputs-1'})
    report.update(status='blocked', inference_executed=False, scored_records=0,
                  evidence_class='source_annotation_agreement_only', records=[], input_exclusions=[])
    report['limitations'].extend([
        'Preservation agrees with supplied source annotations; it does not validate their biological correctness.',
        'Annotation/model eligibility approvals are recorded human decisions, never inferred by software.',
        'FASTA DNA has no documented feature table; agreement rates are null when the denominator is zero.',
    ])
    with Measure() as measured:
        try:
            config = read_json(inputs_path)
            report['input_config_sha256'] = digest_file(inputs_path)
            entries = _validate_entries(config, inputs_path)
            report['requested_records'] = len(entries)
            if not entries:
                report.update(status='no_inputs', reason='No reviewed local records supplied; candidate review/admission remains pending')
            else:
                _, sources, provenance = reviewed_sources(inputs_path)
                report['provenance'] = provenance
                if any(s.review['annotation_status'] != 'approved' for s in sources):
                    raise InputError('Pilot scoring requires approved T02 annotation review before source parsing')
                aliases = [(s.review['checksum_sha256'].lower(), s.sequence_id) for s in sources]
                if len(aliases) != len(set(aliases)):
                    raise InputError('Duplicate checksummed source-record alias in pilot inputs')
                docs, audits, proteins, successful = [], {}, {}, []
                output.mkdir(parents=True, exist_ok=True)
                for source in sources:
                    item = {'input_id': source.input_id, 'construct_id': source.review['record_id'],
                            'primary_archetype': source.review['archetype']}
                    try:
                        checksum = digest_file(source.path)
                        if checksum.lower() != source.review['checksum_sha256'].lower():
                            raise InputError('Source changed after review gate validation')
                        record = read_source(source)
                        if digest_file(source.path) != checksum:
                            raise InputError('Source changed while parsing; no outputs admitted')
                        doc, audit, extracted = normalize_record(record, source.review, checksum)
                        item.update(audit_annotation_preservation(record, doc, audit, extracted))
                        if item['status'] != 'agrees_with_source':
                            raise InputError('Annotation preservation audit failed; source requires manual inspection')
                        # Reuse the accepted T03 structural validator before publishing an artifact.
                        check_path = output / 'schema-validation.json'
                        write_json(check_path, doc)
                        issues = check_schema_example(check_path)
                        check_path.unlink()
                        if issues:
                            raise InputError('Normalized record failed accepted T03 structural validation')
                        docs.append(doc)
                        audits[doc['construct']['construct_id']] = json_safe(audit)
                        proteins.update({source.input_id + '/' + uid: seq for uid,seq in extracted.items()})
                        successful.append(item)
                    except Exception as exc:
                        item.update(status='failed', **failure(exc))
                        report['input_exclusions'].append({'input_id': source.input_id,
                                                            'reason': item['reason'], 'error_type': item['error_type']})
                    report['records'].append(item)
                grouped = defaultdict(list)
                for item in successful:
                    grouped[item['primary_archetype']].append(item)
                experiment = {'experiment': 'A_annotation_preservation', 'status':
                    ('evaluated_source_agreement' if any(i['source_annotation_observable'] for i in successful)
                     else 'no_observable_annotations' if successful else 'blocked'),
                    'claim': 'Agreement with documented source features only; no biological correctness/function inference',
                    'overall': _summarize(successful),
                    'by_archetype': {key: _summarize(group) for key,group in sorted(grouped.items())},
                    'failures_by_archetype': dict(Counter(i['primary_archetype'] for i in report['records'] if i['status'] == 'failed'))}
                report.update(scored_records=sum(i['source_annotation_observable'] for i in successful),
                              parsed_records=len(successful), experiment_a=experiment)
                if successful:
                    write_json(output / 'normalized.json', docs)
                    write_json(output / 'annotation-audit.json', audits)
                    write_json(output / 'experiment-a.json', experiment)
                    with (output / 'units.fasta').open('w', encoding='utf-8') as handle:
                        for iid, sequence in sorted(proteins.items()):
                            handle.write('>' + iid + '\n' + sequence + '\n')
                    report['artifacts'] = ['normalized.json', 'annotation-audit.json', 'experiment-a.json', 'units.fasta']
                report['status'] = 'ok' if len(successful) == len(sources) else 'partial' if successful else 'blocked'
        except Exception as exc:
            report.update(status='blocked', **failure(exc))
    report['resources_total'] = measured.data
    report['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    write_json(output / 'report.json', report)
    return report
