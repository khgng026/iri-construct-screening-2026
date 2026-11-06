"""E1 source-annotation quality and appropriate abstention; no model inference.

Quality strata depend on reviewed source annotations, never retrieval outcomes.
The executable reuses the existing pilot's T02 gate and preservation audit.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from .annotations import annotated_translation, feature_location
from .bench_inputs import read_source, reviewed_sources
from .pilot import audit_annotation_preservation, prepare_pilot
from .runtime import (InputError, Measure, digest_file, failure, local_output,
                      read_json, run_metadata, write_json)
from .schema_check import FEATURE_TYPES


def _fraction(numerator, denominator):
    return {'numerator': numerator, 'denominator': denominator,
            'value': numerator / denominator if denominator else None}


def classify_quality(source_quality: dict) -> str:
    """Classify source metadata only, independent of parser/model success.

    Exact multipart and reverse-strand annotations are not quality defects.
    A source-only or no-CDS record is annotation-limited for protein extraction.
    """
    for name in ('source_features', 'source_cds', 'source_issue_features'):
        if type(source_quality.get(name)) is not int or source_quality[name] < 0:
            raise InputError('Quality classification needs nonnegative source annotation counts')
    if source_quality['source_cds'] > source_quality['source_features']:
        raise InputError('Source CDS count exceeds source feature count')
    if source_quality['source_cds'] == 0:
        return 'annotation_limited'
    return 'partially_annotated' if source_quality['source_issue_features'] else 'well_annotated'


def _source_summary(record):
    """Inspect originals; reuse the accepted source-translation validation policy."""
    from Bio.SeqFeature import ExactPosition
    counts = Counter(source_features=len(record.features), source_cds=0,
                     observable_features=0, source_issue_features=0,
                     missing_location_features=0, fuzzy_location_features=0,
                     remote_location_features=0, unknown_strand_features=0,
                     multipart_features=0, reverse_strand_features=0,
                     mixed_strand_features=0, invalid_location_features=0,
                     missing_translation_cds=0, conflicting_qualifier_cds=0,
                     translation_conflict_cds=0, exceptional_cds=0,
                     validated_source_translation_cds=0, source_noncoding_features=0,
                     source_unknown_features=0, source_ambiguous_features=0)
    entries = []
    for source in record.features:
        parts = [] if source.location is None else source.location.parts
        flags = set()
        if not parts:
            flags.add('missing_location')
        if any(type(part.start) is not ExactPosition or type(part.end) is not ExactPosition for part in parts):
            flags.add('fuzzy_location')
        if any(part.ref is not None or part.ref_db is not None for part in parts):
            flags.add('remote_location')
        if any(part.strand not in (-1, 1) for part in parts):
            flags.add('unknown_strand')
        if len({part.strand for part in parts}) > 1:
            flags.add('mixed_strand')
        counts['multipart_features'] += len(parts) > 1
        counts['reverse_strand_features'] += any(part.strand == -1 for part in parts)
        for key in flags:
            counts[key + '_features'] += 1
        try:
            _, _, issue = feature_location(source, len(record))
        except InputError:
            issue = 'Invalid original source location'
            flags.add('invalid_location')
            counts['invalid_location_features'] += 1
        q, expected, status, reason = source.qualifiers, None, 'known', None
        if source.type == 'CDS':
            counts['source_cds'] += 1
            counts['missing_translation_cds'] += not q.get('translation')
            conflict = any(len(q.get(key, [])) > 1 for key in
                           ('translation', 'product', 'gene', 'codon_start', 'transl_table'))
            counts['conflicting_qualifier_cds'] += conflict
            counts['exceptional_cds'] += any(key in q for key in ('pseudo', 'pseudogene', 'transl_except', 'exception'))
            expected, status, reason = annotated_translation(record, source, issue)
            counts['validated_source_translation_cds'] += expected is not None
            if expected is None:
                flags.add('unsupported_cds_annotation')
                counts['translation_conflict_cds'] += bool(q.get('translation')) and not conflict and not issue
        elif issue:
            status, reason = 'ambiguous', issue
        observable = source.type != 'source'
        counts['observable_features'] += observable
        counts['source_noncoding_features'] += source.type != 'CDS'
        counts['source_unknown_features'] += observable and status == 'unknown'
        counts['source_ambiguous_features'] += observable and status == 'ambiguous'
        counts['source_issue_features'] += observable and bool(flags or issue)
        signature = {'type': source.type, 'location': str(source.location),
                     'operator': getattr(source.location, 'operator', None),
                     'source_feature_id': source.id, 'qualifiers': source.qualifiers}
        suffix = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:20]
        entries.append({'suffix': suffix, 'type': source.type, 'observable': observable,
                        'expected': expected, 'status': status, 'reason': reason,
                        'needs_uncertainty': bool(flags or issue), 'original': source})
    return dict(counts), entries


def audit_quality(record, document: dict, audit: dict, extracted: dict, *,
                  max_residues=1022, source_sha256=None) -> dict:
    """Audit one already admitted source and its actual normalized/export artifacts.

    Low-level helper does not approve records. Production entry point is
    run_annotation_quality, which enforces all T02 gates through prepare_pilot.
    Report contains IDs, reasons and counts; raw qualifiers/sequences stay local.
    """
    if type(max_residues) is not int or not 1 <= max_residues <= 1022:
        raise InputError('E1 max_residues must be in 1..1022')
    source_counts, originals = _source_summary(record)
    cid = document['construct']['construct_id']
    try:
        preservation = audit_annotation_preservation(record, document, audit, extracted)
    except (InputError, KeyError, TypeError):
        preservation = {'status': 'failed', 'failed_checks': ['invalid_source_or_normalized_structure']}
    features = {f['feature_id']: f for f in document['features']}
    expected_ids = {cid + '-f-' + f['suffix'] for f in originals}
    invented = sum(f['feature_id'] not in expected_ids for f in document['features'])
    invented += len(document['features']) - len(features)
    counts = dict(source_counts)
    counts.update(retained_observable_features=0, preserved_observable_qualifiers=0,
                  usable_model_input_cds=0, source_aligned_exported_cds=0,
                  expected_model_abstention_cds=0, justified_model_abstention_cds=0,
                  withheld_annotation_cds=0, unsupported_length_cds=0,
                  uncertainty_required_features=0, uncertainty_preserved_features=0,
                  invented_feature_count=invented, unsupported_feature_label_count=0,
                  unsupported_protein_export_count=0)
    recognized_exports, reasons, unit_states, extra_failures = set(), Counter(), [], []
    for original in originals:
        fid = cid + '-f-' + original['suffix']
        feature = features.get(fid, {})
        linked = [u for u in document['units'] if u.get('feature_ids') == [fid]]
        unit = linked[0] if len(linked) == 1 else {}
        uid = unit.get('unit_id')
        emitted = extracted.get(uid)
        saved = audit.get('features', {}).get(fid, {})
        if original['observable']:
            counts['retained_observable_features'] += bool(feature)
            counts['preserved_observable_qualifiers'] += saved.get('qualifiers') == original['original'].qualifiers
        expected_label = original['type'] if original['type'] in FEATURE_TYPES else 'other'
        if feature and feature.get('feature_label') != expected_label:
            counts['unsupported_feature_label_count'] += 1
        if original['needs_uncertainty'] and original['observable']:
            counts['uncertainty_required_features'] += 1
            counts['uncertainty_preserved_features'] += feature.get('annotation_status') in ('unknown', 'ambiguous')
        if original['type'] != 'CDS':
            if emitted is not None:
                recognized_exports.add(uid)
                counts['unsupported_protein_export_count'] += 1
            unit_states.append({'feature_id': fid, 'unit_id': uid,
                                'source_annotation_status': original['status'],
                                'model_input_status': 'metadata_only',
                                'reason': 'Noncoding/source metadata has no protein-model input',
                                'abstention_justified': None})
            continue
        expected = original['expected']
        valid_export = emitted is not None and expected is not None and emitted == expected
        if emitted is not None:
            recognized_exports.add(uid)
            counts['unsupported_protein_export_count'] += not valid_export
        counts['source_aligned_exported_cds'] += valid_export
        over_cap = expected is not None and len(expected) > max_residues
        should_abstain = expected is None or over_cap
        counts['expected_model_abstention_cds'] += should_abstain
        counts['withheld_annotation_cds'] += expected is None and emitted is None
        counts['unsupported_length_cds'] += over_cap
        usable = valid_export and not over_cap and feature.get('annotation_status') == 'known'
        counts['usable_model_input_cds'] += usable
        justified = ((expected is None and emitted is None and unit.get('model_input_status') == 'needs_review'
                      and feature.get('annotation_status') in ('unknown', 'ambiguous')) or
                     (over_cap and unit.get('model_input_status') != 'ready'))
        counts['justified_model_abstention_cds'] += should_abstain and justified
        reason = ('outside_declared_residue_cap' if over_cap else original['reason'] if expected is None
                  else 'source_aligned_candidate' if usable else 'valid_source_translation_not_exported_correctly')
        if should_abstain:
            reasons[reason] += 1
        elif not usable:
            extra_failures.append('valid_source_translation_lost')
        unit_states.append({'feature_id': fid, 'unit_id': uid,
                            'source_annotation_status': original['status'],
                            'model_input_status': 'supported_candidate' if usable else 'abstain',
                            'reason': reason, 'abstention_justified': bool(justified) if should_abstain else None})
    counts['unsupported_protein_export_count'] += len(set(extracted) - recognized_exports)
    counts['erroneous_invention_count'] = (counts['invented_feature_count'] + counts['unsupported_feature_label_count']
                                          + counts['unsupported_protein_export_count'])
    if counts['erroneous_invention_count']:
        extra_failures.append('unsupported_annotation_or_protein_invented')
    if counts['uncertainty_preserved_features'] != counts['uncertainty_required_features']:
        extra_failures.append('source_uncertainty_not_preserved')
    if counts['justified_model_abstention_cds'] != counts['expected_model_abstention_cds']:
        extra_failures.append('unjustified_or_missing_abstention')
    source_binding = None if source_sha256 is None else audit.get('source_sha256') == source_sha256
    if source_binding is False:
        extra_failures.append('source_checksum_binding_failed')
    failed = sorted(set(preservation.get('failed_checks', []) + extra_failures))
    observation = ('no_observable_annotations' if not counts['observable_features'] else
                   'no_observable_cds' if not counts['source_cds'] else 'observable_source_annotations')
    record_abstention = counts['usable_model_input_cds'] == 0
    return {'construct_id': cid, 'primary_archetype': document['construct']['primary_archetype'],
            'status': 'failed' if failed else observation if observation == 'no_observable_annotations' else 'audited',
            'quality_stratum': classify_quality(source_counts), 'annotation_observation_status': observation,
            'source_checksum_binding': source_binding, 'counts': counts, 'metrics': quality_metrics(counts),
            'record_has_usable_model_input': not record_abstention,
            'record_abstention_justified': (not failed and counts['expected_model_abstention_cds'] == counts['source_cds']) if record_abstention else None,
            'withheld_reasons': dict(sorted(reasons.items())), 'unit_assessments': unit_states,
            'source_preservation_audit_status': preservation['status'],
            'source_preservation_audit': preservation,
            'source_preservation_scope': 'Original features including bookkeeping; observable E1 rates exclude type source',
            'failed_checks': failed}


def quality_metrics(counts):
    return {
        'observable_feature_preservation': _fraction(counts['retained_observable_features'], counts['observable_features']),
        'observable_qualifier_preservation': _fraction(counts['preserved_observable_qualifiers'], counts['observable_features']),
        'confirmed_source_aligned_translation_rate': _fraction(counts['source_aligned_exported_cds'], counts['source_cds']),
        'useful_model_input_rate': _fraction(counts['usable_model_input_cds'], counts['source_cds']),
        'justified_model_abstention_rate': _fraction(counts['justified_model_abstention_cds'], counts['expected_model_abstention_cds']),
        'source_uncertainty_preservation': _fraction(counts['uncertainty_preserved_features'], counts['uncertainty_required_features']),
        'missing_translation_rate': _fraction(counts['missing_translation_cds'], counts['source_cds']),
        'conflicting_qualifier_cds_rate': _fraction(counts['conflicting_qualifier_cds'], counts['source_cds']),
        'source_ambiguous_feature_rate': _fraction(counts['source_ambiguous_features'], counts['observable_features']),
        'source_unknown_feature_rate': _fraction(counts['source_unknown_features'], counts['observable_features']),
    }


def summarize_quality(records):
    """Aggregate counts, keeping every numerator and eligible denominator visible."""
    counts = dict(Counter({key: sum(row['counts'][key] for row in records)
                           for key in records[0]['counts']})) if records else {}
    return {'records': len(records), 'quality_strata': dict(Counter(r['quality_stratum'] for r in records)),
            'counts': counts, 'metrics': quality_metrics(counts) if counts else None,
            'record_usable_model_input_rate': _fraction(sum(r['record_has_usable_model_input'] for r in records), len(records)),
            'record_justified_abstention_rate': _fraction(sum(r['record_abstention_justified'] is True for r in records),
                sum(r['record_abstention_justified'] is not None for r in records)),
            'withheld_reasons': dict(sum((Counter(r['withheld_reasons']) for r in records), Counter()))}


def run_annotation_quality(inputs_path, output_path, *, max_residues=1022):
    """Run E1 only after the existing T08 admission and source-audit pipeline.

    No inputs/pending gates yield no quality measurements. The caller's shared
    v1.1 runner owns protocol freezing; this function records its admitted sources.
    """
    if type(max_residues) is not int or not 1 <= max_residues <= 1022:
        raise InputError('E1 max_residues must be in 1..1022')
    output = local_output(Path(output_path))
    if output.exists() and any(output.iterdir()):
        raise InputError('E1 output directory must be new or empty')
    output.mkdir(parents=True, exist_ok=True)
    report = run_metadata('E1_annotation_quality', {'source_contract': 't05-inputs-1', 'max_residues': max_residues})
    report.update(status='blocked', inference_executed=False, scientific_scores_produced=False,
                  quality_records=[], failed_records=[])
    with Measure() as measured:
        try:
            evidence = read_json(Path(inputs_path)).get('evidence_kind', 'unclassified')
            if evidence not in ('admitted_research', 'fictional_software_fixture', 'unclassified'):
                raise InputError('Unknown explicit evidence_kind')
            report['evidence_kind'] = evidence
            prepared = prepare_pilot(Path(inputs_path), output / 'prepared')
            report['preparation_status'] = prepared['status']
            report['input_config_sha256'] = prepared.get('input_config_sha256')
            report['provenance'] = prepared.get('provenance')
            report['failed_records'] = [r for r in prepared.get('records', []) if r.get('status') == 'failed']
            if prepared['status'] not in ('ok', 'partial') or not prepared.get('parsed_records'):
                report.update(status=prepared['status'], reason=prepared.get('reason', 'No admitted, parsed source annotations'))
            else:
                _, sources, _ = reviewed_sources(Path(inputs_path))
                docs = json.loads((output / 'prepared/normalized.json').read_text(encoding='utf-8'))
                audits = read_json(output / 'prepared/annotation-audit.json')
                from Bio import SeqIO
                with (output / 'prepared/units.fasta').open(encoding='utf-8') as handle:
                    exported = {r.id: str(r.seq) for r in SeqIO.parse(handle, 'fasta')}
                by_parent = {d['construct']['construct_id']: d for d in docs}
                for source in sources:
                    cid = source.review['record_id']
                    if cid not in by_parent:
                        continue
                    checksum = digest_file(source.path)
                    if checksum.lower() != source.review['checksum_sha256'].lower():
                        raise InputError('E1 source changed after admission')
                    record = read_source(source)
                    if digest_file(source.path) != checksum:
                        raise InputError('E1 source changed during re-audit')
                    doc = by_parent[cid]
                    extracted = {u['unit_id']: exported[source.input_id + '/' + u['unit_id']]
                                 for u in doc['units'] if source.input_id + '/' + u['unit_id'] in exported}
                    report['quality_records'].append(audit_quality(record, doc, audits[cid], extracted,
                        source_sha256=checksum, max_residues=max_residues))
                rows = report['quality_records']
                successes = [r for r in rows if r['status'] != 'failed']
                report['audit_failures'] = [r for r in rows if r['status'] == 'failed']
                report['overall'] = summarize_quality(successes)
                report['by_archetype'] = {a: summarize_quality([r for r in successes if r['primary_archetype'] == a])
                                          for a in sorted({r['primary_archetype'] for r in rows})}
                report['by_quality_stratum'] = {a: summarize_quality([r for r in successes if r['quality_stratum'] == a])
                                                for a in sorted({r['quality_stratum'] for r in rows})}
                report['failures_by_archetype'] = dict(Counter(r['primary_archetype']
                    for r in report['failed_records'] + report['audit_failures']))
                report['parser_failure_summary'] = {'attempted_records': len(rows) + len(report['failed_records']),
                    'failed_records': len(report['failed_records']) + len(report['audit_failures']),
                    'erroneous_invention_count': sum(r['counts']['erroneous_invention_count'] for r in rows)}
                failures = report['failed_records'] or report['audit_failures']
                observable = any(r['counts']['observable_features'] for r in successes)
                report.update(status='partial' if failures else 'ok' if observable else 'no_observable_annotations',
                              annotation_measurements_produced=observable and bool(successes),
                              scientific_scores_produced=observable and bool(successes) and evidence == 'admitted_research')
                write_json(output / 'quality-audit.json', {'records': rows, 'overall': report['overall'],
                    'by_archetype': report['by_archetype'], 'by_quality_stratum': report['by_quality_stratum']})
        except Exception as exc:
            report.update(status='blocked', scientific_scores_produced=False, **failure(exc))
    report['resources_total'] = measured.data
    report['limitations'].extend([
        'Source-only bookkeeping features are excluded from observable preservation denominators.',
        'Quality strata describe reviewed annotation completeness, not biological correctness or model retrieval.',
        'Usable means source-aligned canonical translation within the declared cap, not independently validated function.',
        'Long annotated exports remain local metadata; E1 marks them incompatible and runs no pretrained model.',
    ])
    write_json(output / 'report.json', report)
    return report
