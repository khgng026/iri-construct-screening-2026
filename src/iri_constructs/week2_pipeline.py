"""Preparatory ingestion and deterministic source QA; never scientific admission.

Only source-annotated translations are exported. Sequence-bearing source copies,
original annotations and protein exports stay in ignored output directories.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re

from .annotations import json_safe, normalize_record
from .bench_inputs import detect_record_format, parse_sequence_record, validate_sequence
from .runtime import InputError, digest_file, failure, local_output, read_json, require_text, write_json
from .schema_check import ARCHETYPES, check_schema_example


def canonical_sha256(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def _checked_path(entry: dict, key: str, checksum_key: str, base: Path) -> Path:
    path = (base / require_text(entry, key)).resolve()
    checksum = require_text(entry, checksum_key)
    if not re.fullmatch('[0-9a-fA-F]{64}', checksum):
        raise InputError('Expected an exact SHA-256 checksum in source manifest')
    try:
        actual = digest_file(path)
    except OSError as exc:
        raise InputError('Source file is missing or unreadable') from exc
    if checksum.lower() != actual:
        raise InputError('Source checksum does not match preparation manifest')
    return path


def _igem_record(entry: dict, path: Path, base: Path):
    """Conservative adapter for published Registry part/feature API responses.

    No location convention is guessed. Unknown coordinate semantics leave null
    normalized locations and the exact original locations in the audit sidecar.
    No protein is translated from DNA to compensate for missing /translation.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    part = read_json(path)
    annotations_path = _checked_path(entry, 'annotations_path', 'annotations_checksum_sha256', base)
    annotations = read_json(annotations_path)
    sequence_id = require_text(entry, 'sequence_id')
    if part.get('name') != sequence_id:
        raise InputError('Registry identity does not match the preparation manifest')
    if part.get('status') != 'published':
        raise InputError('Registry part is not published; preparation scope withheld')
    sequence = part.get('sequence')
    if isinstance(sequence, dict):
        sequence = sequence.get('sequence') or sequence.get('elements')
    if not isinstance(sequence, str):
        raise InputError('Registry response has no complete sequence string')
    sequence = sequence.upper()
    validate_sequence(sequence, 'dna')
    features = annotations.get('data')
    total = annotations.get('total')
    if (not isinstance(features, list) or type(total) is not int or total != len(features)
            or type(annotations.get('page')) is not int or annotations['page'] != 1):
        raise InputError('Registry feature response is incomplete or has unsupported pagination')
    record = SeqRecord(Seq(sequence), id=sequence_id)
    topology = part.get('topology')
    topology = topology.get('label', '') if isinstance(topology, dict) else ''
    record.annotations = {'molecule_type': 'DNA', 'topology': topology.lower(),
                          'registry_sequence_source': part.get('source'),
                          'registry_coordinate_convention': entry.get('coordinate_convention', 'unverified')}
    one_based = entry.get('coordinate_convention') == 'one_based_inclusive'
    if one_based and not entry.get('coordinate_evidence'):
        raise InputError('Registry coordinate normalization requires source evidence')
    for item in features:
        if not isinstance(item, dict):
            raise InputError('Registry feature must be an object')
        role = item.get('role', {})
        role = role if isinstance(role, dict) else {}
        feature_type = role.get('label') or 'unknown'
        qualifiers = {'registry_original_feature': [json.dumps(item, sort_keys=True)],
                      'registry_coordinate_convention': [entry.get('coordinate_convention', 'unverified')]}
        if isinstance(item.get('label'), str):
            qualifiers['label'] = [item['label']]
        location = None
        if one_based:
            locations = item.get('locations')
            strand = {'forward': 1, 'reverse': -1}.get(item.get('strand'))
            # Joined Registry annotations need a documented ordering convention;
            # retain raw evidence instead of inventing multipart semantics.
            if isinstance(locations, list) and len(locations) == 1:
                start, end = locations[0].get('start'), locations[0].get('end')
                if (type(start) is not int or type(end) is not int
                        or not 1 <= start <= end <= len(sequence)):
                    raise InputError('Registry feature coordinates are outside the parent sequence')
                location = SimpleLocation(start - 1, end, strand=strand)
        record.features.append(SeqFeature(location, type=feature_type,
                                          id=str(item.get('uuid') or '<unknown id>'), qualifiers=qualifiers))
    if digest_file(annotations_path).lower() != entry['annotations_checksum_sha256'].lower():
        raise InputError('Registry annotation checksum changed during parsing')
    return record


def load_preparatory_record(entry: dict, base: Path):
    path = _checked_path(entry, 'path', 'checksum_sha256', base)
    fmt = entry.get('format', 'auto')
    molecule = entry.get('molecule_type', 'dna')
    if molecule != 'dna':
        raise InputError('Week 2 construct QA requires DNA; protein references use their separate review pipeline')
    if fmt == 'igem-json':
        record = _igem_record(entry, path, base)
    else:
        detected = detect_record_format(path)
        if fmt == 'auto':
            fmt = detected
        record = parse_sequence_record(path, fmt, require_text(entry, 'sequence_id'),
                                       'dna', entry['checksum_sha256'])
        if fmt == 'genbank':
            # Exact sequence_id is always mandatory; accession metadata, when
            # supplied, must independently agree and may not silently override it.
            accession, version = entry.get('accession'), entry.get('accession_version')
            if accession and version:
                expected = accession if accession.endswith('.' + str(version)) else accession + '.' + str(version)
                if expected != record.id:
                    raise InputError('GenBank accession/version metadata differs from source record')
    if digest_file(path).lower() != entry['checksum_sha256'].lower():
        raise InputError('Source checksum changed during preparation parsing')
    for key, actual in [('sequence_length_nt', len(record)), ('feature_count', len(record.features))]:
        if key in entry and entry[key] is not None and entry[key] != actual:
            raise InputError('Source length or feature-count metadata does not match original record')
    return record, fmt


def _classification(exc: Exception) -> str:
    text = str(exc) if isinstance(exc, InputError) else ''
    for needle, category in [('Duplicate', 'duplicate_identity'), ('duplicate', 'duplicate_identity'),
                             ('checksum', 'source_checksum_mismatch'), ('format', 'unsupported_format'),
                             ('pagination', 'incomplete_annotations'), ('identity', 'identity_mismatch'),
                             ('sequence_id', 'identity_mismatch'), ('accession/version', 'identity_mismatch'),
                             ('coordinates', 'invalid_coordinates'), ('alphabet', 'invalid_sequence')]:
        if needle in text:
            return category
    return 'invalid_or_unsupported_source' if isinstance(exc, InputError) else 'processing_failure'


def _strict_checks(document: dict, audit: dict, extracted: dict, record) -> list[str]:
    errors = []
    features = {f['feature_id']: f for f in document['features']}
    units = {u['unit_id']: u for u in document['units']}
    if len(features) != len(record.features) or set(features) != set(audit['features']):
        errors.append('Feature count or audit linkage changed')
    for unit in units.values():
        linked = unit['feature_ids']
        if len(linked) != 1 or linked[0] not in features:
            errors.append('Unit does not have one traceable originating feature')
            continue
        if unit['source_segments'] != features[linked[0]]['segments']:
            errors.append('Unit and feature coordinate segments differ')
    for uid, protein in extracted.items():
        unit = units.get(uid)
        if unit is None or unit['molecule_type'] != 'protein':
            errors.append('Extracted protein has no protein-unit provenance')
            continue
        original = audit['features'][unit['feature_ids'][0]]
        translation = original['qualifiers'].get('translation', [])
        if original['type'] != 'CDS' or len(translation) != 1:
            errors.append('Protein is not an original annotated CDS translation')
        elif ''.join(translation[0].split()).upper().rstrip('*') != protein:
            errors.append('Protein differs from the original annotated translation')
    return errors


def run_preparatory_qa(manifest_path: Path, output: Path) -> dict:
    """Prepare local normalized objects and a source-free summary for review."""
    manifest_path = Path(manifest_path).resolve()
    manifest = read_json(manifest_path)
    if manifest.get('schema_version') != 'week2-source-inputs-1':
        raise InputError('Expected schema_version week2-source-inputs-1')
    entries = manifest.get('records')
    if not isinstance(entries, list) or not entries:
        raise InputError('Preparation manifest requires a nonempty records list')
    ids = [entry.get('source_id') if isinstance(entry, dict) else None for entry in entries]
    if any(not isinstance(i, str) or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.:-]*', i) for i in ids):
        raise InputError('Preparation source_id must be a stable identifier without path separators')
    if len(set(ids)) != len(ids):
        raise InputError('Duplicate preparation source_id')
    output = local_output(output)
    output.mkdir(parents=True, exist_ok=False)
    documents, audits, exports, unit_manifest, feature_manifest, rows = [], {}, {}, [], [], []
    parent_sequences, protein_sequences, lineages = defaultdict(list), defaultdict(list), defaultdict(list)
    aliases = set()
    for entry in entries:
        source_id = entry['source_id']
        row = {'source_id': source_id, 'status': 'started', 'scientific_admission': False,
               'scored_experiment_executed': False}
        try:
            proposed_lineages = entry.get('related_lineage_proposals', [])
            if (not isinstance(proposed_lineages, list)
                    or any(not isinstance(g, str) or not g.strip() for g in proposed_lineages)):
                raise InputError('Proposed lineage groups require a list of nonempty stable labels')
            if entry.get('lineage_group'):
                if not isinstance(entry['lineage_group'], str):
                    raise InputError('Legacy lineage group requires a stable text label')
                proposed_lineages = proposed_lineages + [entry['lineage_group']]
            proposed_lineages = sorted(set(proposed_lineages))
            alias = (entry.get('checksum_sha256', '').lower(), entry.get('sequence_id'))
            if alias in aliases:
                raise InputError('Duplicate checksummed source-record alias')
            record, fmt = load_preparatory_record(entry, manifest_path.parent)
            archetype = entry.get('archetype') or entry.get('archetype_proposed') or entry.get('proposed_archetype')
            if archetype not in ARCHETYPES:
                raise InputError('Source requires an explicit supported archetype proposal')
            review = {'record_id': source_id, 'source': entry.get('source', 'public-source-candidate'),
                      'source_url': entry.get('source_url', ''), 'accession': entry.get('accession', record.id),
                      'accession_version': str(entry.get('accession_version', 'source-snapshot')),
                      'archetype': archetype,
                      'notes': 'Preparatory source QA only; no source/scoring admission. See original source manifest.'}
            doc, audit, extracted = normalize_record(record, review, entry['checksum_sha256'].lower())
            audit = json_safe(audit)
            errors = _strict_checks(doc, audit, extracted, record)
            # Repeat from the same record to verify content determinism, including
            # exact source qualifiers, source-feature IDs and segment order.
            doc2, audit2, extracted2 = normalize_record(record, review, entry['checksum_sha256'].lower())
            if (doc, audit, extracted) != (doc2, json_safe(audit2), extracted2):
                errors.append('Repeated normalization changed content')
            path = output / 'normalized' / f'{source_id}.json'
            write_json(path, doc)
            errors.extend(check_schema_example(path))
            if errors:
                raise InputError('Normalized representation invariant failed: ' + '; '.join(errors))
            aliases.add(alias)
            documents.append(doc)
            audits[source_id] = audit
            for feature in doc['features']:
                original = audit['features'][feature['feature_id']]
                feature_manifest.append({'feature_id': feature['feature_id'], 'source_id': source_id,
                                         'source_feature_id': original['source_feature_id'],
                                         'source_feature_type': original['type'],
                                         'feature_label': feature['feature_label'],
                                         'location': feature['location'], 'segments': feature['segments'],
                                         'annotation_status': feature['annotation_status'],
                                         'source_sha256': entry['checksum_sha256'].lower(),
                                         'original_annotation_sha256': canonical_sha256(original)})
            sequence_hash = hashlib.sha256(str(record.seq).upper().encode('ascii')).hexdigest()
            parent_sequences[sequence_hash].append(source_id)
            for group in proposed_lineages:
                lineages[group].append(source_id)
            for uid, protein in sorted(extracted.items()):
                iid = source_id + '/' + uid
                exports[iid] = protein
                unit = next(u for u in doc['units'] if u['unit_id'] == uid)
                digest = hashlib.sha256(protein.encode('ascii')).hexdigest()
                protein_sequences[digest].append(iid)
                unit_manifest.append({'input_id': iid, 'source_id': source_id, 'unit_id': uid,
                                      'feature_ids': unit['feature_ids'], 'source_segments': unit['source_segments'],
                                      'source_sha256': entry['checksum_sha256'].lower(),
                                      'sequence_sha256': digest, 'length_aa': len(protein),
                                      'selected_interface_eligibility': ('eligible_sequence_preflight'
                                                                          if len(protein) <= 1022 else 'blocked_length_cap'),
                                      'scientific_admission': False})
            row.update(status='prepared_for_review', format=fmt, archetype_proposed=archetype,
                       related_lineage_proposals=proposed_lineages,
                       source_record_id=record.id, source_type=entry.get('source_type', 'unspecified'),
                       source_sha256=entry['checksum_sha256'].lower(), sequence_sha256=sequence_hash,
                       source_feature_count=len(record.features), source_cds_count=sum(f.type == 'CDS' for f in record.features),
                       validated_annotated_protein_count=len(extracted),
                       annotation_states=dict(Counter(f['annotation_status'] for f in doc['features'])),
                       unsupported_unit_reasons=dict(Counter(u['representation_note'] for u in doc['units'] if u['unit_id'] not in extracted)),
                       source_content_sha256=canonical_sha256({'normalized': doc, 'audit': audit,
                                                              'protein_sha256': {k: hashlib.sha256(v.encode()).hexdigest() for k, v in extracted.items()}}),
                       deterministic=True, schema_errors=[],
                       coordinate_convention=entry.get('coordinate_convention', 'GenBank: 1-based inclusive source; normalized 0-based half-open' if fmt == 'genbank' else 'not_applicable_or_unverified'))
        except Exception as exc:
            row.update(status='blocked', failure_class=_classification(exc), **failure(exc))
        rows.append(row)
    write_json(output / 'normalized.json', documents)
    write_json(output / 'annotation-audit.json', audits)
    write_json(output / 'unit-provenance.json', {'schema_version': 'week2-unit-provenance-1',
                                              'evidence_kind': 'preparatory_source_QA',
                                              'model_interface': {'molecule_type': 'protein', 'alphabet': 'canonical_amino_acids',
                                                                  'max_residues': 1022,
                                                                  'status': 'sequence_preflight_only_not_inference_or_scientific_admission'},
                                              'units': unit_manifest})
    write_json(output / 'feature-provenance.json', {'schema_version': 'week2-feature-provenance-1',
                                                 'evidence_kind': 'preparatory_source_QA', 'features': feature_manifest})
    with (output / 'units.fasta').open('w', encoding='utf-8') as handle:
        for iid, protein in sorted(exports.items()):
            handle.write('>' + iid + '\n' + protein + '\n')
    prepared = [r for r in rows if r['status'] == 'prepared_for_review']
    summary = {'schema_version': 'week2-source-QA-1', 'evidence_kind': 'preparatory_source_QA',
               'status': 'prepared_for_review' if len(prepared) == len(rows) else 'partial_preparation',
               'scientific_admission': False, 'scientific_scores_produced': False,
               'inference_executed': False, 'input_manifest_sha256': digest_file(manifest_path),
               'requested_records': len(rows), 'prepared_records': len(prepared),
               'blocked_records': len(rows) - len(prepared),
               'coverage_by_archetype_proposed': dict(Counter(r['archetype_proposed'] for r in prepared)),
               'exported_annotated_proteins': len(unit_manifest), 'records': rows,
               'protein_interface_coverage': dict(Counter(u['selected_interface_eligibility'] for u in unit_manifest)),
               'duplicate_parent_sequence_groups': [{'sequence_sha256': h, 'source_ids': sorted(v)}
                                                     for h, v in sorted(parent_sequences.items()) if len(v) > 1],
               'duplicate_protein_sequence_groups': [{'sequence_sha256': h, 'input_ids': sorted(v)}
                                                      for h, v in sorted(protein_sequences.items()) if len(v) > 1],
               'declared_lineage_groups': dict(sorted(lineages.items())),
               'schema': {'version': '0.1', 'status': 'technical_validation_only_pending_scientific_admission',
                          'source_schema_sha256': digest_file(Path(__file__).parents[2] / 'docs/week1/T03_archetypes_and_schema.md')},
               'limitations': ['Source agreement and deterministic preparation are not biological validation or experimental E1 scores.',
                               'Candidate archetypes and lineage groups remain proposals; no family/reference panel is approved here.',
                               'Missing source translations remain withheld; no ORFs or missing sequences are inferred.',
                               'Original annotations and protein sequences remain in ignored local artifacts.']}
    summary['content_sha256'] = canonical_sha256(summary)
    write_json(output / 'report.json', summary)
    return summary
