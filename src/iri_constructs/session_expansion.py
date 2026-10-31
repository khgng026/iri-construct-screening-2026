"""Read-only, independent admission validation for frozen expanded catalogues.

No downloads, model imports or scoring. Returned sequences are limited to the
union of achieved, frozen endpoints; an acquired overshoot is provenance only.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import re
from urllib.parse import parse_qs, urlparse

from .session_store import canonical, digest, file_hash

CANONICAL = frozenset('ACDEFGHIKLMNPQRSTVWY')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def checked_path(root, relative, sha=None):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    require(path.is_relative_to(root), 'Evidence path escapes its owning root')
    require(path.is_file(), 'Missing source evidence: ' + str(path))
    if sha is not None:
        require(file_hash(path) == sha, 'Evidence checksum mismatch: ' + str(path))
    return path


def source_families(text):
    return sorted(set(re.findall(r'Belongs to the [^.{}]*? family(?=\.|\s*\{|$)', text)))


def _labels(value):
    result = {}
    def visit(item):
        if isinstance(item, dict):
            if 'name' in item and 'label' in item:
                result[item['name']] = item['label']
            for child in item.values(): visit(child)
        elif isinstance(item, list):
            for child in item: visit(child)
    visit(value)
    return result


def _eligible(metadata):
    require(metadata['reviewed'].lower() == 'reviewed' and metadata['organism_id'] == '559292',
            'Admission requires reviewed S288c source identity')
    require('-' not in metadata['accession'] and not metadata['fragment'].strip(),
            'Isoform or fragment cannot enter this catalogue')
    require(metadata['version'].isdigit() and metadata['sequence_version'].isdigit(), 'Missing source versions')
    require(metadata['length'].isdigit() and 1 <= int(metadata['length']) <= 1022, 'Invalid source length')
    require(metadata['protein_existence'].lower() != 'uncertain', 'Uncertain protein existence')


def _response(session, relative, fields, labels, expected_sha=None):
    path = checked_path(session, relative, expected_sha)
    provenance = read(path.with_suffix(path.suffix + '.provenance.json'))
    require(file_hash(path) == provenance['sha256'] and path.stat().st_size == provenance['bytes'],
            'Response provenance does not match saved bytes')
    url = urlparse(provenance['url'])
    require(url.scheme == 'https' and url.hostname == 'rest.uniprot.org' and provenance['status'] == 200,
            'Response is not a successful official public source')
    params = parse_qs(url.query)
    require(params.get('fields') == [','.join(fields)] and params.get('format') == ['tsv'],
            'Response requested fields differ from frozen policy')
    require(provenance['retrieved_utc'], 'Missing retrieval time')
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(reader.fieldnames == [labels[f] for f in fields], 'Response headers differ from source field configuration')
        rows = []
        for row in reader:
            require(None not in row and set(row) == set(reader.fieldnames), 'Malformed source TSV row')
            rows.append({f: row[labels[f]] for f in fields})
    require(len({r['accession'] for r in rows}) == len(rows), 'Duplicate response accession')
    return provenance, rows


def _fasta(path):
    records = {}; accession = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            accession = line[1:]
            require(accession and accession not in records and not any(c.isspace() for c in accession), 'Malformed FASTA identifier')
            records[accession] = ''
        else:
            require(accession is not None and line and set(line) <= CANONICAL, 'Noncanonical FASTA sequence')
            records[accession] += line
    require(records and all(records.values()), 'Empty FASTA record')
    return records


def load_expansion(repo_root, session, manifest_path):
    """Validate original source bytes, reviews, aliases and deterministic prefixes.

    Returns the agreed API plus a small ``validation`` evidence summary. This is
    a source consistency check; recorded role reviews are not institutional or
    independent experimental confirmation of function.
    """
    from .t05_composition import load_reviewed
    repo_root, session = Path(repo_root).resolve(), Path(session).resolve()
    manifest_path = Path(manifest_path).resolve()
    require(manifest_path.is_relative_to(session), 'Expanded manifest must belong to this session')
    manifest = read(manifest_path)
    require(manifest['version'] == 'iri-expanded-admitted-v1' and manifest['status'] == 'frozen_before_expanded_scoring',
            'Unsupported or unfrozen expansion')
    design_path = checked_path(session, 'inputs/expanded_catalogue_design_v1.json', manifest['design_sha256'])
    design = read(design_path)
    require(design['version'] == 'iri-expanded-catalogue-v1' and design['status'] == 'frozen_before_catalogue_acquisition_or_scoring',
            'Expansion design is not frozen')
    checked_path(repo_root, 'scripts/session_expand_catalogue.py', manifest['acquisition_script_sha256'])
    index = read(checked_path(session, 'expansion/metadata_index.json', manifest['metadata_index_sha256']))
    milestone = manifest['current_milestone']
    require(milestone in design['nested_unique_sequence_targets'], 'Unplanned expansion milestone')
    review = read(checked_path(session, f'expansion/admission_review_{milestone}_v1.json', manifest['review_sha256']))
    progress = read(checked_path(session, manifest['sequence_acquisition_progress'], manifest['sequence_acquisition_progress_sha256']))
    require(index['design_sha256'] == review['design_sha256'] == progress['design_sha256'] == manifest['design_sha256'],
            'Source/design identity differs across stages')
    require(review['metadata_index_sha256'] == manifest['metadata_index_sha256'] and
            progress['review_sha256'] == manifest['review_sha256'], 'Review/acquisition identity mismatch')
    require(review['status'] == 'frozen_before_new_sequence_acquisition' and progress['requested_milestone'] == milestone,
            'Sequence acquisition was not bound to this frozen review')
    require(index['metadata_only'] and not index['sequences_requested'] and index['remaining_url'] is None and
            index['status'] == 'complete_metadata_pending_review', 'Metadata acquisition is incomplete or includes sequences')
    require('sequence' not in design['metadata_fields'], 'Metadata stage contains sequences')
    labels = _labels(read(checked_path(session, 'expansion/evidence/result_fields.json', design['official_field_configuration_sha256'])))
    license_path = checked_path(session, 'expansion/evidence/license.json', design['license_evidence_sha256'])
    require('Creative Commons Attribution 4.0' in license_path.read_text(), 'Licence evidence lacks expected reuse terms')
    fields = design['metadata_fields']; indexed = index['records']; source_cache = {}
    rebuilt = {}; releases = set()
    for page in index['pages']:
        relative = f"expansion/metadata/page-{page['page']:03d}.tsv"
        provenance, rows = _response(session, relative, fields, labels, page['response_sha256'])
        headers = {k.lower(): v for k, v in provenance['headers'].items()}
        require(int(headers['x-total-results']) == index['source_total_results'] and
                headers['x-uniprot-release'] == page['release'], 'Source pagination release/count mismatch')
        require(len(rows) == page['records'], 'Source page row count mismatch')
        releases.add(page['release']); source_cache[relative] = (provenance, rows)
        for row_number, metadata in enumerate(rows, 2):
            accession = metadata['accession']
            require(accession not in rebuilt, 'Repeated accession across metadata pages')
            expected = indexed[accession]
            require(expected['metadata'] == metadata and expected['metadata_row_sha256'] == digest(canonical(metadata)) and
                    expected['source_page'] == relative and expected['source_page_sha256'] == provenance['sha256'] and
                    expected['source_url'] == provenance['url'] and expected['source_retrieved_utc'] == provenance['retrieved_utc'] and
                    expected['source_row_number_1based'] == row_number, 'Indexed row differs from original response')
            rebuilt[accession] = metadata
    require(set(rebuilt) == set(indexed) and len(rebuilt) == index['record_count'] == index['source_total_results'],
            'Metadata index does not exactly cover source responses')
    require(sorted(releases) == index['source_releases'] == manifest['source_releases'] and len(releases) == 1,
            'Acquisition mixes source releases')
    decisions = {}
    for evidence in review['decision_files']:
        path = checked_path(session, evidence['path'], evidence['sha256'])
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            require(reader.fieldnames == ['accession', 'decision', 'ordinary_role_basis', 'review_note'], 'Changed individual review fields')
            for row in reader:
                accession = row['accession']
                require(accession not in decisions and row['ordinary_role_basis'] and row['review_note'], 'Missing or duplicate individual review')
                decisions[accession] = (row, evidence['sha256'])
    require(set(decisions) == set(review['records']), 'Review records differ from explicit decisions')
    for accession, decision in review['records'].items():
        row, decision_hash = decisions[accession]; source = indexed[accession]; metadata = source['metadata']
        require(row['decision'] in ('approved', 'pending', 'uncertain', 'rejected'), 'Unknown review decision')
        require(decision['review_status'] == decision['safety_review'] == row['decision'] and decision['licensing_status'] == 'approved',
                'Individual admission flags differ')
        require(decision['accession'] == accession and decision['rationale'] == row['ordinary_role_basis'] and
                decision['review_note'] == row['review_note'] and decision['manual_decision_source_sha256'] == decision_hash,
                'Review rationale/decision source mismatch')
        require(all(decision[k] == source[k] for k in ('metadata_row_sha256', 'source_page', 'source_page_sha256', 'source_url', 'source_retrieved_utc')),
                'Review is bound to different metadata')
        for basis in ('source_annotation_basis', 'source_identity_basis'):
            require(decision[basis] and all(metadata[k] == v for k, v in decision[basis].items()), 'Review annotation/identity basis mismatch')
        require(decision['license_evidence_sha256'] == design['license_evidence_sha256'] and
                decision['license_url'] == design['license_url'] and decision['reuse_scope'] == design['license'], 'Review licensing evidence mismatch')
        if row['decision'] == 'approved':
            _eligible(metadata)
            require(decision['annotation_status'] == 'known', 'Approved record lacks known annotation')
    counts = {status: sum(r['review_status'] == status for r in review['records'].values())
              for status in ('approved', 'pending', 'uncertain', 'rejected')}
    require(counts == review['counts'] == manifest['review_decision_counts'], 'Review denominator mismatch')
    require(len(indexed) - len(decisions) == manifest['unreviewed_metadata_records'] == review['unreviewed_metadata_records'],
            'Unreviewed denominator mismatch')
    for prior in review['prior_frozen_reviews']:
        previous = read(checked_path(session, prior['path'], prior['sha256']))
        require(previous['design_sha256'] == review['design_sha256'] and previous['metadata_index_sha256'] == review['metadata_index_sha256'],
                'Prior review design/source mismatch')
        stable_fields = ('review_status', 'safety_review', 'licensing_status', 'metadata_row_sha256',
                         'rationale', 'review_note', 'source_annotation_basis', 'source_identity_basis',
                         'license_evidence_sha256')
        require(all(all(review['records'][a][k] == r[k] for k in stable_fields)
                    for a, r in previous['records'].items()), 'Previously frozen review changed')
        previous_max = max((digest((design['selection_seed'] + '|' + a).encode()), a) for a in previous['records'])
        require(all((digest((design['selection_seed'] + '|' + a).encode()), a) > previous_max
                    for a in set(review['records']) - set(previous['records'])), 'New review precedes an earlier frozen prefix')
    base_path = checked_path(repo_root, 'data/manifests/t05_expanded_catalogue.json', manifest['source_catalogue_manifest_sha256'])
    require(file_hash(base_path) == design['source_catalogue_manifest_sha256'], 'Expansion base differs from design')
    base_manifest = read(base_path); base_sequences = load_reviewed(repo_root, base_manifest)
    require(len(base_sequences) == 105, 'Original 105-accession source changed')
    families = read(checked_path(session, 'inputs/family_labels_v1.json', manifest['original_family_labels_sha256']))
    challenge = read(checked_path(session, 'inputs/challenge_manifest.json', design['source_challenge_manifest_sha256']))
    targets = sorted({a for q in challenge['queries'] for a in q['components']})
    require(len(targets) == 16 and set(targets) <= set(base_sequences), 'Original target membership changed')
    by_hash = {}; aliases = {}; sequences = {}; expected_catalogue = {}
    for record in sorted(base_manifest['records'], key=lambda r: r['accession']):
        accession = record['accession']; sequence = base_sequences[accession]; sha = digest(sequence.encode())
        rep = by_hash.setdefault(sha, accession); aliases[accession] = rep
        if rep not in expected_catalogue:
            expected_catalogue[rep] = {'accession': rep, 'sequence_sha256': sha, 'length': len(sequence),
                'origin': 'original_reviewed_105', 'aliases': [], 'original_source': record,
                'explicit_family_statements': sorted({e['family_statement'] for e in families['source_evidence'][accession]['annotation_evidence']})}
            sequences[rep] = sequence
        expected_catalogue[rep]['aliases'].append(accession)
    base_order = sorted(expected_catalogue)
    require(len(base_order) == 101 and base_order == manifest['base_representatives'], 'Original unique-sequence representatives changed')
    approved = sorted((a for a, d in review['records'].items() if d['review_status'] == 'approved' and a not in base_sequences),
                      key=lambda a: (digest((design['selection_seed'] + '|' + a).encode()), a))[:design['maximum_new_sequence_records']]
    processed = progress['processed_approved_records']; accepted = progress['accepted']; rejected = progress['rejected']
    require(isinstance(processed, int) and 0 <= processed <= len(approved), 'Invalid acquisition prefix size')
    require(set(accepted).isdisjoint(rejected) and set(accepted) | set(rejected) == set(approved[:processed]),
            'Acquisition is not exactly the deterministic approved prefix')
    new_order = []; duplicates = {}; responses = {}
    for accession in approved[:processed]:
        if accession in rejected: continue
        record = accepted[accession]; relative = record['source_response']
        if relative not in responses:
            provenance, rows = _response(session, relative, fields + ['sequence'], labels, record['source_response_sha256'])
            require(all(r['accession'] in approved[:processed] for r in rows), 'Sequence response includes an unapproved/unprocessed accession')
            responses[relative] = (provenance, {r['accession']: r for r in rows})
        provenance, rows = responses[relative]; row = rows[accession]; sequence = row['sequence']
        metadata = {f: row[f] for f in fields}; sha = digest(sequence.encode())
        require(metadata == indexed[accession]['metadata'] == record['metadata'] and
                record['metadata_row_sha256'] == indexed[accession]['metadata_row_sha256'] and record['review'] == review['records'][accession],
                'Acquired source differs from the specifically approved metadata/review')
        require(sequence == record['sequence'] and sha == record['sequence_sha256'] and provenance['sha256'] == record['source_response_sha256'],
                'Acquired sequence or provenance changed')
        require(sequence and set(sequence) <= CANONICAL and len(sequence) == int(metadata['length']) <= 1022, 'Noncanonical source sequence/length')
        rep = by_hash.setdefault(sha, accession); aliases[accession] = rep
        if rep != accession:
            duplicates[accession] = rep; expected_catalogue[rep]['aliases'].append(accession); continue
        new_order.append(accession); sequences[accession] = sequence
        expected_catalogue[accession] = {'accession': accession, 'sequence_sha256': sha, 'length': len(sequence),
            'origin': 'new_reviewed_metadata', 'aliases': [accession], 'metadata': metadata, 'review': review['records'][accession],
            'metadata_row_sha256': record['metadata_row_sha256'], 'source_response': relative,
            'source_response_sha256': record['source_response_sha256'], 'source_provenance': provenance,
            'explicit_family_statements': source_families(metadata['cc_similarity'])}
    require(manifest['catalogue'] == expected_catalogue and manifest['aliases_to_representative'] == aliases and
            manifest['new_representatives_hash_order'] == new_order and manifest['new_exact_sequence_duplicates'] == duplicates,
            'Frozen catalogue/aliases differ from independently reconstructed sources')
    require(len(sequences) == manifest['total_unique_sequences'] == progress['unique_sequences_including_base'], 'Unique-group denominator mismatch')
    require(manifest['acquired_new_records'] == processed and manifest['reviewed_approved_new_records'] == len(approved) and
            manifest['unacquired_approved_records'] == len(approved)-processed and manifest['rejected_after_acquisition'] == rejected,
            'Acquisition denominator mismatch')
    fasta = _fasta(checked_path(session, manifest['canonical_fasta'], manifest['canonical_fasta_sha256']))
    require(fasta == sequences, 'Canonical FASTA differs from exact source strings')
    target_aliases = {a: aliases[a] for a in targets}
    require(manifest['target_aliases_to_representative'] == target_aliases, 'Target aliases changed')
    endpoints = design['nested_unique_sequence_targets']; full_order = base_order + new_order
    require(manifest['requested_endpoints'] == endpoints and manifest['later_endpoints_pending_review'] == [n for n in endpoints if n > milestone],
            'Unplanned endpoint labels')
    require(set(manifest['nested_catalogues']) == {str(n) for n in endpoints if n <= milestone}, 'Missing/extra frozen prefix')
    subsets = {'101': base_order}
    for size, panel in manifest['nested_catalogues'].items():
        n = int(size); included = full_order[:n]
        require(panel == {'requested_unique_sequences': n, 'achieved_unique_sequences': len(included),
                         'endpoint_reached': len(full_order) >= n, 'accessions': included}, 'Frozen prefix differs from predetermined selection')
        if panel['endpoint_reached']:
            require(set(target_aliases.values()) <= set(included), 'An exact target is absent from a reported endpoint')
            subsets[size] = sorted(included)
    for prior in manifest['earlier_frozen_milestones']:
        previous = read(checked_path(session, prior['path'], prior['sha256']))
        require(previous['design_sha256'] == manifest['design_sha256'] and previous['metadata_index_sha256'] == manifest['metadata_index_sha256'],
                'Earlier frozen milestone source changed')
        for size, panel in previous['nested_catalogues'].items():
            current = manifest['nested_catalogues'][size]['accessions']; before = panel['accessions']
            require(current[:len(before)] == before and (not panel['endpoint_reached'] or current == before), 'Earlier frozen prefix changed')
            require(all(manifest['catalogue'][a]['sequence_sha256'] == previous['catalogue'][a]['sequence_sha256'] for a in before),
                    'Earlier frozen sequence changed')
    target_families = {a: sorted({e['family_statement'] for e in families['source_evidence'][a]['annotation_evidence']}) for a in targets}
    composition = {size: {a: [r for r in panel['accessions'] if set(expected_catalogue[r]['explicit_family_statements']) & set(statements)]
                         for a, statements in target_families.items()} for size, panel in manifest['nested_catalogues'].items()}
    require(manifest['target_explicit_families'] == target_families and manifest['target_family_composition_by_subset'] == composition,
            'Descriptive source-family labels/composition changed')
    returned = set().union(*(set(panel) for panel in subsets.values()))
    return {'manifest': manifest, 'manifest_sha256': file_hash(manifest_path),
            'sequences': {a: sequences[a] for a in sorted(returned)}, 'subsets': subsets, 'target_aliases': target_aliases,
            'validation': {'metadata_records': len(indexed), 'individually_reviewed_records': len(decisions),
                'approved_new_records': len(approved), 'acquired_new_records': processed,
                'all_acquired_unique_groups': len(sequences), 'returned_unique_groups': len(returned),
                'acquired_outside_frozen_endpoints': sorted(set(sequences) - returned),
                'subset_counts': {k: len(v) for k, v in subsets.items()}, 'target_alias_count': len(target_aliases),
                'source_release': sorted(releases), 'no_scoring_or_function_validation': True}}
