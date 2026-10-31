#!/usr/bin/env python3
"""Inventory every paired family-retrieval outcome under a prescoring design."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash, atomic_json
from iri_constructs.session_benchmark import load_suite
from session_experiment import check_deadline, public_query
from session_catalogue_analysis import latest_completion, source_arms, validate_family_freeze, subset_denominators
from session_independent_audit import verify_unit, validate_ranking


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def open_source(path, suite):
    path = path.resolve()
    identity = read(path / 'identity.json')
    identity_hash = digest(canonical(identity))
    require(identity['mode'] in ('model', 'alignment_normalization'), 'Unexpected comparison source')
    require(identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'], 'Source challenge differs')
    require(identity['catalogue_accessions'] == sorted(suite['sequences']), 'Source catalogue differs')
    require(identity['catalogue_sequence_sha256'] == {a: digest(s.encode()) for a, s in suite['sequences'].items()},
            'Source catalogue sequence hashes differ')
    ids = identity['query_ids']
    require(len(ids) == len(set(ids)) and set(ids) == {q['query_id'] for q in suite['queries']}, 'Source coverage differs')
    expected = {'query-' + q for q in ids} | {'summary'}
    expected |= {'catalogue-' + a for a in suite['sequences']} if identity['mode'] == 'model' else {'reference-selfscores'}
    completion_path = latest_completion(path)
    completion = read(completion_path)
    require(completion['identity_hash'] == identity_hash and set(completion['units']) == expected
            and completion['expected_count'] == len(expected), 'Source completion coverage or identity differs')
    paths = {key: verify_unit(path, key, entry, identity_hash) for key, entry in completion['units'].items()}
    return {'path': path, 'identity': identity, 'units': paths,
            'evidence': {'path': str(path), 'identity_sha256': identity_hash,
                         'completion_sha256': file_hash(completion_path)}}


def surviving_ranking(rows, included, full_accessions):
    validate_ranking(rows, full_accessions)
    require(set(included) <= set(full_accessions), 'Mask includes unknown accession')
    return [{'accession': row['accession'], 'score': row['score'], 'rank': i + 1}
            for i, row in enumerate(row for row in rows if row['accession'] in included)]


def best_accepted(rows, accepted, k):
    """Absent labels are null, not misses. Input is a complete surviving ranking."""
    require(type(k) is int and k > 0, 'Invalid fixed hit budget')
    accessions = [r['accession'] for r in rows]
    require(len(accessions) == len(set(accessions)), 'Duplicate surviving accessions')
    require([r['rank'] for r in rows] == list(range(1, len(rows) + 1)), 'Surviving ranks not contiguous')
    surviving = sorted(set(accepted) & set(accessions))
    if not surviving:
        return {'scorable': False, 'accepted_survivors': [], 'best_accession': None,
                'rank': None, 'score': None, 'hit': None}
    best = next(row for row in rows if row['accession'] in surviving)
    return {'scorable': True, 'accepted_survivors': surviving, 'best_accession': best['accession'],
            'rank': best['rank'], 'score': best['score'], 'hit': best['rank'] <= k}


def outcome(left, right):
    require(left['scorable'] == right['scorable'] and left['accepted_survivors'] == right['accepted_survivors'],
            'Paired methods use different surviving labels')
    if not left['scorable']:
        return 'unscorable'
    return {(True, True): 'both_hit', (True, False): 'model_only',
            (False, True): 'alignment_only', (False, False): 'both_miss'}[(left['hit'], right['hit'])]


def summarize_rows(rows):
    counts = Counter(row['outcome'] for row in rows)
    scorable = len(rows) - counts['unscorable']
    values = {key: counts[key] for key in ('both_hit', 'model_only', 'alignment_only', 'both_miss', 'unscorable')}
    return {'original_components': len(rows), 'scorable_components': scorable, **values,
            'model_recall_at_5': (counts['both_hit'] + counts['model_only']) / scorable if scorable else None,
            'alignment_recall_at_5': (counts['both_hit'] + counts['alignment_only']) / scorable if scorable else None,
            'distinct_groups': len({row['group_id'] for row in rows}),
            'distinct_targets': len({row['target'] for row in rows})}


def case_selection(rows):
    return {direction: sorted((r for r in rows if r['outcome'] == direction),
                              key=lambda r: (r['query_id'], r['target']))[:5]
            for direction in ('model_only', 'alignment_only')}


def summarize_queries(rows):
    by_query = defaultdict(list)
    for row in rows:
        by_query[row['query_id']].append(row)
    groups = [group for group in by_query.values() if all(row['outcome'] != 'unscorable' for row in group)]
    counts = Counter(outcome({'scorable': True, 'accepted_survivors': [], 'hit': all(r['model']['hit'] for r in group)},
                             {'scorable': True, 'accepted_survivors': [], 'hit': all(r['alignment']['hit'] for r in group)})
                     for group in groups)
    return {'original_queries': len(by_query), 'all_components_scorable_queries': len(groups),
            'excluded_unscorable_queries': len(by_query) - len(groups),
            'all_components_outcome_counts': {key: counts[key] for key in ('both_hit', 'model_only', 'alignment_only', 'both_miss')},
            'model_all_components_recall_at_5': (counts['both_hit'] + counts['model_only']) / len(groups) if groups else None,
            'alignment_all_components_recall_at_5': (counts['both_hit'] + counts['alignment_only']) / len(groups) if groups else None}


def run(args):
    started = time.monotonic()
    session = args.session.resolve()
    check_deadline(session)
    design_path = session / 'inputs/family_discordance_design_v1.json'
    design = read(design_path)
    require(design['version'] == 'iri-family-discordance-v1' and design['k'] == 5 and design['arm'] == 'whole',
            'Unsupported frozen disagreement design')
    require(design['status'] == 'frozen_before_paired_discordance_counts', 'Design is not frozen')
    require(design['policies'] == ['exact_removed', 'near90_removed'] and design['reference_method'] == 'swnorm',
            'Frozen comparison policies or reference method differ')
    suite = load_suite(ROOT, session / 'inputs')
    labels_path = session / 'inputs/family_labels_v1.json'
    labels = read(labels_path)
    label_evidence = validate_family_freeze(labels, suite, labels_path, session / 'audit/independent_family_label_review.json')
    assignments = [value.split('=', 1) for value in args.model]
    require(all(len(value) == 2 for value in assignments) and len({name for name, _ in assignments}) == len(assignments),
            'Duplicate or invalid model source assignment')
    require({name for name, _ in assignments} == set(design['models']), 'Model roster differs from frozen design')
    sources = {name: open_source(Path(path), suite) for name, path in assignments}
    reference = open_source(args.reference_job, suite)
    require(reference['identity']['mode'] == 'alignment_normalization', 'Reference must be normalized alignment')
    normalization_path = session / 'audit/alignment_normalization_design_v1.json'
    normalization = read(normalization_path)
    require(normalization['version'] == 'alignment-reference-selfscore-normalization-v1'
            and normalization['status'] == 'frozen_before_normalized_scoring'
            and normalization['secondary_score'] == 'S(query, reference) / S(reference, reference)'
            and normalization['source_challenge_sha256'] == suite['verification']['challenge_manifest_sha256']
            and reference['identity']['design_sha256'] == file_hash(normalization_path)
            and reference['identity']['design'] == normalization, 'Reference normalization design differs')
    require(all(source['identity']['mode'] == 'model' for source in sources.values()), 'Whole-query sources must be models')
    model_ids = {'8m': 'facebook/esm2_t6_8M_UR50D', '650m': 'facebook/esm2_t33_650M_UR50D',
                 'esmc300m': 'biohub/ESMC-300M'}
    require(all(source['identity']['config']['model_id'] == model_ids[name] for name, source in sources.items()),
            'Roster name is assigned to the wrong model')
    queries = [q for q in suite['queries'] if q['component_count'] > 1]
    identity = {'mode': 'family_discordance_v1', 'design_sha256': file_hash(design_path),
                'challenge_sha256': suite['verification']['challenge_manifest_sha256'],
                'label_evidence': label_evidence, 'sources': {name: source['evidence'] for name, source in sources.items()},
                'reference': reference['evidence'], 'query_ids': [q['query_id'] for q in queries],
                'policies': design['policies'], 'k': 5,
                'code': {name: file_hash(ROOT / name) for name in ['scripts/session_family_discordance.py',
                         'scripts/session_catalogue_analysis.py', 'scripts/session_independent_audit.py',
                         'src/iri_constructs/session_store.py']},
                'interpretation': design['interpretation']}
    store = UnitStore(session / 'jobs' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    expected = []
    for query in queries:
        check_deadline(session)
        key = 'query-' + query['query_id']
        expected.append(key)
        records = {name: read(source['units'][key] / 'record.json') for name, source in sources.items()}
        raw_reference = read(reference['units'][key] / 'record.json')
        require(all(record['query'] == public_query(query) for record in [*records.values(), raw_reference]),
                'Source query metadata differs')
        reference_rows = source_arms(reference['identity'], raw_reference)['normalized_alignment'][0]
        rows = []
        for policy in design['policies']:
            removed = labels['removal_policies'][policy]['common_removed_accessions_by_group'][query['group_id']]
            included = set(suite['sequences']) - set(removed)
            reference_ranking = surviving_ranking(reference_rows, included, sorted(suite['sequences']))
            for name, record in records.items():
                model_ranking = surviving_ranking(record['rankings']['whole'], included, sorted(suite['sequences']))
                for index, span in enumerate(query['spans']):
                    target = span['accession']
                    model = best_accepted(model_ranking, labels['acceptable_matches'].get(target, []), 5)
                    alignment = best_accepted(reference_ranking, labels['acceptable_matches'].get(target, []), 5)
                    rows.append({'query_id': query['query_id'], 'group_id': query['group_id'], 'suite': query['suite'],
                                 'target': target, 'component_count': query['component_count'], 'component_index': index,
                                 'position': 'N' if index == 0 else ('C' if index == query['component_count'] - 1 else 'internal'),
                                 'fraction': span['fraction'], 'start0': span['start0'], 'end0': span['end0'],
                                 'query_length_aa': query['length_aa'], 'components_in_order': query['components'],
                                 'policy': policy, 'model_name': name, 'model': model, 'alignment': alignment,
                                 'outcome': outcome(model, alignment), 'catalogue_count': len(included),
                                 'model_top5': model_ranking[:5], 'alignment_top5': reference_ranking[:5]})
        record = {'query': public_query(query), 'observations': rows}
        old = store.read(key)
        if old is None:
            store.publish(key, record)
        else:
            require(old == record, 'Recomputed paired outcomes differ from saved record')
    summaries = {}
    for name in sorted(sources):
        summaries[name] = {}
        for policy in design['policies']:
            rows = [row for key in expected for row in store.read(key)['observations']
                    if row['model_name'] == name and row['policy'] == policy]
            counts = summarize_rows(rows)
            frozen = subset_denominators(labels, policy, {q['query_id'] for q in queries})['combined']
            require(counts['original_components'] == frozen['original_component_observations'] and
                    counts['scorable_components'] == frozen['scorable_component_observations'], 'Frozen component denominator drift')
            query_summary = summarize_queries(rows)
            require(query_summary['all_components_scorable_queries'] == frozen['all_components_scorable_queries'],
                    'Frozen full-query denominator drift')
            strata = {}
            for field in ('component_count', 'position', 'target'):
                strata[field] = {str(value): summarize_rows([row for row in rows if row[field] == value])
                                 for value in sorted({row[field] for row in rows})}
            summaries[name][policy] = {'components': counts, 'by_stratum': strata, 'queries': query_summary,
                'queries_by_component_count': {str(count): summarize_queries([row for row in rows if row['component_count'] == count])
                                              for count in sorted({row['component_count'] for row in rows})},
                'cases': case_selection(rows)}
    summary = {'summaries': summaries, 'design': design, 'all_outcomes_saved': True,
               'elapsed_seconds_excluded_from_semantic_summary': True,
               'complementarity_note': 'A model-only recovery is descriptive complementarity. Reviewing both top5 lists uses up to10 unique hits, not a demonstrated total5 fusion improvement.'}
    old = store.read('summary')
    if old is None:
        store.publish('summary', summary)
    else:
        require(old == summary, 'Recomputed disagreement summary differs')
    expected.append('summary')
    store.complete(expected, {'elapsed_seconds': time.monotonic() - started})
    atomic_json(store.root / 'summary_export.json', summary)
    print(json.dumps({'status': 'complete', 'job': str(store.root)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--model', action='append', required=True, help='Frozen roster name=completed-job-path')
    parser.add_argument('--reference-job', type=Path, required=True)
    parser.add_argument('--name', required=True)
    run(parser.parse_args())
