#!/usr/bin/env python3
"""Known-boundary retrieval controls at a fixed total returned-accession budget.

This exploratory analysis reuses completed component rankings. It performs no
inference and does not claim automatic component discovery.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from iri_constructs.session_benchmark import evaluate_query, load_suite
from iri_constructs.session_segmentation import fuse_rankings
from session_catalogue_analysis import SummaryAccumulator
from session_experiment import check_deadline, public_query
from session_segment_experiment import open_store, completed_read, validate_source_suite


def require(value, message):
    if not value:
        raise ValueError(message)


def independent_fusion_check(rankings, fused, method):
    """Reconstruct scores/order independently of the production fusion helper."""
    scores = {}
    for rows in rankings:
        require([r['rank'] for r in rows] == list(range(1, len(rows) + 1)), 'Nonconsecutive source ranks')
        require(len({r['accession'] for r in rows}) == len(rows), 'Duplicate source accession')
        for row in rows:
            value = 1.0 / (60 + row['rank']) if method == 'rrf' else row['cosine_similarity']
            require(math.isfinite(value), 'Nonfinite source score')
            scores.setdefault(row['accession'], []).append(value)
    require(all(len(v) == len(rankings) for v in scores.values()), 'Component catalogues differ')
    expected = {a: math.fsum(v) if method == 'rrf' else max(v) for a, v in scores.items()}
    order = sorted(expected, key=lambda a: (-expected[a], a))
    require([r['accession'] for r in fused] == order, 'Independent fusion ordering differs')
    require(all(row['fusion_score'] == expected[row['accession']] for row in fused), 'Fusion score differs')


def run(args):
    session = args.session.resolve()
    design = json.loads(args.design.read_text())
    require(design['version'] == 'iri-oracle-total-budget-v1', 'Unknown frozen design')
    require(design['status'] == 'frozen_before_oracle_fusion_scoring', 'Design is not frozen')
    require(design['budgets'] == [1, 5, 10] and design['known_boundaries'] is True, 'Budget or information policy differs')
    require(design['arms'] == ['contextual_component', 'isolated_component'], 'Oracle arms changed')
    require(design['fusions'] == ['max_score', 'rrf'] and design['rrf_constant'] == 60, 'Fusion design changed')
    suite = load_suite(ROOT, session / 'inputs')
    require(design['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'], 'Challenge changed')
    source = open_store(args.source_job, session)
    require(source.identity['mode'] == 'model', 'Requires a completed representation job')
    validate_source_suite(source.identity, suite, design['challenge_sha256'])
    queries = {q['query_id']: q for q in suite['queries']}
    accessions = sorted(suite['sequences'])
    ids = source.identity['query_ids']
    require(set(ids) <= set(queries), 'Source contains unknown query')
    identity = {
        'mode': 'oracle_total_budget_v1', 'design_sha256': file_hash(args.design),
        'source_path': str(source.root), 'source_identity': source.identity_hash,
        'source_completion_sha256': source.completion_sha256,
        'query_ids': ids, 'challenge_sha256': design['challenge_sha256'],
        'coverage': {'full_suite': set(ids) == set(queries), 'observed_queries': len(ids),
                     'expected_suite_queries': len(queries)},
        'catalogue_accessions': accessions,
        'code': {name: file_hash(ROOT / name) for name in [
            'scripts/session_oracle_fusion.py', 'scripts/session_catalogue_analysis.py',
            'scripts/session_segment_experiment.py', 'src/iri_constructs/session_segmentation.py',
            'src/iri_constructs/session_benchmark.py', 'src/iri_constructs/session_store.py']},
        'information_budget': 'Known complete component boundaries; k total unique accessions per query',
        'interpretation': 'Exploratory fixed-budget oracle control; no automatic discovery or causal decomposition claim',
    }
    store = UnitStore(session / 'jobs' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    accumulators = {arm + '__' + fusion: {scope: SummaryAccumulator() for scope in ('clean', 'combined')}
                    for arm in design['arms'] for fusion in design['fusions']}
    expected = []
    start = time.monotonic()
    for qid in ids:
        check_deadline(session)
        query = queries[qid]
        parent = completed_read(source, 'query-' + qid)
        require(parent['query'] == public_query(query), 'Source query metadata differs')
        record = {'query': public_query(query), 'rankings': {}, 'evaluations': {},
                  'known_component_count': query['component_count'], 'boundaries_supplied': True}
        for arm in design['arms']:
            source_rankings = parent['rankings'][arm]
            require(set(source_rankings) == set(query['components']), 'Component rankings missing or extra')
            rankings = [source_rankings[a] for a in query['components']]
            for fusion in design['fusions']:
                name = arm + '__' + fusion
                result = fuse_rankings(rankings, fusion, catalogue_accessions=accessions)
                independent_fusion_check(rankings, result['ranking'], fusion)
                evaluation = evaluate_query(result['ranking'], query, suite['equivalence'])
                # Verify exact fixed-budget component/query hits without the evaluator.
                ranks = {r['accession']: r['rank'] for r in result['ranking']}
                for k in (1, 5, 10):
                    hits = [ranks[a] <= k for a in query['components']]
                    require(sum(c['exact'][f'hit_at_{k}'] for c in evaluation['components']) == sum(hits), 'Exact hit denominator mismatch')
                    require(evaluation['query_metrics']['exact'][f'all_components_hit_at_{k}'] == all(hits), 'All-hit mismatch')
                    require(evaluation['query_metrics']['exact'][f'any_component_hit_at_{k}'] == any(hits), 'Any-hit mismatch')
                record['rankings'][name] = result['ranking']
                record['evaluations'][name] = evaluation
                scope = 'clean' if query['component_count'] == 1 else 'combined'
                accumulators[name][scope].add(evaluation)
        key = 'query-' + qid
        old = store.read(key)
        if old is None:
            store.publish(key, record)
        else:
            require(old == record, 'Recomputed fixed-budget result differs')
        expected.append(key)
    summary = {'summaries': {arm: {scope: acc.result() for scope, acc in by_scope.items()}
                             for arm, by_scope in accumulators.items()}, 'design': design}
    old = store.read('summary')
    if old is None:
        store.publish('summary', summary)
    else:
        require(old == summary, 'Recomputed oracle summary differs')
    expected.append('summary')
    store.complete(expected, {'elapsed_seconds': time.monotonic() - start})
    print(json.dumps({'status': 'complete', 'job': str(store.root), 'queries': len(ids)}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--session', type=Path, required=True)
    p.add_argument('--design', type=Path, required=True)
    p.add_argument('--source-job', type=Path, required=True)
    p.add_argument('--name', required=True)
    run(p.parse_args())
