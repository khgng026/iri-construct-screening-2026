#!/usr/bin/env python3
"""Encode admitted catalogue competitors, then rerank unchanged saved queries.

The original 105-accession experiment remains separate. Nested comparisons use
one representative per exact sequence, including a common 101-sequence baseline.
"""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from iri_constructs.session_benchmark import load_suite, cosine_ranking, evaluate_query, evaluate_oracle_query, filter_ranking
from session_experiment import check_deadline, mean_residues, model_identity, public_query, peak_mib
from session_segment_experiment import open_store, completed_read, validate_source_suite
from session_catalogue_analysis import SummaryAccumulator, latest_completion
from session_independent_audit import verify_unit, validate_evaluation, validate_summary, validate_ranking


def require(value, message):
    if not value:
        raise ValueError(message)


def source_hashes():
    return {name: file_hash(ROOT / name) for name in [
        'scripts/session_expanded_model.py', 'src/iri_constructs/session_expansion.py',
        'src/iri_constructs/session_models.py', 'src/iri_constructs/session_benchmark.py',
        'src/iri_constructs/session_store.py', 'scripts/session_experiment.py',
        'scripts/session_catalogue_analysis.py', 'scripts/session_segment_experiment.py',
        'scripts/session_independent_audit.py']}


def setup(args):
    from iri_constructs.session_expansion import load_expansion
    session = args.session.resolve()
    check_deadline(session)
    expansion = load_expansion(ROOT, session, args.manifest.resolve())
    suite = load_suite(ROOT, session / 'inputs')
    design_path = session / 'inputs/expanded_retrieval_design_v1.json'
    design = json.loads(design_path.read_text())
    require(design['version'] == 'iri-expanded-retrieval-v1' and design['status'] == 'frozen_before_expanded_scoring', 'Expanded retrieval design is not frozen')
    require(design['budgets'] == [1, 5, 10] and design['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'], 'Expanded retrieval budgets or challenge differ')
    require(design['source_selection_design_sha256'] == expansion['manifest']['design_sha256'], 'Expanded source selection design differs')
    expansion['retrieval_design_sha256'] = file_hash(design_path)
    parent = open_store(args.parent_job.resolve(), session)
    require(parent.identity['mode'] == 'model', 'Requires a completed representation job')
    require(set(parent.identity['query_ids']) == {q['query_id'] for q in suite['queries']}, 'Primary expansion requires all frozen queries')
    validate_source_suite(parent.identity, suite, suite['verification']['challenge_manifest_sha256'])
    require(parent.identity['code_hashes']['scripts/session_experiment.py'] == file_hash(ROOT / 'scripts/session_experiment.py'),
            'Saved query pooling implementation differs')
    require(parent.identity['code_hashes']['src/iri_constructs/session_models.py'] == file_hash(ROOT / 'src/iri_constructs/session_models.py'),
            'Saved query model implementation differs')
    accessions = sorted(set().union(*(set(v) for v in expansion['subsets'].values())))
    require('101' in expansion['subsets'] and len(expansion['subsets']['101']) == 101, 'Common unique-sequence baseline missing')
    require(set(accessions) <= set(expansion['sequences']), 'Nested catalogue lacks source sequences')
    aliases = expansion['target_aliases']
    for q in suite['queries']:
        require(len({aliases[a] for a in q['components']}) == q['component_count'], 'Distinct component labels collapse in expanded catalogue')
        require(all(aliases[a] in subset for a in q['components'] for subset in expansion['subsets'].values()), 'A target is absent from a frozen nested catalogue')
        for a in q['components']:
            require(expansion['sequences'][aliases[a]] == suite['sequences'][a], 'Target representative sequence differs')
    return session, expansion, suite, parent, accessions


def validate_vector(record, arrays, sequence, dimension=None):
    import numpy as np
    require(set(arrays) == {'mean'}, 'Unexpected catalogue arrays')
    vector = arrays['mean']
    require(record['sequence_sha256'] == digest(sequence.encode()) and record['length'] == len(sequence), 'Catalogue sequence metadata differs')
    require(vector.ndim == 1 and vector.dtype == np.float32 and np.isfinite(vector).all(), 'Invalid catalogue vector')
    require(dimension is None or vector.shape == (dimension,), 'Catalogue dimensions differ')
    return vector


def validate_base_reproduction(parent_rows, expanded_rows, accessions):
    expected = filter_ranking(parent_rows, accessions)
    require([r['accession'] for r in expected] == [r['accession'] for r in expanded_rows],
            '101-reference reranking changed native parent rank order')
    error = max(abs(a['score'] - b['score']) for a, b in zip(expected, expanded_rows))
    require(error <= 1e-12, '101-reference reranking changed source cosine beyond tolerance')
    return error


def encode(args):
    import numpy as np
    from iri_constructs.session_models import load_encoder
    session, expansion, suite, parent, accessions = setup(args)
    identity = {'mode': 'expanded_catalogue_vectors_v1', 'expansion_manifest_sha256': expansion['manifest_sha256'],
                'retrieval_design_sha256': expansion['retrieval_design_sha256'],
                'parent_identity': parent.identity_hash, 'parent_completion_sha256': parent.completion_sha256,
                'method': parent.identity['method'], 'config': parent.identity['config'],
                'reference_accessions': accessions, 'reference_sequence_sha256': {a: digest(expansion['sequences'][a].encode()) for a in accessions},
                'code': source_hashes(), 'pooling': 'same float64 residue mean then stored float32 as parent'}
    store = UnitStore(session / 'jobs' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    cache_identity = {'mode': 'expanded_reference_mean_cache_v1', 'method': parent.identity['method'],
                      'pooling': identity['pooling'], 'code': identity['code']}
    cache = UnitStore(session / 'caches' / ('expanded-means-' + digest(canonical(cache_identity))[:16]), cache_identity, session / 'job_ledger.jsonl')
    base = {}
    dimension = None
    for accession in parent.identity['catalogue_accessions']:
        record, arrays = completed_read(parent, 'catalogue-' + accession, arrays=True)
        vector = validate_vector(record, arrays, suite['sequences'][accession], dimension)
        dimension = len(vector)
        sha = record['sequence_sha256']
        require(sha not in base or np.array_equal(base[sha], vector), 'Duplicate source sequences have different vectors')
        base[sha] = vector
    encoder = None
    start = time.monotonic()
    try:
        for i, accession in enumerate(accessions):
            check_deadline(session)
            sequence = expansion['sequences'][accession]
            sha = digest(sequence.encode())
            key = 'catalogue-' + accession
            old = store.read(key, arrays=True)
            if old is not None:
                validate_vector(*old, sequence, dimension)
                continue
            saved = cache.read(sha, arrays=True)
            if saved is not None:
                record, arrays = saved
                vector = validate_vector(record, arrays, sequence, dimension)
                origin = 'verified_cross_milestone_cache'
            elif sha in base:
                vector = base[sha]
                origin = 'verified_original_catalogue'
            else:
                if encoder is None:
                    encoder = load_encoder(parent.identity['config'])
                    require(model_identity(encoder) == parent.identity['method'], 'Current model differs from frozen query model')
                t = time.monotonic()
                states = encoder.encode(sequence)
                vector = mean_residues(states)
                del states
                record = {'sequence_sha256': sha, 'length': len(sequence), 'encode': encoder.last_encode_metadata,
                          'seconds': time.monotonic() - t, 'peak_process_rss_mib': peak_mib()}
                validate_vector(record, {'mean': vector}, sequence, dimension)
                cache.publish(sha, record, {'mean': vector})
                origin = 'new_forward'
            store.publish(key, {'accession': accession, 'sequence_sha256': sha, 'length': len(sequence), 'origin': origin}, {'mean': vector})
            if i % 50 == 0:
                print(json.dumps({'stage': 'expanded_encode', 'references': i + 1, 'total': len(accessions)}), flush=True)
        store.complete(['catalogue-' + a for a in accessions], {'elapsed_seconds': time.monotonic() - start})
        print(json.dumps({'status': 'complete', 'job': str(store.root)}), flush=True)
    finally:
        if encoder is not None:
            encoder.close()


def open_encoding(path, session, expansion, parent, accessions):
    identity = json.loads((path / 'identity.json').read_text())
    identity_hash = digest(canonical(identity))
    require(identity['mode'] == 'expanded_catalogue_vectors_v1' and identity['expansion_manifest_sha256'] == expansion['manifest_sha256'], 'Encoding job refers to different expansion')
    require(identity['retrieval_design_sha256'] == expansion['retrieval_design_sha256'], 'Encoding retrieval design differs')
    require(identity['parent_identity'] == parent.identity_hash and identity['method'] == parent.identity['method'], 'Encoding/query representation mismatch')
    require(identity['parent_completion_sha256'] == parent.completion_sha256, 'Encoding parent completion differs')
    require(identity['reference_accessions'] == accessions and identity['code'] == source_hashes(), 'Encoding reference set or implementation differs')
    require(identity['reference_sequence_sha256'] == {a: digest(expansion['sequences'][a].encode()) for a in accessions}, 'Encoding reference sequence identities differ')
    descriptor = latest_completion(path)
    completion = json.loads(descriptor.read_text())
    expected = {'catalogue-' + a for a in accessions}
    require(completion['identity_hash'] == identity_hash and set(completion['units']) == expected and completion['expected_count'] == len(expected), 'Incomplete catalogue vectors')
    vectors = {}
    dimension = None
    import numpy as np
    for key, entry in completion['units'].items():
        unit = verify_unit(path, key, entry, identity_hash)
        record = json.loads((unit / 'record.json').read_text())
        accession = key.removeprefix('catalogue-')
        require(record['accession'] == accession, 'Encoding accession mismatch')
        with np.load(unit / 'arrays.npz', allow_pickle=False) as z:
            arrays = {k: z[k].copy() for k in z.files}
        vector = validate_vector(record, arrays, expansion['sequences'][accession], dimension)
        dimension = len(vector)
        vectors[accession] = vector
    return identity_hash, file_hash(descriptor), vectors


def rank(args):
    import numpy as np
    session, expansion, suite, parent, accessions = setup(args)
    encoding_hash, completion_hash, vectors = open_encoding(args.encoding_job.resolve(), session, expansion, parent, accessions)
    identity = {'mode': 'expanded_model_retrieval_v1', 'expansion_manifest_sha256': expansion['manifest_sha256'],
                'retrieval_design_sha256': expansion['retrieval_design_sha256'],
                'encoding_identity': encoding_hash, 'encoding_completion_sha256': completion_hash,
                'parent_identity': parent.identity_hash, 'parent_completion_sha256': parent.completion_sha256,
                'query_ids': parent.identity['query_ids'], 'subsets': expansion['subsets'],
                'target_aliases': expansion['target_aliases'], 'code': source_hashes(),
                'scope': 'Exploratory exact-component catalogue competition; size and composition change together; no training-set holdout'}
    store = UnitStore(session / 'jobs' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    queries = {q['query_id']: q for q in suite['queries']}
    arms = ['whole', 'isolated_mixture', 'contextual_component', 'isolated_component']
    summaries = {size: {arm: {scope: SummaryAccumulator() for scope in ('clean', 'combined')} for arm in arms}
                 for size in expansion['subsets']}
    independently_checked = {size: {arm: {scope: [] for scope in ('clean', 'combined')} for arm in arms}
                             for size in expansion['subsets']}
    matrices = {size: np.stack([vectors[a] for a in subset]) for size, subset in expansion['subsets'].items()}
    expected = []
    for qid in parent.identity['query_ids']:
        check_deadline(session)
        original = queries[qid]
        source, arrays = completed_read(parent, 'query-' + qid, arrays=True)
        require(source['query'] == public_query(original), 'Saved query metadata changed')
        query = public_query(original)
        query['components'] = [expansion['target_aliases'][a] for a in original['components']]
        query['spans'] = [{**span, 'accession': expansion['target_aliases'][span['accession']]} for span in original['spans']]
        record = {'query': public_query(original), 'representative_label_query': query, 'rankings': {}, 'evaluations': {},
                  'base101_cosine_max_abs': {}}
        for size, subset in expansion['subsets'].items():
            refs = matrices[size]
            rankings = {arm: cosine_ranking(arrays[arm], refs, subset, dtype=np.float64) for arm in ('whole', 'isolated_mixture')}
            for arm, prefix in [('contextual_component', 'contextual_'), ('isolated_component', 'isolated_')]:
                rankings[arm] = {expansion['target_aliases'][a]: cosine_ranking(arrays[prefix + a], refs, subset, dtype=np.float64)
                                 for a in original['components']}
            if size == '101':
                for arm in arms:
                    if isinstance(rankings[arm], dict):
                        errors = [validate_base_reproduction(source['rankings'][arm][a], rankings[arm][expansion['target_aliases'][a]], subset)
                                  for a in original['components']]
                        record['base101_cosine_max_abs'][arm] = max(errors)
                    else:
                        record['base101_cosine_max_abs'][arm] = validate_base_reproduction(source['rankings'][arm], rankings[arm], subset)
            evaluations = {}
            for arm, rows in rankings.items():
                evaluator = evaluate_oracle_query if isinstance(rows, dict) else evaluate_query
                evaluation = evaluator(rows, query, {})
                for native_rows in rows.values() if isinstance(rows, dict) else [rows]:
                    validate_ranking(native_rows, subset)
                independent = validate_evaluation(evaluation, query, rows, {a: [a] for a in subset})
                evaluations[arm] = evaluation
                scope = 'clean' if query['component_count'] == 1 else 'combined'
                summaries[size][arm][scope].add(evaluation)
                independently_checked[size][arm][scope].append(independent)
            record['rankings'][size] = rankings
            record['evaluations'][size] = evaluations
        key = 'query-' + qid
        old = store.read(key)
        if old is None:
            store.publish(key, record)
        else:
            require(old == record, 'Recomputed expanded rankings differ')
        expected.append(key)
    result = {'summaries': {size: {arm: {scope: acc.result() for scope, acc in scopes.items()} for arm, scopes in by_arm.items()}
                            for size, by_arm in summaries.items()}}
    for size, by_arm in independently_checked.items():
        for arm, scopes in by_arm.items():
            for scope, independent in scopes.items():
                summary = result['summaries'][size][arm][scope]
                validate_summary(summary, independent)
                require(summary['query_count'] == (16 if scope == 'clean' else 1130), 'Expanded query denominator changed')
                require(summary['component_count'] == (16 if scope == 'clean' else 3520), 'Expanded component denominator changed')
    result['independent_acceptance'] = {'all_rankings_evaluations_summaries_checked': True,
                                        'clean_queries_components': [16, 16], 'combined_queries_components': [1130, 3520]}
    old = store.read('summary')
    if old is None:
        store.publish('summary', result)
    else:
        require(old == result, 'Recomputed expanded summary differs')
    store.complete(expected + ['summary'])
    print(json.dumps({'status': 'complete', 'job': str(store.root)}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['encode', 'rank'])
    p.add_argument('--session', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--parent-job', type=Path, required=True)
    p.add_argument('--encoding-job', type=Path)
    p.add_argument('--name', required=True)
    args = p.parse_args()
    if args.stage == 'rank' and args.encoding_job is None:
        p.error('rank requires --encoding-job')
    (encode if args.stage == 'encode' else rank)(args)
