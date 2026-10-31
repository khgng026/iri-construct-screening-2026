#!/usr/bin/env python3
"""Incremental raw local alignment for immutable expanded catalogue prefixes.

Only new, individually admitted references are scored. Original 105-accession
rows are verified and reduced to 101 unique representatives. Reference-selfscore
normalization remains the previously frozen outcome-informed exploratory arm.
Run only in the lead's exclusive resource-guarded compute slot.
"""
from __future__ import annotations
import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from iri_constructs.session_expansion import load_expansion

RETRIEVAL_DESIGN_SHA = '8cf5f2a2184d24465324a589afb5016cddd1421e78055b63c9759f3a08fd1cfc'
ARMS = ('local_alignment_raw', 'local_alignment_reference_selfscore')


def require(value, message):
    if not value: raise ValueError(message)


def read(path): return json.loads(Path(path).read_text())


def public_query(query): return {k: v for k, v in query.items() if k != 'sequence'}


def representative_query(query, aliases, sequences):
    """Evaluation-only label translation, after sequence-only prediction."""
    result = public_query(query)
    result['components'] = [aliases[a] for a in query['components']]
    result['spans'] = [{**span, 'accession': aliases[span['accession']]} for span in query['spans']]
    require(len(set(result['components'])) == query['component_count'], 'Distinct true components collapse to one representative')
    for span in result['spans']:
        require(query['sequence'][span['start0']:span['end0']] == sequences[span['accession']], 'Mapped target sequence differs')
    return result


def rerank(rows, accessions):
    from iri_constructs.session_benchmark import filter_ranking
    require(len(rows) == len(accessions) and {r['accession'] for r in rows} == set(accessions), 'Missing or repeated raw reference')
    return filter_ranking(sorted(rows, key=lambda r: (-r['score'], r['accession'])), accessions)


def validate_raw_rows(rows, sequences, query_length, selfscores):
    from session_independent_audit import validate_ranking
    validate_ranking(rows, sorted(sequences))
    for row in rows:
        a = row['accession']; score = row['score']
        require(math.isfinite(score) and 0 <= score <= selfscores[a] and score * 2 == int(score * 2), 'Invalid raw score or diagonal upper bound')
        require(row['query_length'] == query_length and row['reference_length'] == len(sequences[a]), 'Alignment row sequence length mismatch')
        require(0 <= row['identity_fraction'] <= 1 and 0 <= row['query_coverage'] <= 1 and 0 <= row['reference_coverage'] <= 1,
                'Invalid alignment identity/coverage')
        if score:
            require(row['query_span'] is not None and row['reference_span'] is not None, 'Positive score lacks geometry')
            for span, length in ((row['query_span'], query_length), (row['reference_span'], len(sequences[a]))):
                require(len(span) == 2 and all(isinstance(x, int) for x in span) and 0 <= span[0] < span[1] <= length,
                        'Alignment span outside source sequence')
        else:
            require(row['query_span'] is None and row['reference_span'] is None, 'Zero local score has a spurious span')


def validate_base_reproduction(source_rows, reproduced_rows, base_accessions):
    from iri_constructs.session_benchmark import filter_ranking
    expected = filter_ranking(source_rows, base_accessions)
    require(expected == reproduced_rows, '101-reference raw baseline changed score, geometry, ordering or tie metadata')


def latest_completion(job):
    candidates = [(read(p)['completed_utc'], str(p), p) for p in job.glob('completion-*.json')]
    require(candidates, 'Raw source has no completion descriptor')
    return max(candidates)[2]


def cache_rows(cache, sequence_hashes, query_length, selfscores, sequences):
    """Read only owned immutable shards; exact sequence hashes allow later reuse."""
    result = {}; provenance = []
    units = cache.root / 'units'
    if not units.exists(): return result, provenance
    for key_path in sorted(units.iterdir()):
        require(key_path.is_dir() and key_path.name.startswith('refs-'), 'Unexpected alignment cache unit')
        path = cache.find(key_path.name)
        if path is None: continue  # Corrupt attempts remain, then get recomputed.
        record = read(path / 'record.json')
        hashes = record['reference_sequence_sha256']
        require(record['query_sequence_sha256'] == cache.identity['query_sequence_sha256'], 'Cache query hash mismatch')
        require(key_path.name == 'refs-' + digest(canonical(hashes)), 'Cache shard key/hash mismatch')
        require(len(record['ranking']) == len(hashes) and {r['accession'] for r in record['ranking']} == set(hashes), 'Cache shard is incomplete')
        available = {a: sequences[a] for a in hashes if a in sequence_hashes and hashes[a] == sequence_hashes[a]}
        if not available: continue
        rows = rerank([r for r in record['ranking'] if r['accession'] in available], sorted(available))
        validate_raw_rows(rows, available, query_length, selfscores)
        for row in rows:
            a = row['accession']
            # Ranking/tie fields depend on the shard; geometry/raw score do not.
            semantic = {k: v for k, v in row.items() if k not in TIE_FIELDS}
            require(a not in result or result[a] == semantic, 'Cached optimal score/geometry is inconsistent')
            result[a] = semantic
        provenance.append({'key': key_path.name, 'descriptor_sha256': file_hash(path / 'complete.json')})
    return result, provenance


TIE_FIELDS = {'rank', 'exact_score_tie_count', 'tie_best_rank', 'tie_worst_rank', 'near_tie_count',
              'near_tie_count_1e6', 'near_tie_best_rank', 'near_tie_worst_rank', 'near_tie_tolerance'}


def run(args):
    from iri_constructs.session_alignment import alignment_metadata, rank_alignment
    from iri_constructs.session_benchmark import load_suite, filter_ranking, evaluate_query
    from session_alignment_analysis import normalized_ranking, reference_selfscores
    from session_catalogue_analysis import SummaryAccumulator
    from session_independent_audit import verify_unit, validate_ranking, validate_evaluation, validate_summary
    session, rawjob = args.session.resolve(), args.raw_job.resolve()
    deadline = dt.datetime.fromisoformat(read(session / 'launch.json')['hard_stop_utc'])
    def check_deadline():
        require(dt.datetime.now(dt.timezone.utc) < deadline, 'Session hard stop reached')
    check_deadline(); started = time.monotonic()
    expansion = load_expansion(ROOT, session, args.manifest.resolve())
    suite = load_suite(ROOT, session / 'inputs'); queries = {q['query_id']: q for q in suite['queries']}
    design_path = session / 'inputs/expanded_retrieval_design_v1.json'
    require(file_hash(design_path) == RETRIEVAL_DESIGN_SHA, 'Frozen expanded retrieval design changed')
    design = read(design_path)
    require(design['source_selection_design_sha256'] == expansion['manifest']['design_sha256'] and
            design['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'] and design['budgets'] == [1, 5, 10],
            'Expansion design/source/query identity mismatch')
    norm_path = session / 'audit/alignment_normalization_design_v1.json'; norm_design = read(norm_path)
    require(norm_design['status'] == 'frozen_before_normalized_scoring' and
            norm_design['version'] == 'alignment-reference-selfscore-normalization-v1' and
            norm_design['secondary_score'] == 'S(query, reference) / S(reference, reference)' and
            norm_design['source_challenge_sha256'] == design['challenge_sha256'], 'Existing exploratory normalization design changed')
    raw_identity = read(rawjob / 'identity.json'); raw_hash = digest(canonical(raw_identity))
    require(raw_identity['mode'] == 'alignment' and raw_identity['challenge_sha256'] == design['challenge_sha256'], 'Wrong raw source arm/challenge')
    query_ids = raw_identity['query_ids']
    require(len(query_ids) == len(set(query_ids)) == 1146 and set(query_ids) == set(queries), 'Requires the same completed full 1146-query source')
    require(raw_identity['catalogue_accessions'] == sorted(suite['sequences']) and
            raw_identity['catalogue_sequence_sha256'] == {a: digest(s.encode()) for a, s in suite['sequences'].items()},
            'Cached original reference sequences differ')
    method = alignment_metadata()
    require(raw_identity['method'] == method and raw_identity['code_hashes']['src/iri_constructs/session_alignment.py'] ==
            file_hash(ROOT / 'src/iri_constructs/session_alignment.py'), 'Cached scoring method/code differs')
    completion_path = latest_completion(rawjob); completion = read(completion_path)
    source_expected = {'query-' + q for q in query_ids} | {'summary'}
    require(completion['identity_hash'] == raw_hash and set(completion['units']) == source_expected and
            completion['expected_count'] == len(source_expected), 'Raw completion does not cover the full source')
    source_paths = {key: verify_unit(rawjob, key, entry, raw_hash) for key, entry in completion['units'].items()}
    sequences, subsets = expansion['sequences'], expansion['subsets']
    accessions = sorted(sequences); base = subsets['101']; aliases = expansion['target_aliases']
    require(all(sequences[a] == suite['sequences'][a] for a in base), 'Cached base representative sequences differ')
    for query in queries.values():
        mapped = representative_query(query, aliases, sequences)
        require(all(set(mapped['components']) <= set(panel) for panel in subsets.values()), 'Target missing from a reported endpoint')
    scores, _, matrix_hash = reference_selfscores(sequences, check_examples=False)
    require(matrix_hash == norm_design['matrix_sha256'] == method['matrix_sha256'], 'Normalized and raw matrix identities differ')
    hashes = {a: digest(s.encode()) for a, s in sequences.items()}
    new_order = [a for a in expansion['manifest']['new_representatives_hash_order'] if a in sequences and a not in base]
    require(set(new_order) == set(sequences)-set(base), 'New reference partition differs')
    require(1 <= args.shard_size <= 32, 'Bounded alignment shard size must be 1..32')
    shards = [new_order[i:i + args.shard_size] for i in range(0, len(new_order), args.shard_size)]
    code = {name: file_hash(ROOT / name) for name in ('scripts/session_expanded_alignment.py',
        'src/iri_constructs/session_expansion.py', 'src/iri_constructs/session_alignment.py', 'src/iri_constructs/session_store.py',
        'src/iri_constructs/session_benchmark.py', 'scripts/session_alignment_analysis.py',
        'scripts/session_catalogue_analysis.py', 'scripts/session_independent_audit.py')}
    identity = {'mode': 'expanded_alignment_retrieval_v1', 'expansion_manifest_sha256': expansion['manifest_sha256'],
        'retrieval_design_sha256': RETRIEVAL_DESIGN_SHA, 'normalization_design_sha256': file_hash(norm_path),
        'source_raw_identity': raw_hash, 'source_raw_completion_sha256': file_hash(completion_path),
        'challenge_sha256': design['challenge_sha256'], 'query_ids': query_ids,
        'query_sequence_sha256': {qid: digest(queries[qid]['sequence'].encode()) for qid in query_ids},
        'method': method, 'subsets': subsets, 'target_aliases': aliases, 'reference_sequence_sha256': hashes,
        'reference_selfscores': scores, 'new_reference_shards': shards, 'code': code,
        'scope': 'Outcome-informed exploratory catalogue competition; size and composition change together; raw primary and fixed reference-selfscore secondary; K total unique representatives per query'}
    jobid = args.name + '-' + digest(canonical(identity))[:12]
    store = UnitStore(session / 'jobs' / jobid, identity, session / 'job_ledger.jsonl')
    cache_method = {'mode': 'expanded_alignment_sequence_cache_v1', 'method': method,
                    'alignment_code_sha256': code['src/iri_constructs/session_alignment.py']}
    cache_root = session / 'caches' / ('expanded-alignment-' + digest(canonical(cache_method))[:16])
    expected = ['reference-selfscores', 'summary']
    score_record = {'scores': scores, 'matrix_sha256': matrix_hash, 'derivation': 'sum canonical BLOSUM62 diagonal; previous independent full-source score audit; no new self-alignment'}
    old = store.read('reference-selfscores')
    if old is None: store.publish('reference-selfscores', score_record)
    else: require(old == score_record, 'Saved selfscores differ')
    summaries = {size: {arm: {scope: SummaryAccumulator() for scope in ('clean', 'combined')} for arm in ARMS} for size in subsets}
    independent = {size: {arm: {scope: [] for scope in ('clean', 'combined')} for arm in ARMS} for size in subsets}
    total_new = 0; total_cache = 0
    store.event(event='job_started', state='running', job_id=jobid, pid=os.getpid())
    try:
        for index, qid in enumerate(query_ids):
            check_deadline(); query = queries[qid]; sequence = query['sequence']; qsha = digest(sequence.encode())
            source_path = source_paths['query-' + qid]; source_record = read(source_path / 'record.json')
            require(source_record['query'] == public_query(query), 'Raw source query identity/sequence hash differs')
            source_rows = source_record['rankings']['local_alignment']
            validate_ranking(source_rows, sorted(suite['sequences']))
            base_rows = filter_ranking(source_rows, base)
            validate_raw_rows(base_rows, {a: sequences[a] for a in base}, len(sequence), scores)
            cache = UnitStore(cache_root / qsha, {**cache_method, 'query_sequence_sha256': qsha, 'query_length': len(sequence)})
            saved_rows, cache_provenance = cache_rows(cache, hashes, len(sequence), scores, sequences)
            rows = list(base_rows); shard_provenance = []
            for shard_index, shard in enumerate(shards):
                check_deadline(); key = f'raw-{qid}-{shard_index:03d}'; expected.append(key)
                shard_sequences = {a: sequences[a] for a in shard}; shard_hashes = {a: hashes[a] for a in shard}
                record = store.read(key)
                if record is None:
                    missing = [a for a in shard if a not in saved_rows]
                    if missing:
                        predicted = rank_alignment(sequence, {a: sequences[a] for a in missing})
                        require(predicted['method'] == method, 'Scoring method changed during a shard')
                        validate_raw_rows(predicted['ranking'], {a: sequences[a] for a in missing}, len(sequence), scores)
                        missing_hashes = {a: hashes[a] for a in missing}; cache_key = 'refs-' + digest(canonical(missing_hashes))
                        cache_record = {'query_sequence_sha256': qsha, 'reference_sequence_sha256': missing_hashes,
                                        'ranking': predicted['ranking']}
                        cached = cache.read(cache_key)
                        if cached is None: cache.publish(cache_key, cache_record)
                        else: require(cached == cache_record, 'Immutable sequence cache differs')
                        for row in predicted['ranking']:
                            saved_rows[row['accession']] = {k: v for k, v in row.items() if k not in TIE_FIELDS}
                        total_new += len(missing)
                    total_cache += len(shard) - len(missing)
                    shard_rows = rerank([saved_rows[a] for a in shard], shard)
                    record = {'query_sequence_sha256': qsha, 'reference_sequence_sha256': shard_hashes,
                              'ranking': shard_rows, 'newly_computed_pairs': len(missing),
                              'reused_pairs': len(shard)-len(missing), 'existing_cache_descriptors': cache_provenance}
                    store.publish(key, record)
                require(record['query_sequence_sha256'] == qsha and record['reference_sequence_sha256'] == shard_hashes,
                        'Saved shard query/reference identities differ')
                validate_raw_rows(record['ranking'], shard_sequences, len(sequence), scores)
                rows.extend(record['ranking'])
                path = store.find(key)
                shard_provenance.append({'key': key, 'descriptor_sha256': file_hash(path / 'complete.json')})
            raw = rerank(rows, accessions)
            validate_raw_rows(raw, sequences, len(sequence), scores)
            validate_base_reproduction(source_rows, filter_ranking(raw, base), base)
            normalized = normalized_ranking(raw, scores)
            mapped = representative_query(query, aliases, sequences)
            record = {'query': public_query(query), 'representative_label_query': mapped,
                      'rankings': {}, 'evaluations': {}, 'source_raw_unit_descriptor_sha256': file_hash(source_path / 'complete.json'),
                      'raw_shard_descriptors': shard_provenance, 'base101_exact_reproduction': True}
            scope = 'clean' if query['component_count'] == 1 else 'combined'
            for size, subset in subsets.items():
                equivalence = {a: [a] for a in subset}
                rankings = {ARMS[0]: filter_ranking(raw, subset), ARMS[1]: filter_ranking(normalized, subset)}
                evaluations = {}
                for arm, ranking in rankings.items():
                    validate_ranking(ranking, subset)
                    evaluation = evaluate_query(ranking, mapped, equivalence)
                    independent[size][arm][scope].append(validate_evaluation(evaluation, mapped, ranking, equivalence))
                    summaries[size][arm][scope].add(evaluation); evaluations[arm] = evaluation
                record['rankings'][size] = rankings; record['evaluations'][size] = evaluations
            key = 'query-' + qid; expected.append(key)
            old = store.read(key)
            if old is None: store.publish(key, record)
            else: require(old == record, 'Recomputed expanded alignment query differs')
            if (index + 1) % 50 == 0 or index + 1 == len(query_ids):
                print(json.dumps({'job': jobid, 'queries_completed': index + 1, 'expected': len(query_ids),
                                  'new_pairs_this_invocation': total_new, 'cache_pairs_this_invocation': total_cache}), flush=True)
        result = {'summaries': {}, 'all_query_metrics_and_summary_denominators_independently_checked': True,
                  'original105_unchanged': True, 'base101_exact_reproduction_query_count': len(query_ids)}
        for size, arms in summaries.items():
            result['summaries'][size] = {}
            for arm, scopes in arms.items():
                result['summaries'][size][arm] = {}
                for scope, accumulator in scopes.items():
                    summary = accumulator.result()
                    validate_summary(summary, independent[size][arm][scope])
                    count = 16 if scope == 'clean' else 1130
                    components = 16 if scope == 'clean' else 3520
                    require(summary['query_count'] == count and summary['component_count'] == components, 'Expanded summary denominator mismatch')
                    result['summaries'][size][arm][scope] = summary
        old = store.read('summary')
        if old is None: store.publish('summary', result)
        else: require(old == result, 'Recomputed expanded alignment summary differs')
        store.complete(expected, {'elapsed_seconds': time.monotonic()-started,
                                 'new_pairs_this_invocation': total_new, 'cache_pairs_this_invocation': total_cache})
        print(json.dumps({'status': 'complete', 'job': str(store.root), 'elapsed_seconds': time.monotonic()-started}), flush=True)
    except BaseException as exc:
        store.event(event='job_failed', state='failed', job_id=jobid, exception=f'{type(exc).__name__}: {exc}')
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--session', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--raw-job', type=Path, required=True)
    p.add_argument('--name', default='expanded-alignment')
    p.add_argument('--shard-size', type=int, default=32)
    run(p.parse_args())
