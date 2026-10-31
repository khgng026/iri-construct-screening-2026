#!/usr/bin/env python3
"""Exploratory equal raw-component means from pinned saved vectors only."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
import session_extension_audit as audit
from iri_constructs.session_benchmark import load_suite, cosine_ranking
from iri_constructs.session_store import UnitStore, canonical, digest, safe_name

require, read, sha = audit.require, audit.read, audit.sha
DESIGN_SHA = 'b62aadb8f0085103984174d2c4674a9a1d8af5bf5a1eca80e66afcca24f7c3a0'
MODELS = ('esm2_8m', 'esm2_150m', 'esm2_650m', 'esmc_300m', 'esmc_600m')
MATCHED = {'equal_isolated': 'isolated_mixture', 'equal_contextual': 'whole'}
ARMS = ('whole', 'isolated_mixture', 'equal_isolated', 'equal_contextual')
KS = (1, 5, 10)
FRACTIONS = ('lt_0.25', '0.25_to_lt_0.40', '0.40_to_lt_0.60', '0.60_to_lt_0.75', 'ge_0.75')
VECTOR_TOL = 1e-6
SCORE_TOL = 1e-12


def code_hashes():
    return {name: sha(ROOT / name) for name in (
        'scripts/session_weighting_controls.py', 'tests/test_session_weighting_controls.py',
        'scripts/session_extension_audit.py', 'scripts/session_independent_audit.py',
        'src/iri_constructs/session_store.py', 'src/iri_constructs/session_benchmark.py',
        'src/iri_constructs/t05_composition.py', 'src/iri_constructs/t05_data.py')}


def peak_mib():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (1024 ** 2 if sys.platform == 'darwin' else 1024)


class Deadline(audit.Deadline):
    def check(self):
        super().check()
        require(peak_mib() <= 512, 'Weighting controls exceeded512 MiB')


def aggregate(vectors, lengths=None):
    """Raw vectors, sorted summation, float64 reduction, one float32 cast."""
    import numpy as np
    names = sorted(vectors)
    require(names and (lengths is None or set(lengths) == set(names)), 'Invalid aggregation components')
    width = vectors[names[0]].shape
    require(len(width) == 1 and width[0] > 0 and all(v.shape == width and v.dtype == np.float32 and
            np.isfinite(v).all() for v in vectors.values()), 'Invalid raw component vectors')
    weights = {a: 1 if lengths is None else lengths[a] for a in names}
    require(all(type(w) is int and w > 0 for w in weights.values()), 'Invalid component lengths')
    total = np.zeros(width, dtype=np.float64)
    for a in names:
        total += weights[a] * vectors[a].astype(np.float64)
    result = total / sum(weights.values())
    return result.astype(np.float32), result


def read_arrays(snapshot, key):
    import numpy as np
    path = audit.verify_unit(snapshot.path, key, snapshot.entries[key], snapshot.identity_hash)
    with np.load(path / 'arrays.npz', allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}


def baseline_reproduction(saved, recomputed, accessions):
    audit.validate_ranking(saved, accessions); audit.validate_ranking(recomputed, accessions)
    original = {r['accession']: r for r in saved}; new = {r['accession']: r for r in recomputed}
    maximum = max(abs(audit.score(original[a]) - audit.score(new[a])) for a in accessions)
    changes = [{'accession': a, 'saved_rank': original[a]['rank'], 'recomputed_rank': new[a]['rank']}
               for a in accessions if original[a]['rank'] != new[a]['rank']]
    inversions = []
    for i, left in enumerate(saved):
        for right in saved[i+1:]:
            a, b = left['accession'], right['accession']
            if new[a]['rank'] > new[b]['rank']:
                gaps = (abs(audit.score(left)-audit.score(right)), abs(audit.score(new[a])-audit.score(new[b])))
                inversions.append({'accessions': [a, b], 'saved_gap': gaps[0], 'recomputed_gap': gaps[1],
                                   'within_fixed_tolerance': max(gaps) <= 2*SCORE_TOL})
    return {'status': 'passed' if maximum <= SCORE_TOL and all(p['within_fixed_tolerance'] for p in inversions) else 'failed',
            'maximum_absolute_score_error': maximum, 'changed_candidate_ranks': changes, 'inverted_pairs': inversions,
            'saved_baseline_ranks_preserved': True}


def evaluate(query, ranking):
    rows = {r['accession']: r for r in ranking}
    require(set(query['components']) <= set(rows), 'Exact target absent from full catalogue')
    require(len(query['components']) == len(query['spans']) == query['component_count'] and
            [s['accession'] for s in query['spans']] == query['components'], 'Component geometry differs')
    components = []
    for index, span in enumerate(query['spans']):
        a = span['accession']; row = rows[a]
        position = 'single' if query['component_count'] == 1 else ('N' if index == 0 else 'C' if index == query['component_count']-1 else 'internal')
        require(span['fraction'] == (span['end0']-span['start0'])/query['length_aa'], 'Original fraction differs')
        components.append({'accession': a, 'ordinal_position': index+1, 'position': position,
            'component_count': query['component_count'], 'fraction': span['fraction'],
            'fraction_bin': audit.bin_name(span['fraction']), 'length_aa': span['end0']-span['start0'],
            'start0': span['start0'], 'end0': span['end0'], 'rank': row['rank'], 'score': audit.score(row),
            'hits': {str(k): int(row['rank'] <= k) for k in KS}})
    return {'components': components, 'query': {str(k): {
        'any': int(any(c['hits'][str(k)] for c in components)),
        'all': int(all(c['hits'][str(k)] for c in components))} for k in KS},
        'retrieval_budget': 'k total unique catalogue accessions per original query', 'label_rule': 'exact_accession'}


def outcome_pair(before, after):
    require(before in (0, 1) and after in (0, 1), 'Nonbinary paired outcome')
    return ('both_hit' if before else 'gain') if after else ('loss' if before else 'neither')


def paired(before, after):
    require([c['accession'] for c in before['components']] == [c['accession'] for c in after['components']], 'Paired target order differs')
    return {'components': [{**{key: b[key] for key in ('accession', 'ordinal_position', 'position', 'component_count',
                'fraction', 'fraction_bin', 'length_aa', 'start0', 'end0')},
            'baseline_rank': b['rank'], 'control_rank': a['rank'], 'rank_delta_control_minus_baseline': a['rank']-b['rank'],
            'score_delta_control_minus_baseline': a['score']-b['score'],
            'outcomes': {str(k): outcome_pair(b['hits'][str(k)], a['hits'][str(k)]) for k in KS}}
            for b, a in zip(before['components'], after['components'])],
        'query': {str(k): {name: outcome_pair(before['query'][str(k)][name], after['query'][str(k)][name])
                           for name in ('any', 'all')} for k in KS}}


def derive(query, record, arrays, references):
    import numpy as np
    names = sorted(references); matrix = np.stack([references[a] for a in names])
    width = matrix.shape[1]; ids = sorted(query['components'])
    require(record['query'] == query, 'Saved query metadata differs')
    keys = {'whole', 'isolated_mixture'} | {kind+'_'+a for kind in ('isolated', 'contextual') for a in ids}
    require(set(arrays) == keys and all(v.shape == (width,) and v.dtype == np.float32 and np.isfinite(v).all()
                                      for v in arrays.values()), 'Saved query vector schema differs')
    require(set(record['rankings']) == {'whole', 'isolated_mixture', 'isolated_component', 'contextual_component'}, 'Saved ranking schema differs')
    isolated = {a: arrays['isolated_'+a] for a in ids}; contextual = {a: arrays['contextual_'+a] for a in ids}
    require(all(np.array_equal(isolated[a], references[a]) for a in ids), 'Isolated source vector differs from catalogue')
    lengths = {s['accession']: s['end0']-s['start0'] for s in query['spans']}
    require(sum(lengths.values()) == query['length_aa'], 'Component lengths differ from whole length')
    weighted_isolated, _ = aggregate(isolated, lengths)
    _, weighted_context = aggregate(contextual, lengths)
    new = {'equal_isolated': aggregate(isolated)[0], 'equal_contextual': aggregate(contextual)[0]}
    errors = {'isolated_mixture': float(np.max(np.abs(weighted_isolated.astype(np.float64)-arrays['isolated_mixture']))),
              'whole': float(np.max(np.abs(weighted_context-arrays['whole'].astype(np.float64))))}
    checks = {'weighted_vector_max_abs': errors, 'baseline_cosine': {}}
    rankings = {arm: record['rankings'][arm] for arm in ('whole', 'isolated_mixture')}
    for arm in ('whole', 'isolated_mixture'):
        regenerated = cosine_ranking(arrays[arm], matrix, names, dtype=np.float64)
        checks['baseline_cosine'][arm] = baseline_reproduction(rankings[arm], regenerated, names)
    for arm, vector in new.items():
        rankings[arm] = cosine_ranking(vector, matrix, names, dtype=np.float64)
        audit.validate_ranking(rankings[arm], names)
    checks['clean_vector_max_abs'] = None
    checks['clean_cosine'] = None
    if len(ids) == 1:
        checks['clean_vector_max_abs'] = {arm: float(np.max(np.abs(new[arm].astype(np.float64)-arrays[baseline])))
                                         for arm, baseline in MATCHED.items()}
        checks['clean_cosine'] = {arm: baseline_reproduction(rankings[baseline], rankings[arm], names)
                                  for arm, baseline in MATCHED.items()}
    checks['status'] = 'passed' if (all(e <= VECTOR_TOL for e in errors.values()) and
        all(c['status'] == 'passed' for c in checks['baseline_cosine'].values()) and
        all(e <= VECTOR_TOL for e in (checks['clean_vector_max_abs'] or {}).values()) and
        all(c['status'] == 'passed' for c in (checks['clean_cosine'] or {}).values())) else 'failed'
    evaluations = {arm: evaluate(query, rankings[arm]) for arm in ARMS}
    return {'query': query, 'checks': checks, 'new_rankings': {arm: rankings[arm] for arm in MATCHED},
        'baseline_top10': {arm: rankings[arm][:10] for arm in ('whole', 'isolated_mixture')},
        'evaluations': evaluations,
        'pairs': {arm: {'matching_baseline': baseline, **paired(evaluations[baseline], evaluations[arm])}
                  for arm, baseline in MATCHED.items()}}, new


def strata(component):
    return {'by_fraction': component['fraction_bin'], 'by_component_count': str(component['component_count']),
            'by_position': component['position'], 'by_ordinal_position': str(component['ordinal_position']),
            'by_target': component['accession']}


def rate(hits, n):
    return {'hits': hits, 'denominator': n, 'recall': hits/n if n else None}


class Totals:
    def __init__(self, targets, paired_mode=False):
        self.paired_mode = paired_mode; self.nq = 0; self.nc = 0
        self.component = defaultdict(Counter); self.query = defaultdict(Counter)
        self.groups = defaultdict(lambda: defaultdict(Counter))
        self.buckets = {key: {name: defaultdict(Counter) for name in values} for key, values in {
            'by_fraction': FRACTIONS, 'by_component_count': ('1','2','3','4'),
            'by_position': ('single','N','internal','C'), 'by_ordinal_position': ('1','2','3','4'),
            'by_target': targets}.items()}
        self.query_count = {str(k): defaultdict(Counter) for k in (1,2,3,4)}

    def add(self, query, evaluation):
        self.nq += 1; self.nc += len(evaluation['components'])
        for c in evaluation['components']:
            for k in KS:
                value = c['outcomes'][str(k)] if self.paired_mode else ('hit' if c['hits'][str(k)] else 'miss')
                self.component[k][value] += 1; self.groups[query['group_id']][k][value] += 1
                for key, name in strata(c).items(): self.buckets[key][name][k][value] += 1
        for k in KS:
            for kind in ('any', 'all'):
                value = evaluation['query'][str(k)][kind]
                if not self.paired_mode: value = 'hit' if value else 'miss'
                self.query[(k, kind)][value] += 1
                self.query_count[str(query['component_count'])][(k, kind)][value] += 1

    def counts(self, counter):
        n = sum(counter.values())
        if not self.paired_mode: return rate(counter['hit'], n)
        outcomes = {key: counter[key] for key in ('both_hit', 'gain', 'loss', 'neither')}
        return {'denominator': n, **outcomes, 'net_gain': outcomes['gain']-outcomes['loss'],
                'net_gain_rate': (outcomes['gain']-outcomes['loss'])/n if n else None}

    def result(self):
        require(all(sum(self.component[k].values()) == self.nc and
                    all(sum(self.query[(k, kind)].values()) == self.nq for kind in ('all','any')) for k in KS),
                'Aggregate denominators differ')
        for dimension in self.buckets.values():
            require(all(sum(sum(v[k].values()) for v in dimension.values()) == self.nc for k in KS), 'Stratum coverage differs')
        result = {'queries': self.nq, 'components': self.nc,
            'by_k': {str(k): {'component': self.counts(self.component[k]),
                             **{kind: self.counts(self.query[(k, kind)]) for kind in ('any','all')}} for k in KS},
            'component_strata': {key: {name: {str(k): self.counts(counts[k]) for k in KS}
                                       for name, counts in sorted(values.items())} for key, values in self.buckets.items()},
            'query_by_component_count': {name: {str(k): {kind: self.counts(counts[(k, kind)]) for kind in ('any','all')}
                                               for k in KS} for name, counts in self.query_count.items()}}
        if not self.paired_mode:
            for k in KS:
                for name, groups in [('target_macro', self.buckets['by_target']), ('group_macro', self.groups)]:
                    supported = [counter[k]['hit']/sum(counter[k].values()) for counter in groups.values() if sum(counter[k].values())]
                    result['by_k'][str(k)][name] = {'support': len(supported), 'recall': math.fsum(supported)/len(supported) if supported else None}
        return result


def persist(store, key, record, arrays=None):
    import numpy as np
    old = store.read(key, arrays=arrays is not None)
    if old is None: store.publish(key, record, arrays)
    elif arrays is None: require(old == record, 'Resumed record differs: '+key)
    else:
        require(old[0] == record and set(old[1]) == set(arrays) and
                all(old[1][a].dtype == arrays[a].dtype and np.array_equal(old[1][a], arrays[a]) for a in arrays),
                'Resumed record/vector differs: '+key)


def validate_source(snapshot, pin, suite):
    identity = snapshot.identity; ids = {q['query_id'] for q in suite['queries']}; names = sorted(suite['sequences'])
    require(snapshot.identity_hash == pin['identity_sha256'] and len(snapshot.units) == 1252, 'Frozen source identity/coverage differs')
    require(identity['mode'] == 'model' and len(identity['query_ids']) == len(set(identity['query_ids'])) == 1146 and
            set(identity['query_ids']) == ids and identity['catalogue_accessions'] == names and len(names) == 105,
            'Exact full common panel/catalogue required')
    require(identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'] and
            identity['catalogue_sequence_sha256'] == {a: digest(seq.encode()) for a, seq in suite['sequences'].items()}, 'Canonical source sequences differ')
    config = identity['config']
    require(config['model_id'] == pin['model_id'] and config['device'] == 'cpu' and config['dtype'] == 'float32' and
            config['layer'] == -1, 'Full-depth CPUFP32 representation differs')
    for relative, expected in identity['code_hashes'].items():
        require(sha(audit.safe_file(ROOT, relative)) == expected, 'Source scientific dependency changed')


def run(args):
    require(args.execute_saved_analysis, 'Explicit lead saved-analysis slot required')
    session = args.session.resolve(); require(session.is_relative_to(ROOT / 'results'), 'Session escaped results')
    design_path = session / 'inputs/equal_component_weight_design_v1.json'
    require(sha(__file__) == args.expected_code_sha256 and sha(design_path) == args.expected_design_sha256 == DESIGN_SHA,
            'Frozen code/design invocation differs')
    code = code_hashes()
    require(digest(canonical(code)) == args.expected_dependencies_sha256, 'Frozen consumer dependency dictionary differs')
    design = read(design_path)
    require([p['label'] for p in design['sources']] == list(MODELS), 'Complete five-model roster required')
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    import numpy as np
    started = time.monotonic(); deadline = Deadline(900, session)
    suite = load_suite(ROOT, session / 'inputs')
    queries = [{k: v for k, v in q.items() if k != 'sequence'} for q in suite['queries']]
    require(suite['verification']['challenge_manifest_sha256'] == design['challenge_sha256'] and
            len(queries) == 1146 and dict(Counter(q['suite'] for q in queries)) == design['population']['suite_queries'] and
            sum(q['component_count'] for q in queries if q['component_count'] > 1) == 3520, 'Frozen query support changed')
    snapshots = {}
    for pin in design['sources']:
        snapshot = audit.Snapshot(session / pin['job'], session, deadline, pin['completion_sha256'])
        validate_source(snapshot, pin, suite)
        snapshot.entries = read(snapshot.completion_path)['units']; snapshots[pin['label']] = snapshot
    identity = {'mode': 'equal_component_weight_controls_v1', 'design_sha256': DESIGN_SHA, 'design': design,
                'sources': {label: s.evidence() for label, s in snapshots.items()}, 'code_hashes': code,
                'numeric_runtime': {'numpy': np.__version__, 'python': sys.version, 'threads': 1},
                'query_ids': sorted(q['query_id'] for q in queries), 'catalogue_accessions': sorted(suite['sequences'])}
    store = UnitStore(session / 'analysis' / (safe_name(args.name)+'-'+digest(canonical(identity))[:12]), identity)
    expected = []; model_summaries = {}
    for label, snapshot in snapshots.items():
        refs = {}
        for a in sorted(suite['sequences']):
            deadline.check(); source = snapshot.record('catalogue-'+a); arrays = read_arrays(snapshot, 'catalogue-'+a)
            require(source['accession'] == a and source['sequence_sha256'] == digest(suite['sequences'][a].encode()) and
                    source['length'] == len(suite['sequences'][a]) and set(arrays) == {'mean'}, 'Reference vector metadata differs')
            v = arrays['mean']; require(v.ndim == 1 and v.size and v.dtype == np.float32 and np.isfinite(v).all(), 'Invalid reference mean')
            refs[a] = v
        require(len({v.shape for v in refs.values()}) == 1, 'Reference widths differ')
        totals = {(scope, arm): Totals(suite['manifest']['targets']) for scope in ('clean','combined') for arm in ARMS}
        pairs = {(scope, arm): Totals(suite['manifest']['targets'], True) for scope in ('clean','combined') for arm in MATCHED}
        failures = []; maxima = {'whole': 0., 'isolated_mixture': 0.}
        for query in sorted(queries, key=lambda q: q['query_id']):
            deadline.check(); qkey = 'query-'+query['query_id']
            record, arrays = derive(query, snapshot.record(qkey), read_arrays(snapshot, qkey), refs)
            record['model_label'] = label; record['source_query_descriptor_sha256'] = snapshot.entries[qkey]['descriptor_sha256']
            record['source_identity_sha256'] = snapshot.identity_hash
            key = label+'--'+query['query_id']; persist(store, key, record, arrays); expected.append(key)
            if record['checks']['status'] != 'passed': failures.append(query['query_id'])
            for arm, value in record['checks']['weighted_vector_max_abs'].items(): maxima[arm] = max(maxima[arm], value)
            scope = 'clean' if query['component_count'] == 1 else 'combined'
            for arm in ARMS: totals[(scope, arm)].add(query, record['evaluations'][arm])
            for arm in MATCHED: pairs[(scope, arm)].add(query, record['pairs'][arm])
        summary = {'status': 'failed_sanity' if failures else 'passed', 'failed_query_ids': failures,
            'weighted_reconstruction_max_abs': maxima,
            'arms': {scope: {arm: totals[(scope, arm)].result() for arm in ARMS} for scope in ('clean','combined')},
            'pairs': {scope: {arm: {'matching_baseline': MATCHED[arm], **pairs[(scope, arm)].result()}
                               for arm in MATCHED} for scope in ('clean','combined')}}
        for scope, nq, nc in (('clean',16,16),('combined',1130,3520)):
            require(all(totals[(scope, arm)].nq == nq and totals[(scope, arm)].nc == nc for arm in ARMS) and
                    all(pairs[(scope, arm)].nq == nq and pairs[(scope, arm)].nc == nc for arm in MATCHED), 'Model denominator differs')
        key = label+'--summary'; persist(store, key, summary); expected.append(key); model_summaries[label] = summary
        del refs, totals, pairs
        print(json.dumps({'model': label, 'status': summary['status'], 'queries': 1146}), flush=True)
    for snapshot in snapshots.values():
        require(audit.canonical_hash(read(snapshot.path/'identity.json')) == snapshot.identity_hash and
                sha(snapshot.completion_path) == snapshot.completion_sha256, 'Source changed during analysis')
    require(code_hashes() == code, 'Implementation changed during analysis')
    summary = {'status': 'passed' if all(v['status'] == 'passed' for v in model_summaries.values()) else 'failed_sanity',
               'models': model_summaries, 'interpretation': design['interpretation'], 'design_sha256': DESIGN_SHA,
               'query_units': 5730, 'total_expected_units': 5736, 'no_inference': True, 'budgets': list(KS)}
    deadline.check(); persist(store, 'summary', summary); expected.append('summary')
    require(len(expected) == len(set(expected)) == 5736, 'Missing weighting-control output unit')
    store.complete(expected, {'analysis_status': summary['status'], 'elapsed_seconds': time.monotonic()-started, 'peak_rss_mib': peak_mib()})
    print(json.dumps({'job': str(store.root), 'status': summary['status']}), flush=True)
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True); parser.add_argument('--name', default='equal-component-weight-v1')
    parser.add_argument('--expected-code-sha256', required=True); parser.add_argument('--expected-design-sha256', required=True)
    parser.add_argument('--expected-dependencies-sha256', required=True)
    parser.add_argument('--execute-saved-analysis', action='store_true')
    sys.exit(run(parser.parse_args()))
