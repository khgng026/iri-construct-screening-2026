#!/usr/bin/env python3
"""Separate scalar numerical and full-metric audit of saved weighting controls.

Shares checksum/source loading helpers, not producer calculations or metrics.
The author also wrote the producer; numerical independence is not authorship
independence. Fixed first-five/count selection is frozen before outcome access.
"""
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
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'scripts'), str(ROOT/'src')]
import session_extension_audit as integrity
from iri_constructs.session_store import atomic_json

require, read, sha, hash_json = integrity.require, integrity.read, integrity.sha, integrity.canonical_hash
DESIGN_SHA = 'b62aadb8f0085103984174d2c4674a9a1d8af5bf5a1eca80e66afcca24f7c3a0'
PLAN_SHA = '5b9554d0e6ed920744a31ed354553bfb76c791d09375f0a37af0e54b1b7c017a'
MODELS = ('esm2_8m','esm2_150m','esm2_650m','esmc_300m','esmc_600m')
ARMS = ('whole','isolated_mixture','equal_isolated','equal_contextual')
MATCHED = {'equal_isolated': 'isolated_mixture', 'equal_contextual': 'whole'}
KS = (1,5,10)
VECTOR_TOL, SCORE_TOL = 1e-6, 1e-12
FRACTION_BINS = ('lt_0.25','0.25_to_lt_0.40','0.40_to_lt_0.60','0.60_to_lt_0.75','ge_0.75')
PAIR_LABELS = {(0,0): 'neither', (0,1): 'gain', (1,0): 'loss', (1,1): 'both_hit'}


def code_hashes():
    return {p: sha(ROOT/p) for p in ('scripts/session_weighting_audit.py','tests/test_session_weighting_audit.py',
        'scripts/session_extension_audit.py','scripts/session_independent_audit.py','src/iri_constructs/session_store.py')}


def same(actual, expected, path='root', tolerance=1e-12):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and set(actual) == set(expected), path+': keys differ')
        for k in expected: same(actual[k], expected[k], path+'/'+str(k), tolerance)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), path+': length differs')
        for i, value in enumerate(expected): same(actual[i], value, path+'/'+str(i), tolerance)
    elif isinstance(expected, float):
        require(type(actual) in (int,float) and math.isfinite(actual) and abs(actual-expected) <= tolerance,
                path+': numerical value differs: '+repr((actual,expected)))
    else:
        require(type(actual) is type(expected) and actual == expected, path+': value/type differs: '+repr((actual,expected)))


def select_queries(queries):
    require(len(queries) == len({q['query_id'] for q in queries}), 'Duplicate sample population ID')
    result = []
    for count in (1,2,3,4):
        candidates = sorted(q['query_id'] for q in queries if q['component_count'] == count)
        require(len(candidates) >= 5, 'Sample stratum has fewer than five queries')
        result.extend(candidates[:5])
    require(len(result) == len(set(result)) == 20, 'Frozen sample is not20 unique queries')
    return result


def rank_scores(scores):
    require(scores and all(math.isfinite(v) and abs(v) <= 1+SCORE_TOL for v in scores.values()), 'Invalid cosine values')
    counts = Counter(scores.values())
    return [{'accession': a, 'score': scores[a], 'rank': i+1, 'exact_score_tie_count': counts[scores[a]]}
            for i,a in enumerate(sorted(scores, key=lambda a: (-scores[a],a)))]


def ranking_map(rows, names):
    require(len(rows) == len(names) and [r['rank'] for r in rows] == list(range(1,len(names)+1)), 'Incomplete ranking/ranks')
    require(len({r['accession'] for r in rows}) == len(names) and {r['accession'] for r in rows} == set(names), 'Ranking roster differs')
    scores = {r['accession']: r['score'] for r in rows}; rebuilt = rank_scores(scores)
    require([r['accession'] for r in rebuilt] == [r['accession'] for r in rows], 'Score/accession ordering differs')
    for row in rows:
        require(row['exact_score_tie_count'] == sum(v == row['score'] for v in scores.values()), 'Exact ties differ')
        if 'cosine_similarity' in row: same(row['cosine_similarity'], row['score'], 'Cosine alias')
    return {r['accession']: r for r in rows}


def compare_rankings(saved, independent, names):
    left, right = ranking_map(saved,names), ranking_map(independent,names)
    maximum = max(abs(left[a]['score']-right[a]['score']) for a in names)
    inversions = []
    for index, a in enumerate(r['accession'] for r in saved):
        for later in saved[index+1:]:
            b = later['accession']
            if right[a]['rank'] > right[b]['rank']:
                oldgap, newgap = abs(left[a]['score']-left[b]['score']), abs(right[a]['score']-right[b]['score'])
                inversions.append({'accessions': [a,b], 'saved_score_gap': oldgap, 'scalar_score_gap': newgap,
                                   'near_tie': max(oldgap,newgap) <= 2*SCORE_TOL})
    changed = [{'accession': a, 'saved_rank': left[a]['rank'], 'scalar_rank': right[a]['rank']}
               for a in names if left[a]['rank'] != right[a]['rank']]
    return {'passed': maximum <= SCORE_TOL and all(x['near_tie'] for x in inversions),
            'maximum_absolute_score_error': maximum, 'rank_changes': changed, 'inversions': inversions}


def fraction_bin(value):
    require(0 < value <= 1, 'Invalid original residue fraction')
    return FRACTION_BINS[sum(value >= edge for edge in (.25,.4,.6,.75))]


def evaluation(query, rows, names):
    ranked = ranking_map(rows,names); components = []
    for i, span in enumerate(query['spans']):
        target = query['components'][i]; require(target == span['accession'], 'Span target mismatch')
        rank = ranked[target]['rank']
        components.append({'accession': target, 'ordinal_position': i+1,
            'position': 'single' if len(query['components']) == 1 else ('N' if i == 0 else 'C' if i+1 == len(query['components']) else 'internal'),
            'component_count': len(query['components']), 'fraction': span['fraction'], 'fraction_bin': fraction_bin(span['fraction']),
            'length_aa': span['end0']-span['start0'], 'start0': span['start0'], 'end0': span['end0'],
            'rank': rank, 'score': ranked[target]['score'], 'hits': {str(k): int(rank <= k) for k in KS}})
    return {'components': components, 'query': {str(k): {'any': int(sum(c['rank'] <= k for c in components) > 0),
                'all': int(sum(c['rank'] <= k for c in components) == len(components))} for k in KS},
            'retrieval_budget': 'k total unique catalogue accessions per original query', 'label_rule': 'exact_accession'}


def pair_evaluation(baseline, equal):
    components = []
    for before, after in zip(baseline['components'], equal['components']):
        require(before['accession'] == after['accession'], 'Paired target mismatch')
        components.append({k: before[k] for k in ('accession','ordinal_position','position','component_count','fraction',
            'fraction_bin','length_aa','start0','end0')} | {'baseline_rank': before['rank'], 'control_rank': after['rank'],
            'rank_delta_control_minus_baseline': after['rank']-before['rank'],
            'score_delta_control_minus_baseline': after['score']-before['score'],
            'outcomes': {str(k): PAIR_LABELS[(before['hits'][str(k)],after['hits'][str(k)])] for k in KS}})
    return {'components': components, 'query': {str(k): {kind: PAIR_LABELS[(baseline['query'][str(k)][kind],equal['query'][str(k)][kind])]
                                                       for kind in ('any','all')} for k in KS}}


def tally(values, paired=False):
    values = list(values); n = len(values)
    if not paired:
        require(all(v in (0,1) for v in values), 'Nonbinary hit')
        return {'hits': sum(values), 'denominator': n, 'recall': sum(values)/n if n else None}
    counts = {key: values.count(key) for key in ('both_hit','gain','loss','neither')}
    require(sum(counts.values()) == n, 'Invalid pair outcome')
    net = counts['gain']-counts['loss']
    return {'denominator': n, **counts, 'net_gain': net, 'net_gain_rate': net/n if n else None}


def summarize(observations, targets, paired=False):
    """Batch selections/sums independent of the producer's streaming counters."""
    cells = [(q,c) for q,e in observations for c in e['components']]
    dimensions = {'by_fraction': ('fraction_bin',FRACTION_BINS), 'by_component_count': ('component_count',(1,2,3,4)),
        'by_position': ('position',('single','N','internal','C')), 'by_ordinal_position': ('ordinal_position',(1,2,3,4)),
        'by_target': ('accession',targets)}
    def hit(c,k): return c['outcomes' if paired else 'hits'][str(k)]
    result = {'queries': len(observations), 'components': len(cells), 'by_k': {},
              'component_strata': {}, 'query_by_component_count': {}}
    for k in KS:
        row = {'component': tally((hit(c,k) for _,c in cells),paired)}
        row.update({kind: tally((e['query'][str(k)][kind] for _,e in observations),paired) for kind in ('any','all')})
        if not paired:
            for name, cell_key in [('target_macro',lambda q,c:c['accession']), ('group_macro',lambda q,c:q['group_id'])]:
                groups = defaultdict(list)
                for q,c in cells: groups[cell_key(q,c)].append(hit(c,k))
                means = [math.fsum(values)/len(values) for values in groups.values()]
                row[name] = {'support': len(means), 'recall': math.fsum(means)/len(means) if means else None}
        result['by_k'][str(k)] = row
    for dimension,(field,names) in dimensions.items():
        result['component_strata'][dimension] = {str(name): {str(k): tally((hit(c,k) for _,c in cells if c[field] == name),paired)
                                                                 for k in KS} for name in names}
    for count in (1,2,3,4):
        result['query_by_component_count'][str(count)] = {str(k): {kind: tally((e['query'][str(k)][kind]
            for q,e in observations if q['component_count'] == count),paired) for kind in ('any','all')} for k in KS}
    return result


def scalar_mean(vectors, weights=None):
    import numpy as np
    ids = sorted(vectors); require(ids, 'Empty component set')
    width = len(vectors[ids[0]])
    require(all(len(v) == width and v.dtype == np.float32 and v.ndim == 1 and np.isfinite(v).all() for v in vectors.values()), 'Invalid scalar input')
    weights = {a: 1 for a in ids} if weights is None else weights
    require(set(weights) == set(ids) and all(type(w) is int and w > 0 for w in weights.values()), 'Invalid scalar weights')
    total = sum(weights.values())
    unrounded = [math.fsum(weights[a]*float(vectors[a][j]) for a in ids)/total for j in range(width)]
    return np.array(unrounded,dtype=np.float32), unrounded


def scalar_cosines(vector, references):
    q = [float(v) for v in vector]; qnorm = math.sqrt(math.fsum(v*v for v in q))
    require(qnorm > 0, 'Zero query norm')
    scores = {}
    for a, original in references.items():
        r = [float(v) for v in original]; require(len(q) == len(r), 'Scalar vector width differs')
        rnorm = math.sqrt(math.fsum(v*v for v in r)); require(rnorm > 0, 'Zero reference norm')
        scores[a] = math.fsum(x*y for x,y in zip(q,r))/(qnorm*rnorm)
    return rank_scores(scores)


def replay(query, parent_arrays, output_arrays, source_rows, derived, references):
    import numpy as np
    ids = query['components']; names = sorted(references)
    require(set(output_arrays) == set(MATCHED), 'Output vectors omitted or added')
    expected_keys = {'whole','isolated_mixture'} | {kind+'_'+a for kind in ('isolated','contextual') for a in ids}
    require(set(parent_arrays) == expected_keys, 'Parent vector schema differs')
    width = len(references[names[0]])
    require(all(v.shape == (width,) and v.dtype == np.float32 and np.isfinite(v).all()
                for v in [*parent_arrays.values(),*output_arrays.values()]), 'Saved vector shape/dtype differs')
    require(all(np.array_equal(parent_arrays['isolated_'+a],references[a]) for a in ids), 'Isolated reference differs')
    weights = {s['accession']: s['end0']-s['start0'] for s in query['spans']}
    iso = {a: parent_arrays['isolated_'+a] for a in ids}; ctx = {a: parent_arrays['contextual_'+a] for a in ids}
    new = {'equal_isolated': scalar_mean(iso)[0], 'equal_contextual': scalar_mean(ctx)[0]}
    wiso = scalar_mean(iso,weights)[0]; wctx = scalar_mean(ctx,weights)[1]
    errors = {'isolated_mixture': max(abs(float(x)-float(y)) for x,y in zip(wiso,parent_arrays['isolated_mixture'])),
              'whole': max(abs(x-float(y)) for x,y in zip(wctx,parent_arrays['whole']))}
    vector_checks = {arm: {'exact_float32_equal': bool(np.array_equal(v,output_arrays[arm])),
                          'max_abs': max(abs(float(x)-float(y)) for x,y in zip(v,output_arrays[arm]))}
                     for arm,v in new.items()}
    replay_vectors = {'whole': parent_arrays['whole'], 'isolated_mixture': parent_arrays['isolated_mixture'], **new}
    rankings, checks = {}, {}
    for arm, vector in replay_vectors.items():
        ranks = scalar_cosines(vector,references); rankings[arm] = ranks
        checks[arm] = compare_rankings(source_rows[arm] if arm in source_rows else derived['new_rankings'][arm],ranks,names)
    return {'query_id': query['query_id'], 'passed': all(v <= VECTOR_TOL for v in errors.values()) and
                all(v['exact_float32_equal'] for v in vector_checks.values()) and all(v['passed'] for v in checks.values()),
            'weighted_vector_max_abs': errors, 'equal_vector_checks': vector_checks, 'ranking_checks': checks,
            'scalar_rankings': rankings, 'candidate_score_comparisons': 4*len(names)}


def peak_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform == 'darwin' else 1024)


class Deadline(integrity.Deadline):
    def check(self):
        super().check(); require(peak_mib() <= 512, 'Audit exceeded512 MiB')


class DerivedSnapshot:
    def __init__(self,path,session,completion_hash,expected,deadline):
        self.path = path.resolve(); require(self.path.is_relative_to(session/'analysis'), 'Derived source escaped analysis')
        self.identity = read(self.path/'identity.json'); self.identity_hash = hash_json(self.identity)
        require(self.path.name.endswith('-'+self.identity_hash[:12]), 'Derived job identity mismatch')
        matches = [p for p in self.path.glob('completion-*.json') if sha(p) == completion_hash]
        require(len(matches) == 1, 'Exact completed derived snapshot missing')
        self.completion_path = matches[0]; self.completion_sha256 = completion_hash
        completion = read(matches[0]); self.entries = completion['units']; self.units = {}
        require(completion['identity_hash'] == self.identity_hash and completion['expected_count'] == len(expected) and
                set(self.entries) == expected and completion['metadata']['analysis_status'] == 'passed', 'Incomplete or failed derived snapshot')
        for key,entry in self.entries.items():
            deadline.check(); self.units[key] = integrity.verify_unit(self.path,key,entry,self.identity_hash)
    def record(self,key):
        return read(integrity.verify_unit(self.path,key,self.entries[key],self.identity_hash)/'record.json')


def arrays(snapshot,key):
    import numpy as np
    unit = integrity.verify_unit(snapshot.path,key,snapshot.entries[key],snapshot.identity_hash)
    with np.load(unit/'arrays.npz',allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}


def validate_checks(checks):
    require(checks['status'] == 'passed', 'Failed query sanity cannot be accepted')
    require(set(checks['weighted_vector_max_abs']) == {'whole','isolated_mixture'} and
            all(math.isfinite(v) and 0 <= v <= VECTOR_TOL for v in checks['weighted_vector_max_abs'].values()), 'Weighted sanity differs')
    for arm in ('whole','isolated_mixture'):
        row = checks['baseline_cosine'][arm]
        require(row['status'] == 'passed' and math.isfinite(row['maximum_absolute_score_error']) and
                0 <= row['maximum_absolute_score_error'] <= SCORE_TOL and row['saved_baseline_ranks_preserved'] is True and
                all(p['within_fixed_tolerance'] and max(p['saved_gap'],p['recomputed_gap']) <= 2*SCORE_TOL for p in row['inverted_pairs']),
                'Saved baseline numerical sanity failed')


def clean_checks(output,baselines,names):
    """Reconstruct the new acceptance gate from saved scores, not its flags."""
    result = {}
    require(set(output['checks']['clean_cosine']) == set(MATCHED) and
            set(output['checks']['clean_vector_max_abs']) == set(MATCHED), 'Clean check coverage differs')
    for arm,baseline in MATCHED.items():
        value = compare_rankings(baselines[baseline],output['new_rankings'][arm],names); result[arm] = value
        expected = {'status':'passed' if value['passed'] else 'failed',
            'maximum_absolute_score_error':value['maximum_absolute_score_error'],
            'changed_candidate_ranks':[{'accession':r['accession'],'saved_rank':r['saved_rank'],'recomputed_rank':r['scalar_rank']}
                                       for r in value['rank_changes']],
            'inverted_pairs':[{'accessions':r['accessions'],'saved_gap':r['saved_score_gap'],'recomputed_gap':r['scalar_score_gap'],
                               'within_fixed_tolerance':r['near_tie']} for r in value['inversions']],
            'saved_baseline_ranks_preserved':True}
        same(output['checks']['clean_cosine'][arm],expected,'Clean numerical check')
        error = output['checks']['clean_vector_max_abs'][arm]
        require(math.isfinite(error) and 0 <= error <= VECTOR_TOL,'Clean vector tolerance failed')
    return result


def audit_run(args,report):
    import numpy as np
    session = args.session.resolve(); require(session.is_relative_to(ROOT/'results'), 'Session escaped results')
    design_path = session/'inputs/equal_component_weight_design_v1.json'; plan_path = session/'audit/independent_weighting_audit_plan_v1.json'
    require(sha(design_path) == DESIGN_SHA and sha(plan_path) == PLAN_SHA and sha(__file__) == args.expected_code_sha256,
            'Frozen design/audit plan/code changed')
    codes = code_hashes(); require(hash_json(codes) == args.expected_dependencies_sha256, 'Audit dependencies changed')
    report.update(code_hashes=codes, design_sha256=DESIGN_SHA, plan_sha256=PLAN_SHA, plan=read(plan_path))
    deadline = Deadline(900,session); ctx = integrity.context(session,deadline)
    queries = ctx['queries']; names = sorted(ctx['hashes']); design = read(design_path)
    selected = select_queries(list(queries.values())); report['selected_query_ids'] = selected
    report['selection_sha256'] = hash_json(selected)
    require([p['label'] for p in design['sources']] == list(MODELS) and len(queries) == 1146 and len(names) == 105,
            'Exact five-model/panel/catalogue required')
    expected = {label+'--'+qid for label in MODELS for qid in queries} | {label+'--summary' for label in MODELS} | {'summary'}
    derived = DerivedSnapshot(args.job,session,args.completion_sha256,expected,deadline); identity = derived.identity
    require(identity['mode'] == 'equal_component_weight_controls_v1' and identity['design_sha256'] == DESIGN_SHA and
            identity['design'] == design and identity['query_ids'] == sorted(queries) and identity['catalogue_accessions'] == names,
            'Derived identity/design/panel differs')
    require(set(identity['sources']) == set(MODELS),'Derived source roster differs')
    require(identity['code_hashes']['scripts/session_weighting_controls.py'] == args.expected_producer_sha256 and
            hash_json(identity['code_hashes']) == args.expected_producer_dependencies_sha256, 'Producer source freeze differs')
    for name,expected_hash in identity['code_hashes'].items():
        require(sha(integrity.safe_file(ROOT,name)) == expected_hash, 'Producer dependency changed')
    report['derived'] = {'job': str(derived.path), 'identity_sha256': derived.identity_hash, 'completion_sha256': derived.completion_sha256,
                         'verified_units': len(expected)}
    report['models'] = {}; final_model_summaries = {}
    for pin in design['sources']:
        label = pin['label']; model_report = {'status':'running','replays':[],'query_units_checked':0}; report['models'][label] = model_report
        parent = integrity.Snapshot(session/pin['job'],session,deadline,pin['completion_sha256'])
        parent.entries = read(parent.completion_path)['units']; pident = parent.identity
        require(parent.identity_hash == pin['identity_sha256'] and set(pident['query_ids']) == set(queries) and
                len(pident['query_ids']) == 1146 and pident['catalogue_accessions'] == names and
                pident['catalogue_sequence_sha256'] == ctx['hashes'] and pident['challenge_sha256'] == design['challenge_sha256'], 'Parent identity/panel differs')
        same(identity['sources'][label],parent.evidence(),'Parent source evidence')
        config = pident['config']; require(config['model_id'] == pin['model_id'] and config['device'] == 'cpu' and
            config['dtype'] == 'float32' and config['layer'] == -1, 'Parent full-depth CPUFP32 representation differs')
        references = {}
        for a in names:
            deadline.check(); record = parent.record('catalogue-'+a); data = arrays(parent,'catalogue-'+a)
            require(record['accession'] == a and record['sequence_sha256'] == ctx['hashes'][a] and record['length'] == ctx['lengths'][a] and
                    set(data) == {'mean'}, 'Reference metadata differs')
            references[a] = data['mean']
            require(references[a].ndim == 1 and references[a].size > 0 and references[a].dtype == np.float32 and
                    np.isfinite(references[a]).all(), 'Invalid reference vector')
        require(len({v.shape for v in references.values()}) == 1, 'Reference dimensions differ')
        observations = {arm:[] for arm in ARMS}; paired_observations = {arm:[] for arm in MATCHED}
        maxima = {'whole':0.,'isolated_mixture':0.}; model_report['clean_ranking_checks'] = []
        for qid,query in sorted(queries.items()):
            deadline.check(); key = label+'--'+qid; output = derived.record(key); original = parent.record('query-'+qid)
            same(output['query'],query,'Derived canonical query'); same(original['query'],query,'Parent canonical query')
            require(output['model_label'] == label and output['source_identity_sha256'] == parent.identity_hash and
                    output['source_query_descriptor_sha256'] == parent.entries['query-'+qid]['descriptor_sha256'], 'Query source binding differs')
            validate_checks(output['checks'])
            require(set(output['new_rankings']) == set(MATCHED) and set(output['evaluations']) == set(ARMS) and
                    set(output['pairs']) == set(MATCHED), 'Output arm coverage differs')
            calculated = {}
            baselines = {arm: original['rankings'][arm] for arm in ('whole','isolated_mixture')}
            for arm in ARMS:
                ranked = baselines[arm] if arm in baselines else output['new_rankings'][arm]
                calculated[arm] = evaluation(query,ranked,names)
                same(output['evaluations'][arm],calculated[arm],label+'/'+qid+'/'+arm)
                observations[arm].append((query,calculated[arm]))
            same(output['baseline_top10'],{arm:rows[:10] for arm,rows in baselines.items()},'Preserved baseline returns')
            for arm,baseline in MATCHED.items():
                rebuilt = pair_evaluation(calculated[baseline],calculated[arm])
                same(output['pairs'][arm],{'matching_baseline':baseline,**rebuilt},'Paired outcomes')
                paired_observations[arm].append((query,rebuilt))
            if query['component_count'] == 1:
                clean = {'query_id':qid,'arms':clean_checks(output,baselines,names)}
                model_report['clean_ranking_checks'].append(clean)
                require(all(c['passed'] for c in clean['arms'].values()), 'Clean equal-control score/ranking differs')
            else:
                require(output['checks']['clean_cosine'] is None and output['checks']['clean_vector_max_abs'] is None,
                        'Combined query has non-null clean checks')
            for arm,value in output['checks']['weighted_vector_max_abs'].items(): maxima[arm] = max(maxima[arm],value)
            if qid in selected:
                sample = replay(query,arrays(parent,'query-'+qid),arrays(derived,key),baselines,output,references)
                sample['source_query_descriptor_sha256'] = parent.entries['query-'+qid]['descriptor_sha256']
                sample['derived_query_descriptor_sha256'] = derived.entries[key]['descriptor_sha256']
                model_report['replays'].append(sample)
                same(output['checks']['weighted_vector_max_abs'],sample['weighted_vector_max_abs'],'Scalar weighted reconstruction')
                require(sample['passed'], 'Independent scalar replay failed')
            model_report['query_units_checked'] += 1
        require({s['query_id'] for s in model_report['replays']} == set(selected) and len(model_report['replays']) == 20,
                'Frozen replay sample coverage differs')
        model_summary = {'status':'passed','failed_query_ids':[],'weighted_reconstruction_max_abs':maxima,'arms':{},'pairs':{}}
        for scope,nq,nc in (('clean',16,16),('combined',1130,3520)):
            eligible = lambda q: (q['component_count'] == 1) == (scope == 'clean')
            model_summary['arms'][scope] = {}; model_summary['pairs'][scope] = {}
            for arm in ARMS:
                selected_rows = [(q,e) for q,e in observations[arm] if eligible(q)]
                result = summarize(selected_rows,ctx['manifest']['targets'])
                require(result['queries'] == nq and result['components'] == nc, 'Frozen summary denominator differs')
                model_summary['arms'][scope][arm] = result
            for arm,baseline in MATCHED.items():
                selected_rows = [(q,e) for q,e in paired_observations[arm] if eligible(q)]
                model_summary['pairs'][scope][arm] = {'matching_baseline':baseline,**summarize(selected_rows,ctx['manifest']['targets'],True)}
        same(derived.record(label+'--summary'),model_summary,'All model summary fields')
        final_model_summaries[label] = model_summary
        require(hash_json(read(parent.path/'identity.json')) == parent.identity_hash and sha(parent.completion_path) == parent.completion_sha256,
                'Parent changed during audit')
        model_report.update(status='passed',source=parent.evidence(),all_model_summary_fields_checked=True,
                            exact_component_evaluations=3536*4,paired_component_evaluations=3536*2)
        print(json.dumps({'model':label,'audit':'passed','queries':1146,'scalar_replays':20}),flush=True)
    expected_final = {'status':'passed','models':final_model_summaries,'interpretation':design['interpretation'],
        'design_sha256':DESIGN_SHA,'query_units':5730,'total_expected_units':5736,'no_inference':True,'budgets':list(KS)}
    same(derived.record('summary'),expected_final,'Final complete summary')
    require(hash_json(read(derived.path/'identity.json')) == derived.identity_hash and sha(derived.completion_path) == derived.completion_sha256 and
            code_hashes() == codes,'Snapshot or audit implementation changed')
    deadline.check(); report.update(status='passed',all_final_summary_fields_checked=True,
        query_units_checked=5730,arm_evaluations_checked=22920,pair_evaluations_checked=11460,scalar_model_query_replays=100,
        candidate_cosine_comparisons=42000,numerical_scope='20/1146 queries per model; all saved ranking evaluations and summaries checked')


def main(args):
    require(args.execute_saved_audit, 'Explicit lead saved-audit slot required')
    session = args.session.resolve(); output = args.output.resolve()
    require(output.is_relative_to(session/'audit') and not output.exists(), 'Audit output must be new and owned')
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    started = time.monotonic(); report = {'status':'running','errors':[]}
    try:
        audit_run(args,report)
    except Exception as exc:
        report.update(status='failed'); report['errors'].append({'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()})
    finally:
        report.update(elapsed_seconds=time.monotonic()-started,peak_rss_mib=peak_mib())
        atomic_json(output,report)
    print(json.dumps({'status':report['status'],'output':str(output),'sha256':sha(output)}),flush=True)
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session',type=Path,required=True); parser.add_argument('--job',type=Path,required=True)
    parser.add_argument('--completion-sha256',required=True); parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--expected-code-sha256',required=True); parser.add_argument('--expected-dependencies-sha256',required=True)
    parser.add_argument('--expected-producer-sha256',required=True); parser.add_argument('--expected-producer-dependencies-sha256',required=True)
    parser.add_argument('--execute-saved-audit',action='store_true')
    sys.exit(main(parser.parse_args()))
