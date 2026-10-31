#!/usr/bin/env python3
"""Resumable catalogue/mask analysis of completed immutable retrieval jobs.

No sequence inference or alignment is performed. Allocate the lead's analysis
slot before running: complete source records are read and per-query derived
metrics are saved, with bounded-memory streaming summaries.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
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
from iri_constructs.session_benchmark import (
    KS, catalogue_masks, evaluate_oracle_query, evaluate_query, filter_ranking, load_suite,
)
from session_independent_audit import verify_unit, validate_ranking


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def latest_completion(job):
    candidates = [(read(p)['completed_utc'], p) for p in job.glob('completion-*.json')]
    require(candidates, 'Source job has no completion descriptor')
    return max(candidates)[1]


def source_arms(identity, record):
    """Expose native full-catalogue arms without mixing methods or masks."""
    mode = identity['mode']
    if mode == 'model':
        expected = {'whole', 'isolated_mixture', 'contextual_component', 'isolated_component'}
        require(set(record['rankings']) == expected, 'Model source has unexpected/missing arms')
        return {name: (record['rankings'][name], '/rankings/'+name) for name in sorted(expected)}
    if mode == 'alignment':
        require(set(record['rankings']) == {'local_alignment'}, 'Raw alignment arms changed')
        return {'local_alignment': (record['rankings']['local_alignment'], '/rankings/local_alignment')}
    if mode == 'alignment_normalization':
        require('accession105' in record['rankings'], 'Normalized alignment lacks full catalogue ranking')
        return {'normalized_alignment': (record['rankings']['accession105'], '/rankings/accession105')}
    raise ValueError('Only completed model/raw-alignment/normalized-alignment source jobs are supported')


def validate_family_freeze(labels, suite, label_path, independent_review_path):
    require(labels['version'] == 'iri-family-labels-v1' and labels['status'] == 'frozen_before_family_scoring',
            'Unknown or unfrozen family labels')
    require(labels['source_challenge_manifest_sha256'] == suite['verification']['challenge_manifest_sha256'],
            'Family labels refer to another challenge suite')
    require(labels['source_catalogue_manifest_sha256'] == suite['manifest']['source_files_sha256']['data/manifests/t05_expanded_catalogue.json'],
            'Family labels refer to another catalogue')
    audit = read(independent_review_path)
    require(audit['status'] == 'passed' and audit['family_label_sha256'] == file_hash(label_path),
            'Independent family-label audit does not verify these exact bytes')
    require(set(labels['acceptable_matches']) == set(suite['manifest']['targets']), 'Family target coverage mismatch')
    for target, accessions in labels['acceptable_matches'].items():
        require(target in accessions and len(accessions) == len(set(accessions))
                and set(accessions) <= set(suite['sequences']), 'Invalid acceptable-family accession set')
    for accession, evidence in labels['source_evidence'].items():
        record = suite['catalogue'][accession]['metadata']
        require(evidence['raw_sha256'] == record['checksum_sha256']
                and evidence['sequence_sha256'] == record['sequence_sha256'], 'Family evidence/source hashes differ')
    for policy, value in labels['removal_policies'].items():
        require(policy in ('exact_removed', 'near90_removed'), 'Unexpected removal policy')
        require(audit['masks_and_denominators_recomputed'][policy] == value['denominators'],
                'Frozen family denominators disagree with independent review')
        for query in suite['queries']:
            rebuilt = sorted(set().union(*(set(value['removed_accessions_by_target'][a]) for a in query['components'])))
            require(rebuilt == value['common_removed_accessions_by_group'][query['group_id']],
                    'Common family mask does not remove every component consistently')
    return {'label_sha256': file_hash(label_path), 'independent_review_sha256': file_hash(independent_review_path),
            'independent_review_path': str(independent_review_path)}


class ComponentAccumulator:
    def __init__(self):
        self.n = self.scorable = 0
        self.rank_sum = self.rr_sum = 0.
        self.hits = {k: 0 for k in KS}
        self.parents, self.groups = set(), set()
        self.targets = {}

    def add(self, evaluation, component, kind):
        self.n += 1
        self.parents.add(evaluation['query_id']); self.groups.add(evaluation['group_id'])
        match = component[kind]
        if not match['scorable']:
            return
        self.scorable += 1
        self.rank_sum += match['rank']; self.rr_sum += match['reciprocal_rank']
        target = self.targets.setdefault(component['accession'], {'n': 0, **{k: 0 for k in KS}})
        target['n'] += 1
        for k in KS:
            self.hits[k] += match[f'hit_at_{k}']; target[k] += match[f'hit_at_{k}']

    def result(self):
        n = self.scorable
        return {'original_components': self.n, 'scorable_components': n, 'unscorable_components': self.n-n,
                'parent_queries': len(self.parents), 'distinct_groups': len(self.groups),
                'scorable_targets': len(self.targets), 'mean_rank': self.rank_sum/n if n else None,
                'mrr': self.rr_sum/n if n else None,
                **{f'hits_at_{k}': self.hits[k] for k in KS},
                **{f'recall_at_{k}': self.hits[k]/n if n else None for k in KS},
                **{f'per_target_macro_recall_at_{k}':
                   math.fsum(v[k]/v['n'] for v in self.targets.values())/len(self.targets)
                   if self.targets else None for k in KS}}


class QueryAccumulator:
    def __init__(self):
        self.n = self.scorable = 0
        self.counts = {f'{metric}_{k}': 0 for metric in ('any_component_hit', 'all_components_hit') for k in KS}

    def add(self, evaluation, kind):
        self.n += 1
        row = evaluation['query_metrics'][kind]
        if not row['all_components_scorable']:
            return
        self.scorable += 1
        for k in KS:
            for metric in ('any_component_hit', 'all_components_hit'):
                self.counts[f'{metric}_{k}'] += row[f'{metric}_at_{k}']

    def result(self):
        result = {'original_queries': self.n, 'all_components_scorable_queries': self.scorable,
                  'excluded_unscorable_queries': self.n-self.scorable}
        for k in KS:
            for metric in ('any_component_hit', 'all_components_hit'):
                value = self.counts[f'{metric}_{k}']
                result[f'{metric}_count_at_{k}'] = value
                result[f'{metric}_at_{k}'] = value/self.scorable if self.scorable else None
        return result


class SummaryAccumulator:
    """Same scientific fields as benchmark.summarize without retaining records."""
    def __init__(self):
        self.ids, self.groups = set(), set()
        self.components = 0
        self.budget = self.kinds = None
        self.data = {}

    def add(self, evaluation):
        require(evaluation['query_id'] not in self.ids, 'Duplicate query evaluation in streaming summary')
        kinds = tuple(evaluation['query_metrics'])
        if self.kinds is None:
            self.kinds, self.budget = kinds, evaluation['retrieval_budget']
            self.data = {kind: {'total': ComponentAccumulator(), 'queries': QueryAccumulator(),
                               'strata': {name: {} for name in ('by_component_count','by_fraction','by_position','by_target')},
                               'count_queries': {}} for kind in kinds}
        require(kinds == self.kinds and evaluation['retrieval_budget'] == self.budget,
                'Mixed label policies or retrieval budgets in summary')
        self.ids.add(evaluation['query_id']); self.groups.add(evaluation['group_id'])
        self.components += evaluation['component_count']
        for kind, data in self.data.items():
            data['queries'].add(evaluation, kind)
            data['count_queries'].setdefault(str(evaluation['component_count']), QueryAccumulator()).add(evaluation, kind)
            for component in evaluation['components']:
                data['total'].add(evaluation, component, kind)
                values = {'by_component_count': str(evaluation['component_count']),
                          'by_fraction': component['fraction_bin'], 'by_position': component['position'],
                          'by_target': component['accession']}
                for name, key in values.items():
                    data['strata'][name].setdefault(key, ComponentAccumulator()).add(evaluation, component, kind)

    def result(self):
        if not self.ids:
            return {'query_count': 0, 'component_count': 0, 'kinds': {}}
        result = {'query_count': len(self.ids), 'distinct_group_count': len(self.groups),
                  'component_count': self.components, 'retrieval_budget': self.budget, 'kinds': {}}
        for kind, data in self.data.items():
            value = data['total'].result(); value['query_metrics'] = data['queries'].result()
            for name, groups in data['strata'].items():
                value[name] = {key: group.result() for key, group in sorted(groups.items())}
            for key, queries in data['count_queries'].items():
                value['by_component_count'][key]['query_metrics'] = queries.result()
            result['kinds'][kind] = value
        return result


def condition_masks(suite, labels, query):
    masks = catalogue_masks(suite)
    result = {name: {'included_accessions': accessions, 'label_kind': 'exact_and_duplicate',
                     'source': 'frozen_challenge.catalogue_subsets'} for name, accessions in masks['masks'].items()}
    for policy, value in labels['removal_policies'].items():
        removed = value['common_removed_accessions_by_group'][query['group_id']]
        result['family_'+policy] = {'included_accessions': sorted(set(suite['sequences'])-set(removed)),
                                   'removed_accessions': removed, 'label_kind': 'labelled_family',
                                   'source': 'family_labels.removal_policies.'+policy,
                                   'policy': policy}
    return result


def _masked(rankings, included):
    if isinstance(rankings, dict):
        return {accession: filter_ranking(rows, included) for accession, rows in rankings.items()}
    return filter_ranking(rankings, included)


def _removal_diagnostics(raw, masked, query):
    """Nearest surviving candidates only; absent exact accuracy is undefined."""
    by_component = {}
    for accession in query['components']:
        original = raw[accession] if isinstance(raw,dict) else raw
        surviving = masked[accession] if isinstance(masked,dict) else masked
        exact = next(row for row in original if row['accession'] == accession)
        ranks = {row['accession']: row['rank'] for row in original}
        top = [{'accession': row['accession'], 'masked_rank': row['rank'],
                'source_rank': ranks[row['accession']], 'score': row['score']} for row in surviving[:10]]
        by_component[accession] = {'removed_exact_source_rank': exact['rank'], 'removed_exact_score': exact['score'],
                                  'surviving_top10': top,
                                  'best_survivor_minus_removed_exact_score': top[0]['score']-exact['score'] if top else None}
    return {'exact_accession_accuracy_defined': False,
            'interpretation': 'Exact references are absent; surviving identities/scores are diagnostics, not functional correctness.',
            'components': by_component}


def derive_record(source_identity, source_record, suite, labels):
    """Pure derivation used by runner and hand-built numerical tests."""
    query = source_record['query']; accessions = sorted(suite['sequences'])
    arms = source_arms(source_identity, source_record)
    for name, (rankings, _) in arms.items():
        if isinstance(rankings, dict):
            require(set(rankings) == set(query['components']), 'Oracle arm omits a component')
            for rows in rankings.values(): validate_ranking(rows, accessions)
        else:
            validate_ranking(rankings, accessions)
    evaluations, recipes, diagnostics = {}, {}, {}
    for condition, mask in condition_masks(suite, labels, query).items():
        included = mask['included_accessions']
        require(set(included) <= set(accessions), 'Mask refers outside source catalogue')
        evaluations[condition] = {}
        recipes[condition] = {**mask, 'included_accessions_sha256': digest(canonical(included)),
                              'ordered_ranking_recipe': 'Read immutable source ranking pointer; retain included accessions in existing score/accession order; renumber from 1; recompute exact/near ties at tolerance 1e-6.'}
        for arm, (raw, pointer) in arms.items():
            masked = _masked(raw, included)
            evaluator = evaluate_oracle_query if isinstance(masked,dict) else evaluate_query
            family = labels['acceptable_matches'] if mask['label_kind'] == 'labelled_family' else None
            evaluation = evaluator(masked, query, suite['equivalence'], acceptable_matches=family)
            if family is not None:
                expected = next(q for q in labels['removal_policies'][mask['policy']]['query_scoring_status']
                                if q['query_id'] == query['query_id'])
                metrics = evaluation['query_metrics']['labelled']
                require(metrics['scorable_components'] == expected['scorable_components']
                        and metrics['unscorable_components'] == expected['unscorable_components']
                        and metrics['all_components_scorable'] == expected['all_components_scorable'],
                        'Derived family denominator differs from prescoring freeze')
                require(evaluation['query_metrics']['exact']['scorable_components'] == 0
                        and evaluation['query_metrics']['duplicate_equivalent']['scorable_components'] == 0,
                        'Common mask left an exact target or duplicate scorable')
            evaluations[condition][arm] = evaluation
            if condition == 'family_exact_removed':
                diagnostics[arm] = _removal_diagnostics(raw, masked, query)
    return {'query': query, 'evaluations': evaluations, 'ranking_recipes': recipes,
            'source_ranking_pointers': {name: pointer for name, (_,pointer) in arms.items()},
            'exact_reference_removal_diagnostics': diagnostics,
            'full_rankings_omitted_from_derived_output': 'Immutable complete source ranking plus exact mask/recipe is authoritative and regenerates all ordered accession/score lists without re-inference.'}


def denominator_check(summary, expected):
    require(summary['query_count'] == expected['original_queries']
            and summary['component_count'] == expected['original_component_observations'],
            'Family summary original denominator mismatch')
    if not expected['original_queries']:
        require(summary['kinds'] == {}, 'Empty family summary has metrics')
        return
    labelled = summary['kinds']['labelled']; queries = labelled['query_metrics']
    require(labelled['scorable_components'] == expected['scorable_component_observations']
            and labelled['unscorable_components'] == expected['unscorable_component_observations']
            and queries['all_components_scorable_queries'] == expected['all_components_scorable_queries']
            and queries['excluded_unscorable_queries'] == expected['excluded_from_all_components_metric'],
            'Family summary surviving denominator mismatch')
    if not expected['all_components_scorable_queries']:
        require(all(queries[f'all_components_hit_at_{k}'] is None for k in KS),
                'All-components metric must be undefined when no fully scorable query exists')


def validate_derived_coverage(record,query,suite,labels,mode):
    """Do not let a resumed cache omit an arm, condition, or parent query."""
    require(record['query']=={k:v for k,v in query.items() if k!='sequence'},
            'Derived cached query metadata changed')
    expected_masks=condition_masks(suite,labels,query)
    expected_arms={'model':{'whole','isolated_mixture','contextual_component','isolated_component'},
                   'alignment':{'local_alignment'},'alignment_normalization':{'normalized_alignment'}}[mode]
    require(set(record['evaluations'])==set(expected_masks)==set(record['ranking_recipes']),
            'Derived cache has missing/extra catalogue conditions')
    require(set(record['source_ranking_pointers'])==expected_arms,'Derived cache has missing native ranking pointers')
    for condition,mask in expected_masks.items():
        require(set(record['evaluations'][condition])==expected_arms,'Derived cache has missing/extra native arms')
        recipe=record['ranking_recipes'][condition]
        require(recipe['included_accessions']==mask['included_accessions']
                and recipe['included_accessions_sha256']==digest(canonical(mask['included_accessions'])),
                'Derived cached catalogue mask differs')
        for arm,evaluation in record['evaluations'][condition].items():
            require(evaluation['query_id']==query['query_id'] and evaluation['group_id']==query['group_id']
                    and evaluation['component_count']==len(query['components'])
                    and [c['accession'] for c in evaluation['components']]==query['components']
                    and evaluation['catalogue_count']==len(mask['included_accessions']),
                    'Derived evaluation parent, component or catalogue coverage differs')
            oracle=arm in ('contextual_component','isolated_component')
            require(evaluation['retrieval_budget']==('k per known component' if oracle else 'k total per original query'),
                    'Derived cache changed retrieval budget')


def subset_denominators(labels, policy, query_ids):
    rows = [row for row in labels['removal_policies'][policy]['query_scoring_status'] if row['query_id'] in query_ids]
    require(len(rows) == len(query_ids), 'Frozen denominator panel does not cover source queries')
    def counts(values):
        return {'original_queries': len(values), 'original_component_observations': sum(x['component_count'] for x in values),
                'scorable_component_observations': sum(x['scorable_components'] for x in values),
                'unscorable_component_observations': sum(x['unscorable_components'] for x in values),
                'all_components_scorable_queries': sum(x['all_components_scorable'] for x in values),
                'excluded_from_all_components_metric': sum(not x['all_components_scorable'] for x in values)}
    return {'all': counts(rows), 'combined': counts([x for x in rows if x['component_count']>1]),
            'by_suite': {suite: counts([x for x in rows if x['suite']==suite]) for suite in sorted({x['suite'] for x in rows})}}


def run(args):
    started = time.monotonic(); session = args.session.resolve(); source = args.source_job.resolve()
    deadline = dt.datetime.fromisoformat(read(session/'launch.json')['hard_stop_utc'])
    def check_deadline():
        require(dt.datetime.now(dt.timezone.utc) < deadline, 'Session hard stop reached')
    check_deadline()
    suite = load_suite(ROOT, session/'inputs')
    label_path = args.family_labels.resolve() if args.family_labels else session/'inputs/family_labels_v1.json'
    labels = read(label_path)
    label_audit = validate_family_freeze(labels, suite, label_path, session/'audit/independent_family_label_review.json')
    source_identity = read(source/'identity.json'); source_hash = digest(canonical(source_identity))
    require(source_identity['mode'] in ('model','alignment','alignment_normalization'), 'Unsupported source job mode')
    require(source_identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'], 'Source challenge mismatch')
    accessions = sorted(suite['sequences'])
    require(source_identity['catalogue_accessions'] == accessions, 'Source catalogue mismatch')
    require(source_identity['catalogue_sequence_sha256'] == {a:digest(suite['sequences'][a].encode()) for a in accessions},
            'Source catalogue sequence hashes changed')
    query_ids = source_identity['query_ids']; queries = {q['query_id']:q for q in suite['queries']}
    require(query_ids and len(query_ids) == len(set(query_ids)) and set(query_ids) <= set(queries), 'Invalid source query coverage')
    if not args.allow_partial: require(set(query_ids) == set(queries), 'Full catalogue analysis requires completed full source suite')
    completion_path = latest_completion(source); completion = read(completion_path)
    require(completion['identity_hash'] == source_hash and completion['expected_count'] == len(completion['units']),
            'Source completion identity/count mismatch')
    expected_source = {'query-'+qid for qid in query_ids} | {'summary'}
    if source_identity['mode'] == 'model': expected_source |= {'catalogue-'+a for a in accessions}
    if source_identity['mode'] == 'alignment_normalization': expected_source.add('reference-selfscores')
    require(set(completion['units']) == expected_source, 'Source completion has omitted or unexpected units')
    source_paths = {key:verify_unit(source,key,index,source_hash) for key,index in completion['units'].items()}
    masks = catalogue_masks(suite)
    identity = {'version':1, 'mode':'catalogue_and_family_analysis', 'source_job':str(source),
                'source_identity_hash':source_hash, 'source_completion_path':str(completion_path),
                'source_completion_sha256':file_hash(completion_path), 'source_mode':source_identity['mode'],
                'source_model_configuration':source_identity.get('config'),
                'challenge_sha256':suite['verification']['challenge_manifest_sha256'], 'query_ids':query_ids,
                'catalogue_accessions':accessions, 'catalogue_sequence_sha256':source_identity['catalogue_sequence_sha256'],
                'family_labels_path':str(label_path), 'family_label_validation':label_audit,
                'catalogue_masks':masks, 'partial_source':set(query_ids)!=set(queries),
                'ranking_recipe':'Subset immutable complete raw rankings, keeping original score/accession order; renumber and recompute tie metadata. No representation, score or precision changes.',
                'code_hashes':{str(p.relative_to(ROOT)):file_hash(p) for p in [Path(__file__),ROOT/'scripts/session_independent_audit.py',
                                      ROOT/'src/iri_constructs/session_benchmark.py',ROOT/'src/iri_constructs/session_store.py']},
                'interpretation':'Catalogue size sensitivity and exploratory source-family retrieval; common masks; no pretrained holdout, whole-protein function claim or calibrated detection.'}
    jobid = args.name+'-'+digest(canonical(identity))[:12]
    store = UnitStore(session/'jobs'/jobid,identity,session/'job_ledger.jsonl')
    store.event(event='job_started',state='running',job_id=jobid,pid=os.getpid(),query_count=len(query_ids))
    expected = []
    summaries = {}
    try:
        for index,qid in enumerate(query_ids):
            check_deadline(); key='query-'+qid; expected.append(key)
            old = store.read(key)
            if old is None:
                path=source_paths[key]; raw=read(path/'record.json')
                require(raw['query']=={k:v for k,v in queries[qid].items() if k!='sequence'},'Source query metadata differs from frozen input')
                derived=derive_record(source_identity,raw,suite,labels)
                derived['source_record']={'unit_path':str(path),'record_sha256':file_hash(path/'record.json'),
                                          'descriptor_sha256':file_hash(path/'complete.json'),
                                          'source_identity_hash':source_hash,'source_completion_sha256':file_hash(completion_path)}
                store.publish(key,derived)
                # Summaries consume the committed bytes, not an unsaved object.
                old=store.read(key)
            validate_derived_coverage(old,queries[qid],suite,labels,source_identity['mode'])
            for condition,arms in old['evaluations'].items():
                for arm,evaluation in arms.items():
                    bundle=summaries.setdefault((condition,arm),{'all':SummaryAccumulator(),'combined':SummaryAccumulator(),'by_suite':{}})
                    bundle['all'].add(evaluation)
                    if evaluation['component_count']>1:bundle['combined'].add(evaluation)
                    bundle['by_suite'].setdefault(evaluation['suite'],SummaryAccumulator()).add(evaluation)
            if (index+1)%100==0 or index+1==len(query_ids):
                print(json.dumps({'job_id':jobid,'queries_completed':index+1,'queries_expected':len(query_ids)}),flush=True)
        output={}
        frozen_counts={policy:subset_denominators(labels,policy,set(query_ids)) for policy in labels['removal_policies']}
        for (condition,arm),bundle in summaries.items():
            summary={'all':bundle['all'].result(),'combined':bundle['combined'].result(),
                     'by_suite':{name:value.result() for name,value in sorted(bundle['by_suite'].items())}}
            require(summary['all']['query_count']==len(query_ids)
                    and summary['all']['component_count']==sum(len(queries[qid]['components']) for qid in query_ids),
                    'Catalogue summary omitted a source query or component')
            if condition.startswith('family_'):
                counts=frozen_counts[condition.removeprefix('family_')]
                for scope in ('all','combined'):denominator_check(summary[scope],counts[scope])
                for name,value in summary['by_suite'].items():denominator_check(value,counts['by_suite'][name])
            output.setdefault(condition,{})[arm]=summary
        report={'summary':output,'query_count':len(query_ids),'catalogue_masks':masks,
                'family_frozen_denominators_for_source_subset':frozen_counts,'all_family_denominators_match_freeze':True,
                'unique_catalogue_conditions':len(masks['masks']), 'prepared_subset_design_count':len(suite['manifest']['catalogue_subsets']),
                'full101_shared_endpoint_scored_once':True,
                'omitted_full_rankings_regeneration':'Each committed query has checksum-pinned source_record, source_ranking_pointers and exact ranking_recipes.'}
        expected.append('summary'); old=store.read('summary')
        if old is None:store.publish('summary',report)
        else:require(old==report,'Recomputed catalogue summary differs from saved output')
        complete=store.complete(expected,{'query_count':len(query_ids),'elapsed_seconds':time.monotonic()-started,
                                          'no_reinference':True,'no_new_alignments':True})
        print(json.dumps({'status':'complete','job':str(store.root),'units':complete['expected_count'],
                          'elapsed_seconds':time.monotonic()-started}),flush=True)
    except BaseException as exc:
        store.event(event='job_failed',state='failed',job_id=jobid,exception=f'{type(exc).__name__}: {exc}')
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--session',type=Path,required=True);p.add_argument('--source-job',type=Path,required=True)
    p.add_argument('--family-labels',type=Path);p.add_argument('--name',default='catalogue-family-analysis')
    p.add_argument('--allow-partial',action='store_true',help='Explicitly label completed trial/subset coverage')
    run(p.parse_args())
