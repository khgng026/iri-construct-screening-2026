"""Small hand-built retrieval/report fixtures; no biological/model computation."""
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_segment_compare import (Snapshot, expected_source_keys, validate_panel, validate_parent,
    return_volume, union_observation, union_tables, ReportBook, resource_row, summarize_diagnostics,
    summarize_volumes, matched_baseline, cached_family_rankings, summarize_baseline,
    derive_family_baseline, add_family_core_baselines)
from iri_constructs.session_benchmark import rank_scores, evaluate_query, evaluate_oracle_query
from iri_constructs.session_segmentation import fixed_windows, interval_metrics
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash, atomic_json


def query(names=('A', 'B'), qid='pair', group='AB'):
    n = len(names)
    return {'query_id': qid, 'group_id': group, 'suite': {1:'clean', 2:'pairs', 3:'triples', 4:'quadruples'}[n],
            'component_count': n, 'components': list(names), 'length_aa': 10*n,
            'spans': [{'accession': a, 'start0': i*10, 'end0': (i+1)*10, 'fraction': 1/n} for i, a in enumerate(names)]}


def ranking(order):
    return rank_scores(list(range(len(order), 0, -1)), list(order), 'cosine_similarity')


def meta(policy='original', arm='windows'):
    return {'job_id': 'segment8m', 'mode': 'automatic_segment_retrieval_v1', 'arm': arm,
            'policy': policy, 'role': 'selected_segmentation' if policy == 'original' else 'exploratory_family_transfer',
            'configuration_id': 'windows', 'segmenter': 'fixed_windows', 'label': 'fixture windows'}


def proposal(q):
    prediction = fixed_windows(q['length_aa'], 10, 5)
    encoded = sum(s['end0']-s['start0'] for s in prediction['segments'])
    return {'prediction': prediction, 'boundary_evaluation': interval_metrics(prediction['segments'], q['spans']),
            'logical_candidate_scan_residues': 0, 'logical_encoded_segment_residues': encoded,
            'logical_total_queried_residues': encoded, 'seconds': .2}


def family_baseline_fixture():
    queries = [query(), query(('A','B','C','D'),'quad','ABCD'), query(('A',),'clean','A')]
    suite = {'sequences':dict.fromkeys('ABCDEFGH','fixture'),
             'equivalence':{'A':['A','H'],'B':['B'],'C':['C'],'D':['D']}}
    labels = {'acceptable_matches':{'A':['A','E','H'],'B':['B','F'],'C':['C','G'],'D':['D']},
              'removal_policies':{}}
    for policy, target_masks, scorable in [
        ('exact_removed',{'A':['A','H'],'B':['B'],'C':['C'],'D':['D']},[2,3,1]),
        ('near90_removed',{'A':['A','E','H'],'B':['B'],'C':['C'],'D':['D']},[1,2,0])]:
        labels['removal_policies'][policy] = {'removed_accessions_by_target':target_masks,
            'common_removed_accessions_by_group':{q['group_id']:sorted(set().union(*(set(target_masks[a]) for a in q['components']))) for q in queries},
            'query_scoring_status':[{'query_id':q['query_id'],'suite':q['suite'],'group_id':q['group_id'],
                'component_count':len(q['components']),'scorable_components':n,
                'unscorable_components':len(q['components'])-n,'all_components_scorable':n==len(q['components'])}
                for q,n in zip(queries,scorable)]}
    sources = []
    for i, mode in enumerate(['model']*4+['alignment','alignment_normalization']):
        records = {}; rows = ranking('ABCDEFGH')
        for q in queries:
            arms = ({'whole':rows,'isolated_mixture':rows,
                     'contextual_component':{a:rows for a in q['components']},
                     'isolated_component':{a:rows for a in q['components']}} if mode=='model' else
                    {'local_alignment':rows} if mode=='alignment' else {'accession105':rows})
            records['query-'+q['query_id']] = {'query':q,'rankings':arms}
        identity = {'mode':mode,'query_ids':[q['query_id'] for q in queries],
                    'method':{'model_id':f'model{i}'},'config':{'model_id':f'model{i}'}}
        source = SimpleNamespace(root=Path(f'/fixture/source{i}'), identity=identity, identity_hash=f'id{i}',
            record=lambda key, records=records:records[key],
            reference=lambda key,pointer:{'unit_key':key,'json_pointer':pointer,'completion_sha256':'fixture'})
        sources.append(source)
    return sources,queries,suite,labels


class SegmentCompareTest(unittest.TestCase):
    def test_panel_requires_every_frozen_query_and_rejects_extra_segmentation(self):
        self.assertEqual(validate_panel(['clean', 'pair', 'dev'], ['pair', 'clean'])['source_queries_outside_comparison'], 1)
        with self.assertRaisesRegex(ValueError, 'missing held-out'):
            validate_panel(['clean'], ['pair', 'clean'])
        with self.assertRaisesRegex(ValueError, 'outside frozen'):
            validate_panel(['clean', 'pair', 'dev'], ['pair', 'clean'], exact=True)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            validate_panel(['pair', 'pair'], ['pair'])

    def test_snapshot_only_uses_pinned_units_and_detects_tampered_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            identity = {'mode': 'alignment', 'query_ids': ['pair']}
            store = UnitStore(root, identity)
            store.publish('query-pair', {'value': 1})
            store.publish('summary', {})
            store.complete(['query-pair', 'summary'])
            snap = Snapshot(root)
            store.publish('unrelated-loose-unit', {'value': 99})
            self.assertEqual(snap.record('query-pair'), {'value': 1})
            with self.assertRaisesRegex(ValueError, 'outside pinned'):
                snap.record('unrelated-loose-unit')
            with (snap.paths['query-pair'] / 'record.json').open('ab') as stream: stream.write(b' ')
            with self.assertRaises(ValueError): Snapshot(root)

    def test_semantically_incomplete_completion_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            store = UnitStore(root, {'mode': 'alignment', 'query_ids': ['pair']})
            store.publish('query-pair', {})
            store.complete(['query-pair'])
            with self.assertRaisesRegex(ValueError, 'exact expected unit coverage'): Snapshot(root)
        expected = expected_source_keys({'mode': 'alignment_normalization', 'query_ids': ['pair']})
        self.assertEqual(expected, {'query-pair', 'summary', 'reference-selfscores'})

    def test_method_identity_and_completed_parent_both_required(self):
        method = {'model_id': '8m', 'attention': 'eager', 'tokenizer': 'pinned'}
        parent = SimpleNamespace(identity={'mode':'model', 'method':method}, identity_hash='parent', completion_sha256='completion')
        child = SimpleNamespace(identity={'model':method, 'parent_model_identity':'parent', 'parent_completion_sha256':'completion'})
        validate_parent(child, parent, development_method=method)
        bad = deepcopy(child)
        bad.identity['model']['attention'] = 'another-backend'
        with self.assertRaisesRegex(ValueError, 'representations differ'): validate_parent(bad, parent)
        with self.assertRaisesRegex(ValueError, 'frozen development'): validate_parent(child, parent, development_method={'model_id':'8m'})
        bad = deepcopy(child); bad.identity['parent_completion_sha256'] = 'new-completion'
        with self.assertRaisesRegex(ValueError, 'snapshot differs'): validate_parent(bad, parent)

    def test_same_total_vs_segment_union_volume_and_extraneous_labels(self):
        q = query(); lists = [ranking('CABD'), ranking('CBAD')]
        total = return_volume(lists[0], q, 1)
        larger = return_volume(lists, q, 1, budget='k per predicted segment')
        self.assertEqual(total['returned_hit_occurrences'], 1)
        self.assertEqual((larger['returned_hit_occurrences'], larger['unique_returned_hits'], larger['duplicate_returned_hits']), (2, 1, 1))
        self.assertEqual(larger['outside_exact_labels'], 1)
        wide = return_volume(lists, q, 5, budget='k per predicted segment')
        self.assertEqual((wide['returned_hit_occurrences'], wide['unique_returned_hits']), (8, 4))
        self.assertEqual(wide['outside_exact_labels'], 2)

    def test_oracle_hits_remain_component_specific_while_volume_is_union(self):
        q = query(); rankings = {'A':ranking('BAC'), 'B':ranking('ABC')}
        evaluation = evaluate_oracle_query(rankings, q, {})
        self.assertEqual(evaluation['query_metrics']['exact']['component_hits_at_1'], 0)
        volume = return_volume(rankings, q, 1, budget=evaluation['retrieval_budget'])
        self.assertEqual(volume['unique_returned_hits'], 2)
        self.assertEqual(volume['returned_hit_occurrences'], 2)
        self.assertEqual(volume['retrieval_budget'], 'k per known component')

    def test_union_component_and_query_counts_strata_and_macro(self):
        q = query(); lists = [ranking('ACBD'), ranking('BCAD')]
        evaluation = evaluate_query(lists[0], q, {})
        observation = union_observation(evaluation, lists, q)
        components, queries = union_tables([observation], meta())
        row = next(r for r in queries if r['scope']=='combined' and r['k']==1 and r['stratum_dimension']=='overall')
        self.assertEqual(row['all_components_hit_count'], 1)
        self.assertEqual(row['any_component_hit_count'], 1)
        component = next(r for r in components if r['scope']=='combined' and r['k']==1 and r['stratum_dimension']=='overall')
        self.assertEqual(component['component_hits'], 2)
        self.assertEqual(component['per_target_macro_component_recall'], 1.)
        self.assertEqual({r['stratum_dimension'] for r in components}, {'overall','component_count','fraction','position','target'})
        self.assertFalse(any('rank' in key for row in components for key in row))

    def test_near90_quadruple_query_rates_null_despite_some_component_hits(self):
        q = query(('A','B','C','D'), 'quad', 'ABCD')
        labels = {'A':['E'], 'B':['F'], 'C':[], 'D':[]}
        lists = [ranking('EFG'), ranking('FEG')]
        evaluation = evaluate_query(lists[0], q, {}, acceptable_matches=labels)
        observation = union_observation(evaluation, lists, q, acceptable=labels)
        components, queries = union_tables([observation], meta('near90_removed'))
        row = next(r for r in queries if r['scope']=='suite:quadruples' and r['k']==1 and r['stratum_dimension']=='overall')
        self.assertEqual(row['all_components_scorable_queries'], 0)
        self.assertEqual(row['excluded_unscorable_queries'], 1)
        self.assertIsNone(row['all_components_hit'])
        self.assertIsNone(row['any_component_hit'])
        c = next(r for r in components if r['scope']=='combined' and r['k']==1 and r['stratum_dimension']=='overall')
        self.assertEqual((c['scorable_components'], c['unscorable_components'], c['component_hits']), (2,2,2))
        self.assertEqual(c['component_recall'], 1.)

    def test_completely_unscorable_family_union_is_not_a_zero_percent_failure(self):
        q = query(); labels = {'A':[], 'B':[]}; rows = ranking('CDE')
        observation = union_observation(evaluate_query(rows, q, {}, acceptable_matches=labels), [rows], q, acceptable=labels)
        crows, qrows = union_tables([observation], meta('near90_removed'))
        self.assertTrue(all(r['component_recall'] is None for r in crows))
        self.assertTrue(all(r['all_components_hit'] is None for r in qrows))

    def test_headlines_exclude_clean_oracles_and_keep_family_separate(self):
        book = ReportBook(); rows = ranking('ABCDE')
        for q in [query(), query(('A',), 'clean', 'A')]:
            book.add(meta(), q, evaluate_query(rows,q,{}), rows, {'source':'fixture'})
        family_labels = {'A':['C'], 'B':[]}
        q = query(); survivors = ranking('CDE')
        book.add(meta('near90_removed'), q, evaluate_query(survivors,q,{},family_labels), survivors, {}, acceptable=family_labels)
        oracle = {a:rows for a in q['components']}
        oracle_meta = {**meta(arm='isolated_component'), 'role':'oracle_control'}
        book.add(oracle_meta,q,evaluate_oracle_query(oracle,q,{}),oracle,{})
        report = book.result()
        self.assertEqual(len(report['tables']['headlines']), 1)
        self.assertEqual(report['tables']['headlines'][0]['original_components'], 2)
        self.assertEqual(report['tables']['headlines'][0]['scope'], 'combined')
        self.assertEqual(len(report['family_tables']['headlines']), 1)
        self.assertIsNone(report['family_tables']['headlines'][0]['all_components_hit_at_5'])
        self.assertTrue(all(r['policy']=='original' for r in report['tables']['component_metrics']))
        self.assertTrue(any(r['scope']=='clean' for r in report['tables']['component_metrics']))
        self.assertTrue(all(r['policy']=='original' for r in report['summaries']))
        self.assertTrue(all(r['policy']=='near90_removed' for r in report['family_summaries']))

    def test_resource_scan_plus_encoding_and_boundary_aggregation(self):
        q = query(); p = proposal(q)
        p['logical_candidate_scan_residues'] = 60
        p['logical_total_queried_residues'] += 60
        r = resource_row(q,p,{'seconds':.3,'cache_hits':2,'new_encodes':1},meta())
        self.assertEqual(r['logical_total_queried_residues'], 90)
        self.assertEqual(r['logical_encoded_segment_residues'], 30)
        book = ReportBook(); rows = ranking('ABCD')
        record = {'ranking':rows,'evaluation':evaluate_query(rows,q,{}),'seconds':.3,'new_encodes':1}
        book.add_segments(meta(),q,record,p,[rows,rows,rows],{})
        report = book.result()
        boundary = next(r for r in report['tables']['boundary_metrics'] if r['stratum_dimension']=='overall')
        self.assertEqual(boundary['mean_best_interval_iou'], 1.)
        self.assertEqual(boundary['exact_boundary_match_count'], 2)
        resource = next(r for r in report['tables']['resource_metrics'] if r['stratum_dimension']=='overall')
        self.assertEqual(resource['logical_total_queried_residues'], 90)
        p['logical_total_queried_residues'] = 60
        with self.assertRaisesRegex(ValueError, 'scanning AND encoding'): resource_row(q,p,record,meta())

    def test_volume_aggregates_count_unique_per_query_not_globally(self):
        rows = []
        for qid in ['q1','q2']:
            q = query(qid=qid)
            rows.append({**meta(), 'scope':'combined','query_id':qid,'component_count':2,
                         **return_volume(ranking('ABC'),q,1)})
        result = next(r for r in summarize_volumes(rows) if r['stratum_dimension']=='overall')
        self.assertEqual(result['unique_returned_hits'], 2)
        self.assertEqual(result['mean_unique_returned_hits'], 1.)

    def test_matched_baseline_uses_same_query_and_equal_total_budget(self):
        q = query(); whole = evaluate_query(ranking('CDEAB'),q,{})
        segment = evaluate_query(ranking('ABCDE'),q,{})
        row = matched_baseline(meta(),q,segment,whole)
        self.assertEqual(row['component_hit_difference_at_1'],1)
        self.assertEqual(row['all_hit_difference_at_5'],0)
        self.assertEqual(row['retrieval_budget'],'k total per original query')
        bad = deepcopy(whole); bad['query_id'] = 'another'
        with self.assertRaisesRegex(ValueError,'Paired baseline'): matched_baseline(meta(),q,segment,bad)

    def test_small_end_to_end_book_serialization_has_separate_budget_tables(self):
        book = ReportBook()
        for q in [query(), query(('A',), 'clean', 'A')]:
            rows = ranking('ABCDE'); p = proposal(q)
            record = {'ranking':rows,'evaluation':evaluate_query(rows,q,{}),'seconds':.1,'cache_hits':2,'new_encodes':0}
            segment_lists = [rows] * p['prediction']['segment_count']
            book.add_segments(meta(),q,record,p,segment_lists,{'completion_sha256':'fixture'})
            book.tables['matched_whole_baseline'].append(matched_baseline(meta(),q,record['evaluation'],evaluate_query(ranking('CDEAB'),q,{})))
        family_q = query(('A','B','C','D'), 'quad', 'ABCD')
        labels = {'A':['E'],'B':['F'],'C':[],'D':[]}
        rows = ranking('EFG'); p = proposal(family_q)
        record = {'ranking':rows,'evaluation':evaluate_query(rows,family_q,{},labels),'seconds':.1,'cache_counts':{'segment_cache':3}}
        book.add_segments(meta('near90_removed'),family_q,record,p,[rows]*p['prediction']['segment_count'],{},acceptable=labels)
        report = book.result()
        # Full output must be finite canonical JSON, including all null denominators.
        self.assertIn(b'null', canonical(report))
        self.assertEqual(len(report['tables']['union_per_query']),2)
        self.assertEqual(len(report['family_tables']['union_per_query']),1)
        row = next(r for r in report['tables']['matched_whole_baseline_metrics'] if r['scope']=='combined' and r['stratum_dimension']=='overall')
        self.assertEqual(row['component_recall_difference_at_1'],.5)
        null_quad = [r for r in report['family_tables']['union_query_metrics'] if r['scope']=='suite:quadruples']
        self.assertTrue(null_quad)
        self.assertTrue(all(r['all_components_hit'] is None for r in null_quad))

    def test_family_cache_exact_sequence_representation_and_descriptor(self):
        with tempfile.TemporaryDirectory() as tmp:
            session = Path(tmp); method = {'model_id':'fixture','attention':'eager'}
            identity = {'model':method,'pooling':'fixture','catalogue_sequence_sha256':{'A':'hash'},
                        'code_hashes':{'fixture':'hash'},'catalogue_accessions':list('ABC')}
            ci = {'mode':'family_transfer_segment_cache_v1', **{k:identity[k] for k in (
                'model','pooling','catalogue_sequence_sha256','code_hashes')}}
            cache_hash = digest(canonical(ci)); cache = UnitStore(session/'caches'/('family-segments-'+cache_hash[:16]),ci)
            sequence = 'XX'; sha = digest(sequence.encode())
            unit = cache.publish(sha,{'sequence_sha256':sha,'length':2,'representation_identity_hash':digest(canonical(method)),
                                     'full_catalogue_ranking':ranking('ABC')})
            evidence = {'sequence_sha256':sha,'length':2,'cache_identity_hash':cache_hash,'unit_path':str(unit),
                        'descriptor_sha256':file_hash(unit/'complete.json'),'ranking_pointer':'/full_catalogue_ranking'}
            record = {'segment_ranking_references':[evidence],'prediction':{'segments':[{'start0':0,'end0':2}]}}
            result = cached_family_rankings(record,{'sequence':sequence},identity,session,{})
            self.assertEqual(result,[ranking('ABC')])
            with self.assertRaisesRegex(ValueError,'reference differs'):
                cached_family_rankings(record,{'sequence':'XY'},identity,session,{})
            bad = deepcopy(record); bad['segment_ranking_references'][0]['descriptor_sha256'] = 'bad'
            with self.assertRaises(ValueError): cached_family_rankings(bad,{'sequence':sequence},identity,session,{})

    def test_all_six_family_primary_baselines_have_both_policies_same_panel(self):
        sources,queries,suite,labels = family_baseline_fixture()
        book = ReportBook(); paired = add_family_core_baselines(book,sources,queries,suite,labels)
        report = book.result(); headline = report['family_tables']['headlines']
        self.assertEqual(len(headline),12)
        self.assertEqual({(r['job_id'],r['policy']) for r in headline},
                         {(s.root.name,p) for s in sources for p in labels['removal_policies']})
        self.assertTrue(all(r['original_components']==6 and r['original_queries']==2 for r in headline))
        self.assertTrue(all(r['retrieval_budget']=='k total per original query' for r in headline))
        self.assertTrue(all(r['scope']=='combined' and r['role']=='family_retrieval_baseline' for r in headline))
        self.assertEqual(len(paired),4*3*2)
        normalized = [r for r in headline if r['mode']=='alignment_normalization']
        self.assertTrue(all(r['baseline_method_role']=='exploratory_normalization' for r in normalized))
        self.assertFalse(report['tables']['headlines'])
        clean = [r for r in report['family_tables']['query_metrics'] if r['scope']=='clean' and r['label_kind']=='labelled' and r['stratum_dimension']=='overall']
        self.assertEqual(len(clean),12)
        self.assertTrue(all(r['original_queries']==1 for r in clean))
        quad = [r for r in report['family_tables']['query_metrics'] if r['scope']=='suite:quadruples' and r['label_kind']=='labelled' and r['policy']=='near90_removed']
        self.assertTrue(quad)
        self.assertTrue(all(r['all_components_scorable_queries']==0 and r['all_components_hit_at_5'] is None and r['any_component_hit_at_5'] is None for r in quad))

    def test_family_baseline_common_mask_removes_all_targets_and_duplicates(self):
        sources,queries,suite,labels = family_baseline_fixture(); source=sources[0]; q=queries[0]
        derived = derive_family_baseline(source.identity,source.record('query-pair'),q,suite,labels,'exact_removed')
        self.assertEqual(derived['common_mask']['removed_accessions'],['A','B','H'])
        self.assertEqual([r['accession'] for r in derived['ranking']],list('CDEFG'))
        self.assertEqual([r['rank'] for r in derived['ranking']],list(range(1,6)))
        self.assertEqual(derived['evaluation']['query_metrics']['exact']['scorable_components'],0)
        bad=deepcopy(labels);bad['removal_policies']['exact_removed']['common_removed_accessions_by_group']['AB']=['A','H']
        with self.assertRaisesRegex(ValueError,'all components together'):
            derive_family_baseline(source.identity,source.record('query-pair'),q,suite,bad,'exact_removed')
        source.identity['query_ids']=['pair','clean']
        with self.assertRaisesRegex(ValueError,'missing held-out'):
            add_family_core_baselines(ReportBook(),[source],queries,suite,labels)

    def test_family_matched_whole_contrasts_preserve_null_and_scorable_denominators(self):
        sources,queries,suite,labels = family_baseline_fixture(); source=sources[0]
        observations=[]
        for q in queries:
            derived=derive_family_baseline(source.identity,source.record('query-'+q['query_id']),q,suite,labels,'near90_removed')
            a=derived['evaluation']; b=evaluate_query(list(derived['ranking']),q,suite['equivalence'],labels['acceptable_matches'])
            observations.append(matched_baseline(meta('near90_removed'),q,a,b,kind='labelled'))
        self.assertTrue(all(r['all_hit_difference_at_5'] is None for r in observations))
        combined=next(r for r in summarize_baseline(observations) if r['scope']=='combined' and r['stratum_dimension']=='overall')
        self.assertEqual((combined['components'],combined['scorable_components'],combined['unscorable_components']),(6,3,3))
        self.assertEqual(combined['component_recall_difference_at_5'],0.)
        self.assertIsNone(combined['all_hit_rate_difference_at_5'])
        self.assertIsNone(combined['any_hit_rate_difference_at_5'])
        clean=next(r for r in summarize_baseline(observations) if r['scope']=='clean' and r['stratum_dimension']=='overall')
        self.assertIsNone(clean['component_recall_difference_at_5'])
        bad=deepcopy(a);bad['components'][0]['labelled']['accepted_surviving']=['different']
        with self.assertRaisesRegex(ValueError,'surviving labels'):
            matched_baseline(meta('near90_removed'),q,a,bad,kind='labelled')


if __name__ == '__main__':
    unittest.main()
