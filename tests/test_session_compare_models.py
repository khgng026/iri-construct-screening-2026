"""Small numerical fixtures for fair comparison; no real scientific scoring."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from session_compare_models import (analyse_jobs,deterministic_cases,order_movements,
                                     paired_method_comparisons,csv_bytes)
from iri_constructs.session_benchmark import rank_scores


def query(names,qid,group):
    n=len(names)
    return {'query_id':qid,'group_id':group,'suite':{1:'clean',2:'pairs',3:'triples'}[n],
            'components':names,'length_aa':10*n,
            'spans':[{'accession':a,'start0':i*10,'end0':(i+1)*10,'fraction':1/n} for i,a in enumerate(names)]}


def fixture():
    accessions=list('ABCDEFGH');equivalence={a:[a] for a in accessions}
    queries={'pair':query(['A','B'],'pair','AB'),'triple':query(['A','B','C'],'triple','ABC'),
             'clean':query(['A'],'clean','clean_A')}
    good=rank_scores(list(range(8,0,-1)),accessions)
    shifted=rank_scores([1,8,7,6,5,4,3,2],accessions)
    model_records=[]
    for q in queries.values():
        model_records.append({'query':q,'rankings':{'whole':good,'isolated_mixture':good,
             'contextual_component':{a:good for a in q['components']},
             'isolated_component':{a:good for a in q['components']}}})
    raw_records=[{'query':queries[q],'rankings':{'local_alignment':shifted}} for q in ['pair','clean']]
    jobs=[{'job_id':'model','label':'model','mode':'model','config':{},'query_ids':list(queries),'records':model_records},
          {'job_id':'raw','label':'raw','mode':'alignment','config':{},'query_ids':['pair','clean'],'records':raw_records}]
    return jobs,equivalence,queries


class CompareModelsTest(unittest.TestCase):
    def test_common_subset_and_clean_exclusion_in_headlines(self):
        jobs,eq,queries=fixture();report=analyse_jobs(jobs,eq,set(queries),queries)
        self.assertEqual(report['common_query_ids'],['clean','pair'])
        rows=[r for r in report['tables']['headlines'] if r['evaluation_population']=='common_queries']
        self.assertEqual(len(rows),2)
        for row in rows:
            self.assertEqual(row['scope'],'combined')
            self.assertEqual(row['original_components'],2)
            self.assertEqual(row['original_queries'],1)
        rates={r['job_id']:r['recall_at_5'] for r in rows}
        self.assertEqual(rates,{'model':1.,'raw':.5})
        coverage={r['job_id']:r for r in report['tables']['coverage']}
        self.assertEqual(coverage['model']['job_only_queries'],1)
        self.assertEqual(coverage['raw']['missing_from_full_suite'],1)

    def test_control_budget_labels_and_mechanism_remain_separate(self):
        jobs,eq,queries=fixture();report=analyse_jobs(jobs,eq,set(queries),queries)
        rows=[r for r in report['tables']['component_metrics'] if r['arm']=='contextual_component']
        self.assertTrue(rows)
        self.assertTrue(all(r['retrieval_budget']=='k per known component' for r in rows))
        mechanism=report['tables']['mechanism_contrasts']
        self.assertEqual({r['contrast'] for r in mechanism},{'joined_context_under_whole_pooling',
            'joined_context_with_known_component_pooling','mixing_isolated_components','mixing_contextual_components'})
        self.assertTrue(all(r['mean_rank_change']==0 for r in mechanism))

    def test_paired_comparison_uses_ranks_not_scores(self):
        jobs,eq,queries=fixture();report=analyse_jobs(jobs,eq,set(queries),queries)
        row=next(r for r in report['tables']['paired_method_comparisons'] if r['scope']=='combined')
        self.assertEqual(row['component_observations'],2)
        self.assertEqual(row['right_rank_better'],1)
        self.assertEqual(row['right_rank_worse'],1)
        self.assertEqual(row['mean_right_minus_left_rank'],3.)
        self.assertEqual(row['right_minus_left_recall_at_5'],-.5)
        self.assertFalse(any('cosine' in key for key in row))

    def test_cases_are_deterministic_and_group_capped(self):
        jobs,eq,queries=fixture();report=analyse_jobs(jobs,eq,set(queries),queries)
        discordant=[r for r in report['cases'] if r['case_type']=='between_methods']
        self.assertEqual(len(discordant),1)
        self.assertEqual(discordant[0]['left']['accession'],'A')
        self.assertEqual(abs(discordant[0]['right_minus_left_rank']),7)
        self.assertNotIn('score',discordant[0]['left'])

    def test_order_movement_is_within_same_group_and_target(self):
        base={'job_id':'j','label':'j','arm':'whole','retrieval_budget':'k total per original query',
              'group_id':'AB','accession':'A','component_count':2,'fraction':.3}
        rows=[{**base,'query_id':'AB','rank':3,'position':'N'},
              {**base,'query_id':'BA','rank':1,'position':'C'}]
        result=order_movements(rows)
        self.assertEqual(len(result),1)
        self.assertEqual(result[0]['rank_range'],2)
        self.assertEqual(result[0]['C_minus_N_mean_rank'],-2)
        self.assertEqual(result,order_movements(rows[::-1]))

    def test_missing_source_record_cannot_silently_reduce_comparison(self):
        jobs,eq,queries=fixture();jobs[0]['records'].pop()
        with self.assertRaisesRegex(ValueError,'omitted'):
            analyse_jobs(jobs,eq,set(queries),queries)

    def test_normalized_alignment_is_labelled_exploratory(self):
        jobs,eq,queries=fixture()
        normalized=deepcopy(jobs[1]);normalized.update(job_id='norm',label='normalized',mode='alignment_normalization')
        for record in normalized['records']:
            record['rankings']={'accession105':record['rankings']['local_alignment']}
        report=analyse_jobs([jobs[0],normalized],eq,set(queries),queries)
        normalized_rows=[r for r in report['tables']['headlines'] if r['job_id']=='norm']
        self.assertTrue(all(r['comparison_role']=='exploratory_normalization' for r in normalized_rows))
        self.assertTrue(all(r['contains_exploratory_method'] for r in report['tables']['paired_method_comparisons']))

    def test_csv_nested_metadata_is_stable_json(self):
        result=csv_bytes([{'job':'j','counts':{'b':2,'a':1},'missing':None}]).decode()
        self.assertIn('"{""a"":1,""b"":2}"',result)


if __name__=='__main__':unittest.main()
