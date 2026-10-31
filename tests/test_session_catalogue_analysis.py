"""Synthetic numerical/metadata checks only; no real-result analysis or inference."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_catalogue_analysis import (
    SummaryAccumulator, denominator_check, derive_record, source_arms, subset_denominators, validate_derived_coverage,
)
from iri_constructs.session_benchmark import evaluate_query, rank_scores, summarize


def query(components, qid='q', group='g'):
    return {'query_id':qid,'group_id':group,'suite':'pairs' if len(components)==2 else 'clean',
            'components':components,'spans':[{'accession':a,'start0':i*10,'end0':(i+1)*10,
                                             'fraction':1/len(components)} for i,a in enumerate(components)]}


def fixture():
    q=query(['A','B'])
    suite={'sequences':{'A':'reference-1','A2':'reference-1','B':'reference-2','C':'reference-3','D':'reference-4'},
           'equivalence':{'A':['A','A2'],'A2':['A','A2'],'B':['B'],'C':['C'],'D':['D']},
           'manifest':{'catalogue_subsets':[
               {'subset_id':'unique4_rep0','entries':[{'representative_accession':a} for a in ['A','B','C','D']]},
               {'subset_id':'unique4_rep1','entries':[{'representative_accession':a} for a in ['A','B','C','D']]}]}}
    status={'query_id':'q','group_id':'g','suite':'pairs','component_count':2,
            'scorable_components':1,'unscorable_components':1,'all_components_scorable':False}
    labels={'acceptable_matches':{'A':['A','A2','C'],'B':['B']},
            'removal_policies':{name:{'common_removed_accessions_by_group':{'g':['A','A2','B']},
                                      'query_scoring_status':[status]} for name in ['exact_removed','near90_removed']}}
    ranking=rank_scores([5,4,3,2,1],['A','A2','B','C','D'])
    record={'query':q,'rankings':{'local_alignment':ranking}}
    return suite,labels,record


class CatalogueAnalysisTest(unittest.TestCase):
    def assert_nested_close(self, a, b):
        if isinstance(a,dict):
            self.assertEqual(set(a),set(b))
            for key in a:self.assert_nested_close(a[key],b[key])
        elif isinstance(a,float):self.assertAlmostEqual(a,b,places=12)
        else:self.assertEqual(a,b)

    def test_streaming_summary_matches_native_with_unscorable_and_strata(self):
        rankings=rank_scores([5,4,3,2,1],['A','A2','B','C','D'])
        rows=[evaluate_query(rankings,query(['A','B'],'q1','g1'),{'A':['A','A2']},
                             acceptable_matches={'A':['C'],'B':[]}),
              evaluate_query(rankings,query(['B'],'q2','g2'),{},acceptable_matches={'B':['D']}),
              evaluate_query(rankings,query(['A','B'],'q3','g1'),{},acceptable_matches={'A':['C'],'B':[]})]
        accumulator=SummaryAccumulator()
        for row in rows:accumulator.add(row)
        self.assert_nested_close(accumulator.result(),summarize(rows))
        with self.assertRaises(ValueError):accumulator.add(rows[0])

    def test_family_masks_make_absent_exact_unscorable_and_keep_family_labels(self):
        suite,labels,record=fixture()
        result=derive_record({'mode':'alignment'},record,suite,labels)
        evaluation=result['evaluations']['family_exact_removed']['local_alignment']
        self.assertEqual(evaluation['query_metrics']['exact']['scorable_components'],0)
        self.assertEqual(evaluation['query_metrics']['duplicate_equivalent']['scorable_components'],0)
        self.assertEqual(evaluation['query_metrics']['labelled']['scorable_components'],1)
        self.assertEqual(evaluation['components'][0]['labelled']['rank'],1)
        self.assertFalse(evaluation['components'][1]['labelled']['scorable'])
        self.assertIsNone(evaluation['query_metrics']['labelled']['all_components_hit_at_5'])
        self.assertNotIn('rankings',result)
        self.assertEqual(result['source_ranking_pointers']['local_alignment'],'/rankings/local_alignment')
        self.assertEqual(result['ranking_recipes']['family_exact_removed']['included_accessions'],['C','D'])
        diagnostics=result['exact_reference_removal_diagnostics']['local_alignment']
        self.assertFalse(diagnostics['exact_accession_accuracy_defined'])
        self.assertEqual(diagnostics['components']['A']['surviving_top10'][0]['source_rank'],4)

    def test_duplicate_full_endpoint_is_not_scored_again(self):
        suite,labels,record=fixture()
        result=derive_record({'mode':'alignment'},record,suite,labels)
        self.assertIn('unique4_rep0',result['evaluations'])
        self.assertNotIn('unique4_rep1',result['evaluations'])

    def test_oracle_lists_keep_larger_budget_and_component_ownership(self):
        suite,labels,record=fixture();ranking=record['rankings']['local_alignment']
        record['rankings']={'whole':ranking,'isolated_mixture':ranking,
                            'isolated_component':{'A':ranking,'B':ranking},
                            'contextual_component':{'A':ranking,'B':ranking}}
        result=derive_record({'mode':'model'},record,suite,labels)
        evaluation=result['evaluations']['family_exact_removed']['contextual_component']
        self.assertEqual(evaluation['retrieval_budget'],'k per known component')
        self.assertEqual(evaluation['returned_hits_ceiling']['1'],2)
        self.assertEqual(evaluation['query_metrics']['labelled']['unscorable_components'],1)

    def test_normalized_source_uses_only_native_full_ranking(self):
        suite,labels,record=fixture();ranking=record['rankings']['local_alignment']
        record['rankings']={'accession105':ranking,'unique25_rep0':ranking[:3]}
        result=derive_record({'mode':'alignment_normalization'},record,suite,labels)
        self.assertEqual(set(result['evaluations']['accession105']),{'normalized_alignment'})
        self.assertEqual(result['source_ranking_pointers']['normalized_alignment'],'/rankings/accession105')
        with self.assertRaises(ValueError):source_arms({'mode':'unrecognized'},record)

    def test_frozen_denominator_mismatch_is_rejected(self):
        suite,labels,record=fixture()
        labels['removal_policies']['exact_removed']['query_scoring_status'][0]['scorable_components']=2
        with self.assertRaisesRegex(ValueError,'denominator'):
            derive_record({'mode':'alignment'},record,suite,labels)

    def test_zero_fully_scorable_queries_report_null_not_zero(self):
        suite,labels,record=fixture()
        result=derive_record({'mode':'alignment'},record,suite,labels)
        evaluation=result['evaluations']['family_exact_removed']['local_alignment']
        acc=SummaryAccumulator();acc.add(evaluation);summary=acc.result()
        counts=subset_denominators(labels,'exact_removed',{'q'})
        denominator_check(summary,counts['all'])
        self.assertIsNone(summary['kinds']['labelled']['query_metrics']['all_components_hit_at_5'])
        broken=deepcopy(summary);broken['kinds']['labelled']['query_metrics']['all_components_hit_at_5']=0
        with self.assertRaises(ValueError):denominator_check(broken,counts['all'])

    def test_missing_ranking_accession_cannot_silently_change_catalogue(self):
        suite,labels,record=fixture();record['rankings']['local_alignment'].pop()
        with self.assertRaises(ValueError):derive_record({'mode':'alignment'},record,suite,labels)

    def test_resumed_cache_must_cover_every_condition_arm_and_query(self):
        suite,labels,record=fixture()
        result=derive_record({'mode':'alignment'},record,suite,labels)
        validate_derived_coverage(result,record['query'],suite,labels,'alignment')
        incomplete=deepcopy(result);incomplete['evaluations'].pop('unique4_rep0')
        with self.assertRaisesRegex(ValueError,'conditions'):
            validate_derived_coverage(incomplete,record['query'],suite,labels,'alignment')
        wrong=deepcopy(result);wrong['evaluations']['accession105']['local_alignment']['query_id']='other'
        with self.assertRaisesRegex(ValueError,'coverage'):
            validate_derived_coverage(wrong,record['query'],suite,labels,'alignment')


if __name__=='__main__':
    unittest.main()
