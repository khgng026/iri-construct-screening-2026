"""Independent hand-checkable weighting auditor fixtures; no producer imports."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'scripts'),str(ROOT/'src')]
import numpy as np
import session_weighting_audit as audit


def query():
    return {'query_id':'toy','group_id':'A__B','components':['A','B'],'component_count':2,'length_aa':4,'suite':'pairs',
            'spans':[{'accession':'A','start0':0,'end0':1,'fraction':.25},
                     {'accession':'B','start0':1,'end0':4,'fraction':.75}]}


def scores(a=.9,b=.7,c=.8): return audit.rank_scores({'A':a,'B':b,'C':c})


class SelectionAndRankingTests(unittest.TestCase):
    def test_fixed_first_five_per_count_independent_of_input_order(self):
        population = [{'query_id':f'count{count}-{i:02d}','component_count':count} for count in (1,2,3,4) for i in range(9)]
        expected = [f'count{count}-{i:02d}' for count in (1,2,3,4) for i in range(5)]
        self.assertEqual(audit.select_queries(list(reversed(population))),expected)
        population.append(population[0])
        with self.assertRaisesRegex(ValueError,'Duplicate'): audit.select_queries(population)

    def test_incomplete_sample_stratum_is_not_replaced(self):
        population = [{'query_id':f'{count}-{i}','component_count':count} for count in (1,2,3,4) for i in range(5)]
        with self.assertRaisesRegex(ValueError,'fewer'): audit.select_queries(population[:-1])

    def test_exact_ties_use_accession_and_tie_counts(self):
        ranked = audit.rank_scores({'B':.5,'A':.5,'C':0.})
        self.assertEqual([r['accession'] for r in ranked],['A','B','C'])
        self.assertEqual([r['exact_score_tie_count'] for r in ranked],[2,2,1])
        audit.ranking_map(ranked,['A','B','C'])
        ranked[0]['exact_score_tie_count'] = 1
        with self.assertRaisesRegex(ValueError,'ties'): audit.ranking_map(ranked,['A','B','C'])

    def test_near_tie_inversions_retained_but_material_inversions_fail(self):
        saved = audit.rank_scores({'A':.5,'B':.5-1e-13})
        new = audit.rank_scores({'A':.5-1e-13,'B':.5})
        check = audit.compare_rankings(saved,new,['A','B'])
        self.assertTrue(check['passed']); self.assertEqual(len(check['inversions']),1)
        self.assertEqual(saved[0]['accession'],'A'); self.assertEqual(new[0]['accession'],'B')
        self.assertFalse(audit.compare_rankings(saved,audit.rank_scores({'A':.4,'B':.6}),['A','B'])['passed'])


class ScalarTests(unittest.TestCase):
    def test_raw_mean_does_not_normalize_components(self):
        raw = {'A':np.array([10,0],np.float32),'B':np.array([0,1],np.float32)}
        stored,unrounded = audit.scalar_mean(raw)
        np.testing.assert_array_equal(stored,[5,.5]); self.assertEqual(unrounded,[5,.5])
        weighted,_ = audit.scalar_mean(raw,{'A':1,'B':3}); np.testing.assert_array_equal(weighted,[2.5,.75])

    def test_scalar_fsum_and_fp32_cast(self):
        raw = {'A':np.array([1e8],np.float32),'B':np.array([1],np.float32),'C':np.array([-1e8],np.float32)}
        stored,unrounded = audit.scalar_mean(raw)
        self.assertEqual(unrounded,[1/3]); self.assertEqual(stored[0],np.float32(1/3))
        self.assertEqual(stored.dtype,np.float32)

    def test_independent_cosine_matches_hand_values_and_rejects_zero_norm(self):
        refs = {'A':np.array([1.,0.],np.float32),'B':np.array([0.,1.],np.float32),'C':np.array([-1.,0.],np.float32)}
        ranked = audit.scalar_cosines(np.array([1.,0.],np.float32),refs)
        self.assertEqual([(r['accession'],r['score']) for r in ranked],[('A',1.),('B',0.),('C',-1.)])
        with self.assertRaisesRegex(ValueError,'Zero'): audit.scalar_cosines(np.zeros(2,np.float32),refs)


class OutcomeTests(unittest.TestCase):
    def test_exact_total_budget_metrics_and_original_geometry(self):
        result = audit.evaluation(query(),scores(),['A','B','C'])
        self.assertEqual(result['query']['1'],{'any':1,'all':0})
        self.assertEqual(result['query']['5'],{'any':1,'all':1})
        self.assertEqual([c['ordinal_position'] for c in result['components']],[1,2])
        self.assertEqual([c['fraction'] for c in result['components']],[.25,.75])
        self.assertEqual([c['fraction_bin'] for c in result['components']],['0.25_to_lt_0.40','ge_0.75'])

    def test_direct_summary_counts_macro_metrics_and_null_strata(self):
        q = query(); evaluation = audit.evaluation(q,scores(),['A','B','C'])
        result = audit.summarize([(q,evaluation)],['A','B','unobserved'])
        self.assertEqual(result['queries'],1); self.assertEqual(result['components'],2)
        self.assertEqual(result['by_k']['1']['component'],{'hits':1,'denominator':2,'recall':.5})
        self.assertEqual(result['by_k']['1']['target_macro'],{'support':2,'recall':.5})
        self.assertEqual(result['by_k']['1']['group_macro'],{'support':1,'recall':.5})
        self.assertIsNone(result['component_strata']['by_target']['unobserved']['1']['recall'])
        self.assertIsNone(result['component_strata']['by_fraction']['lt_0.25']['5']['recall'])
        self.assertIsNone(result['query_by_component_count']['4']['10']['all']['recall'])

    def test_pair_summary_preserves_offsetting_gains_losses_and_all_query_denominator(self):
        q = query(); before = audit.evaluation(q,scores(),['A','B','C'])
        after = audit.evaluation(q,scores(a=.7,b=.9),['A','B','C'])
        pairs = audit.pair_evaluation(before,after)
        self.assertEqual([r['outcomes']['1'] for r in pairs['components']],['loss','gain'])
        self.assertEqual([r['rank_delta_control_minus_baseline'] for r in pairs['components']],[2,-2])
        result = audit.summarize([(q,pairs)],['A','B'],True)
        self.assertEqual(result['by_k']['1']['component'],{'denominator':2,'both_hit':0,'gain':1,'loss':1,'neither':0,'net_gain':0,'net_gain_rate':0.})
        self.assertEqual(result['by_k']['1']['all']['denominator'],1)
        self.assertEqual(result['by_k']['1']['all']['neither'],1)
        self.assertIsNone(result['query_by_component_count']['1']['1']['all']['net_gain_rate'])

    def test_structural_comparison_detects_missing_extra_denominator_and_nonfinite(self):
        expected = {'count':2,'recall':.5,'unsupported':None}
        for changed in ({'count':1,'recall':.5,'unsupported':None},{'count':2,'recall':float('nan'),'unsupported':None},
                        {'count':2,'recall':.5},expected|{'extra':1}):
            with self.assertRaises(ValueError): audit.same(changed,expected)

    def test_failed_sanity_rejected_even_if_coverage_would_complete(self):
        with self.assertRaisesRegex(ValueError,'Failed query'):
            audit.validate_checks({'status':'failed'})

    def test_clean_declared_comparison_must_match_independent_score_evidence(self):
        check = {'status':'passed','maximum_absolute_score_error':0.,'changed_candidate_ranks':[],
                 'inverted_pairs':[],'saved_baseline_ranks_preserved':True}
        output = {'new_rankings':{arm:scores() for arm in audit.MATCHED},
                  'checks':{'clean_vector_max_abs':{arm:0. for arm in audit.MATCHED},
                            'clean_cosine':{arm:deepcopy(check) for arm in audit.MATCHED}}}
        baselines = {arm:scores() for arm in ('whole','isolated_mixture')}
        self.assertTrue(all(x['passed'] for x in audit.clean_checks(output,baselines,['A','B','C']).values()))
        changed = deepcopy(output); changed['checks']['clean_cosine']['equal_contextual']['status'] = 'failed'
        with self.assertRaisesRegex(ValueError,'Clean numerical check'):
            audit.clean_checks(changed,baselines,['A','B','C'])
        changed = deepcopy(output); changed['checks']['clean_cosine']['equal_isolated']['maximum_absolute_score_error'] = .01
        with self.assertRaisesRegex(ValueError,'Clean numerical check'):
            audit.clean_checks(changed,baselines,['A','B','C'])
        changed = deepcopy(output); changed['checks']['clean_cosine'].pop('equal_contextual')
        with self.assertRaisesRegex(ValueError,'coverage'):
            audit.clean_checks(changed,baselines,['A','B','C'])


if __name__ == '__main__': unittest.main()
