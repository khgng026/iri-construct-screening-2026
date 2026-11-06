"""Independent hand-checkable tests of denominators, masks and design checks."""
from copy import deepcopy
import os
from pathlib import Path
import tempfile
import unittest

from iri_constructs.session_benchmark import (
    catalogue_masks, common_reference_removal_mask, cosine_ranking, evaluate_oracle_query,
    evaluate_query, filter_ranking, fraction_bin, load_suite, rank_scores, summarize,
)


def query(components, qid='query', group='group'):
    count = len(components)
    return {'query_id': qid, 'suite': {1: 'clean', 2: 'pairs', 3: 'triples', 4: 'quadruples'}[count],
            'group_id': group, 'components': components,
            'spans': [{'accession': accession, 'start0': i*10, 'end0': (i+1)*10,
                       'fraction': 1/count} for i, accession in enumerate(components)]}


class SessionMetricsTest(unittest.TestCase):
    def test_exact_ties_are_accession_sorted_and_near_ties_retained(self):
        rows = rank_scores([1, 1, 1-5e-7, -1], ['B', 'A', 'C', 'D'])
        self.assertEqual([r['accession'] for r in rows], ['A', 'B', 'C', 'D'])
        e = evaluate_query(rows, query(['B']), {'B': ['B']})
        c = e['components'][0]['exact']
        self.assertEqual((c['rank'], c['tie_best_rank'], c['tie_worst_rank']), (2, 1, 2))
        self.assertEqual(c['near_tie_count'], 3)
        self.assertEqual(c['near_tie_worst_rank'], 3)

    def test_every_component_count_has_correct_ceiling_and_all_hit(self):
        rows = rank_scores([4, 3, 2, 1], list('ABCD'))
        for n in range(1, 5):
            q = query(list('ABCD')[:n])
            e = evaluate_query(rows, q, {})
            self.assertEqual(len(e['components']), n)
            m = e['query_metrics']['exact']
            self.assertEqual(m['component_recall_at_1'], 1/n)
            self.assertEqual(e['exact_recall_ceiling']['1'], 1/n)
            self.assertEqual(m['any_component_hit_at_1'], 1)
            self.assertEqual(m['all_components_hit_at_1'], int(n == 1))
            self.assertEqual(m['all_components_hit_at_5'], 1)

    def test_duplicate_equivalence_is_separate_from_accession(self):
        rows = rank_scores([3, 2, 1], ['A', 'X', 'B'])
        e = evaluate_query(rows, query(['B']), {'B': ['A', 'B']})
        self.assertEqual(e['components'][0]['exact']['rank'], 3)
        self.assertEqual(e['components'][0]['duplicate_equivalent']['rank'], 1)
        self.assertEqual(e['query_metrics']['exact']['component_recall_at_1'], 0)
        self.assertEqual(e['query_metrics']['duplicate_equivalent']['component_recall_at_1'], 1)

    def test_absent_reference_is_unscorable_not_a_failure(self):
        e = evaluate_query(rank_scores([1], ['A']), query(['A', 'B']), {})
        self.assertIsNone(e['components'][1]['exact']['rank'])
        self.assertIsNone(e['components'][1]['exact']['hit_at_1'])
        m = e['query_metrics']['exact']
        self.assertEqual((m['scorable_components'], m['unscorable_components']), (1, 1))
        self.assertEqual(m['component_recall_at_1'], 1)
        self.assertIsNone(m['all_components_hit_at_1'])
        self.assertIsNone(m['any_component_hit_at_1'])
        s = summarize([e])['kinds']['exact']
        self.assertEqual(s['original_components'], 2)
        self.assertEqual(s['scorable_components'], 1)
        self.assertEqual(s['recall_at_1'], 1)
        self.assertEqual(s['query_metrics']['excluded_unscorable_queries'], 1)
        self.assertIsNone(s['query_metrics']['all_components_hit_at_1'])

    def test_common_mask_removes_every_component_duplicate(self):
        q = query(['A', 'B'])
        equivalence = {'A': ['A', 'A2'], 'B': ['B', 'B2']}
        mask = common_reference_removal_mask(q, equivalence, ['A', 'A2', 'B', 'B2', 'C', 'D'])
        self.assertEqual(mask, ['C', 'D'])
        rows = filter_ranking(rank_scores([6, 5, 4, 3, 2, 1], ['A', 'A2', 'B', 'B2', 'C', 'D']), mask)
        e = evaluate_query(rows, q, equivalence, {'A': ['A2', 'C'], 'B': []})
        self.assertEqual([r['rank'] for r in rows], [1, 2])
        self.assertEqual(e['query_metrics']['exact']['scorable_components'], 0)
        self.assertEqual(e['query_metrics']['labelled']['scorable_components'], 1)
        self.assertEqual(e['components'][0]['labelled']['accepted_surviving'], ['C'])

    def test_catalogue_mask_recomputes_saved_tie_intervals(self):
        rows = rank_scores([5, 3, 3, 3-5e-7], ['X', 'A', 'B', 'C'])
        for row in rows:
            row.update(tie_best_rank=2, tie_worst_rank=3, near_tie_count_1e6=3,
                       near_tie_best_rank=2, near_tie_worst_rank=4)
        filtered = filter_ranking(rows, ['B', 'C'])
        self.assertEqual(filtered[0]['rank'], 1)
        self.assertEqual(filtered[0]['tie_best_rank'], 1)
        self.assertEqual(filtered[0]['tie_worst_rank'], 1)
        self.assertEqual(filtered[0]['near_tie_count_1e6'], 2)
        self.assertEqual(filtered[0]['near_tie_best_rank'], 1)
        self.assertEqual(filtered[0]['near_tie_worst_rank'], 2)

    def test_summary_macro_and_micro_denominators_are_distinct(self):
        rows = rank_scores([2, 1], ['A', 'B'])
        evaluations = [evaluate_query(rows, query(['A'], 'a1', 'g1'), {}),
                       evaluate_query(rows, query(['A'], 'a2', 'g1'), {}),
                       evaluate_query(rows, query(['A', 'B'], 'ab', 'g2'), {})]
        s = summarize(evaluations)
        self.assertEqual((s['query_count'], s['component_count'], s['distinct_group_count']), (3, 4, 2))
        exact = s['kinds']['exact']
        self.assertEqual(exact['recall_at_1'], .75)
        self.assertEqual(exact['per_target_macro_recall_at_1'], .5)
        self.assertEqual(exact['query_metrics']['all_components_hit_at_1'], 2/3)
        self.assertEqual(exact['by_component_count']['2']['original_components'], 2)
        self.assertEqual(exact['by_position']['N']['parent_queries'], 1)

    def test_oracle_does_not_borrow_a_hit_from_another_component(self):
        rows = rank_scores([2, 1], ['A', 'B'])
        q = query(['A', 'B'])
        e = evaluate_oracle_query({'A': rows, 'B': rows}, q, {})
        self.assertEqual(e['query_metrics']['exact']['component_recall_at_1'], .5)
        self.assertEqual(e['retrieval_budget'], 'k per known component')
        self.assertEqual(e['returned_hits_ceiling']['1'], 2)
        self.assertEqual(e['exact_recall_ceiling']['1'], 1)
        reversed_rows = rank_scores([1, 2], ['A', 'B'])
        good = evaluate_oracle_query({'A': rows, 'B': reversed_rows}, q, {})
        self.assertEqual(good['query_metrics']['exact']['all_components_hit_at_1'], 1)
        with self.assertRaises(ValueError):
            evaluate_oracle_query({'A': rows, 'B': rows[:1]}, q, {})

    def test_invalid_and_mixed_results_are_rejected(self):
        with self.assertRaises(ValueError):
            rank_scores([float('nan')], ['A'])
        with self.assertRaises(ValueError):
            rank_scores([1, 2], ['A', 'A'])
        rows = rank_scores([2, 1], ['A', 'B'])
        with self.assertRaises(ValueError):
            evaluate_query(rows[::-1], query(['A']), {})
        e = evaluate_query(rows, query(['A']), {})
        with self.assertRaises(ValueError):
            summarize([e, e])
        with self.assertRaises(ValueError):
            filter_ranking(rows, ['X'])

    def test_fraction_boundaries(self):
        expected = {0.1: 'lt_0.25', .25: '0.25_to_lt_0.40', .4: '0.40_to_lt_0.60',
                    .6: '0.60_to_lt_0.75', .75: 'ge_0.75', 1: 'ge_0.75'}
        for value, name in expected.items():
            self.assertEqual(fraction_bin(value), name)
        for value in (0, -1, 1.1, float('nan')):
            with self.assertRaises(ValueError):
                fraction_bin(value)

    def test_hand_computed_cosines(self):
        rows = cosine_ranking([3., 4.], [[3., 4.], [4., -3.], [-3., -4.]], ['A', 'B', 'C'])
        self.assertEqual([r['score'] for r in rows], [1, 0, -1])
        with self.assertRaises(ValueError):
            cosine_ranking([0., 0.], [[1., 0.]], ['A'])


class PreparedInputTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.packet = cls.repo / 'config/reproduction/packet'
        if os.environ.get('IRI_RUN_DATA_TESTS') != '1':
            raise unittest.SkipTest('Install the verified input bundle and set IRI_RUN_DATA_TESTS=1 for frozen-data integration tests')
        cls.suite = load_suite(cls.repo, cls.packet)

    def test_all_reviewed_inputs_and_designs(self):
        v = self.suite['verification']
        self.assertEqual(v['query_count'], 1146)
        self.assertEqual(v['all_component_observations'], 3536)
        self.assertEqual(v['reviewed_records_checked'], 105)
        self.assertEqual(v['unique_sequences'], 101)
        self.assertEqual(v['max_query_length'], 1011)

    def test_nested_catalogues_retain_targets_and_share_full_endpoint(self):
        result = catalogue_masks(self.suite)
        self.assertEqual(len(result['masks']), 8)
        self.assertEqual(result['aliases'], {'unique101_rep1': 'unique101_rep0', 'unique101_rep2': 'unique101_rep0'})
        targets = set(self.suite['manifest']['targets'])
        for replicate in range(3):
            small = set(result['masks'][f'unique25_rep{replicate}'])
            medium = set(result['masks'][f'unique50_rep{replicate}'])
            large = set(result['masks']['unique101_rep0'])
            self.assertTrue(targets <= small < medium < large)

    def test_corrupted_manifest_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for filename in ('challenge_manifest.json', 'preparation_checks.json'):
                (root / filename).write_bytes((self.packet / filename).read_bytes())
            with (root / 'challenge_manifest.json').open('a') as f:
                f.write(' ')
            with self.assertRaisesRegex(ValueError, 'Challenge manifest checksum'):
                load_suite(self.repo, root)


if __name__ == '__main__':
    unittest.main()
