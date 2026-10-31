"""Fixed-total oracle fusion preserves all component inputs and exact tie policy."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_oracle_fusion import independent_fusion_check
from iri_constructs.session_segmentation import fuse_rankings
from iri_constructs.session_benchmark import rank_scores


class OracleFusionTests(unittest.TestCase):
    def setUp(self):
        self.rows = [rank_scores([.9, .8, .1], ['A', 'B', 'C'], 'cosine_similarity'),
                     rank_scores([.1, .8, .9], ['A', 'B', 'C'], 'cosine_similarity')]

    def test_max_score_and_rrf_use_each_component_without_extra_slots(self):
        for method in ('max_score', 'rrf'):
            ranking = fuse_rankings(self.rows, method)['ranking']
            independent_fusion_check(self.rows, ranking, method)
            self.assertEqual([r['accession'] for r in ranking], ['A', 'C', 'B'])
            self.assertEqual([r['rank'] for r in ranking], [1, 2, 3])
        ranking = fuse_rankings(self.rows, 'max_score')['ranking']
        self.assertEqual([r['score'] for r in ranking], [.9, .9, .8])

    def test_changed_scores_order_and_missing_reference_are_rejected(self):
        ranking = fuse_rankings(self.rows, 'max_score')['ranking']
        broken = deepcopy(ranking); broken[0]['fusion_score'] += .01
        with self.assertRaises(ValueError): independent_fusion_check(self.rows, broken, 'max_score')
        with self.assertRaises(ValueError): independent_fusion_check(self.rows, ranking[::-1], 'max_score')
        broken_rows = deepcopy(self.rows); broken_rows[0].pop()
        with self.assertRaises(ValueError): independent_fusion_check(broken_rows, ranking, 'max_score')

    def test_single_component_leaves_native_order_and_values(self):
        ranking = fuse_rankings(self.rows[:1], 'max_score')['ranking']
        independent_fusion_check(self.rows[:1], ranking, 'max_score')
        self.assertEqual([(r['accession'], r['score']) for r in ranking],
                         [(r['accession'], r['score']) for r in self.rows[0]])


if __name__ == '__main__':
    unittest.main()
