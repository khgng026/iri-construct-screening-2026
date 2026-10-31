from pathlib import Path
import sys
import unittest
import copy
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_context_compare import require_common_panel, require_shared_policy, contrast


class CommonProfileTests(unittest.TestCase):
    def test_common_panel_cannot_hide_different_geometry_or_source_policy(self):
        design = {'geometry': [[0, 8], [8, 32]]}
        identity = {'design': design, 'design_sha256': 'pinned', 'quantiles': [0., .05, .25, .5, .75, .95, 1.],
                    'tolerance': 1e-6, 'numerical_policy': 'per-residue comparisons and reductions float64; stored means float32',
                    'storage_policy': 'standalone residue states only; joined matrices and per-residue deltas never persisted',
                    'sequence_sha256': {'A': 'canonical'},
                    'code_hashes': {'scripts/session_context_profile.py': '6c7cf52848981adf1b707cc25948c9ea5850e55fc9fa205c2b73fadd8666e18d'}}
        require_shared_policy([identity, copy.deepcopy(identity)], design, 'pinned', {'A': 'canonical'})
        for key, changed in [('design', {'geometry': [[0, 7], [7, 32]]}), ('sequence_sha256', {'A': 'changed'}),
                             ('quantiles', [0., .5, 1.]), ('numerical_policy', 'float32'), ('code_hashes', {})]:
            other = copy.deepcopy(identity)
            other[key] = changed
            with self.assertRaises(ValueError):
                require_shared_policy([identity, other], design, 'pinned', {'A': 'canonical'})

    def test_missing_or_duplicated_query_cannot_shrink_common_panel(self):
        ids = [str(i) for i in range(32)]
        self.assertEqual(require_common_panel([{'query_ids': ids}, {'query_ids': ids + ['extra']}], ids), sorted(ids))
        for malformed in (ids[:-1], ids + [ids[0]]):
            with self.assertRaises(ValueError):
                require_common_panel([{'query_ids': malformed}], ids)

    def test_group_target_contrast_is_not_weighted_by_repeated_orders(self):
        rows = [{'accession': 'A', 'matched_query_count': 2, 'cosine_far_minus_near': .2, 'relative_l2_far_minus_near': -.1},
                {'accession': 'B', 'matched_query_count': 8, 'cosine_far_minus_near': .4, 'relative_l2_far_minus_near': -.3}]
        result = contrast(rows)
        self.assertAlmostEqual(result['cosine_far_minus_near_mean'], .3)
        self.assertEqual(result['matched_query_component_observations'], 10)
        self.assertIsNone(contrast([])['cosine_far_minus_near_mean'])


if __name__ == '__main__':
    unittest.main()
