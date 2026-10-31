from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_context_profile_audit import direct_bins, direct_distribution, weighted_mean, independent_mean, independent_summaries


class IndependentProfileArithmeticTests(unittest.TestCase):
    def test_only_artificial_join_distances_define_bins(self):
        self.assertEqual(direct_bins(3, False, False), {-1: [0, 1, 2]})
        right = direct_bins(300, False, True)
        self.assertEqual(right[0], list(range(292, 300)))
        self.assertEqual(right[4], list(range(172)))
        internal = direct_bins(300, True, True)
        self.assertEqual(internal[0], list(range(8)) + list(range(292, 300)))
        self.assertEqual(internal[4], list(range(128, 172)))
        self.assertEqual(sum(map(len, internal.values())), 300)

    def test_weighted_reconstruction_and_linear_quantiles(self):
        value = weighted_mean([np.array([1., 3.]), np.array([7., 9.])], [1, 3])
        np.testing.assert_array_equal(value, np.array([5.5, 7.5]))
        np.testing.assert_array_equal(independent_mean(np.array([[1., 3.], [7., 9.]])), np.array([4., 6.]))
        summary = direct_distribution([9., 1., 3.])
        self.assertEqual(summary['mean'], 13 / 3)
        self.assertEqual(summary['quantiles'][3], 3.)
        self.assertAlmostEqual(summary['quantiles'][1], 1.2)
        self.assertAlmostEqual(summary['quantiles'][5], 8.4)

    def test_matched_contrasts_average_queries_within_group_target(self):
        def row(index, count, cosine, relative):
            return {'accession': 'A', 'bin_index': index, 'bin': '0-7' if index == 0 else '128+', 'residue_count': count,
                    'cosine': {'mean': cosine, 'sum': cosine * count},
                    'relative_l2': {'mean': relative, 'sum': relative * count}}
        records = [
            {'query': {'query_id': 'q1', 'group_id': 'same', 'suite': 'pairs'}, 'scope': 'combined',
             'bins': [row(0, 5, .3, 2.), row(4, 2, .5, 1.)]},
            {'query': {'query_id': 'q2', 'group_id': 'same', 'suite': 'pairs'}, 'scope': 'combined',
             'bins': [row(0, 7, .4, 4.), row(4, 3, .5, 1.)]}]
        strata, matched = independent_summaries(records)
        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['query_ids'], ['q1', 'q2'])
        self.assertAlmostEqual(matched[0]['cosine_far_minus_near'], .15)
        self.assertEqual(matched[0]['relative_l2_far_minus_near'], -2.)
        near = next(row for row in strata if row['suite'] == 'all' and row['bin'] == '0-7')
        self.assertEqual(near['residue_observations'], 12)
        self.assertEqual(near['group_target_pairs'], 1)
        self.assertAlmostEqual(near['cosine']['residue_weighted_mean'], 4.3 / 12)
        self.assertAlmostEqual(near['cosine']['component_bin_mean_distribution']['mean'], .35)


if __name__ == '__main__': unittest.main()
