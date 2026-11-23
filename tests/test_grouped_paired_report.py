"""Exercise the group summary on paired counts without model inputs."""
from copy import deepcopy
import unittest

from experiments.paired_preprocessing.grouped_report import GROUP_SCHEMA, grouped_report
from experiments.paired_preprocessing.verify import SCHEMA


def _fixture():
    baseline = [{'construct_id': key, 'eligible_components': n, 'recovered_components': recovered}
                for key, n, recovered in [('a', 4, 0), ('b', 4, 0), ('c', 1, 1)]]
    candidate = [dict(row, recovered_components=recovered) for row, recovered in zip(baseline, (2, 2, 0), strict=True)]
    report = {'schema': SCHEMA, 'baseline': 'minimal', 'annotation_reference': None,
              'protocol': {'query_caps': [2], 'primary_query_cap': 2},
              'runs': {'minimal': {'2': {'per_construct': baseline}},
                       'candidate': {'2': {'per_construct': candidate}}}}
    groups = {'schema': GROUP_SCHEMA, 'constructs': [
        {'id': key, 'group': group, 'partition': 'evaluation'} for key, group in [('a', 'big'), ('b', 'big'), ('c', 'small')]
    ]}
    return report, groups


class GroupedReportTests(unittest.TestCase):
    def test_large_group_can_make_the_pooled_gain_positive(self):
        report, groups = _fixture()
        result = grouped_report(report, groups)['candidate_minus_baseline']['candidate'][2]
        self.assertAlmostEqual(result['pooled_recall_gain'], 1 / 3)
        self.assertEqual(result['equal_group_recall_gain'], -0.25)
        self.assertAlmostEqual(result['largest_group_share_of_components'], 8 / 9)
        self.assertEqual(result['leave_one_group_out_pooled_range'], [-1, 0.5])

    def test_rejects_a_declared_group_crossing_partitions(self):
        report, groups = _fixture()
        groups['constructs'].append({'id': 'd', 'group': 'big', 'partition': 'development'})
        with self.assertRaisesRegex(ValueError, 'group crosses partitions'):
            grouped_report(report, groups)

    def test_rejects_unpaired_denominator_and_missing_group(self):
        report, groups = _fixture()
        wrong = deepcopy(report)
        wrong['runs']['candidate']['2']['per_construct'][0]['eligible_components'] += 1
        with self.assertRaisesRegex(ValueError, 'eligible counts differ'):
            grouped_report(wrong, groups)
        groups['constructs'].pop()
        with self.assertRaisesRegex(ValueError, 'not all scored constructs have a group'):
            grouped_report(report, groups)

    def test_one_group_has_no_omission_range(self):
        report, groups = _fixture()
        for row in groups['constructs']:
            row['group'] = 'shared'
        result = grouped_report(report, groups)['candidate_minus_baseline']['candidate'][2]
        self.assertEqual(result['groups_with_eligible_components'], 1)
        self.assertIsNone(result['leave_one_group_out_pooled_range'])


if __name__ == '__main__':
    unittest.main()
