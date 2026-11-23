"""Ensure a crossed comparison uses paired inputs and computes both effects."""
from copy import deepcopy
import unittest

from experiments.paired_preprocessing.compare_characterizers import METHODS, compare
from experiments.paired_preprocessing.verify import SCHEMA


def _fixture():
    common = {'schema': SCHEMA, 'query_caps': [1], 'primary_query_cap': 1, 'prediction_limit': 1,
              'minimum_component_overlap': 0.8, 'minimum_interesting_gain': 0.05,
              'reference_panel_sha256': 'panel', 'sequence_input_sha256': 'same-sequence'}
    truth = {'schema': SCHEMA, 'constructs': [{
        'id': 'c', 'archetype': 'single_cds', 'length': 90,
        'components': [{'id': 'annotated', 'family': 'a', 'spans': [[0, 90]], 'strand': 1}],
    }]}
    protocols, runs = {}, {}
    for scorer in ('esm2', 'blosum62'):
        protocols[scorer] = dict(common, characterization_id=scorer, characterizer=scorer)
        runs[scorer] = {}
        for method in METHODS:
            correct = scorer == 'esm2' and method == 'gene_caller'
            runs[scorer][method] = {
                'schema': SCHEMA, 'characterization_id': scorer, 'characterizer': scorer, 'method': method,
                'reference_panel_sha256': 'panel', 'input_sha256': 'same-sequence',
                'selected_inputs_sha256': '0' * 64,
                'tiers': [{'query_cap': 1, 'constructs': [{
                    'id': 'c', 'queries': [{'id': 'q', 'strand': 1, 'spans': [[0, 90]], 'input_residues': 30}],
                    'predictions': [{'query_id': 'q', 'family': 'a' if correct else 'b'}],
                    'candidate_count': 1, 'status': 'ok', 'preprocess_seconds': 0, 'characterize_seconds': 0,
                }]}],
            }
    return protocols, runs, truth


class CharacterizerComparisonTests(unittest.TestCase):
    def test_reports_paired_recoveries_and_difference_in_preprocessing_gain(self):
        report = compare(*_fixture())
        result = report['by_query_cap'][1]['gene_caller']
        self.assertEqual(result['scorer_only_recoveries'], {'esm2': 1, 'blosum62': 0})
        self.assertEqual(result['recall_gain_over_minimal'], {'esm2': 1.0, 'blosum62': 0.0})
        self.assertEqual(result['difference_in_gain'], 1.0)

    def test_rejects_different_reference_panel_or_candidate_selection(self):
        protocols, runs, truth = _fixture()
        wrong = deepcopy(protocols)
        wrong['esm2']['reference_panel_sha256'] = 'other-panel'
        with self.assertRaisesRegex(ValueError, 'share the evaluation protocol'):
            compare(wrong, runs, truth)
        runs['blosum62']['minimal']['selected_inputs_sha256'] = '1' * 64
        with self.assertRaisesRegex(ValueError, 'selected inputs differ'):
            compare(protocols, runs, truth)


if __name__ == '__main__':
    unittest.main()
