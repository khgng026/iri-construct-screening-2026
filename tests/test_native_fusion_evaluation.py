"""Check family ranking and region support without model or source records."""
import unittest

from experiments.benign_protein_pilot.native_fusions import _recovery, _top_labels


class NativeFusionEvaluationTests(unittest.TestCase):
    def test_top_labels_uses_the_highest_scoring_input_for_each_label(self):
        spans = [(0, 40), (40, 80)]
        calls = _top_labels(spans, [{'left': 0.6, 'right': 0.7, 'other': 0.1},
                                    {'left': 0.8, 'right': 0.9, 'other': 0.2}])
        self.assertEqual(calls, {'right': (40, 80), 'left': (40, 80)})
        regions = [{'label': 'left', 'span': [0, 40]}, {'label': 'right', 'span': [40, 80]}]
        self.assertEqual(_recovery(calls, regions), {'label_hits': 2, 'region_supported_hits': 1})


if __name__ == '__main__':
    unittest.main()
