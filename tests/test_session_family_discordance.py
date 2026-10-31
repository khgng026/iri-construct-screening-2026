from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_family_discordance import best_accepted, outcome, summarize_rows, case_selection, summarize_queries


class PairedFamilyTests(unittest.TestCase):
    def test_absent_label_and_present_miss_are_distinct(self):
        ranking = [{'accession': str(i), 'rank': i + 1, 'score': 10 - i} for i in range(8)]
        absent = best_accepted(ranking, ['outside'], 5)
        miss = best_accepted(ranking, ['7'], 5)
        self.assertIsNone(absent['hit'])
        self.assertFalse(miss['hit'])
        self.assertEqual(outcome(absent, absent), 'unscorable')
        self.assertEqual(outcome(miss, miss), 'both_miss')

    def test_any_accepted_survivor_and_paired_mask_equality(self):
        ranking = [{'accession': 'A', 'rank': 1, 'score': 1}, {'accession': 'B', 'rank': 2, 'score': .5}]
        selected = best_accepted(ranking, ['B', 'A', 'missing'], 1)
        self.assertEqual(selected['best_accession'], 'A')
        self.assertEqual(selected['accepted_survivors'], ['A', 'B'])
        with self.assertRaises(ValueError):
            outcome(selected, best_accepted(ranking, ['B'], 1))

    def test_cases_are_lexical_not_largest_score_gap_and_denominators_keep_null(self):
        rows = [{'query_id': str(i), 'target': 'T', 'group_id': str(i), 'outcome': 'model_only'} for i in range(8, -1, -1)]
        self.assertEqual([r['query_id'] for r in case_selection(rows)['model_only']], ['0', '1', '2', '3', '4'])
        row = dict(rows[0], outcome='unscorable')
        summary = summarize_rows([row])
        self.assertEqual(summary['original_components'], 1)
        self.assertEqual(summary['scorable_components'], 0)
        self.assertIsNone(summary['model_recall_at_5'])
        query_summary = summarize_queries([row])
        self.assertEqual(query_summary['all_components_scorable_queries'], 0)
        self.assertIsNone(query_summary['model_all_components_recall_at_5'])


if __name__ == '__main__':
    unittest.main()
