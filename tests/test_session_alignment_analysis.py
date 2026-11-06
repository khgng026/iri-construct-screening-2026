import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("alignment_analysis", ROOT / "scripts/session_alignment_analysis.py")
analysis = importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)


class AlignmentAnalysisTests(unittest.TestCase):
    def test_reference_normalization_preserves_primary_and_changes_order(self):
        raw = [{"accession": "long", "rank": 1, "score": 150., "exact_score_tie_count": 1},
               {"accession": "short", "rank": 2, "score": 50., "exact_score_tie_count": 1}]
        result = analysis.normalized_ranking(raw, {"long": 200., "short": 50.})
        self.assertEqual([r["accession"] for r in result], ["short", "long"])
        self.assertEqual([r["score"] for r in result], [1., .75])
        self.assertEqual([r["raw_alignment_score"] for r in result], [50., 150.])
        self.assertEqual(raw[0]["score"], 150.)
        self.assertEqual(raw[0]["rank"], 1)

    def test_normalized_tie_uses_accession_not_raw_score(self):
        raw = [{"accession": "Z", "rank": 1, "score": 200.},
               {"accession": "A", "rank": 2, "score": 50.}]
        result = analysis.normalized_ranking(raw, {"Z": 200., "A": 50.})
        self.assertEqual([r["accession"] for r in result], ["A", "Z"])
        self.assertEqual(result[0]["exact_score_tie_count"], 2)
        self.assertEqual((result[1]["tie_best_rank"], result[1]["tie_worst_rank"]), (1, 2))

    def test_scores_above_upper_bound_are_not_clipped(self):
        with self.assertRaises(ValueError):
            analysis.normalized_ranking([{"accession": "A", "score": 5.}], {"A": 4.})
        with self.assertRaises(ValueError):
            analysis.normalized_ranking([{"accession": "A", "score": 5.}], {"B": 5.})


if __name__ == "__main__": unittest.main()
