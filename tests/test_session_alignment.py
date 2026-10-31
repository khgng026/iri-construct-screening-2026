"""Toy strings test alignment arithmetic; they are not proposed biological inputs."""
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from iri_constructs.session_alignment import (
    AlignmentConfig, align_pair, alignment_metadata, independent_affine_score,
    make_aligner, rank_alignment,
)


class SessionAlignmentTests(unittest.TestCase):
    def test_gap_convention_one_and_four_residues(self):
        for gap, expected in [("A", 78.0), ("AAAA", 76.5)]:
            q, r = "WWWW" + gap + "WWWW", "WWWWWWWW"
            row = align_pair(q, r)
            self.assertEqual(row["score"], expected)
            self.assertEqual(row["score"], independent_affine_score(q, r))
            self.assertEqual(row["gap_runs"], 1)
            self.assertEqual(row["gap_columns"], len(gap))
            self.assertEqual(row["aligned_length"], len(q))
            self.assertEqual(row["identities"], 8)
            self.assertEqual(row["query_span"], [0, len(q)])
            self.assertEqual(row["reference_span"], [0, 8])
            self.assertEqual(row["query_coverage"], 1.0)
            self.assertEqual(row["query_paired_coverage"], 8 / len(q))

    def test_known_local_span_orientation(self):
        row = align_pair("AAAAWWWWAAAA", "WWWW")
        self.assertEqual(row["score"], 44.0)
        self.assertEqual(row["query_span"], [4, 8])
        self.assertEqual(row["reference_span"], [0, 4])
        self.assertEqual(row["query_coverage"], 1 / 3)
        self.assertEqual(row["reference_coverage"], 1.0)
        self.assertEqual(row["identity_fraction"], 1.0)

    def test_no_positive_alignment(self):
        row = align_pair("AAAA", "WWWW")
        self.assertEqual(row["score"], 0.0)
        self.assertIsNone(row["query_span"])
        self.assertEqual(row["coordinates"], [[], []])
        self.assertEqual(row["query_coverage"], 0.0)

    def test_complete_catalogue_and_accession_ties(self):
        result = rank_alignment("ARND", {"Z": "ARND", "A": {"sequence": "ARND"}, "B": "WWWW"})
        rows = result["ranking"]
        self.assertEqual([r["accession"] for r in rows], ["A", "Z", "B"])
        self.assertEqual([r["rank"] for r in rows], [1, 2, 3])
        self.assertEqual(rows[0]["score"], 21.0)
        self.assertEqual(rows[0]["exact_score_tie_count"], 2)
        self.assertEqual((rows[1]["tie_best_rank"], rows[1]["tie_worst_rank"]), (1, 2))
        self.assertIn("length - 1", result["method"]["gap_cost"])

    def test_independent_dp_exhaustive_toy_panel(self):
        # Exhaust all 1..3-symbol words over a fixed two-letter test alphabet.
        panel = ["".join(x) for n in (1, 2, 3) for x in itertools.product("AW", repeat=n)]
        aligner = make_aligner()
        for q in panel:
            for r in panel:
                self.assertEqual(float(aligner.score(q, r)), independent_affine_score(q, r))

    def test_reject_unapproved_symbols_and_configuration_drift(self):
        for q in ("", "ACX", "ac", None):
            with self.assertRaises(ValueError):
                align_pair(q, "AC")
        with self.assertRaises(ValueError):
            rank_alignment("AC", {})
        with self.assertRaises(ValueError):
            AlignmentConfig(gap_open=11)
        metadata = alignment_metadata()
        self.assertEqual(metadata["implementation_version"], "1.88")
        self.assertEqual(len(metadata["matrix_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
