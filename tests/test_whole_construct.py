"""Regression checks for the full-pool benchmark's independent scoring and evaluation."""

from contextlib import contextmanager
from unittest import TestCase
from unittest.mock import patch

from experiments.preprocessing_exploration.hybrid import Candidate
from experiments.whole_construct.availability import component_evidence
from experiments.whole_construct.compare import ARMS, compare
from experiments.whole_construct.evaluate import match_calls, summarize
from experiments.whole_construct.scoring import esm2_calls


class WholeConstructTest(TestCase):
    def setUp(self):
        self.construct = {
            "id": "example", "length": 2000,
            "components": [{"id": "a", "family": "A", "strand": 1, "spans": [[100, 400]]},
                           {"id": "b", "family": "A", "strand": 1, "spans": [[800, 1100]]}],
        }

    def test_held_out_matching_keeps_distinct_copies_and_consolidates_repeated_calls(self):
        calls = [
            {"family": "A", "score": 60, "strand": 1, "spans": [[100, 400]]},
            {"family": "A", "score": 50, "strand": 1, "spans": [[101, 398]]},
            {"family": "A", "score": 45, "strand": -1, "spans": [[800, 1100]]},
            {"family": "A", "score": 40, "strand": 1, "spans": [[800, 1100]]},
            {"family": "B", "score": 30, "strand": 1, "spans": [[1300, 1600]]},
        ]
        self.assertEqual(match_calls(self.construct, calls, 0, 0.8), ({"a", "b"}, 4, 2))
        self.assertEqual(match_calls(self.construct, calls, 41, 0.8), ({"a"}, 2, 2))
        summary = summarize({"constructs": [self.construct]}, {"example": calls}, 0, 0.8)
        self.assertEqual(summary["unverified_partial"], 2)
        self.assertEqual(summary["unmatched_complete"], 0)

    def test_a_family_call_can_be_locally_correct_without_covering_the_locus(self):
        call = {"family": "A", "score": 1, "strand": 1, "spans": [[100, 250]]}
        self.assertEqual(match_calls(self.construct, [call], 0, 0.5), ({"a"}, 1, 2))
        self.assertEqual(match_calls(self.construct, [call], 0, 0.8), (set(), 1, 2))

    def test_missing_cached_score_is_unknown_and_does_not_count_as_a_hit(self):
        old = Candidate("minimal", 0, 1, 100, 400, 100)
        new = Candidate("minimal", 1, 1, 97, 400, 101)
        evidence = component_evidence((1, ((100, 400),)), "A", [old, new], {old: ["B", "A"]}, 0.8)
        self.assertFalse(evidence["top_family_in_cached_covering_candidate"])
        self.assertTrue(evidence["possible_top_family_in_full_pool"])

    def test_protein_arm_scores_the_entire_pool(self):
        candidates = [{"id": str(i), "strand": 1, "spans": [[i * 120, (i + 1) * 120]],
                       "protein": "A" * 40} for i in range(7)]

        @contextmanager
        def fixture(_):
            yield (lambda protein: [("A", len(protein) / 100), ("B", 0.1)]), {"test_double": True}

        with patch("experiments.whole_construct.scoring.candidate_units", return_value=candidates):
            calls, costs = esm2_calls("minimal", {"example": "A" * 900}, [], scorer=fixture)
        self.assertEqual([row["query_id"] for row in calls["example"]], [str(i) for i in range(7)])
        self.assertEqual(calls["example"][-1]["ranking"], [("A", 0.4), ("B", 0.1)])
        self.assertEqual(costs["per_record"]["example"]["residues"], 280)

    def test_calibrated_threshold_is_chosen_without_held_out_labels(self):
        development = self.construct | {"id": "development"}
        low_locus = {"family": "A", "strand": 1, "spans": [[100, 400]]}
        off_locus = {"family": "B", "strand": 1, "spans": [[1600, 1900]]}
        dev = {
            arm: {"input_sha256": "dev", "references_sha256": "panel", "calls": {
                "development": [low_locus | {"score": 0.6}, off_locus | {"score": 0.4}],
            }} for arm in ARMS
        }
        held = {
            arm: {"input_sha256": "held", "references_sha256": "panel", "calls": {
                "example": [low_locus | {"score": 0.9 if arm == "gene_caller" else 0.5}],
            }} for arm in ARMS
        }
        metadata = {"example": {"lineage_group": "single-lineage", "components": [
            {"id": "a", "similarity": {"tier": "low"}}, {"id": "b", "similarity": {"tier": "high"}},
        ]}}
        result = compare({"constructs": [development]}, {"constructs": [self.construct]}, metadata,
                         dev, held, [0], 0.8)[0]
        self.assertTrue(all(row["cutoff"] == 0.6 for row in result["development"].values()))
        self.assertEqual(result["held_out_by_identity"]["low"]["recovered"],
                         {"minimal": 0, "gene_caller": 1, "blastx": 0})
