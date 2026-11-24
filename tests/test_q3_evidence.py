import unittest

from experiments.q3_followup import domain_replay, evidence_audit, open_reporters, rank_replay


def case(query, model, conventional, covered=None, identity=0.2):
    return {
        "query": query,
        "fold": "fold",
        "fine_family": "family",
        "label": "group",
        "max_global_identity": identity,
        "covered_positive_local_identity": covered,
        "correct": {
            **{method: conventional for method in domain_replay.CONVENTIONAL},
            domain_replay.MODEL: model,
        },
    }


class Q3EvidenceTest(unittest.TestCase):
    def test_union_rescue_requires_all_conventional_to_miss(self):
        row = {method: False for method in evidence_audit.METHODS}
        model = evidence_audit.MODELS[0]
        row[model] = True
        row[evidence_audit.CONVENTIONAL[1]] = True
        result = evidence_audit.decision_summary([row, {**row, evidence_audit.CONVENTIONAL[1]: False}])
        self.assertEqual(result["model_correct_on_conventional_misses"][model], 1)
        self.assertEqual(result["paired"][f"{model}_vs_alignment"]["only_model"], 2)

    def test_coverage_subset_excludes_newly_covering_wider_reference(self):
        a = {"cases": [case("1", True, False), case("2", True, False)]}
        b = {"cases": [case("2", True, False), case("1", False, True, covered=0.23)]}
        results = domain_replay.replay(a, b)["bands"]
        self.assertEqual(
            results["under_25_both"]["original_to_wider_vs_conventional_union"]["model_only"],
            {"model_only": 1, "conventional_only": 1},
        )
        subset = results["under_25_and_no_covered_same_label_hit_both"]
        self.assertEqual(subset["wider"]["domains"], 1)
        self.assertEqual(subset["wider"]["paired_150m_local"][domain_replay.UNION]["only_model"], 1)

    def test_refuses_incompatible_domain_pair(self):
        a = {"cases": [case("1", True, True)]}
        with self.assertRaisesRegex(ValueError, "no longer qualifies"):
            domain_replay.replay(a, {"cases": [case("1", True, True, identity=0.25)]})
        with self.assertRaisesRegex(ValueError, "same domains"):
            domain_replay.replay(a, {"cases": [case("2", True, True)]})

    def test_zero_event_planning_bound(self):
        self.assertGreater(evidence_audit.zero_event_upper_bound(58), 0.05)
        self.assertLess(evidence_audit.zero_event_upper_bound(59), 0.05)

    def test_rank_combination_uses_both_ranks_and_conventional_tiebreak(self):
        self.assertEqual(rank_replay.combine_ranks(["a", "b", "c"], ["c", "b", "a"]), "a")
        self.assertEqual(rank_replay.combine_ranks(["a", "b", "c", "d"], ["d", "b", "c", "a"]), "b")
        with self.assertRaisesRegex(ValueError, "same families"):
            rank_replay.combine_ranks(["a", "b"], ["c", "b"])

    def test_reach_keeps_locus_and_exact_peptide_diagnostics_distinct(self):
        target = {"strand": 1, "spans": [(6, 15)], "protein": "MKS"}
        pool = [
            {"strand": 1, "spans": [(0, 15)], "protein": "AAAAA"},
            {"strand": 1, "spans": [(0, 15)], "protein": "AAMKS"},
        ]
        tiers = open_reporters._reach(target, pool, 0.8)["caps"]
        self.assertTrue(tiers["1"]["reached"])
        self.assertFalse(tiers["1"]["contains_documented_peptide"])
        self.assertTrue(tiers["2"]["contains_documented_peptide"])
        self.assertFalse(tiers["2"]["exact"])
