"""Lightweight orchestration/eligibility tests; never import upstream JIT or models."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from iri_constructs.session_store import UnitStore
import session_plmblast_experiment as runner


def query(components):
    return {"query_id": "q", "suite": "quadruples" if len(components) == 4 else "pairs", "group_id": "g", "components": components}


def candidate(accession, scores=()):
    return {"accession": accession, "status": "eligible_alignment" if scores else "no_eligible_alignment",
            "eligible_hits": [{"cli_score": score, "raw_score": score + 0.0001} for score in scores]}


class EligibleRetrievalTests(unittest.TestCase):
    def test_no_hit_is_scorable_miss_without_filling_topk(self):
        rows = [candidate("A"), candidate("B", [.8]), candidate("C")]
        ranking = runner.ranking_from_pairs(rows)
        self.assertEqual([r["accession"] for r in ranking], ["B"])
        result = runner.eligible_metrics(ranking, ["A", "B", "C"], query(["A", "B"]), {})
        a, b = [c["matches"]["exact"] for c in result["components"]]
        self.assertTrue(a["scorable"] and a["abstention"])
        self.assertIsNone(a["rank"])
        self.assertEqual(a["reciprocal_rank"], 0)
        self.assertEqual(a["hit_at_10"], 0)
        self.assertEqual(b["hit_at_1"], 1)
        self.assertEqual(result["returned_counts"], {"1": 1, "5": 1, "10": 1})
        self.assertEqual(result["query_metrics"]["exact"]["all_components_hit_at_10"], 0)

    def test_empty_predictions_retain_every_present_target(self):
        result = runner.eligible_metrics([], ["A", "B"], query(["A", "B"]), {})
        self.assertEqual(result["query_metrics"]["exact"]["scorable_components"], 2)
        self.assertEqual(result["query_metrics"]["exact"]["any_component_hit_at_10"], 0)
        self.assertEqual(result["returned_counts"]["10"], 0)

    def test_missing_catalogue_label_is_unscorable_but_family_survivor_is_miss(self):
        result = runner.eligible_metrics([], ["B", "C"], query(["A", "D"]), {}, {"A": ["A", "C"], "D": ["D"]})
        self.assertFalse(result["components"][0]["matches"]["exact"]["scorable"])
        self.assertTrue(result["components"][0]["matches"]["family"]["scorable"])
        self.assertEqual(result["components"][0]["matches"]["family"]["hit_at_5"], 0)
        self.assertIsNone(result["query_metrics"]["family"]["all_components_hit_at_5"])

    def test_family_any_member_recovery_and_unique_accession_budget(self):
        ranking = runner.ranking_from_pairs([candidate("X", [.8, .9]), candidate("Y", [.8]), candidate("A")])
        result = runner.eligible_metrics(ranking, ["X", "Y", "A"], query(["A", "B"]), {}, {"A": ["X"], "B": ["Y"]})
        self.assertEqual(len(ranking), 2)
        self.assertEqual(result["query_metrics"]["family"]["component_hits_at_1"], 1)
        self.assertEqual(result["query_metrics"]["family"]["all_components_hit_at_5"], 1)

    def test_rounded_tie_uses_accession_not_unrounded_score(self):
        rows = [candidate("Z", [.7]), candidate("A", [.7])]
        rows[0]["eligible_hits"][0]["raw_score"] = .70049
        rows[1]["eligible_hits"][0]["raw_score"] = .70041
        self.assertEqual([r["accession"] for r in runner.ranking_from_pairs(rows)], ["A", "Z"])

    def test_failed_candidate_cannot_be_coerced_to_no_hit(self):
        with self.assertRaises(ValueError):
            runner.ranking_from_pairs([{"accession": "A", "status": "failed", "eligible_hits": []}])
        with self.assertRaises(ValueError):
            runner.ranking_from_pairs([candidate("A"), candidate("A")])
        with self.assertRaises(ValueError):
            runner.ranking_from_pairs([{"accession": "A", "status": "eligible_alignment", "eligible_hits": []}])

    def test_zero_all_scorable_quadruple_denominator_stays_null(self):
        evaluation = runner.eligible_metrics([], ["X"], query(["A", "B", "C", "D"]), {}, {"A": ["X"]})
        summary = runner.summarize_eligible([{"status": "success", "suite": "quadruples", "evaluation": evaluation}])
        family = summary["quadruples"]["kinds"]["family"]
        self.assertEqual(family["scorable_components"], 1)
        self.assertEqual(family["no_hit_components"], 1)
        self.assertEqual(family["fully_scorable_queries"], 0)
        self.assertIsNone(family["all_components_hit_at_5"])
        self.assertEqual(family["component_recall_at_5"], 0)

    def test_technical_failures_have_separate_coverage(self):
        evaluation = runner.eligible_metrics([], ["A"], query(["A"]), {})
        summary = runner.summarize_eligible([{"status": "success", "suite": "clean", "evaluation": evaluation},
                                            {"status": "failed", "suite": "pairs"}])
        self.assertEqual(summary["all"]["expected_queries"], 2)
        self.assertEqual(summary["all"]["evaluated_queries"], 1)
        self.assertEqual(summary["all"]["failed_queries"], 1)
        self.assertEqual(summary["all"]["kinds"]["exact"]["mrr"], 0)


class PathAndSourceTests(unittest.TestCase):
    def test_reference_query_orientation_bounds_and_gaps(self):
        runner.validate_path([[0, 2], [1, 2], [2, 3], [2, 4]], 5, 3)
        for path in ([[0, 2], [0, 2]], [[0, 2], [2, 3]], [[-1, 0]], [[3, 0]], [[0, 5]], [[1, 1], [0, 2]]):
            with self.subTest(path=path), self.assertRaises(ValueError):
                runner.validate_path(path, 5, 3)

    def test_probe_order_is_fixed_cartesian_without_substitution(self):
        self.assertEqual(runner.pair_list(["q67", "q1011"], ["r67", "r642"]),
                         [("q67", "r67"), ("q67", "r642"), ("q1011", "r67"), ("q1011", "r642")])

    def test_completed_sources_reject_missing_units_and_changed_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            session = Path(directory)
            source = UnitStore(session / "jobs/source", {"unit_keys": ["a", "summary"]})
            source.publish("a", {"status": "failed"})
            with self.assertRaises(ValueError):
                runner.open_completed(source.root, session)
            source.publish("summary", {"coverage": "failed retained"})
            source.complete(["a", "summary"])
            opened = runner.open_completed(source.root, session)
            self.assertEqual(runner.completed_unit(opened, "a")["status"], "failed")
            (opened.completed_paths["a"] / "record.json").write_text("{}")
            with self.assertRaises(ValueError):
                runner.open_completed(source.root, session)

    def test_sparse_ranking_duplicates_and_nonfinite_scores_rejected(self):
        for rows in ([dict(accession="A", rank=1, score=float("nan"))],
                     [dict(accession="A", rank=1, score=1), dict(accession="A", rank=2, score=.5)]):
            with self.assertRaises(ValueError):
                runner.eligible_metrics(rows, ["A"], query(["A"]), {})

    def test_source_boundary_still_rejects_external_job(self):
        with tempfile.TemporaryDirectory() as directory:
            session = Path(directory)
            with self.assertRaises(ValueError):
                runner.open_completed(session / "outside_jobs", session)


if __name__ == "__main__":
    unittest.main()
