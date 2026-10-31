import copy
import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("independent_audit", ROOT / "scripts/session_independent_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class IndependentAuditTests(unittest.TestCase):
    def fixture(self):
        # Production outputs are a test fixture only; the auditor never imports
        # their implementations and computes its own expected values.
        from iri_constructs.session_benchmark import rank_scores, evaluate_query, summarize
        names = [f"R{i:03d}" for i in range(12)]
        ranking = rank_scores(list(range(12, 0, -1)), names)
        query = {"query_id": "test", "group_id": "group", "suite": "quadruples",
                 "components": [names[i] for i in (0, 1, 5, 10)], "component_count": 4,
                 "spans": [{"accession": names[j], "start0": i, "end0": i + 1, "fraction": .25}
                           for i, j in enumerate((0, 1, 5, 10))]}
        equivalence = {a: [a] for a in names}
        evaluation = evaluate_query(ranking, query, equivalence)
        independent = audit.validate_evaluation(evaluation, query, ranking, equivalence)
        summary = summarize([evaluation])
        return ranking, query, equivalence, evaluation, independent, summary

    def test_independent_all_component_denominators(self):
        ranking, query, eq, ev, independent, summary = self.fixture()
        self.assertEqual(ev["query_metrics"]["exact"]["component_recall_at_5"], .5)
        self.assertEqual(ev["query_metrics"]["exact"]["any_component_hit_at_5"], 1)
        self.assertEqual(ev["query_metrics"]["exact"]["all_components_hit_at_5"], 0)
        audit.validate_summary(summary, [independent])
        bad = copy.deepcopy(summary)
        bad["kinds"]["exact"]["scorable_components"] = 3
        with self.assertRaises(ValueError): audit.validate_summary(bad, [independent])
        bad = copy.deepcopy(ev)
        bad["components"][2]["exact"]["hit_at_5"] = 1
        with self.assertRaises(ValueError): audit.validate_evaluation(bad, query, ranking, eq)
        if "query_metrics" in summary["kinds"]["exact"]["by_component_count"]["4"]:
            bad = copy.deepcopy(summary)
            bad["kinds"]["exact"]["by_component_count"]["4"]["query_metrics"]["original_queries"] = 2
            with self.assertRaises(ValueError): audit.validate_summary(bad, [independent])

    def test_missing_reordered_and_duplicate_rankings_fail(self):
        rows = [{"accession": "A", "rank": 1, "score": 2., "exact_score_tie_count": 1},
                {"accession": "B", "rank": 2, "score": 1., "exact_score_tie_count": 1}]
        audit.validate_ranking(rows, ["A", "B"])
        for bad in [rows[:1], rows[::-1], [rows[0], {**rows[1], "accession": "A"}]]:
            with self.assertRaises(ValueError): audit.validate_ranking(bad, ["A", "B"])

    def test_tie_hit_and_rank_accounting(self):
        rows = [{"accession": "A", "rank": 1, "score": 1., "exact_score_tie_count": 2},
                {"accession": "B", "rank": 2, "score": 1., "exact_score_tie_count": 2}]
        exact = audit.expected_match(rows, ["B"])
        equivalent = audit.expected_match(rows, ["A", "B"])
        self.assertEqual((exact["rank"], exact["tie_best_rank"], exact["tie_worst_rank"]), (2, 1, 2))
        self.assertEqual(exact["hit_at_1"], 0)
        self.assertEqual(equivalent["hit_at_1"], 1)

    def test_oracle_budget_is_explicitly_larger(self):
        from iri_constructs.session_benchmark import evaluate_oracle_query
        ranking, query, eq, _, _, _ = self.fixture()
        lists = {a: ranking for a in query["components"]}
        ev = evaluate_oracle_query(lists, query, eq)
        audit.validate_evaluation(ev, query, lists, eq)
        self.assertEqual(ev["returned_hits_ceiling"]["5"], 20)
        self.assertEqual(ev["exact_recall_ceiling"]["1"], 1.)
        ev["retrieval_budget"] = "k total per original query"
        with self.assertRaises(ValueError): audit.validate_evaluation(ev, query, lists, eq)


if __name__ == "__main__": unittest.main()
