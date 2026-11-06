"""Small mock/interface fixtures: no real alignments, models or child launches."""
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import session_fragment_replay_audit as audit


def sample(qid="q", policy="exact_removed", suite="pairs", length=80, fingerprint=None, mask="m"):
    return {"policy": policy, "suite": suite, "query_id": qid, "start0": 10, "end0": 10 + length,
            "length_aa": length, "fragment_sequence_sha256": fingerprint or qid,
            "included_accessions_sha256": mask, "expected_attempt": {"best_eligible_score": 100}}


def hit(a, score=60., span=None, coverage=.5, identity=.25):
    return {"accession": a, "score": score, "query_span": [0, 32] if span is None else span,
            "reference_coverage": coverage, "identity_fraction": identity}


PARAMS = {"minimum_span": 32, "minimum_score": 60., "minimum_reference_coverage": .5,
          "minimum_identity_fraction": .25, "maximum_segments": 4, "empty_prediction_fallback": "whole_query"}


def case_and_rows():
    rows = audit.rank_rows([hit("B", 70), hit("A", 70), hit("C", 90, coverage=.1)])
    case = {"case_id": "case-x", "start0": 20, "end0": 100,
            "common_mask": {"included_accessions": ["A", "B", "C"]},
            "expected_attempt": {"start0": 20, "end0": 100, "catalogue_count": 3, "eligible_reference_count": 2,
                                 "best_eligible_accession": "A", "best_eligible_score": 70},
            "expected_accepted_segment": {**audit.independent_candidates(rows, (20, 100), PARAMS)[0], "selection_order": 1}}
    return case, rows


class SelectionTests(unittest.TestCase):
    def test_strata_stress_and_outcome_blind_selection(self):
        pool = [sample(p + s + str(i), p, s, 80 + i) for p in audit.POLICIES for s in audit.SUITES for i in range(3)]
        cases, evidence = audit.select_cases(pool)
        self.assertLessEqual(len(cases), 7)
        self.assertEqual(sum("stratum" in c["selection_roles"] for c in cases), 6)
        self.assertEqual(sum("longest_fragment" in c["selection_roles"] for c in cases), 1)
        for value in pool:
            value["expected_attempt"]["best_eligible_score"] = -999
        replay, _ = audit.select_cases(list(reversed(pool)))
        self.assertEqual([c["case_id"] for c in cases], [c["case_id"] for c in replay])
        self.assertEqual(evidence["fresh_attempts"], 18)

    def test_dedup_uses_sequence_and_full_mask(self):
        left, right = sample("a", fingerprint="same", length=100), sample("b", fingerprint="same", length=100)
        other = sample("c", fingerprint="same", mask="different")
        cases, counts = audit.select_cases([right, other, left])
        self.assertEqual(counts["distinct_fragment_mask_cases"], 2)
        representative = min((left, right), key=audit.selection_hash)
        self.assertIn(representative["query_id"], [c["query_id"] for c in cases])

    def test_empty_strata_are_not_replaced_and_longest_not_backfilled(self):
        cases, evidence = audit.select_cases([sample()])
        self.assertEqual(len(cases), 1)
        self.assertTrue(evidence["longest_already_selected"])
        self.assertEqual(sum(s["status"] == "unavailable_no_replacement" for s in evidence["strata"]), 5)


class ScientificComparisonTests(unittest.TestCase):
    def test_ties_and_near_ties_have_distinct_semantics(self):
        rows = audit.rank_rows([hit("Z", 60), hit("B", 60), hit("A", 60 - 5e-7)])
        self.assertEqual([r["accession"] for r in rows], ["B", "Z", "A"])
        self.assertEqual(rows[0]["exact_score_tie_count"], 2)
        self.assertEqual(rows[0]["near_tie_count_1e6"], 3)
        self.assertEqual(rows[1]["tie_best_rank"], 1)
        self.assertEqual(rows[1]["tie_worst_rank"], 2)

    def test_thresholds_are_inclusive_and_coordinates_translate(self):
        rows = audit.rank_rows([hit("A"), hit("B", 59.5), hit("C", 100, span=[0, 31]),
                                hit("D", 100, coverage=.49), hit("E", 100, identity=.249)])
        candidates = audit.independent_candidates(rows, (50, 130), PARAMS)
        self.assertEqual(len(candidates), 1)
        self.assertEqual((candidates[0]["start0"], candidates[0]["end0"]), (50, 82))

    def test_complete_masked_ranking_and_saved_seed_match(self):
        case, rows = case_and_rows()
        result = audit.compare_case(case, rows, PARAMS)
        self.assertEqual(result["status"], "passed")
        self.assertTrue(result["accepted_seed_geometry_checked"])
        self.assertFalse(result["old_full_ranking_available"])
        self.assertEqual(result["best_eligible_candidate"]["seed_hit"]["rank"], 2)

    def test_changed_mask_rank_score_and_seed_are_rejected(self):
        case, rows = case_and_rows()
        with self.assertRaises(ValueError):
            audit.compare_case(case, rows[:-1], PARAMS)
        bad = deepcopy(rows); bad[0]["rank"] = 9
        with self.assertRaises(ValueError):
            audit.compare_case(case, bad, PARAMS)
        bad_case = deepcopy(case); bad_case["expected_attempt"]["best_eligible_score"] += 1e-9
        self.assertEqual(audit.compare_case(bad_case, rows, PARAMS)["status"], "failed")
        bad_case = deepcopy(case); bad_case["expected_accepted_segment"]["start0"] += 1
        self.assertEqual(audit.compare_case(bad_case, rows, PARAMS)["status"], "failed")

    def test_no_eligible_reference_remains_null(self):
        case, rows = case_and_rows()
        rows = audit.rank_rows([hit("A", 2), hit("B", 1), hit("C", 0)])
        case["expected_attempt"].update(eligible_reference_count=0, best_eligible_accession=None, best_eligible_score=None)
        case["expected_accepted_segment"] = None
        result = audit.compare_case(case, rows, PARAMS)
        self.assertEqual(result["status"], "passed")
        self.assertIsNone(result["best_eligible_candidate"])

    def test_nonfinite_or_almost_equal_score_fails(self):
        with self.assertRaises(ValueError):
            audit.rank_rows([hit("A", float("nan"))])
        with self.assertRaises(ValueError):
            audit.exact(60 + 1e-12, 60)


class OfficialCountInterfaceTests(unittest.TestCase):
    def test_official_counts_reconstruct_gap_metrics(self):
        counts = SimpleNamespace(aligned=4, gaps=2, open_gaps=1, identities=3, substitution_score=24., score=13.5)
        coords = SimpleNamespace(tolist=lambda: [[1, 3, 5, 7], [0, 2, 2, 4]])
        alignment = SimpleNamespace(counts=lambda aligner: counts, coordinates=coords, shape=(2, 6))
        class Alignments:
            score = 13.5
            def __iter__(self): return iter([alignment])
        aligner = SimpleNamespace(align=lambda query, reference: Alignments())
        result = audit.independent_pair(aligner, "A" * 8, "A" * 5)
        self.assertEqual(result["aligned_length"], 6)
        self.assertEqual(result["query_coverage"], .75)
        self.assertEqual(result["reference_coverage"], .8)
        self.assertEqual(result["identity_fraction"], .5)
        self.assertEqual(result["query_paired_coverage"], .5)
        self.assertEqual(result["reconstructed_score"], 13.5)

    def test_zero_score_does_not_request_nonexistent_coordinates(self):
        class Zero:
            score = 0.
            def __iter__(self): raise AssertionError("Zero score should not enumerate alignments")
        result = audit.independent_pair(SimpleNamespace(align=lambda q, r: Zero()), "AA", "WW")
        self.assertIsNone(result["query_span"])
        self.assertEqual(result["coordinates"], [[], []])
        self.assertNotIn("reconstructed_score", result)


class OwnershipAccountingTests(unittest.TestCase):
    def test_missing_guard_is_rejected_before_execution(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(ValueError, "root-owned"):
                audit.guard_admission()

    def test_guard_rejects_wrong_hash_and_accepts_only_recorded_group(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "guard.json"
            record = {"token": "test-token", "guard_source_sha256": "bad", "wrapper_pid": 10,
                      "child_groups": [{"pid": 123, "pgid": 123, "start_new_session": True}]}
            path.write_text(json.dumps(record))
            with patch.dict(os.environ, {"IRI_GUARD_OWNERSHIP_PATH": str(path), "IRI_GUARD_OWNERSHIP_TOKEN": "test-token"}), patch.object(audit.os, "getpgrp", return_value=123):
                with self.assertRaises(ValueError): audit.guard_admission()
                record["guard_source_sha256"] = audit.GUARD_SHA
                path.write_text(json.dumps(record))
                self.assertEqual(audit.guard_admission()["owned_pgid"], 123)

    def test_failed_and_unattempted_cases_do_not_become_passed_coverage(self):
        cases = [{"case_id": str(i)} for i in range(3)]
        records = [{"case_id": str(i), "status": value} for i, value in enumerate(("passed", "failed", "unattempted"))]
        result = audit.accounting(cases, records, {"fresh_attempts": 1931, "distinct_fragment_mask_cases": 1200,
                                                  "original_full_panel_fresh_attempts": 1933,
                                                  "excluded_clean_sanity_fresh_attempts": 2})
        self.assertEqual(result["status"], "incomplete_or_failed")
        self.assertEqual(result["attempted_cases"], 2)
        self.assertEqual(result["unreplayed_distinct_cases"], 1199)
        self.assertEqual(result["original_fresh_attempts"], 1933)
        self.assertEqual(result["eligible_combined_fresh_attempts"], 1931)
        self.assertFalse(result["all_1933_attempts_validated"])
        with self.assertRaises(ValueError): audit.accounting(cases, records[:-1], {})

    def bounded_child(self, *, elapsed, observed_rss):
        class Child:
            pid, returncode, killed, reaped = 333, None, False, False
            def poll(self): return self.returncode
            def kill(self): self.killed = True; self.returncode = -9
            def wait(self, timeout): self.reaped = True; return self.returncode
        child = Child()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            design = root / "inputs/design.json"
            design.parent.mkdir(); design.write_text("{}")
            with patch.object(audit.subprocess, "Popen", return_value=child), \
                 patch.object(audit.os, "getpgid", return_value=os.getpgrp()), \
                 patch.object(audit.time, "monotonic", side_effect=[0, elapsed, elapsed]), \
                 patch.object(audit, "rss_sum", return_value=observed_rss):
                result = audit.execute_case({"case_id": "case-bounded"}, design, SimpleNamespace(root=root / "job"),
                    {"owned_pgid": os.getpgrp()}, SimpleNamespace(check=lambda: None))
            self.assertTrue(child.killed and child.reaped)
            self.assertTrue(result["child_reaped"])
            self.assertEqual(result["status"], "failed")
            return result

    def test_case_timeout_reaps_only_its_owned_child(self):
        result = self.bounded_child(elapsed=61, observed_rss=100)
        self.assertEqual(result["reason"], "case_wall_time_exceeded")

    def test_memory_overrun_reaps_owned_child(self):
        result = self.bounded_child(elapsed=1, observed_rss=513 * 1024**2)
        self.assertEqual(result["reason"], "combined_memory_bound_exceeded")


if __name__ == "__main__":
    unittest.main()
