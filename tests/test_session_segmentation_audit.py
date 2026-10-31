"""Independent audit fixtures; stdlib only, no production metrics or model code."""
from collections import Counter
from copy import deepcopy
import hashlib
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from types import SimpleNamespace
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import session_segmentation_audit as audit


def ranking(values):
    values = sorted(values, key=lambda row: (-row[1], row[0]))
    counts = Counter(v for _, v in values)
    return [{"accession": a, "rank": rank, "score": v, "cosine_similarity": v, "exact_score_tie_count": counts[v]}
            for rank, (a, v) in enumerate(values, 1)]


def query(names, qid="q"):
    return {"query_id": qid, "group_id": "g-" + qid, "suite": "pairs", "component_count": len(names), "components": names,
            "length_aa": 10 * len(names), "spans": [{"accession": a, "start0": i * 10, "end0": (i + 1) * 10, "fraction": 1 / len(names)} for i, a in enumerate(names)]}


class GeometryBudgetTests(unittest.TestCase):
    def test_windows_cover_short_input_and_both_termini(self):
        short = audit.expected_windows(50, {"width": 128, "stride": 64})
        self.assertEqual(short["segments"], [{"start0": 0, "end0": 50}])
        result = audit.expected_windows(300, {"width": 128, "stride": 64})
        self.assertEqual([s["start0"] for s in result["segments"]], [0, 64, 128, 172])
        self.assertEqual(result["segments"][-1]["end0"], 300)

    def test_interval_union_does_not_double_count_overlap(self):
        segments = [{"start0": 0, "end0": 12}, {"start0": 8, "end0": 20}]
        target = [{"accession": "A", "start0": 0, "end0": 20}]
        result = audit.interval_metrics(segments, target, 20)["components"][0]
        self.assertEqual(result["covered_fraction_by_union"], 1)
        self.assertEqual(result["best_interval_iou"], .6)
        self.assertFalse(result["exact_boundary_match"])
        with self.assertRaises(ValueError):
            audit.interval_metrics([{"start0": 0, "end0": 21}], target, 20)

    def test_union_counts_duplicates_and_explicit_larger_budget(self):
        lists = [ranking([("A", 1.), ("B", .5), ("X", .1)]), ranking([("A", 1.), ("X", .5), ("B", .1)])]
        result = audit.union_result(lists, 1, query(["A", "B"]))
        self.assertEqual(result["returned_hit_occurrences"], 2)
        self.assertEqual(result["unique_returned_hits"], 1)
        self.assertEqual(result["duplicate_returned_hits"], 1)
        self.assertEqual(result["evaluation"]["component_recall"], .5)
        self.assertEqual(result["retrieval_budget"], "k per predicted segment")

    def test_family_union_unscorable_query_does_not_become_a_miss(self):
        value = audit.union_result([ranking([("X", .9), ("Y", .8)])], 1, query(["A", "B"]), {"A": ["A", "X"], "B": ["B"]}, ["X", "Y"])
        self.assertEqual(value["evaluation"]["scorable_components"], 1)
        self.assertEqual(value["evaluation"]["component_recall"], 1)
        self.assertIsNone(value["evaluation"]["all_components_hit"])
        acc = audit.UnionRollup(); acc.add("q", value)
        result = acc.result()
        self.assertEqual(result["all_components_scorable_queries"], 0)
        self.assertEqual(result["excluded_all_component_queries"], 1)
        self.assertIsNone(result["all_components_hit"])
        self.assertEqual(result["per_target_macro_component_recall"], 1)
        bad = deepcopy(result); bad["all_components_hit"] = 0
        with self.assertRaises(ValueError):
            audit.same(bad, result)

    def test_mask_before_rrf_changes_scores_and_can_change_order(self):
        lists = [ranking([("X", 1.), ("A", .8), ("B", .7)]), ranking([("B", 1.), ("X", .8), ("A", .7)])]
        masked = [audit.independently_ranked(rows, ["A", "B"]) for rows in lists]
        correct = audit.base.fused_ranking(masked, "rrf", ["A", "B"])
        wrong = audit.independently_ranked(audit.base.fused_ranking(lists, "rrf", ["A", "B", "X"]), ["A", "B"])
        self.assertEqual(correct[0]["exact_score_tie_count"], 2)
        self.assertEqual([row["accession"] for row in correct], ["A", "B"])
        self.assertEqual([row["accession"] for row in wrong], ["B", "A"])
        self.assertNotEqual(correct[0]["score"], wrong[0]["score"])
        with self.assertRaises(ValueError):
            audit.same(wrong, correct)

    def test_boundary_informed_or_wrong_window_proposal_rejected(self):
        q, config = query(["A", "B"]), {"configuration_id": "w", "segmenter": "fixed_windows", "width": 12, "stride": 6, "fusion": "max_score", "rrf_constant": 60}
        prediction = audit.expected_windows(q["length_aa"], config)
        encoded = prediction["total_queried_residues"]
        record = {"query": q, "prediction_config": audit.geometry(config), "prediction": prediction,
                  "logical_candidate_scan_residues": 0, "logical_encoded_segment_residues": encoded, "logical_total_queried_residues": encoded,
                  "boundary_evaluation": audit.interval_metrics(prediction["segments"], q["spans"], q["length_aa"])}
        audit.audit_proposal(record, q, config, "A" * 20, ["A", "B"], None, {}, Counter())
        record["prediction"]["boundaries_supplied"] = True
        with self.assertRaises(ValueError):
            audit.audit_proposal(record, q, config, "A" * 20, ["A", "B"], None, {}, Counter())

    def test_masked_seed_recomputes_near_ties_and_legacy_alias(self):
        rows = ranking([("A", 1.), ("B", 1.), ("C", 1. - 5e-7)])
        for row in rows:
            row["near_tie_count_1e6"] = 3
        result = audit.independently_ranked(rows, ["B", "C"], diagnostic_ties=True)
        self.assertEqual(result[0]["rank"], 1)
        self.assertEqual(result[0]["exact_score_tie_count"], 1)
        self.assertEqual(result[0]["tie_best_rank"], 1)
        self.assertEqual(result[0]["tie_worst_rank"], 1)
        self.assertEqual(result[0]["near_tie_count"], 2)
        self.assertEqual(result[0]["near_tie_count_1e6"], 2)
        self.assertEqual(result[0]["near_tie_worst_rank"], 2)


class SelectionCoverageTests(unittest.TestCase):
    def fixture(self):
        configs = [{"configuration_id": "a", "segmenter": "fixed_windows"}, {"configuration_id": "b", "segmenter": "fixed_windows"}]
        design = {"candidate_configurations": configs, "selection_query_ids": ["q1", "q2", "q3"],
                  "partition": {"queries": {"q1": "dev", "q2": "dev", "q3": "dev"}}}
        ctx = {"queries": {"q1": {"group_id": "large"}, "q2": {"group_id": "large"}, "q3": {"group_id": "small"}}}
        rows = {"a": {"q1": (1., 1, 10), "q2": (1., 1, 10), "q3": (0., 0, 10)},
                "b": {"q1": (.6, 0, 10), "q2": (.6, 0, 10), "q3": (.6, 0, 10)}}
        return design, ctx, rows

    def test_selection_is_group_macro_not_query_mean(self):
        design, ctx, rows = self.fixture()
        result = audit.selection_from_outcomes(rows, design, ctx)
        self.assertEqual(result["selected_configuration_id"], "b")
        self.assertEqual(result["development_scores"][1]["group_macro_component_recall_at_5"], .5)

    def test_selection_rejects_missing_candidate_and_test_leakage(self):
        design, ctx, rows = self.fixture()
        with self.assertRaises(ValueError):
            audit.selection_from_outcomes({"a": rows["a"]}, design, ctx)
        design["partition"]["queries"]["q3"] = "test"
        with self.assertRaises(ValueError):
            audit.selection_from_outcomes(rows, design, ctx)

    def test_selection_tie_breaks_work_then_configuration_id(self):
        design, ctx, rows = self.fixture()
        rows["b"] = deepcopy(rows["a"])
        self.assertEqual(audit.selection_from_outcomes(rows, design, ctx)["selected_configuration_id"], "a")
        rows["b"] = {q: (v[0], v[1], 9) for q, v in rows["b"].items()}
        self.assertEqual(audit.selection_from_outcomes(rows, design, ctx)["selected_configuration_id"], "b")

    def test_complete_unit_schema_includes_all_queries_policies_configs(self):
        identity = {"mode": "family_automatic_retrieval_v1", "query_ids": ["q1", "q2"],
                    "configs": [{"configuration_id": "w"}, {"configuration_id": "a"}], "policies": ["exact_removed", "near90_removed"]}
        keys = audit.expected_keys(identity)
        self.assertEqual(len(keys), 9)
        self.assertIn("q2--near90_removed--w", keys)
        identity["policies"] = ["exact_removed", "exact_removed"]
        with self.assertRaises(ValueError):
            audit.expected_keys(identity)

    def test_common_mask_has_all_targets_and_stable_digest(self):
        ctx = {"hashes": {a: a for a in "ABXY"}, "removals": {"exact_removed": {"A": ["A", "X"], "B": ["B"]}}}
        result = audit.common_mask(ctx, query(["A", "B"]), "exact_removed")
        self.assertEqual(result["included_accessions"], ["Y"])
        self.assertEqual(result["removed_accessions"], ["A", "B", "X"])
        self.assertEqual(result["remaining_catalogue_count"], 1)
        self.assertEqual(result["included_accessions_sha256"], audit.canonical_hash(["Y"]))


class SourceNumericsTests(unittest.TestCase):
    def vector_file(self, path, values, descr="<f4"):
        header = repr({"descr": descr, "fortran_order": False, "shape": (len(values),)}).encode() + b"\n"
        payload = b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header + struct.pack("<" + "f" * len(values), *values)
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("mean.npy", payload)

    def test_saved_vector_reader_and_independent_cosine(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "arrays.npz"
            self.vector_file(path, [3., 4.])
            self.assertEqual(audit.float32_mean(path), (3., 4.))
            self.assertAlmostEqual(audit.cosine(audit.float32_mean(path), [0., 1.]), .8)
            with self.assertRaises(ValueError):
                audit.cosine([0., 0.], [1., 0.])

    def test_saved_vector_reader_rejects_dtype_and_nonfinite_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "arrays.npz"
            self.vector_file(path, [1.], descr="<f2")
            with self.assertRaises(ValueError):
                audit.float32_mean(path)
            self.vector_file(path, [float("nan")])
            with self.assertRaises(ValueError):
                audit.float32_mean(path)

    def trace_fixture(self):
        sequence = "ACDEFGHI"
        source_query = {"query_id": "q", "joined_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}
        raw_rows = ranking([("A", 100.), ("B", 80.)])
        for row in raw_rows:
            row.update(query_span=[0, 8], reference_coverage=1., identity_fraction=1.)
        seed = audit.independently_ranked(raw_rows, ["B"], diagnostic_ties=True)[0]
        params = {"minimum_span": 4, "minimum_reference_coverage": .5, "minimum_identity_fraction": .25,
                  "minimum_score": 60., "maximum_segments": 4, "empty_prediction_fallback": "whole_query"}
        prediction = {"method": "iterative_alignment_candidates", "config": params, "query_length": 8, "boundaries_supplied": False,
            "segments": [{"start0": 0, "end0": 8, "source_interval": [0, 8], "seed_hit": seed, "selection_order": 0}],
            "accepted_alignment_segment_count": 1, "whole_query_fallback_used": False, "segment_count": 1,
            "interval_attempts": [{"start0": 0, "end0": 8, "catalogue_count": 1, "eligible_reference_count": 1,
                                   "best_eligible_accession": "B", "best_eligible_score": 80.}],
            "total_queried_residues": 8, "alignment_comparisons": 1, "intervals_aligned": 1,
            "discarded_intervals": [], "unprocessed_intervals": [], "cap_bound": False, "total_predicted_segment_residues": 8}
        raw = SimpleNamespace(record=lambda key: {"query": source_query, "rankings": {"local_alignment": raw_rows}})
        ctx = {"queries": {"q": source_query}, "hashes": {"A": "a", "B": "b"}}
        return prediction, sequence, {"parameters": params}, raw, ctx

    def test_cached_fragment_is_masked_before_prediction(self):
        prediction, sequence, config, raw, ctx = self.trace_fixture()
        counters = Counter()
        audit.audit_alignment_trace(prediction, sequence, config, ["B"], raw, ctx, counters)
        self.assertEqual(counters["attempts_replayed_from_complete_original_rankings"], 1)
        prediction["segments"][0]["seed_hit"]["accession"] = "A"
        with self.assertRaisesRegex(ValueError, "Removed reference"):
            audit.audit_alignment_trace(prediction, sequence, config, ["B"], raw, ctx, Counter())

    def test_cached_fragment_attempt_tampering_rejected(self):
        prediction, sequence, config, raw, ctx = self.trace_fixture()
        prediction["interval_attempts"][0]["eligible_reference_count"] = 0
        with self.assertRaises(ValueError):
            audit.audit_alignment_trace(prediction, sequence, config, ["B"], raw, ctx, Counter())

    def test_selection_order_is_zero_based_and_shift_rejected(self):
        prediction, sequence, config, raw, ctx = self.trace_fixture()
        audit.audit_alignment_trace(prediction, sequence, config, ["B"], raw, ctx, Counter())
        prediction["segments"][0]["selection_order"] = 1
        with self.assertRaisesRegex(ValueError, "seed order"):
            audit.audit_alignment_trace(prediction, sequence, config, ["B"], raw, ctx, Counter())


if __name__ == "__main__":
    unittest.main()
