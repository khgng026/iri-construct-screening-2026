"""Numeric/interface fixtures only; no protein model, download or inference."""
import copy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
import session_context_profile as profile


class NumericProfileTests(unittest.TestCase):
    def setUp(self):
        import numpy as np
        self.np = np

    def fixture(self, clean=False):
        np = self.np
        isolated = {"A": (np.arange(450).reshape(150, 3) / 100 + 1).astype(np.float32),
                    "B": (np.arange(450).reshape(150, 3) / 200 + 2).astype(np.float32)}
        names = ["A"] if clean else ["A", "B"]
        states = np.concatenate([isolated[a] for a in names]).copy()
        if not clean:
            states[140:160] += np.array([0.3, -0.1, 0.2], dtype=np.float32)
        query = {"query_id": "cleanA" if clean else "A_B", "group_id": "A" if clean else "A__B",
                 "suite": "clean" if clean else "pairs", "component_count": len(names),
                 "components": names, "length_aa": len(states), "sequence": "unused numeric fixture",
                 "spans": [{"accession": a, "start0": i * 150, "end0": (i + 1) * 150}
                           for i, a in enumerate(names)]}
        parent = {"whole": profile.mean_residues(states)}
        for span in query["spans"]:
            a = span["accession"]
            parent["contextual_" + a] = profile.mean_residues(states[span["start0"]:span["end0"]])
            parent["isolated_" + a] = profile.mean_residues(isolated[a])
        return query, states, isolated, parent

    def test_terminal_distance_uses_only_artificial_join(self):
        left = profile.distance_bins(140, True, False)
        right = profile.distance_bins(140, False, True)
        self.assertEqual(left[[0, 7, 8, 31, 32, 63, 64, 127, 128, 139]].tolist(), [0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        self.assertEqual(right.tolist(), left[::-1].tolist())
        self.assertEqual(int(right[0]), 4)  # biological N terminus is far from its only join

    def test_internal_nearest_join_and_clean_separation(self):
        internal = profile.distance_bins(300, True, True)
        self.assertEqual(internal[[0, 299, 149, 150]].tolist(), [0, 0, 4, 4])
        self.assertEqual(profile.distance_bins(140, False, False).tolist(), [-1] * 140)
        with self.assertRaises(ValueError):
            profile.distance_bins(0, True, False)

    def test_hand_calculated_residue_metrics_and_undefined_case(self):
        np = self.np
        cosine, relative = profile.residue_changes([[0, 1], [0, 1]], [[1, 0], [0, 2]])
        np.testing.assert_allclose(cosine, [0, 1])
        np.testing.assert_allclose(relative, [np.sqrt(2), 0.5])
        with self.assertRaises(ValueError):
            profile.residue_changes([[1, 1]], [[0, 0]])

    def test_bin_vectors_reconstruct_parent_and_no_joined_matrix_saved(self):
        query, states, isolated, parent = self.fixture()
        record, arrays = profile.profile_query(query, states, isolated, parent)
        self.assertEqual(len(record["bins"]), 10)
        self.assertEqual(sum(r["residue_count"] for r in record["bins"]), 300)
        self.assertEqual(set(arrays), {"contextual_bin_means", "isolated_bin_means"})
        self.assertEqual(arrays["contextual_bin_means"].shape, (10, 3))
        self.assertLessEqual(record["checks"]["binned_whole_parent_max_abs"], 1e-6)
        changed = dict(parent, contextual_A=parent["contextual_A"] + 0.01)
        with self.assertRaises(ValueError):
            profile.profile_query(query, states, isolated, changed)

    def test_cached_bins_revalidate_geometry_means_and_checks(self):
        query, states, isolated, parent = self.fixture()
        record, arrays = profile.profile_query(query, states, isolated, parent)
        profile.validate_profile_unit(copy.deepcopy(record), arrays, query, isolated, parent)
        wrong = copy.deepcopy(record)
        wrong["bins"][0]["residue_count"] += 1
        with self.assertRaises(ValueError):
            profile.validate_profile_unit(wrong, arrays, query, isolated, parent)
        wrong = copy.deepcopy(record)
        wrong["checks"]["binned_whole_parent_max_abs"] = 0.5
        with self.assertRaises(ValueError):
            profile.validate_profile_unit(wrong, arrays, query, isolated, parent)
        wrong_arrays = {k: v.copy() for k, v in arrays.items()}
        wrong_arrays["isolated_bin_means"][0] += 0.1
        with self.assertRaises(ValueError):
            profile.validate_profile_unit(record, wrong_arrays, query, isolated, parent)

    def test_clean_has_own_bin_and_exact_repeat(self):
        query, states, isolated, parent = self.fixture(clean=True)
        record, arrays = profile.profile_query(query, states, isolated, parent)
        self.assertEqual(record["scope"], "clean")
        self.assertEqual(record["bins"][0]["bin"], "clean_no_artificial_join")
        self.assertEqual(record["checks"]["clean_residue_repeat_max_abs"], 0)
        self.assertEqual(record["bins"][0]["relative_l2"]["mean"], 0)

    def test_summary_group_support_and_matched_contrasts(self):
        q, states, isolated, parent = self.fixture()
        record, _ = profile.profile_query(q, states, isolated, parent)
        q2 = dict(q, query_id="second_order_fixture")
        second, _ = profile.profile_query(q2, states, isolated, parent)
        clean_args = self.fixture(clean=True)
        clean, _ = profile.profile_query(*clean_args)
        summary = profile.summarize_profiles([record, second, clean])
        row = next(r for r in summary["by_scope_suite_bin"] if r["scope"] == "combined" and r["suite"] == "all" and r["bin"] == "0-7")
        self.assertEqual((row["queries"], row["groups"], row["target_count"], row["group_target_pairs"]), (2, 1, 2, 2))
        self.assertEqual(row["residue_observations"], 32)
        matched = summary["matched_far128plus_vs_near0to7_by_group_target"]
        self.assertEqual(len(matched), 2)
        self.assertEqual(matched[0]["matched_query_count"], 2)
        self.assertIn("never pooled", summary["quantile_scope"])

    def test_explicit_parent_and_repeat_selection(self):
        q, *_ = self.fixture()
        suite = {"queries": [q], "manifest": {"repeat_query_ids": [q["query_id"]]}}
        design = {"version": "context-profile-v1", "model_first": "first", "bins_distance_to_nearest_artificial_join_residues": [list(x) for x in profile.BINS]}
        parent = SimpleNamespace(root=Path("first"), identity={"mode": "model", "query_ids": [q["query_id"]]})
        self.assertEqual(profile.select_queries(suite, parent, design, "full"), [q])
        parent.root = Path("another_model")
        with self.assertRaises(ValueError):
            profile.select_queries(suite, parent, design, "full")
        self.assertEqual(profile.select_queries(suite, parent, design, "repeat"), [q])
        parent.identity["query_ids"] = []
        with self.assertRaises(ValueError):
            profile.select_queries(suite, parent, design, "repeat")


class PrecisionTests(unittest.TestCase):
    def fixture(self, parameter_dtype=None, raw_dtype=None):
        import torch
        class ToyModel(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.weight = torch.nn.Parameter(torch.ones(2, dtype=parameter_dtype or torch.float32))
            def forward(self):
                return SimpleNamespace(last_hidden_state=torch.ones((1, 5, 2), dtype=raw_dtype or torch.float32))
        model = ToyModel().eval()
        encoder = SimpleNamespace(torch=torch, model=model, layer=6, num_layers=6,
                    config={"family": "esm2", "device": "cpu", "dtype": "float32"})
        def encode(sequence):
            with torch.no_grad():
                raw = model().last_hidden_state
                return raw[0, 1:-1].float().numpy().copy()  # deliberately mimics potentially hiding a low dtype
        encoder.encode = encode
        return encoder

    def test_raw_fp32_evidence_captured_and_hook_removed(self):
        encoder = self.fixture()
        audit = profile.PrecisionAudit(encoder)
        states, raw = audit.encode("xyz")  # numeric fixture never tokenizes this label
        self.assertEqual(raw["dtype"], "torch.float32")
        profile.validate_precision_record(raw, 3, 2, "cpu")
        audit.close()
        self.assertFalse(encoder.model._forward_hooks)

    def test_parameter_and_precast_raw_dtype_mismatch_rejected(self):
        import torch
        with self.assertRaises(ValueError):
            profile.PrecisionAudit(self.fixture(parameter_dtype=torch.float16))
        encoder = self.fixture(raw_dtype=torch.float16)
        audit = profile.PrecisionAudit(encoder)
        try:
            with self.assertRaisesRegex(ValueError, "Raw forward output"):
                audit.encode("xyz")
        finally:
            audit.close()

    def test_declared_device_and_cached_precision_evidence_rejected(self):
        encoder = self.fixture()
        encoder.config["device"] = "mps"
        with self.assertRaises(ValueError):
            profile.PrecisionAudit(encoder)
        bad = {"dtype": "torch.float16", "device_type": "cpu", "shape": [1, 5, 2],
               "requires_grad": False, "model_training": False, "autocast_enabled": False, "adapter_return_dtype": "float32"}
        with self.assertRaises(ValueError):
            profile.validate_precision_record(bad, 3, 2, "cpu")


if __name__ == "__main__":
    unittest.main()
