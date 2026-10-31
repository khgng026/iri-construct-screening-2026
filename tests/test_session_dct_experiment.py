"""Tensor/geometry/provenance fixtures; no model, checkpoint or C++ execution."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from iri_constructs.session_dct_native import (parse_domains, boundary_metrics, window_lengths,
    ce_evidence, exhaustive_distances, load_module, isolated_working_directory)
from iri_constructs.session_store import UnitStore
import session_dct_experiment as runner


class DomainTests(unittest.TestCase):
    def test_discontinuous_order_partition_and_appended_global(self):
        result = parse_domains(["65-96,1-32", "33-64", "1-96"], 96)
        self.assertEqual(result[0]["segments"], [{"start0": 64, "end0": 96}, {"start0": 0, "end0": 32}])
        self.assertEqual([d["appended_whole_fingerprint"] for d in result], [False, False, True])
        self.assertFalse(parse_domains(["1-96"], 96)[0]["appended_whole_fingerprint"])

    def test_invalid_endpoints_overlap_and_partition_rejected(self):
        for domains in (["0-96"], ["1-97"], ["1-95"], ["1-32,20-96"],
                        ["1-64", "60-96", "1-96"], ["1-32", "34-96", "1-96"], ["96-1"], ["1-96", "1-96"]):
            with self.subTest(domains=domains), self.assertRaises(ValueError):
                parse_domains(domains, 96)

    def test_boundary_iou_uses_sets_and_excludes_appended_whole(self):
        domains = parse_domains(["1-32,65-96", "33-64", "1-96"], 96)
        metrics = boundary_metrics(domains, [{"accession": "numeric_target", "start0": 32, "end0": 64}])
        self.assertEqual(metrics["predicted_domain_count"], 2)
        comparisons = metrics["components"][0]["domain_comparisons"]
        self.assertEqual(comparisons[0]["intersection"], 0)
        self.assertEqual(comparisons[0]["iou"], 0)
        self.assertEqual(comparisons[1]["iou"], 1)

    def test_original_window_boundary_policy(self):
        for length, expected in ((500, [500]), (501, [500, 201]), (800, [500, 500]),
                                 (801, [500, 500, 201]), (1011, [500, 500, 411])):
            self.assertEqual(window_lengths(length), expected)


class NumericFingerprintTests(unittest.TestCase):
    def test_integer_graph_rounding_and_required_short_edges(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "numeric.ce"
            path.write_text("CON   0 5 0.004999,0 9 0.005000,3 8 0.995000\n")
            record, arrays = ce_evidence(path, 10, save_graph=True)
            self.assertEqual(arrays["ce_graph_weights"].tolist(), [0, 1, 100])
            self.assertEqual(arrays["integer_graph"][0, 1:5].tolist(), [100, 100, 100, 0])
            self.assertEqual(record["selected_contact_count"], 3)

    def test_widened_raw_l1_preserves_large_distances_and_domain_pairs(self):
        import numpy as np
        query = np.stack([np.zeros(480, dtype=np.int8), np.full(480, 127, dtype=np.int8)])
        refs = {"A": np.zeros((1, 480), dtype=np.int8), "B": np.full((1, 480), 64, dtype=np.int8)}
        minima, pairs = exhaustive_distances(query, refs)
        self.assertEqual(minima, {"A": 0, "B": 63 * 480})
        self.assertEqual(int(pairs[:, 3].max()), 127 * 480)
        self.assertEqual(pairs.shape, (4, 4))
        ranking = runner.ranking_from_distances(minima)
        self.assertEqual([r["accession"] for r in ranking], ["A", "B"])
        self.assertLess(ranking[1]["unclipped_display_similarity"], 0)

    def test_exact_distance_tie_uses_accession_and_invalid_float_rejected(self):
        import numpy as np
        zeros = np.zeros((1, 480), dtype=np.int16)
        minima, pairs = exhaustive_distances(zeros, {"Z": zeros, "A": zeros})
        self.assertEqual([r["accession"] for r in runner.ranking_from_distances(minima)], ["A", "Z"])
        with self.assertRaises(ValueError):
            exhaustive_distances(zeros.astype(float), {"A": zeros})

    def test_global_control_uses_only_existing_literal_whole(self):
        import numpy as np
        values = {"fingerprints": np.arange(3 * 480).reshape(3, 480)}
        record = {"length": 96, "domains": parse_domains(["1-32,65-96", "33-64", "1-96"], 96)}
        self.assertTrue(np.array_equal(runner.whole_fingerprint(record, values), values["fingerprints"][[2]]))
        missing = {"length": 96, "domains": parse_domains(["1-32,33-96"], 96)}
        self.assertIsNone(runner.whole_fingerprint(missing, values))

    def test_original_zero_quantization_raises_before_int8(self):
        if importlib.util.find_spec("scipy") is None:
            self.skipTest("SciPy is absent in this optional environment; run this original-quantizer assertion in .venv-esm-readiness (SciPy installed).")
        import numpy as np
        path = ROOT / "results/iri-session-20261031-032506/vendor/DCTdomain/src/fingerprint.py"
        if not path.exists():
            self.skipTest("Pinned optional source not installed")
        module = load_module(path, "_test_original_dct_fingerprint")
        fp = module.Fingerprint(embed={15: np.zeros((22, 640)), 21: np.zeros((22, 640))}, domains=["1-22"])
        with self.assertRaises(FloatingPointError), np.errstate(invalid="raise", divide="raise", over="raise"):
            fp.quantize([3, 80, 3, 80])


class ProvenanceTests(unittest.TestCase):
    def test_complete_failure_unit_remains_visible_and_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            session = Path(directory)
            store = UnitStore(session / "source", {"unit_keys": ["query-q", "summary"]})
            store.publish("query-q", {"status": "failed", "exception_type": "FloatingPointError"})
            with self.assertRaises(ValueError):
                runner.open_completed(store.root, session)
            store.publish("summary", {"failed_units": ["query-q"]})
            store.complete(["query-q", "summary"])
            source = runner.open_completed(store.root, session)
            self.assertEqual(runner.completed_unit(source, "query-q")["status"], "failed")
            (source.completed_paths["query-q"] / "record.json").write_text("{}")
            with self.assertRaises(ValueError):
                runner.open_completed(store.root, session)

    def test_failed_discontinuous_fixture_cannot_pass_build_gate(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gate.json"
            path.write_text(json.dumps({"status": "passed", "fixtures": {"continuous": True, "discontinuous": False}}))
            with self.assertRaisesRegex(ValueError, "discontinuous"):
                runner.build_gate(path, {}, ROOT)

    def test_temporary_raw_directory_removed_and_cwd_restored_on_failure(self):
        old = Path.cwd()
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            with self.assertRaises(RuntimeError):
                with isolated_working_directory(base) as temporary:
                    (temporary / "owned.ce").write_text("numeric fixture")
                    raise RuntimeError("fixture")
            self.assertEqual(Path.cwd(), old)
            self.assertEqual(list(base.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
