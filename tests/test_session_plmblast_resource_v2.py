"""Resource amendment regressions: metadata/mocks only, no model/JIT/probe."""
import copy
from contextlib import contextmanager
import datetime as dt
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
import session_plmblast_experiment as runner

SESSION = ROOT / "results/iri-session-20261031-032506"
V1 = SESSION / "inputs/plmblast_design_v1.json"
V2 = SESSION / "inputs/plmblast_design_resource_v2.json"


def admission():
    return {
        "schema": "iri-plmblast-explicit-resource-admission-v2", "status": "passed",
        "decision": "admit_repeat32", "resource_gate": "passed", "decided_by": "root",
        "decision_utc": "2026-10-31T14:00:00+00:00", "priority_work_and_actual_resources_reviewed": True,
        "design_sha256": runner.file_hash(V2), "runner_sha256": runner.file_hash(runner.__file__),
        "probe_job": str(SESSION / "jobs/toy-new-resource-v2"), "probe_identity_sha256": "new-v2-identity",
        "probe_completion_sha256": "new-v2-completion", "original_v1_design_sha256": runner.ORIGINAL_DESIGN_SHA256,
        "original_v1_negative_probe": copy.deepcopy(runner.RESOURCE_AMENDMENT["original_negative_probe"]),
        "resource_reservation": {"session_hard_stop_utc": "2026-10-31T19:25:06+00:00",
                                 "alignment_seconds": 10800, "encoding_seconds": 1800, "cleanup_audit_seconds": 3600},
    }


class ResourceDesignTests(unittest.TestCase):
    def setUp(self):
        self.original, self.amended = runner.read(V1), runner.read(V2)

    def test_frozen_original_and_exact_amendment_accepted(self):
        self.assertEqual(runner.file_hash(V1), runner.ORIGINAL_DESIGN_SHA256)
        runner.validate_resource_design(self.original, self.original)
        runner.validate_resource_design(self.amended, self.original)
        self.assertEqual(self.amended["limits"]["repeat_alignment_seconds"], 10800)
        self.assertEqual(self.amended["status"], runner.RESOURCE_STATUS)

    def test_v1_cannot_be_relabelled_or_cap_changed(self):
        changed = copy.deepcopy(self.original)
        changed["limits"]["repeat_alignment_seconds"] = 10800
        with self.assertRaisesRegex(ValueError, "Original v1 design changed"):
            runner.validate_resource_design(changed, self.original)
        changed = copy.deepcopy(self.amended)
        changed["status"] = self.original["status"]
        with self.assertRaises(ValueError):
            runner.validate_resource_design(changed, self.original)

    def test_changed_examples_safety_scoring_runtime_or_extra_field_rejected(self):
        changes = [
            (("probe_query_ids",), list(reversed(self.amended["probe_query_ids"]))),
            (("repeat_query_ids",), self.amended["repeat_query_ids"][:-1]),
            (("probe_reference_accessions",), ["different"]),
            (("forecast", "safety_factor"), 1.0),
            (("forecast", "alignment_work_proxy"), "N*M"),
            (("extractor", "gap_penalty"), 0.4),
            (("postprocess", "alignment_cutoff"), 0.2),
            (("postprocess", "round_decimals"), 4),
            (("screening", "COS_PER_CUT"), 70),
            (("limits", "repeat_encoding_seconds"), 3600),
            (("limits", "repeat_alignment_seconds"), 10801),
            (("limits", "pair_seconds"), 120),
            (("limits", "maximum_worker_peak_mib"), 4096),
            (("limits", "numba_threads"), 2),
            (("model", "layer"), 29),
            (("extra_unreviewed_annotation",), True),
        ]
        for path, value in changes:
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, "exact cap-only amendment"):
                changed = copy.deepcopy(self.amended)
                target = changed
                for part in path[:-1]:
                    target = target[part]
                target[path[-1]] = value
                runner.validate_resource_design(changed, self.original)

    def test_identical_forecast_v1_fails_5400_v2_passes_10800(self):
        query_ids = sorted(set(self.original["probe_query_ids"] + self.original["repeat_query_ids"]))
        suite = {"queries": [{"query_id": q, "length_aa": 1} for q in query_ids],
                 "sequences": {a: "A" for a in self.original["probe_reference_accessions"]}}
        pairs = [{"query_id": q, "accession": a, "status": "eligible_alignment"}
                 for q, a in runner.pair_list(self.original["probe_query_ids"], self.original["probe_reference_accessions"])]
        timings = {runner.pair_key(p["query_id"], p["accession"]): 50.0 for p in pairs}
        with mock.patch.object(runner, "completed_unit", return_value={"times_seconds": [.001, .001, .001]}):
            old = runner.forecast_repeat(pairs, timings, [{"startup_seconds": 10.0}], None, suite, self.original)
            new = runner.forecast_repeat(pairs, timings, [{"startup_seconds": 10.0}], None, suite, self.amended)
        self.assertEqual(old["alignment_seconds"], 7215.0)
        self.assertFalse(old.pop("within_frozen_budget"))
        self.assertTrue(new.pop("within_frozen_budget"))
        self.assertEqual(old, new)
        self.assertFalse(runner.forecast_within_limits({"alignment_seconds": 10800.01, "encoding_seconds": 1}, self.amended))
        self.assertFalse(runner.forecast_within_limits({"alignment_seconds": 1, "encoding_seconds": 1800.01}, self.amended))

    def test_v2_requires_bytecode_environment_before_startup(self):
        with mock.patch.dict(runner.os.environ, {"PYTHONDONTWRITEBYTECODE": "0"}), mock.patch.object(sys, "dont_write_bytecode", False):
            with self.assertRaisesRegex(ValueError, "before interpreter startup"):
                runner.load_validated_design(SESSION, V2)
        with mock.patch.dict(runner.os.environ, {"PYTHONDONTWRITEBYTECODE": "1"}), mock.patch.object(sys, "dont_write_bytecode", True):
            self.assertEqual(runner.load_validated_design(SESSION, V2), self.amended)


class ExplicitReservationTests(unittest.TestCase):
    def test_full_encoding_and_alignment_reservations(self):
        gate = admission(); stop = gate["resource_reservation"]["session_hard_stop_utc"]
        encoding_start = dt.datetime.fromisoformat("2026-10-31T14:55:06+00:00")
        alignment_start = dt.datetime.fromisoformat("2026-10-31T15:25:06+00:00")
        runner.validate_resource_reservation(gate, stop, "encoding", encoding_start)
        runner.validate_resource_reservation(gate, stop, "alignment", alignment_start)
        for stage, now in (("encoding", encoding_start), ("alignment", alignment_start)):
            with self.subTest(stage=stage), self.assertRaisesRegex(ValueError, "Insufficient remaining"):
                runner.validate_resource_reservation(gate, stop, stage, now + dt.timedelta(microseconds=1))

    def test_root_review_exact_hardstop_and_allowance_mandatory(self):
        now = dt.datetime.fromisoformat("2026-10-31T14:30:00+00:00")
        for field, value in (("decided_by", "automatic"), ("priority_work_and_actual_resources_reviewed", False)):
            changed = admission(); changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                runner.validate_resource_reservation(changed, "2026-10-31T19:25:06+00:00", "encoding", now)
        changed = admission(); changed["resource_reservation"]["cleanup_audit_seconds"] = 3599
        with self.assertRaises(ValueError):
            runner.validate_resource_reservation(changed, "2026-10-31T19:25:06+00:00", "encoding", now)
        with self.assertRaises(ValueError):
            runner.validate_resource_reservation(admission(), "2026-10-31T20:25:06+00:00", "encoding", now)

    def test_late_original_decision_cannot_be_backfilled_at_alignment(self):
        gate = admission(); gate["decision_utc"] = "2026-10-31T15:00:00+00:00"
        with self.assertRaisesRegex(ValueError, "Original resource decision lacked"):
            runner.validate_resource_reservation(gate, "2026-10-31T19:25:06+00:00", "alignment",
                                                 dt.datetime.fromisoformat("2026-10-31T15:05:00+00:00"))

    def test_original_v1_probe_cannot_admit_v2(self):
        gate = admission(); old = runner.RESOURCE_AMENDMENT["original_negative_probe"]
        design = runner.read(V2)
        gate["probe_completion_sha256"] = old["completion_sha256"]
        probe = SimpleNamespace(identity={"mode": "plmblast_alignment_v1", "panel": "probe12",
                                          "design_sha256": runner.ORIGINAL_DESIGN_SHA256},
                                completion_sha256=old["completion_sha256"])
        with mock.patch.object(runner, "read", return_value=gate), mock.patch.object(runner, "open_completed", return_value=probe):
            with self.assertRaisesRegex(ValueError, "another probe/code/design"):
                runner.validate_admission(Path("gate.json"), SESSION, runner.file_hash(V2), {},
                                          design=design, stage="encoding")

    def test_new_probe_admission_checks_reservation_and_pins_completion(self):
        gate = admission(); design = runner.read(V2); code = {"runner": runner.file_hash(runner.__file__)}
        probe = SimpleNamespace(root=Path(gate["probe_job"]), identity_hash=gate["probe_identity_sha256"],
                                completion_sha256=gate["probe_completion_sha256"],
                                identity={"mode": "plmblast_alignment_v1", "panel": "probe12", "source": {},
                                          "design_sha256": runner.file_hash(V2), "code_hashes": code,
                                          "design_version": runner.RESOURCE_VERSION, "method": {"model": "fixed"},
                                          "runtime": {"runtime": "fixed"}})
        summary = {"resource_recommendation": "eligible_for_explicit_repeat32_admission", "failed_pair_count": 0,
                   "resource_gate_failures": [], "expected_pair_count": 12,
                   "forecast": {"alignment_seconds": 8859.0, "encoding_seconds": 620.0, "within_frozen_budget": True}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "explicit-root-gate.json"; path.write_text(json.dumps(gate))
            with mock.patch.object(runner, "open_completed", return_value=probe), \
                 mock.patch.object(runner, "completed_unit", return_value=summary), \
                 mock.patch.object(runner, "code_hashes", return_value=code), \
                 mock.patch.object(runner, "validate_resource_reservation") as reservation:
                result = runner.validate_admission(path, SESSION, runner.file_hash(V2), {}, design=design, stage="encoding")
                self.assertEqual(result["gate_sha256"], runner.file_hash(path))
                self.assertEqual(result["probe"]["completion_sha256"], "new-v2-completion")
                reservation.assert_called_once_with(gate, "2026-10-31T19:25:06+00:00", "encoding")


class SelectedOverlayTests(unittest.TestCase):
    @contextmanager
    def approved_and_sibling(self):
        """Two same-version byte-identical toy trees; only one is admitted."""
        design = runner.read(V2)
        evidence = design["resource_amendment"]
        old = evidence["original_negative_probe"]
        actual_hash = runner.file_hash
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            session = root / "results/session"
            approved = session / "vendor/plmblast-overlay-v1"
            sibling = session / "vendor/plmblast-overlay-sibling"
            for overlay in (approved, sibling):
                (overlay / "numpy").mkdir(parents=True)
                (overlay / "numpy/__init__.py").write_text("__version__ = '1.26.4'\n")
            installation = session / "inputs/plmblast_overlay_installation_v1.json"
            installation.parent.mkdir()
            installation.write_text(json.dumps({"status": "verified_installed", "overlay": str(approved),
                "installed_files": {"numpy/__init__.py": actual_hash(approved / "numpy/__init__.py")}}))
            resolution_path = session / evidence["overlay_resolution"]["path"]
            resolution_path.parent.mkdir(parents=True)
            resolution_path.write_text(json.dumps({"status": "passed", "all_extra_bytes_preserved": True,
                "all_original_installed_files_unchanged": True, "exact_file_membership_restored": True,
                "original_installation_sha256": actual_hash(installation)}))
            # Provenance records are already checksum-verified in this fixture;
            # exercise real installation membership/hash and path binding below.
            pinned = {session / evidence[key]["path"]: evidence[key]["sha256"]
                      for key in ("independent_resource_audit", "overlay_resolution")}
            def hash_with_pinned_provenance(path):
                path = Path(path)
                return pinned[path] if path in pinned else actual_hash(path)
            probe = SimpleNamespace(identity_hash=old["identity_sha256"], completion_sha256=old["completion_sha256"],
                                    identity={"design_sha256": runner.ORIGINAL_DESIGN_SHA256})
            summary = {"resource_recommendation": "not_admitted", "forecast": {"within_frozen_budget": False,
                       "alignment_seconds": old["alignment_forecast_seconds"]}}
            with mock.patch.object(runner, "ROOT", root), \
                 mock.patch.object(runner, "file_hash", side_effect=hash_with_pinned_provenance), \
                 mock.patch.object(runner, "open_completed", return_value=probe), \
                 mock.patch.object(runner, "completed_unit", return_value=summary), \
                 mock.patch.object(runner.importlib, "import_module", side_effect=AssertionError("No numerical imports in path validation")):
                yield SimpleNamespace(session=session, approved=approved, sibling=sibling, design=design,
                                      code_sha=actual_hash(runner.__file__), design_sha=actual_hash(V2))

    def test_only_recorded_overlay_is_accepted_even_with_identical_sibling_bytes(self):
        with self.approved_and_sibling() as value:
            self.assertEqual((value.approved / "numpy/__init__.py").read_bytes(),
                             (value.sibling / "numpy/__init__.py").read_bytes())
            runner.verify_resource_amendment(value.session, value.design, actual_overlay=value.approved)
            with self.assertRaisesRegex(ValueError, "Selected overlay differs"):
                runner.verify_resource_amendment(value.session, value.design, actual_overlay=value.sibling)

    def test_parent_and_worker_reject_sibling_before_numerical_runtime(self):
        with self.approved_and_sibling() as value, \
             mock.patch.object(runner, "load_validated_design", return_value=value.design), \
             mock.patch.object(runner, "load_alignment_runtime") as runtime:
            args = SimpleNamespace(session=value.session, overlay=value.sibling, design=V2,
                                   expected_code_sha256=value.code_sha, expected_design_sha256=value.design_sha)
            with self.assertRaisesRegex(ValueError, "Selected overlay differs"):
                runner.setup(args)
            spec = value.session / "worker-spec.json"
            code = {"runner": value.code_sha}
            spec.write_text(json.dumps({"session": str(value.session), "design_path": str(V2),
                                       "design_sha256": value.design_sha, "code_hashes": code}))
            with mock.patch.object(runner, "code_hashes", return_value=code):
                with self.assertRaisesRegex(ValueError, "Selected overlay differs"):
                    runner.worker(spec, value.sibling)
            runtime.assert_not_called()


if __name__ == "__main__":
    unittest.main()
