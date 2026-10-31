#!/usr/bin/env python3
"""Bounded, explicitly adapted pLM-BLAST comparison using existing ESM2 states.

Encoding and alignment are separate processes/environments. Never installs,
downloads, tunes thresholds, changes upstream algorithms, or admits a full run.
"""
from __future__ import annotations
import argparse
import copy
import datetime as dt
import importlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import select
import signal
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.session_store import UnitStore, atomic_json, canonical, digest, file_hash
from iri_constructs.session_jobs import cleanup_owned_group, OwnedCleanupFailure

KS = (1, 5, 10)
PARAMETERS = dict(enh=False, norm=False, bfactor=2, sigma_factor=2.0,
                  gap_penalty=0.5, min_spanlen=25, window_size=15, filter_results=True)
REVISION = "207c252c96361096d8d5a54faa5705af6e9848be"
ORIGINAL_DESIGN_SHA256 = "98528bfad9d28a6b2c9de7a1c8d50c050cae5c9759bbb3dbd104e12fcfcd83b1"
RESOURCE_VERSION = "iri-plmblast-existing-esm2-resource-v2"
RESOURCE_STATUS = "frozen_after_v1_timing_before_repeat32_benchmark_scoring"
RESOURCE_AMENDMENT = {
    "basis": "Explicit resource allocation after the original v1 timing failure; no repeat32 benchmark scoring observed",
    "only_limit_change": {"path": "limits.repeat_alignment_seconds", "original": 5400, "amended": 10800},
    "original_design": {"path": "inputs/plmblast_design_v1.json", "sha256": ORIGINAL_DESIGN_SHA256},
    "original_negative_probe": {
        "path": "jobs/plmblast-align-probe-v1-77145f0dd4e5",
        "identity_sha256": "77145f0dd4e5672f76f07cba6a25947414cca03f5cfebb794ecf7d2461edc641",
        "completion_sha256": "8ee4cda85ab49638c542ce7c2166f7e4a7047f47027d1cb1a072146ab1e49d06",
        "decision": "not_admitted", "alignment_forecast_seconds": 8859.360704898589},
    "prior_runner": {"git_commit": "9ef47d3", "sha256": "640ad1d05be5abfe1ed0ee8a7a68e690bddf57f1cf75e83b1d6e4995ebb10035"},
    "independent_resource_audit": {
        "path": "audit/independent_plmblast_probe_resource_v2.json",
        "sha256": "84ff7509054023554d347298f10bc53ca4a2eb44b300031c88100231379f8887"},
    "overlay_resolution": {
        "path": "vendor/plmblast-runtime-bytecode-quarantine-v1/report.json",
        "sha256": "efd6c97204a62b5f4856412b396add888350aa3a3bfbcb6cd537214e11fb237a"},
    "command_local_environment": {"PYTHONDONTWRITEBYTECODE": "1"},
    "reservation_seconds": {"alignment": 10800, "encoding": 1800, "cleanup_and_audit": 3600},
    "admission": "Root explicit decision after priority-work/resource review and a new v2 encoding/probe identity; original v1 remains not_admitted",
}
SPAWNING = False
CLEANING = False
PENDING_SIGNAL = None
OWNED_CHILDREN = {}


def interrupted(signum, frame):
    """Record cancellation during spawn; finish ownership before unwinding."""
    global PENDING_SIGNAL
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if SPAWNING or CLEANING:
        PENDING_SIGNAL = signum
        return
    raise SystemExit(128 + signum)


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def code_hashes():
    names = ["scripts/session_plmblast_experiment.py", "scripts/session_context_profile.py",
             "scripts/session_experiment.py", "scripts/session_catalogue_analysis.py",
             "src/iri_constructs/session_store.py", "src/iri_constructs/session_jobs.py", "src/iri_constructs/session_models.py",
             "src/iri_constructs/session_benchmark.py", "src/iri_constructs/t05_data.py",
             "src/iri_constructs/t05_composition.py"]
    return {name: file_hash(ROOT / name) for name in names}


def open_completed(path, session, *, model_probe=False):
    session = Path(session).resolve()
    path = Path(path).resolve()
    require(path.is_relative_to(session / "jobs"), "Source must be an owned session job")
    store = UnitStore(path, read(path / "identity.json"), session / "job_ledger.jsonl")
    candidates = [(read(p)["completed_utc"], p) for p in path.glob("completion-*.json")]
    require(candidates, "Source has no completion descriptor")
    completion_path = max(candidates)[1]
    completion = read(completion_path)
    expected = store.identity["query_ids" if model_probe else "unit_keys"]
    require(expected and len(expected) == len(set(expected)) and
            completion["identity_hash"] == store.identity_hash and
            completion["expected_count"] == len(expected) and set(completion["units"]) == set(expected),
            "Completed source coverage/identity mismatch")
    store.completed_paths = {}
    for key, entry in completion["units"].items():
        unit = (path / entry["path"]).resolve()
        require(unit.is_relative_to(path) and file_hash(unit / "complete.json") == entry["descriptor_sha256"],
                "Source unit descriptor changed")
        require(store.verify(unit)["key"] == key, "Source unit key changed")
        store.completed_paths[key] = unit
    store.completion_sha256 = file_hash(completion_path)
    return store


def completed_unit(store, key, arrays=False):
    path = store.completed_paths[key]
    record = read(path / "record.json")
    if not arrays:
        return record
    import numpy as np
    require("arrays.npz" in store.verify(path)["files"], "Arrays lack a source checksum")
    with np.load(path / "arrays.npz", allow_pickle=False) as archive:
        return record, {name: archive[name].copy() for name in archive.files}


def source_evidence(session, design):
    manifest_path = session / "audit/plmblast_source_acquisition.json"
    require(file_hash(manifest_path) == design["source_manifest_sha256"], "Source acquisition manifest changed")
    manifest = read(manifest_path)
    require(manifest["revision"] == REVISION and manifest["licence"] == "MIT", "Wrong upstream source")
    directory = session / "vendor" / ("pLM-BLAST-" + REVISION)
    records = {}
    for item in manifest["files"]:
        path = (directory / item["path"]).resolve()
        require(path.is_relative_to(directory) and path.stat().st_size == item["bytes"] and
                file_hash(path) == item["sha256"], "Pinned pLM-BLAST file changed")
        records[item["path"]] = item["sha256"]
    return {"revision": REVISION, "manifest_sha256": file_hash(manifest_path), "files": records}


def validate_resource_design(design, original):
    """Allow only the exact resource amendment relative to checksum-pinned v1."""
    require(original["version"] == "iri-plmblast-existing-esm2-v1" and
            original["status"] == "frozen_before_installation_import_timing_or_scoring" and
            original["limits"]["repeat_alignment_seconds"] == 5400 and
            original["limits"]["repeat_encoding_seconds"] == 1800, "Wrong original design")
    if design.get("version") == original["version"]:
        require(canonical(design) == canonical(original), "Original v1 design changed")
        return
    require(design.get("version") == RESOURCE_VERSION, "Unknown pLM-BLAST design version")
    expected = copy.deepcopy(original)
    expected["version"] = RESOURCE_VERSION
    expected["status"] = RESOURCE_STATUS
    expected["limits"]["repeat_alignment_seconds"] = 10800
    require(expected["forecast"]["admission"].count("alignment forecast<=5400s") == 1,
            "Original admission description differs")
    expected["forecast"]["admission"] = expected["forecast"]["admission"].replace(
        "alignment forecast<=5400s", "alignment forecast<=10800s")
    expected["resource_amendment"] = RESOURCE_AMENDMENT
    require(canonical(design) == canonical(expected), "Resource v2 exceeds the exact cap-only amendment")


def load_validated_design(session, design_path):
    original_path = session / "inputs/plmblast_design_v1.json"
    require(file_hash(original_path) == ORIGINAL_DESIGN_SHA256, "Pinned original v1 design bytes changed")
    design = read(design_path)
    validate_resource_design(design, read(original_path))
    if design["version"] == RESOURCE_VERSION:
        require(os.environ.get("PYTHONDONTWRITEBYTECODE") == "1" and sys.dont_write_bytecode,
                "Resource v2 requires command-local PYTHONDONTWRITEBYTECODE=1 before interpreter startup")
    return design


def verify_resource_amendment(session, design, *, actual_overlay=None):
    if design["version"] != RESOURCE_VERSION:
        return
    evidence = design["resource_amendment"]
    original = evidence["original_negative_probe"]
    probe = open_completed(session / original["path"], session)
    require(probe.identity_hash == original["identity_sha256"] and
            probe.completion_sha256 == original["completion_sha256"] and
            probe.identity["design_sha256"] == ORIGINAL_DESIGN_SHA256,
            "Original negative probe identity/completion changed")
    summary = completed_unit(probe, "summary")
    require(summary["resource_recommendation"] == original["decision"] == "not_admitted" and
            summary["forecast"]["within_frozen_budget"] is False and
            summary["forecast"]["alignment_seconds"] == original["alignment_forecast_seconds"],
            "Original v1 negative resource decision changed")
    for key in ("independent_resource_audit", "overlay_resolution"):
        require(file_hash(session / evidence[key]["path"]) == evidence[key]["sha256"],
                "Resource amendment evidence changed: " + key)
    resolution = read(session / evidence["overlay_resolution"]["path"])
    require(resolution["status"] == "passed" and resolution["all_extra_bytes_preserved"] is True and
            resolution["all_original_installed_files_unchanged"] is True and resolution["exact_file_membership_restored"] is True,
            "Overlay resolution did not pass")
    installation_path = session / "inputs/plmblast_overlay_installation_v1.json"
    require(file_hash(installation_path) == resolution["original_installation_sha256"], "Resolved installation record changed")
    installation = read(installation_path)
    overlay = session / "vendor/plmblast-overlay-v1"
    require(installation["status"] == "verified_installed" and Path(installation["overlay"]).resolve() == overlay.resolve(),
            "Resolved overlay differs")
    if actual_overlay is not None:
        require(Path(actual_overlay).resolve(strict=True) == overlay.resolve(strict=True),
                "Selected overlay differs from the pinned installation overlay")
    entries = list(overlay.rglob("*"))
    require(not any(p.is_symlink() for p in entries) and
            {str(p.relative_to(overlay)) for p in entries if p.is_file()} == set(installation["installed_files"]),
            "Resolved overlay exact file membership changed")
    for relative, sha in installation["installed_files"].items():
        require(file_hash(overlay / relative) == sha, "Resolved overlay bytes changed")


def setup(args):
    session = args.session.resolve()
    require(session.is_relative_to(ROOT / "results"), "Session outside owned results")
    require(file_hash(Path(__file__)) == args.expected_code_sha256, "Runner differs from frozen invocation")
    design_path = args.design or session / "inputs/plmblast_design_v1.json"
    require(file_hash(design_path) == args.expected_design_sha256, "Design differs from frozen invocation")
    design = load_validated_design(session, design_path)
    verify_resource_amendment(session, design, actual_overlay=args.overlay)
    require(design["extractor"] == PARAMETERS and design["budgets"] == list(KS) and
            design["postprocess"] == dict(alignment_cutoff=0.3, round_decimals=3, reduce_duplicates=False) and
            design["screening"]["COS_PER_CUT"] == 0 and design["screening"]["workers"] == 1,
            "Frozen scientific settings changed")
    from iri_constructs.session_benchmark import load_suite
    suite = load_suite(ROOT, session / "inputs")
    require(suite["verification"]["challenge_manifest_sha256"] == design["challenge_sha256"] and
            suite["manifest"]["repeat_query_ids"] == design["repeat_query_ids"] and
            len(design["repeat_query_ids"]) == 32 and len(suite["sequences"]) == 105,
            "Frozen suite or repeat panel changed")
    queries = {q["query_id"]: q for q in suite["queries"]}
    require([queries[q]["length_aa"] for q in design["probe_query_ids"]] == design["probe_query_lengths"] and
            [len(suite["sequences"][a]) for a in design["probe_reference_accessions"]] == design["probe_reference_lengths"],
            "Fixed probe lengths differ")
    require(file_hash(ROOT / "src/iri_constructs/session_models.py") == design["model_adapter_sha256"],
            "Frozen model adapter changed")
    return session, design, design_path, suite, source_evidence(session, design)


def verify_model_probe(path, session, design, suite):
    probe = open_completed(path, session, model_probe=True)
    identity = probe.identity
    require(identity["mode"] == "admission_probe_v1" and identity["query_ids"] == design["probe_query_ids"] and
            identity["repeats"] == 3 and identity["probe_code_sha256"] == design["model_probe_code_sha256"],
            "Wrong completed model-probe source")
    model = identity["model"]
    for name in ("model_id", "revision", "layer", "num_layers", "dtype", "device"):
        require(model[name] == design["model"][name], "Model probe differs: " + name)
    require(model["checkpoint_files"]["model.safetensors"]["sha256"] == design["model"]["checkpoint_sha256"] and
            model["adapter_source_sha256"] == design["model_adapter_sha256"] and
            model["representation"] == "last_layer_post_final_norm" and model["batch_size"] == 1,
            "Probe checkpoint/representation differs")
    config_path = session / "inputs/configs/esm2_t30_150M_UR50D-cpu.json"
    require(file_hash(config_path) == design["model_config_sha256"] and read(config_path) == identity["config"],
            "Probe configuration differs from frozen config")
    lookup = {q["query_id"]: q for q in suite["queries"]}
    for qid in design["probe_query_ids"]:
        query = lookup[qid]
        record = completed_unit(probe, qid)
        require(identity["query_hashes"][qid] == digest(query["sequence"].encode()) and
                record["query_id"] == qid and record["length"] == query["length_aa"] and
                record["within_config_max_abs"] <= 1e-6 and record["actual_parameter_dtypes"] == ["torch.float32"] and
                len(record["actual_forward_outputs"]) == 3 and
                all(row == {"dtype": "torch.float32", "device": "cpu"} for row in record["actual_forward_outputs"]),
                "Probe hashes or raw precision evidence differ")
    return probe


def check_states(array, sequence, dimension=640):
    import numpy as np
    require(array.dtype == np.float32 and array.shape == (len(sequence), dimension) and
            np.isfinite(array).all() and (np.linalg.norm(array.astype(np.float64), axis=1) > 0).all(),
            "Residue states have invalid shape/dtype/finite values/norm")


def validate_state_record(record, sequence):
    require(record['status'] == 'success' and record['sequence_sha256'] == digest(sequence.encode()) and
            record['length'] == len(sequence), 'Saved residue-state source/status/length differs')


def close_states(left, right):
    import numpy as np
    require(left.shape == right.shape and np.isfinite(left).all() and np.isfinite(right).all(), "State shapes differ")
    error = float(np.max(np.abs(left.astype(np.float64) - right.astype(np.float64))))
    require(error <= 1e-6, "Fresh states differ from completed model probe")
    return error


def pair_list(query_ids, accessions):
    return [(qid, accession) for qid in query_ids for accession in accessions]


def pair_key(qid, accession):
    return "pair-" + qid + "--" + accession


def encoding_panel(design, suite, probe):
    return (design["probe_query_ids"], design["probe_reference_accessions"]) if probe else (
        design["repeat_query_ids"], sorted(suite["sequences"]))


def evidence_of(store):
    return {"path": str(store.root), "identity_sha256": store.identity_hash,
            "completion_sha256": store.completion_sha256}


def new_store(session, name, identity):
    require(Path(name).name == name and name not in {"", ".", ".."}, "Unsafe job name")
    return UnitStore(session / "jobs" / (name + "-" + digest(canonical(identity))[:12]),
                     identity, session / "job_ledger.jsonl")


def session_seconds(session):
    deadline = dt.datetime.fromisoformat(read(session / "launch.json")["hard_stop_utc"])
    require(deadline.tzinfo is not None, "Session deadline lacks timezone")
    return (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds()


def forecast_within_limits(forecast, design):
    alignment, encoding = forecast["alignment_seconds"], forecast["encoding_seconds"]
    require(all(math.isfinite(value) and value > 0 for value in (alignment, encoding)), "Invalid resource forecast")
    return (alignment <= design["limits"]["repeat_alignment_seconds"] and
            encoding <= design["limits"]["repeat_encoding_seconds"])


def validate_resource_reservation(gate, hard_stop, stage, now=None):
    """Keep the full stage budgets plus one hour for cleanup/audit available."""
    require(stage in ("encoding", "alignment"), "Unknown resource-reservation stage")
    require(gate["schema"] == "iri-plmblast-explicit-resource-admission-v2" and gate["decided_by"] == "root" and
            gate["priority_work_and_actual_resources_reviewed"] is True, "Resource v2 needs explicit root review")
    reservation = gate["resource_reservation"]
    expected = {"session_hard_stop_utc": hard_stop, "alignment_seconds": 10800,
                "encoding_seconds": 1800, "cleanup_audit_seconds": 3600}
    require(canonical(reservation) == canonical(expected), "Resource v2 reservation changed")
    deadline = dt.datetime.fromisoformat(hard_stop)
    decided = dt.datetime.fromisoformat(gate["decision_utc"])
    now = dt.datetime.now(dt.timezone.utc) if now is None else now
    require(deadline.tzinfo is not None and decided.tzinfo is not None and now.tzinfo is not None and decided <= now,
            "Resource reservation timestamps are invalid")
    require((deadline - decided).total_seconds() >= 16200, "Original resource decision lacked the full 16200 s reservation")
    required = 16200 if stage == "encoding" else 14400
    require((deadline - now).total_seconds() >= required, "Insufficient remaining session time for resource v2 " + stage)


def validate_admission(path, session, design_hash, sources, *, design, stage):
    require(path is not None, "Repeat32 requires explicit passed pair-probe resource admission")
    gate = read(path)
    require(gate["status"] == "passed" and gate["decision"] == "admit_repeat32" and gate["resource_gate"] == "passed",
            "Repeat32 admission not granted")
    probe = open_completed(Path(gate["probe_job"]), session)
    require(probe.identity["mode"] == "plmblast_alignment_v1" and probe.identity["panel"] == "probe12" and
            probe.completion_sha256 == gate["probe_completion_sha256"] and
            probe.identity["design_sha256"] == design_hash and probe.identity["source"] == sources and
            probe.identity["code_hashes"] == code_hashes(), "Admission refers to another probe/code/design")
    summary = completed_unit(probe, "summary")
    require(summary["resource_recommendation"] == "eligible_for_explicit_repeat32_admission" and
            summary["failed_pair_count"] == 0 and not summary["resource_gate_failures"] and summary["expected_pair_count"] == 12 and
            summary["forecast"]["within_frozen_budget"] is True and forecast_within_limits(summary["forecast"], design),
            "Probe/resource gates failed")
    if design["version"] == RESOURCE_VERSION:
        require(gate["design_sha256"] == design_hash and gate["runner_sha256"] == file_hash(__file__) and
                gate["probe_identity_sha256"] == probe.identity_hash and probe.identity["design_version"] == RESOURCE_VERSION and
                gate["original_v1_design_sha256"] == ORIGINAL_DESIGN_SHA256 and
                canonical(gate["original_v1_negative_probe"]) == canonical(RESOURCE_AMENDMENT["original_negative_probe"]),
                "Resource v2 admission provenance differs")
        validate_resource_reservation(gate, read(session / "launch.json")["hard_stop_utc"], stage)
    return {"gate_sha256": file_hash(path), "probe": evidence_of(probe),
            "method": probe.identity["method"], "runtime": probe.identity["runtime"],
            "forecast": summary["forecast"]}


def encode(args, session, design, design_path, suite, sources):
    import numpy as np
    from session_experiment import mean_residues, model_identity, peak_mib
    from session_context_profile import PrecisionAudit
    probe_mode = args.mode == "encode-probe"
    original = verify_model_probe(args.model_probe_job, session, design, suite)
    admission = None if probe_mode else validate_admission(args.admission, session, file_hash(design_path), sources,
                                                          design=design, stage="encoding")
    require(admission is None or admission["method"] == original.identity["model"], "Admission model differs")
    query_ids, accessions = encoding_panel(design, suite, probe_mode)
    lookup = {q["query_id"]: q for q in suite["queries"]}
    units = ["reference-" + a for a in accessions] + ["query-" + q for q in query_ids]
    identity = {"mode": "plmblast_encoding_v1", "panel": "probe12" if probe_mode else "repeat32",
                "design_version": design["version"], "design_status": design["status"],
                "design_sha256": file_hash(design_path), "source": sources, "code_hashes": code_hashes(),
                "method": original.identity["model"], "model_probe": evidence_of(original),
                "query_ids": query_ids, "reference_accessions": accessions, "unit_keys": units + ["summary"],
                "sequence_hashes": {"reference-" + a: digest(suite["sequences"][a].encode()) for a in accessions} |
                                   {"query-" + q: digest(lookup[q]["sequence"].encode()) for q in query_ids},
                "admission": admission, "adaptation": "Existing150M last-layer FP32 residue states; no ProtT5 checkpoint"}
    store = new_store(session, args.name, identity)
    encoder = audit = None
    deadline = time.monotonic() + min(session_seconds(session), design["limits"]["repeat_encoding_seconds"])
    try:
        if not probe_mode:
            from iri_constructs.session_models import load_encoder
            encoder = load_encoder(original.identity["config"])
            require(model_identity(encoder) == original.identity["model"], "Encoding model/runtime differs from completed probe")
            audit = PrecisionAudit(encoder)
        for key in units:
            require(time.monotonic() < deadline, "Encoding/session deadline reached")
            is_ref = key.startswith("reference-")
            identifier = key.split("-", 1)[1]
            sequence = suite["sequences"][identifier] if is_ref else lookup[identifier]["sequence"]
            saved = store.read(key, arrays=True)
            if saved is not None:
                record, arrays = saved
                validate_state_record(record, sequence)
                check_states(arrays["residues"], sequence)
                continue
            start = time.monotonic()
            source_id = "clean_" + identifier if is_ref else identifier
            comparison = None
            if probe_mode:
                parent, arrays = completed_unit(original, source_id, arrays=True)
                states = arrays["residues"]
                close_states(mean_residues(states), arrays["mean"])
                precision = {"source_probe_unit": source_id, "actual_parameter_dtypes": parent["actual_parameter_dtypes"],
                             "actual_forward_outputs": parent["actual_forward_outputs"], "reused_without_forward": True}
            else:
                states, raw = audit.encode(sequence)
                precision = {"raw_forward": raw, "parameters": audit.metadata}
                if source_id in original.identity["query_ids"]:
                    _, old = completed_unit(original, source_id, arrays=True)
                    comparison = close_states(states, old["residues"])
            check_states(states, sequence)
            store.publish(key, {"status": "success", "sequence_sha256": digest(sequence.encode()), "length": len(sequence),
                               "precision": precision, "probe_max_abs": comparison, "seconds": time.monotonic() - start,
                               "peak_process_rss_mib": peak_mib()}, {"residues": states, "mean": mean_residues(states)})
        if store.read("summary") is None:
            store.publish("summary", {"status": "success", "reference_count": len(accessions), "query_count": len(query_ids),
                                      "peak_process_rss_mib": peak_mib(), "model_probe": evidence_of(original)})
        store.complete(identity["unit_keys"])
    finally:
        if audit is not None:
            audit.close()
        if encoder is not None:
            encoder.close()
    return store


def validate_path(indices, query_length, reference_length):
    require(isinstance(indices, list) and len(indices) >= 1, "Empty native alignment path")
    previous = None
    for pair in indices:
        require(isinstance(pair, list) and len(pair) == 2 and all(type(v) is int for v in pair), "Invalid path coordinates")
        y, x = pair
        require(0 <= y < reference_length and 0 <= x < query_length, "Native path outside sequence bounds")
        if previous is not None:
            delta = (y - previous[0], x - previous[1])
            require(delta in ((1, 0), (0, 1), (1, 1)), "Native path is not a monotone adjacent alignment")
        previous = pair


def normalize_native_result(frame, formatted, query_length, reference_length):
    """Preserve native raw paths; select only entries retained by original output formatting."""
    import numpy as np
    if frame is None or len(frame) == 0:
        require(formatted is None or len(formatted) == 0, "Formatted hits without native paths")
        return {"status": "no_eligible_alignment", "native_paths": [], "eligible_hits": []}
    require(frame.index.is_unique, "Native path identifiers are not unique")
    native, by_id = [], {}
    for index, row in frame.iterrows():
        indices = np.asarray(row["indices"])
        require(np.issubdtype(indices.dtype, np.integer) and indices.ndim == 2 and indices.shape[1] == 2,
                "Native alignment path dtype/shape differs")
        path = indices.astype(np.int64).tolist()
        validate_path(path, query_length, reference_length)
        score = float(row["score"])
        require(math.isfinite(score) and -1.00001 <= score <= 1.00001, "Invalid unenhanced native score")
        value = {"path_key": str(index), "raw_score": score, "indices_reference_query": path,
                 "query_start0": path[0][1], "query_end0_exclusive": path[-1][1] + 1,
                 "reference_start0": path[0][0], "reference_end0_exclusive": path[-1][0] + 1}
        native.append(value)
        by_id[index] = value
    hits = []
    for index, row in formatted.iterrows():
        require(index in by_id, "Formatted output introduced a path")
        raw = by_id[index]
        require(raw["raw_score"] >= float(np.float32(0.3)), "Original cutoff returned an ineligible path")
        require([int(row[x]) for x in ("qstart", "qend", "tstart", "tend")] ==
                [raw["query_start0"], raw["query_end0_exclusive"] - 1,
                 raw["reference_start0"], raw["reference_end0_exclusive"] - 1], "Formatted coordinate orientation differs")
        score = float(row["score"])
        require(math.isfinite(score) and score == float(np.asarray([raw["raw_score"]], dtype=np.float32).round(3)[0]),
                "CLI three-decimal score differs from original raw FP32 score")
        hits.append({**raw, "cli_score": score})
    expected = set(frame.index[frame["score"] >= 0.3])
    require(set(formatted.index) == expected, "Original cutoff/postprocessing coverage differs")
    return {"status": "eligible_alignment" if hits else "no_eligible_alignment", "native_paths": native, "eligible_hits": hits}


def compare_pair(module, upstream, X, Y, query_sequence, reference_sequence, *, direct=False):
    import pandas as pd
    frame = module.full_compare(X, Y, 0, 0) if not direct else module.embedding_to_span(X, Y)
    if direct and len(frame):
        frame["queryid"], frame["dbid"] = 0, 0
        frame = upstream.postprocess.filter_result_dataframe(frame)
    if frame is None or len(frame) == 0:
        return normalize_native_result(frame, None, len(query_sequence), len(reference_sequence))
    query_frame = frame.copy()
    query_frame["id"], query_frame["sequence"] = "q", query_sequence
    database = pd.DataFrame([{"id": "r", "sequence": reference_sequence}])
    formatted = upstream.postprocess.prepare_output(query_frame, database, alignment_cutoff=0.3, verbose=False)
    return normalize_native_result(frame, formatted, len(query_sequence), len(reference_sequence))


def ranking_from_pairs(pair_records):
    require(len({row["accession"] for row in pair_records}) == len(pair_records), "Duplicate candidate audit")
    require(all(row["status"] in ("eligible_alignment", "no_eligible_alignment") for row in pair_records),
            "Failed candidates cannot be silently ranked as no-hits")
    ranked = []
    for row in pair_records:
        require((row["status"] == "eligible_alignment") == bool(row["eligible_hits"]), "Eligibility/status mismatch")
        if row["eligible_hits"]:
            best = max(hit["cli_score"] for hit in row["eligible_hits"])
            require(math.isfinite(best), "Nonfinite candidate score")
            ranked.append({"accession": row["accession"], "score": best,
                           "best_raw_eligible_score": max(hit["raw_score"] for hit in row["eligible_hits"])})
    ranked.sort(key=lambda row: (-row["score"], row["accession"]))
    return [{**row, "rank": index} for index, row in enumerate(ranked, 1)]


def eligible_metrics(ranking, catalogue, query, equivalence, families=None):
    """Membership defines scorable labels, independently of sparse returned hits."""
    catalogue = set(catalogue)
    require(query["components"] and len(set(query["components"])) == len(query["components"]), "Empty or duplicated components")
    require([r["rank"] for r in ranking] == list(range(1, len(ranking) + 1)) and
            len({r["accession"] for r in ranking}) == len(ranking) and
            {r["accession"] for r in ranking} <= catalogue and
            all(math.isfinite(r["score"]) for r in ranking) and
            ranking == sorted(ranking, key=lambda r: (-r["score"], r["accession"])), "Invalid eligible-only ranking")
    ranks = {r["accession"]: r["rank"] for r in ranking}
    components = []
    for accession in query["components"]:
        labels = {"exact": [accession], "duplicate_equivalent": equivalence.get(accession, [accession])}
        if families is not None:
            labels["family"] = families.get(accession, [])
        matches = {}
        for kind, accepted in labels.items():
            surviving = sorted(set(accepted) & catalogue)
            hit_ranks = [ranks[a] for a in surviving if a in ranks]
            rank = min(hit_ranks) if hit_ranks else None
            matches[kind] = {"scorable": bool(surviving), "accepted_surviving": surviving,
                             "rank": rank, "abstention": bool(surviving) and rank is None,
                             "reciprocal_rank": (1 / rank if rank is not None else 0.0) if surviving else None,
                             **{f"hit_at_{k}": int(rank is not None and rank <= k) if surviving else None for k in KS}}
        components.append({"accession": accession, "matches": matches})
    query_metrics = {}
    for kind in components[0]["matches"]:
        rows = [c["matches"][kind] for c in components]
        scorable = [row for row in rows if row["scorable"]]
        all_scorable = len(scorable) == len(rows)
        query_metrics[kind] = {"scorable_components": len(scorable), "unscorable_components": len(rows) - len(scorable),
                              "all_components_scorable": all_scorable}
        for k in KS:
            hits = sum(row[f"hit_at_{k}"] for row in scorable)
            query_metrics[kind].update({f"component_hits_at_{k}": hits,
                f"any_component_hit_at_{k}": int(hits > 0) if all_scorable else None,
                f"all_components_hit_at_{k}": int(hits == len(rows)) if all_scorable else None})
    return {"query_id": query["query_id"], "suite": query.get("suite"), "group_id": query.get("group_id"),
            "component_count": len(components), "catalogue_count": len(catalogue), "eligible_accession_count": len(ranking),
            "components": components, "query_metrics": query_metrics,
            "retrieval_budget": "K total genuine eligible unique accessions; no no-hit tail",
            "returned_counts": {str(k): min(k, len(ranking)) for k in KS}}


def summarize_eligible(records):
    outcomes = {}
    for scope in ("all", "clean", "combined", "pairs", "triples", "quadruples"):
        subset = [r for r in records if scope == "all" or r["suite"] == scope or
                  (scope == "combined" and r["suite"] != "clean")]
        valid = [r["evaluation"] for r in subset if r["status"] == "success"]
        value = {"expected_queries": len(subset), "evaluated_queries": len(valid),
                 "failed_queries": len(subset) - len(valid), "metrics_scope": "successful complete candidate audits only; failures reported separately",
                 "kinds": {}}
        kinds = set().union(*(set(v["query_metrics"]) for v in valid))
        for kind in sorted(kinds):
            all_rows = [c["matches"][kind] for v in valid for c in v["components"]]
            scorable = [r for r in all_rows if r["scorable"]]
            complete = [v["query_metrics"][kind] for v in valid if v["query_metrics"][kind]["all_components_scorable"]]
            metrics = {"scorable_components": len(scorable), "unscorable_components": len(all_rows) - len(scorable),
                       "no_hit_components": sum(r["abstention"] for r in scorable),
                       "fully_scorable_queries": len(complete),
                       "mrr": math.fsum(r["reciprocal_rank"] for r in scorable) / len(scorable) if scorable else None}
            for k in KS:
                metrics[f"component_recall_at_{k}"] = sum(r[f"hit_at_{k}"] for r in scorable) / len(scorable) if scorable else None
                for name in ("any", "all"):
                    key = f"{name}_component{'s' if name == 'all' else ''}_hit_at_{k}"
                    metrics[key] = sum(r[key] for r in complete) / len(complete) if complete else None
            value["kinds"][kind] = metrics
        outcomes[scope] = value
    return outcomes


def numeric_gate(upstream, sequence):
    """Numeric arrays only; formatter uses an already approved existing sequence."""
    import numpy as np
    import pandas as pd
    length = len(sequence)
    module = upstream.Extractor(**PARAMETERS)
    x = np.eye(length, dtype=np.float32)
    ordinary = compare_pair(module, upstream, x, x, sequence, sequence)
    direct = compare_pair(upstream.Extractor(**PARAMETERS), upstream, x, x, sequence, sequence, direct=True)
    require(ordinary == direct and ordinary["status"] == "eligible_alignment" and
            max(h["cli_score"] for h in ordinary["eligible_hits"]) == 1.0,
            "Original full_compare/direct pipeline identity fixture failed")
    nohit = compare_pair(module, upstream, x, -x, sequence, sequence)
    require(nohit["status"] == "no_eligible_alignment", "Native no-hit fixture failed")
    a = np.asarray([[1, 0], [0, 1], [1, 1], [2, 1]], dtype=np.float32)
    b = np.asarray([[2, 3], [3, 1], [1, 2]], dtype=np.float32)
    cosine = upstream.embedding_local_similarity(a, b)
    expected = (a.astype(float) / np.linalg.norm(a.astype(float), axis=1)[:, None]) @ (
        b.astype(float) / np.linalg.norm(b.astype(float), axis=1)[:, None]).T
    require(cosine.shape == (3, 4) and np.max(np.abs(cosine - expected.T)) <= 2e-6,
            "Original reference-row/query-column cosine fixture failed")
    indices = np.stack([np.arange(40), np.arange(40)], axis=1).astype(np.int32)
    raw = pd.DataFrame([dict(pathid=i, spanid=0, span_start=0, span_end=40, len=40,
                            indices=indices, queryid=0, dbid=0) for i in range(4)])
    raw["score"] = np.asarray([0.29999, 0.30001, 0.70041, 0.70049], dtype=np.float32)
    fmt = raw.copy()
    fmt["id"], fmt["sequence"] = "q", sequence
    formatted = upstream.postprocess.prepare_output(fmt, pd.DataFrame([dict(id="r", sequence=sequence)]),
                                                   alignment_cutoff=0.3, verbose=False)
    checked = normalize_native_result(raw, formatted, length, length)
    require({h["path_key"] for h in checked["eligible_hits"]} == {"1", "2", "3"} and
            len({h["cli_score"] for h in checked["eligible_hits"] if h["path_key"] in {"2", "3"}}) == 1,
            "Original cutoff-before-rounding fixture failed")
    return {"status": "passed", "identity_pipeline_equivalence": True, "native_no_hit": True,
            "orientation_cosine_max_abs": float(np.max(np.abs(cosine - expected.T))),
            "cutoff_before_rounding": True, "rounding_tie": True,
            "coordinate_validation": True, "fixture_inputs": "numeric arrays plus existing approved clean sequence for formatting"}


def load_alignment_runtime(session, design, overlay):
    require(sys.version_info[:2] == tuple(design["alignment_environment"]["python_major_minor"]), "Wrong alignment Python")
    root = session / "vendor" / ("pLM-BLAST-" + REVISION)
    require("alntools" not in sys.modules, "Upstream package already imported")
    sys.dont_write_bytecode = True
    cache = session / "tmp/plmblast-numba-cache" / digest(canonical(design["alignment_environment"]))[:16]
    cache.mkdir(parents=True, exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = str(cache)
    os.environ["NUMBA_NUM_THREADS"] = "1"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "4"
    sys.path.insert(0, str(root))
    modules = {}
    for package in ("numpy", "pandas", "scipy", "numba", "llvmlite", "pytz", "tzdata", "dateutil", "six", "torch", "Bio", "tqdm", "psutil"):
        module = importlib.import_module(package)
        modules[package] = module
        if package in design["alignment_environment"]:
            require(module.__version__ == design["alignment_environment"][package] and
                    Path(module.__file__).resolve().is_relative_to(overlay), "Wrong isolated dependency version/path: " + package)
    modules["torch"].set_num_threads(4)
    modules["torch"].set_num_interop_threads(1)
    modules["numba"].set_num_threads(1)
    upstream = importlib.import_module("alntools")
    require(Path(upstream.__file__).resolve() == root / "alntools/__init__.py", "Upstream import escaped pinned source")
    source_evidence(session, design)
    records = {}
    for name, module in modules.items():
        path = Path(module.__file__).resolve()
        distribution = importlib.metadata.distribution({"Bio": "biopython", "dateutil": "python-dateutil"}.get(name, name))
        metadata = Path(distribution._path) / "METADATA"
        records[name] = {"version": distribution.version, "module_path": str(path), "module_sha256": file_hash(path),
                         "distribution_metadata_path": str(metadata), "distribution_metadata_sha256": file_hash(metadata)}
    runtime = {"python": sys.version, "executable": sys.executable, "overlay": str(overlay), "packages": records,
               "overlay_binary_sha256": {str(path.relative_to(overlay)): file_hash(path)
                                         for path in sorted(overlay.rglob("*")) if path.is_file() and path.suffix in {".so", ".dylib"}},
               "torch_threads": modules["torch"].get_num_threads(), "numba_threads": modules["numba"].get_num_threads(),
               "numerical_policy": "Original FP32 Numba functions and original fastmath/cache settings; sequential alignment"}
    return upstream, runtime


def worker(spec_path, overlay):
    """One persistent heavy child. Parent owns deadlines and process-group termination."""
    started = time.monotonic()
    spec = read(spec_path)
    session = Path(spec["session"])
    require(file_hash(spec["design_path"]) == spec["design_sha256"] and code_hashes() == spec["code_hashes"], "Worker code/design changed")
    design = load_validated_design(session, spec["design_path"])
    verify_resource_amendment(session, design, actual_overlay=overlay)
    upstream, runtime = load_alignment_runtime(session, design, overlay)
    from iri_constructs.session_benchmark import load_suite
    from session_experiment import peak_mib
    suite = load_suite(ROOT, session / "inputs")
    queries = {q["query_id"]: q for q in suite["queries"]}
    encoding = open_completed(Path(spec["encoding_job"]), session)
    require(evidence_of(encoding) == spec["encoding_source"], "Worker encoding source changed")
    fixtures = numeric_gate(upstream, suite["sequences"][design["probe_reference_accessions"][0]])
    print(json.dumps({"event": "ready", "runtime": runtime, "fixtures": fixtures,
                      "startup_seconds": time.monotonic() - started, "peak_process_rss_mib": peak_mib()}), flush=True)
    binding = json.loads(sys.stdin.readline())
    store = UnitStore(Path(binding["job"]), read(Path(binding["job"]) / "identity.json"), session / "job_ledger.jsonl")
    require(store.identity_hash == binding["identity_sha256"] and store.identity["runtime"] == runtime and
            store.identity["encoding_source"] == evidence_of(encoding), "Worker binding differs")
    module = upstream.Extractor(**design["extractor"])
    current_qid, current_states = None, None
    for line in sys.stdin:
        request = json.loads(line)
        if request.get("command") == "stop":
            break
        qid, accession = request["query_id"], request["accession"]
        key = pair_key(qid, accession)
        require(key in store.identity["unit_keys"] and qid in encoding.identity["query_ids"] and
                accession in encoding.identity["reference_accessions"], "Unapproved worker pair")
        started = time.monotonic()
        try:
            if qid != current_qid:
                record, arrays = completed_unit(encoding, "query-" + qid, arrays=True)
                validate_state_record(record, queries[qid]['sequence'])
                current_states = arrays["residues"]
                check_states(current_states, queries[qid]["sequence"])
                current_qid = qid
            record, arrays = completed_unit(encoding, "reference-" + accession, arrays=True)
            validate_state_record(record, suite['sequences'][accession])
            reference = arrays["residues"]
            check_states(reference, suite["sequences"][accession])
            result = compare_pair(module, upstream, current_states, reference, queries[qid]["sequence"], suite["sequences"][accession])
        except Exception as error:
            result = {"status": "failed", "exception_type": type(error).__name__, "stage": "native_pair_extraction_or_validation"}
        result.update(query_id=qid, accession=accession, native_seconds=time.monotonic() - started,
                      query_sequence_sha256=digest(queries[qid]["sequence"].encode()),
                      reference_sequence_sha256=digest(suite["sequences"][accession].encode()),
                      peak_process_rss_mib=peak_mib())
        store.publish(key, result)
        print(json.dumps({"event": "pair_complete", "key": key, "status": result["status"]}), flush=True)


def stop_child(process):
    global CLEANING
    if process is None:
        return
    if getattr(process, "_iri_cleanup_complete", False):
        return process._iri_cleanup_result
    ownership = getattr(process, "_iri_ownership", None)
    if ownership is None or ownership['pid'] != process.pid or ownership['pgid'] != process.pid or not ownership['start_new_session']:
        raise OwnedCleanupFailure("Child lacks valid recorded session ownership")
    evidence = {"status": "running", "pid": process.pid, "pgid": ownership["pgid"],
                "token": ownership["token"], "started_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
    path = Path(ownership["directory"]) / ("cleanup-" + uuid.uuid4().hex + ".json")
    error = None
    pending_before = PENDING_SIGNAL
    CLEANING = True
    try:
        try:
            recorded = read(Path(ownership["directory"]) / "ownership.json")
            if recorded != ownership:
                raise ValueError("Owned child process-group record differs")
            evidence['ownership_record_verified'] = True
        except (OSError, ValueError) as exc:
            # The Popen object and in-memory new-session record still establish
            # ownership when durable publication fails. Clean that exact group,
            # preserve the evidence error, and halt admission afterwards.
            error = exc
            evidence.update(ownership_record_verified=False, ownership_error=type(exc).__name__)
        try:
            os.killpg(ownership["pgid"], signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            pass
        # Group lifetime is independent of leader lifetime. This also checks a
        # resistant descendant when the direct child has already been reaped.
        evidence["group"] = cleanup_owned_group(ownership["pgid"], grace_seconds=2)
        process.wait(timeout=2)
        evidence.update(status="passed" if error is None else "failed", returncode=process.returncode)
    except BaseException as exc:
        error = exc
        evidence.update(status="failed", exception_type=type(exc).__name__, exception=str(exc))
    finally:
        stream_errors = []
        for stream in (process.stdin, process.stdout):
            if stream is not None:
                try:
                    stream.close()
                except BaseException as exc:
                    stream_errors.append(type(exc).__name__)
        evidence["stream_close_errors"] = stream_errors
        if stream_errors:
            evidence["status"] = "failed"
        evidence["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        try:
            atomic_json(path, evidence)
        finally:
            CLEANING = False
    if error is not None or stream_errors:
        raise OwnedCleanupFailure("Owned pLM-BLAST cleanup failed; evidence: " + str(path)) from error
    process._iri_cleanup_complete = True
    process._iri_cleanup_result = evidence
    OWNED_CHILDREN.pop(process.pid, None)
    if pending_before is None and PENDING_SIGNAL is not None:
        raise SystemExit(128 + PENDING_SIGNAL)
    return evidence


def stop_owned_children():
    """Recover ownership even if a signal interrupts a caller's assignment."""
    failures = []
    for child in list(OWNED_CHILDREN.values()):
        try:
            stop_child(child)
        except BaseException as exc:
            failures.append(exc)
    if failures:
        raise OwnedCleanupFailure("Recorded pLM-BLAST child cleanup failed") from failures[0]


def spawn_owned_child(command, directory, environment):
    """Create one new session and durably record it before raising cancellation."""
    global SPAWNING
    process = None
    directory = Path(directory).resolve()
    try:
        SPAWNING = True
        try:
            with (directory / "stderr.txt").open("xb") as stderr:
                process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                           stderr=stderr, text=True, bufsize=1, env=environment,
                                           cwd=directory, start_new_session=True)
                OWNED_CHILDREN[process.pid] = process
                ownership = {"pid": process.pid, "pgid": process.pid, "start_new_session": True,
                             "directory": str(directory), "token": uuid.uuid4().hex,
                             "command_sha256": digest(canonical(command)),
                             "spawned_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
                process._iri_ownership = ownership
                atomic_json(directory / "ownership.json", ownership)
        finally:
            SPAWNING = False
        if PENDING_SIGNAL is not None:
            raise SystemExit(128 + PENDING_SIGNAL)
        return process
    except BaseException:
        if process is not None:
            # Includes the signal-before-assignment-return case: the local
            # Popen object is owned here even though the caller has none yet.
            stop_child(process)
        raise


def child_message(process, timeout):
    require(timeout > 0, "No remaining child time budget")
    deadline = time.monotonic() + timeout
    data = bytearray()
    while b"\n" not in data:
        remaining = deadline - time.monotonic()
        if remaining <= 0 or not select.select([process.stdout], [], [], max(0, remaining))[0]:
            raise TimeoutError("Owned alignment child exceeded its frozen deadline")
        chunk = os.read(process.stdout.fileno(), 65536)
        require(chunk, "Child protocol failed or process exited")
        data.extend(chunk)
        require(len(data) <= 100000, "Child protocol response exceeded bound")
    require(data.endswith(b"\n") and data.count(b"\n") == 1, "Unexpected extra child protocol output")
    return json.loads(data.decode())


def start_child(args, session, design, design_path, encoding, total_deadline):
    attempt = session / "tmp/plmblast-workers" / uuid.uuid4().hex
    attempt.mkdir(parents=True, exist_ok=False)
    spec = {"session": str(session), "design_path": str(design_path.resolve()), "design_sha256": file_hash(design_path),
            "code_hashes": code_hashes(), "encoding_job": str(encoding.root), "encoding_source": evidence_of(encoding)}
    atomic_json(attempt / "spec.json", spec)
    process = None
    try:
        environment = dict(os.environ, TMPDIR=str(attempt), LC_ALL="C")
        process = spawn_owned_child([sys.executable, str(Path(__file__).resolve()), "--mode", "worker",
                                     "--worker-spec", str(attempt / "spec.json"), "--overlay", str(args.overlay.resolve())],
                                    attempt, environment)
        ready = child_message(process, min(design["limits"]["worker_startup_seconds"], total_deadline - time.monotonic()))
        require(ready["event"] == "ready" and ready["fixtures"]["status"] == "passed" and
                ready["peak_process_rss_mib"] <= design["limits"]["maximum_worker_peak_mib"], "Startup fixture/resource gate failed")
        atomic_json(attempt / "ready.json", ready)
        return process, ready, attempt
    except BaseException as error:
        try:
            try:
                stop_child(process)
            finally:
                stop_owned_children()
        finally:
            atomic_json(attempt / "failure.json", {"status": "failed", "stage": "worker_startup",
                                                  "exception_type": type(error).__name__})
        raise


def forecast_repeat(pair_rows, timings, startups, model_probe, suite, design):
    require(len(pair_rows) == 12 and all(r["status"] in ("eligible_alignment", "no_eligible_alignment") for r in pair_rows),
            "Forecast requires all12 successful original pairs")
    lengths = {q["query_id"]: q["length_aa"] for q in suite["queries"]}
    rates = [timings[pair_key(r["query_id"], r["accession"])] /
             (lengths[r["query_id"]] * len(suite["sequences"][r["accession"]]) *
              (lengths[r["query_id"]] + len(suite["sequences"][r["accession"]]))) for r in pair_rows]
    require(all(math.isfinite(t) and t > 0 for t in rates), "Invalid observed pair timing")
    work = sum(lengths[q] * len(seq) * (lengths[q] + len(seq))
               for q in design["repeat_query_ids"] for seq in suite["sequences"].values())
    encoding_rate = max(max(completed_unit(model_probe, q)["times_seconds"]) / lengths[q] ** 2 for q in design["probe_query_ids"])
    encoding_work = sum(len(seq) ** 2 for seq in suite["sequences"].values()) + sum(lengths[q] ** 2 for q in design["repeat_query_ids"])
    factor = design["forecast"]["safety_factor"]
    result = {"alignment_seconds": factor * (max(s["startup_seconds"] for s in startups) + max(rates) * work),
              "encoding_seconds": factor * encoding_rate * encoding_work, "safety_factor": factor,
              "maximum_observed_pair_seconds_per_work": max(rates), "repeat_work": work,
              "basis": "All12 prespecified pairs; no retrieval metric, pair substitution or threshold tuning", "is_guarantee": False}
    result["within_frozen_budget"] = forecast_within_limits(result, design)
    return result


def summary_startups(store, session, existing_summary, memory_limit):
    """Revalidate immutable startup provenance, retaining the published snapshot.

    A resumed publication collects all previous worker attempts. An already
    published summary describes its original fixed attempts; a later verification
    startup is retained separately and cannot silently rewrite those observations.
    """
    evidence = (existing_summary['startup_evidence'] if existing_summary is not None else
                [read(p) for p in sorted(store.root.glob('startup-*.json'))])
    require(evidence and len({r['ready_path'] for r in evidence}) == len(evidence), 'Missing/duplicate startup evidence')
    for record in evidence:
        path = Path(record['ready_path']).resolve()
        require(path.is_relative_to((session / 'tmp/plmblast-workers').resolve()) and
                file_hash(path) == record['ready_sha256'] and read(path) == record['ready'], 'Saved startup evidence changed')
        require(read(store.root / ('startup-' + path.parent.name + '.json')) == record, 'Startup record is not owned by this job')
        ready = record['ready']
        require(ready['event'] == 'ready' and ready['runtime'] == store.identity['runtime'] and
                ready['fixtures']['status'] == 'passed' and math.isfinite(ready['startup_seconds']) and
                ready['startup_seconds'] > 0 and ready['peak_process_rss_mib'] <= memory_limit,
                'Saved startup fixture/runtime/resource evidence differs')
    return evidence


def alignment(args, session, design, design_path, suite, sources):
    from session_experiment import public_query
    from session_catalogue_analysis import validate_family_freeze
    probe_mode = args.mode == "align-probe"
    require(args.overlay is not None and args.overlay.is_dir(), "Alignment requires an explicitly isolated dependency overlay")
    encoding = open_completed(args.encoding_job, session)
    query_ids, accessions = encoding_panel(design, suite, probe_mode)
    require(encoding.identity["mode"] == "plmblast_encoding_v1" and
            encoding.identity["panel"] == ("probe12" if probe_mode else "repeat32") and
            encoding.identity["query_ids"] == query_ids and encoding.identity["reference_accessions"] == accessions and
            encoding.identity["design_sha256"] == file_hash(design_path) and encoding.identity["source"] == sources and
            encoding.identity["code_hashes"] == code_hashes(), "Encoding source identity/panel differs")
    query_lookup = {q['query_id']: q for q in suite['queries']}
    expected_sequence_hashes = ({'reference-' + a: digest(suite['sequences'][a].encode()) for a in accessions} |
                                {'query-' + q: digest(query_lookup[q]['sequence'].encode()) for q in query_ids})
    require(encoding.identity['sequence_hashes'] == expected_sequence_hashes and
            set(encoding.identity['unit_keys']) == set(expected_sequence_hashes) | {'summary'},
            'Encoding identity does not match canonical query/reference strings and coverage')
    for key in expected_sequence_hashes:
        sequence = (suite['sequences'][key.removeprefix('reference-')] if key.startswith('reference-') else
                    query_lookup[key.removeprefix('query-')]['sequence'])
        validate_state_record(completed_unit(encoding, key), sequence)
    original = verify_model_probe(Path(encoding.identity["model_probe"]["path"]), session, design, suite)
    require(evidence_of(original) == encoding.identity["model_probe"] and encoding.identity["method"] == original.identity["model"],
            "Completed model source changed")
    admission = None if probe_mode else validate_admission(args.admission, session, file_hash(design_path), sources,
                                                          design=design, stage="alignment")
    if admission is not None:
        require(encoding.identity["admission"] == admission, "Repeat encoding used another admission")
    label_path = session / "inputs/family_labels_v1.json"
    labels = read(label_path)
    family_evidence = validate_family_freeze(labels, suite, label_path, session / "audit/independent_family_label_review.json")
    pairs = pair_list(query_ids, accessions)
    units = [pair_key(*pair) for pair in pairs] + ["timing-" + pair_key(*pair) for pair in pairs] + ["evaluation-" + q for q in query_ids] + ["startup", "summary"]
    deadline = time.monotonic() + min(session_seconds(session), design["limits"]["probe_total_seconds" if probe_mode else "repeat_alignment_seconds"])
    process = None
    store = None
    startups, startup_files, timings = [], [], {}
    try:
        process, ready, attempt = start_child(args, session, design, design_path, encoding, deadline)
        if admission is not None:
            require(ready["runtime"] == admission["runtime"], "Repeat alignment runtime differs from timing probe")
        identity = {"mode": "plmblast_alignment_v1", "panel": "probe12" if probe_mode else "repeat32",
                    "design_version": design["version"], "design_status": design["status"],
                    "design_sha256": file_hash(design_path), "source": sources, "code_hashes": code_hashes(),
                    "runtime": ready["runtime"], "encoding_source": evidence_of(encoding),
                    "method": encoding.identity["method"], "query_ids": query_ids, "reference_accessions": accessions,
                    "unit_keys": units, "family_evidence": family_evidence, "admission": admission,
                    "information_budget": "No supplied component boundaries; fixed total eligible unique-accession K",
                    "adaptation": design["deviations"]}
        store = new_store(session, args.name, identity)
        def bind(child, evidence, directory):
            require(evidence["runtime"] == identity["runtime"], "Restarted child runtime changed")
            child.stdin.write(json.dumps({"job": str(store.root), "identity_sha256": store.identity_hash}) + "\n")
            child.stdin.flush()
            record = {"ready": evidence, "ready_path": str(directory / "ready.json"), "ready_sha256": file_hash(directory / "ready.json")}
            atomic_json(store.root / ("startup-" + directory.name + ".json"), record)
            startups.append(evidence)
            startup_files.append(record)
        bind(process, ready, attempt)
        if store.read("startup") is None:
            store.publish("startup", ready)
        else:
            initial_startup = store.read("startup")
            require(initial_startup["fixtures"] == ready["fixtures"] and initial_startup["runtime"] == ready["runtime"],
                    "Startup fixture/runtime changed on resume")
            startups.append(initial_startup)
        for qid, accession in pairs:
            key = pair_key(qid, accession)
            if store.read(key) is not None:
                timing = store.read("timing-" + key)
                if timing is None:
                    # Never repeat a scored pair to replace missing timing after
                    # a crash. Its resource gate is explicitly incomplete.
                    timing = {"status": "failed", "seconds": None, "exception_type": "MissingCompletedWallTiming"}
                    store.publish("timing-" + key, timing)
                timings[key] = timing["seconds"]
                continue
            started = time.monotonic()
            gate_status, exception_type, attempted = "passed", None, False
            try:
                require(time.monotonic() < deadline, "Frozen alignment/session deadline reached")
                if process is None:
                    process, ready, attempt = start_child(args, session, design, design_path, encoding, deadline)
                    bind(process, ready, attempt)
                process.stdin.write(json.dumps({"query_id": qid, "accession": accession}) + "\n")
                process.stdin.flush()
                attempted = True
                message = child_message(process, min(design["limits"]["pair_seconds"], deadline - time.monotonic()))
                require(message["event"] == "pair_complete" and message["key"] == key, "Unexpected pair completion message")
                result = store.read(key)
                require(result is not None and result["status"] == message["status"], "Missing complete pair unit")
                require(result["peak_process_rss_mib"] <= design["limits"]["maximum_worker_peak_mib"], "Pair memory gate exceeded")
            except Exception as error:
                gate_status, exception_type = "failed", type(error).__name__
                stop_child(process)
                process = None
                if store.read(key) is None:
                    store.publish(key, {"status": "failed", "query_id": qid, "accession": accession,
                                        "exception_type": type(error).__name__, "attempted": attempted,
                                        "stage": "parent_deadline_process_or_resource_gate"})
            timings[key] = max(time.monotonic() - started, 1e-9)
            store.publish("timing-" + key, {"seconds": timings[key], "status": gate_status,
                                          "exception_type": exception_type, "attempted": attempted})
        stop_child(process)
        process = None
        lookup = {q["query_id"]: q for q in suite["queries"]}
        evaluations, pair_rows = [], []
        gate_failures = [pair_key(*pair) for pair in pairs if store.read("timing-" + pair_key(*pair))["status"] != "passed"]
        for qid in query_ids:
            rows = [store.read(pair_key(qid, a)) for a in accessions]
            pair_rows.extend(rows)
            failed = [r["accession"] for r in rows if r["status"] == "failed" or pair_key(qid, r["accession"]) in gate_failures]
            query = lookup[qid]
            record = {"query": public_query(query), "suite": query["suite"], "status": "failed" if failed else "success",
                      "failed_candidates": failed, "expected_candidates": len(accessions),
                      "candidate_statuses": {r["accession"]: r["status"] for r in rows},
                      "scope": "timing probe partial catalogue; not benchmark performance" if probe_mode else "fixed repeat32 full105 catalogue"}
            if not failed:
                ranking = ranking_from_pairs(rows)
                record.update(ranking=ranking, evaluation=eligible_metrics(ranking, accessions, query, suite["equivalence"], labels["acceptable_matches"]))
            if store.read("evaluation-" + qid) is None:
                store.publish("evaluation-" + qid, record)
            else:
                require(store.read("evaluation-" + qid) == record, "Recomputed eligible metrics differ")
            evaluations.append(record)
        failed_count = sum(row["status"] == "failed" for row in pair_rows)
        existing_summary = store.read('summary')
        startup_files = summary_startups(store, session, existing_summary, design['limits']['maximum_worker_peak_mib'])
        startups = [entry['ready'] for entry in startup_files]
        forecast = None
        if probe_mode and not failed_count and not gate_failures:
            forecast = forecast_repeat(pair_rows, timings, startups, original, suite, design)
        summary = {"expected_pair_count": len(pairs), "failed_pair_count": failed_count,
                   "failed_or_unadmitted_pair_count": len({pair_key(r["query_id"], r["accession"]) for r in pair_rows if r["status"] == "failed"} | set(gate_failures)),
                   "resource_gate_failures": gate_failures, "no_eligible_alignment_pair_count": sum(r["status"] == "no_eligible_alignment" for r in pair_rows),
                   "forecast": forecast, "startup_evidence": startup_files,
                   "pair_timings": timings, "summaries": None if probe_mode else summarize_eligible(evaluations),
                   "resource_recommendation": "eligible_for_explicit_repeat32_admission" if forecast and forecast["within_frozen_budget"] else "not_admitted",
                   "interpretation": "Exploratory adapted existing-method result; no function/screening validation; native no-hits are scorable misses"}
        if existing_summary is None:
            store.publish("summary", summary)
        else:
            require(existing_summary == summary, 'Recomputed pLM-BLAST summary differs from saved scientific/resource evidence')
        store.complete(units)
        return store
    finally:
        try:
            stop_child(process)
        finally:
            stop_owned_children()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=["encode-probe", "encode-repeat", "align-probe", "align-repeat", "worker"])
    parser.add_argument("--session", type=Path)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--expected-code-sha256")
    parser.add_argument("--expected-design-sha256")
    parser.add_argument("--name")
    parser.add_argument("--model-probe-job", type=Path)
    parser.add_argument("--encoding-job", type=Path)
    parser.add_argument("--admission", type=Path)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--worker-spec", type=Path)
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    if args.overlay is not None:
        require("numpy" not in sys.modules and "pandas" not in sys.modules, "Overlay must precede numerical imports")
        args.overlay = args.overlay.resolve(strict=True)
        sys.path.insert(0, str(args.overlay))
    if args.mode == "worker":
        require(args.worker_spec is not None and args.overlay is not None, "Worker requires explicit owned spec and overlay")
        worker(args.worker_spec, args.overlay)
        return
    require(args.session is not None and args.name and args.expected_code_sha256 and args.expected_design_sha256,
            "Session, name, frozen code and design checksums are mandatory")
    require(not args.mode.startswith("encode") or args.overlay is None, "Model encoding must use its original environment, not alignment overlay")
    session, design, design_path, suite, sources = setup(args)
    store = (encode if args.mode.startswith("encode") else alignment)(args, session, design, design_path, suite, sources)
    print(json.dumps({"status": "complete", "job": str(store.root), "mode": args.mode}), flush=True)


if __name__ == "__main__":
    main()
