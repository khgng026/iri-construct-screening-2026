#!/usr/bin/env python3
"""Exploratory residue changes versus artificial-join distance, offline/resumable.

Only 16 standalone residue matrices persist. Joined states are reduced one query
at a time. This is descriptive, uses supplied artificial boundaries, and is not
an additive explanation of rank loss or an automatic boundary detector.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
from iri_constructs.session_store import UnitStore, atomic_json, canonical, digest, file_hash, utcnow
from iri_constructs.session_benchmark import load_suite
from session_experiment import check_deadline, mean_residues, model_identity, peak_mib, public_query
from session_segment_experiment import open_store, completed_read, validate_source_suite

BINS = ((0, 8), (8, 32), (32, 64), (64, 128), (128, None))
BIN_NAMES = ("0-7", "8-31", "32-63", "64-127", "128+")
QUANTILES = (0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0)
TOLERANCE = 1e-6


def require(condition, message):
    if not condition:
        raise ValueError(message)


def code_identity():
    names = ["scripts/session_context_profile.py", "scripts/session_experiment.py",
             "scripts/session_segment_experiment.py", "src/iri_constructs/session_models.py",
             "src/iri_constructs/session_benchmark.py", "src/iri_constructs/session_store.py"]
    return {name: file_hash(ROOT / name) for name in names}


def close_vectors(actual, expected, label):
    import numpy as np
    a, b = np.asarray(actual), np.asarray(expected)
    require(a.shape == b.shape and a.size > 0 and np.isfinite(a).all() and np.isfinite(b).all(),
            label + ": invalid shape or value")
    error = float(np.max(np.abs(a.astype(np.float64) - b.astype(np.float64))))
    require(error <= TOLERANCE, label + f": max_abs={error} exceeds {TOLERANCE}")
    return error


def distance_bins(length, left_join, right_join):
    """Zero is a residue next to a join; biological termini are not joins."""
    import numpy as np
    require(type(length) is int and length > 0, "Invalid component length")
    require(type(left_join) is bool and type(right_join) is bool, "Join flags must be booleans")
    if not left_join and not right_join:
        return np.full(length, -1, dtype=np.int8)
    offsets = np.arange(length)
    distance = (np.minimum(offsets, length - 1 - offsets) if left_join and right_join
                else offsets if left_join else length - 1 - offsets)
    return np.searchsorted([8, 32, 64, 128], distance, side="right").astype(np.int8)


def scalar_summary(values):
    import numpy as np
    values = np.asarray(values, dtype=np.float64)
    require(values.ndim == 1 and values.size > 0 and np.isfinite(values).all(), "Invalid scalar observations")
    return {"count": len(values), "sum": float(values.sum(dtype=np.float64)),
            "mean": float(values.mean()), "quantile_probabilities": list(QUANTILES),
            "quantiles": np.quantile(values, QUANTILES, method="linear").tolist()}


def residue_changes(contextual, isolated):
    import numpy as np
    a, b = np.asarray(contextual, dtype=np.float64), np.asarray(isolated, dtype=np.float64)
    require(a.ndim == 2 and a.shape == b.shape and min(a.shape) > 0 and
            np.isfinite(a).all() and np.isfinite(b).all(), "Invalid paired residue matrices")
    an, bn = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    require((an > 0).all() and (bn > 0).all(), "Cosine/relative L2 undefined for a zero-norm residue")
    cosine = np.einsum("ij,ij->i", a, b) / (an * bn)
    require(np.max(np.abs(cosine)) <= 1 + 1e-12, "Invalid cosine arithmetic")
    return np.clip(cosine, -1, 1), np.linalg.norm(a - b, axis=1) / bn


class PrecisionAudit:
    """Observe actual tensors before the frozen adapter casts returned arrays."""
    def __init__(self, encoder):
        self.encoder = encoder
        self.torch = encoder.torch
        require(encoder.config["dtype"] == "float32", "Profile requires declared FP32")
        require(encoder.layer == encoder.num_layers, "Profile requires the parent's final representation")
        family = encoder.config["family"]
        if family == "esm3":
            require(encoder.config.get("representation_variant") == "official_embeddings_prenorm",
                    "ESM3 profile only supports the official primary final representation")
        self.raw_field = "embeddings" if family == "esm3" else "last_hidden_state"
        self.device = encoder.config["device"]
        self.raw = []
        params, buffers = defaultdict(int), defaultdict(int)
        names = []
        for name, tensor in encoder.model.named_parameters():
            if tensor.is_floating_point():
                require(tensor.dtype == self.torch.float32 and tensor.device.type == self.device,
                        "Floating parameter dtype/device differs: " + name)
                params[f"{tensor.dtype}@{tensor.device}"] += tensor.numel()
                names.append(name)
        require(names, "Model has no floating parameters")
        for name, tensor in encoder.model.named_buffers():
            if tensor.is_floating_point():
                buffers[f"{tensor.dtype}@{tensor.device}"] += tensor.numel()
        self.metadata = {"floating_parameter_counts_by_dtype_device": dict(params),
                         "floating_parameter_name_sha256": digest(canonical(names)),
                         "floating_buffer_counts_by_dtype_device": dict(buffers),
                         "raw_output_field": self.raw_field, "declared_dtype": "float32",
                         "declared_device": self.device,
                         "limitation": "Checks storage/output precision and disabled autocast, not every internal arithmetic operation."}
        self.handle = encoder.model.register_forward_hook(self._observe)

    def _observe(self, module, args, output):
        tensor = getattr(output, self.raw_field, None)
        require(isinstance(tensor, self.torch.Tensor), "Missing raw model representation")
        require(tensor.dtype == self.torch.float32 and tensor.device.type == self.device,
                "Raw forward output differs from declared dtype/device before adapter conversion")
        require(not tensor.requires_grad and not module.training, "Forward unexpectedly uses gradients/training")
        self.raw.append({"dtype": str(tensor.dtype), "device": str(tensor.device),
                         "device_type": tensor.device.type, "shape": list(tensor.shape),
                         "requires_grad": tensor.requires_grad, "model_training": module.training})

    def encode(self, sequence):
        import numpy as np
        self.raw.clear()
        require(not self.torch.is_autocast_enabled(self.device), "Autocast must be disabled")
        result = self.encoder.encode(sequence)
        require(len(self.raw) == 1, "Expected exactly one observed model forward")
        raw = self.raw[0]
        require(result.ndim == 2, "Adapter returned a non-matrix representation")
        require(raw["shape"] == [1, len(sequence) + 2, result.shape[1]], "Raw representation token shape differs")
        require(result.dtype == np.float32 and result.shape[0] == len(sequence), "Invalid adapter result")
        # Preserve small evidence only, never a second reference to raw tensors.
        return result, {**raw, "autocast_enabled": False, "adapter_return_dtype": str(result.dtype)}

    def close(self):
        self.handle.remove()


def profile_query(query, states, isolated_states, parent_arrays):
    """Produce exact per-component/bin summaries; discard residue deltas on return."""
    import numpy as np
    require(states.ndim == 2 and states.shape[0] == query["length_aa"] and states.dtype == np.float32,
            "Joined residue shape/dtype mismatch")
    rows, contextual_means, isolated_means, checks = [], [], [], {}
    ncomponents = len(query["components"])
    checks["whole_parent_max_abs"] = close_vectors(mean_residues(states), parent_arrays["whole"], "Parent whole mean")
    for index, span in enumerate(query["spans"]):
        accession = span["accession"]
        a = states[span["start0"]:span["end0"]]
        b = isolated_states[accession]
        require(a.shape == b.shape, "Component/standalone residue alignment differs")
        bins = distance_bins(len(a), index > 0, index < ncomponents - 1)
        cosine, relative = residue_changes(a, b)
        for bin_index in sorted(set(bins.tolist())):
            mask = bins == bin_index
            row = {"row_index": len(rows), "accession": accession, "component_index": index,
                   "bin_index": int(bin_index), "bin": BIN_NAMES[bin_index] if bin_index >= 0 else "clean_no_artificial_join",
                   "residue_count": int(mask.sum()), "component_length": len(a),
                   "left_artificial_join": index > 0, "right_artificial_join": index < ncomponents - 1,
                   "cosine": scalar_summary(cosine[mask]), "relative_l2": scalar_summary(relative[mask])}
            rows.append(row)
            contextual_means.append(mean_residues(a[mask]))
            isolated_means.append(mean_residues(b[mask]))
        checks[accession] = {"contextual_parent_max_abs": close_vectors(mean_residues(a), parent_arrays["contextual_" + accession], "Parent contextual mean"),
                             "isolated_parent_max_abs": close_vectors(mean_residues(b), parent_arrays["isolated_" + accession], "Parent isolated mean")}
        component_rows = [r for r in rows if r["component_index"] == index]
        for kind, vectors, fresh in (("contextual", contextual_means, a), ("isolated", isolated_means, b)):
            reconstructed = sum(r["residue_count"] * vectors[r["row_index"]].astype(np.float64)
                                for r in component_rows) / len(a)
            checks[accession][kind + "_bin_fresh_mean_max_abs"] = close_vectors(
                reconstructed, mean_residues(fresh), kind + " bins versus fresh component mean")
        if ncomponents == 1:
            checks["clean_residue_repeat_max_abs"] = close_vectors(a, b, "Clean repeated residue states")
    record = {"query": public_query(query), "scope": "clean" if ncomponents == 1 else "combined",
              "bins": rows, "checks": checks, "tolerance": TOLERANCE}
    arrays = {"contextual_bin_means": np.stack(contextual_means), "isolated_bin_means": np.stack(isolated_means)}
    validate_profile_unit(record, arrays, query, isolated_states, parent_arrays)
    return record, arrays


def validate_profile_unit(record, arrays, query, isolated_states, parent_arrays):
    """Also used on resume: rederive bin geometry and parent reconstruction."""
    import numpy as np
    def check_value(container, name, value):
        if name in container:
            require(abs(container[name] - value) <= 1e-12, "Saved reconstruction check differs")
        else:
            container[name] = value
    require(record["query"] == public_query(query), "Profile query identity differs")
    rows = record["bins"]
    require(len(rows) > 0 and [r["row_index"] for r in rows] == list(range(len(rows))), "Invalid bin row indices")
    expected_scope = "clean" if len(query["components"]) == 1 else "combined"
    require(record["scope"] == expected_scope and record["tolerance"] == TOLERANCE, "Profile scope/tolerance mismatch")
    require(set(arrays) == {"contextual_bin_means", "isolated_bin_means"}, "Unexpected profile arrays")
    dimension = next(iter(isolated_states.values())).shape[1]
    for array in arrays.values():
        require(array.shape == (len(rows), dimension) and array.dtype == np.float32 and np.isfinite(array).all(),
                "Invalid bin mean arrays")
    expected_row_keys = []
    reconstructed_components = []
    for index, span in enumerate(query["spans"]):
        accession = span["accession"]
        length = span["end0"] - span["start0"]
        for name in ("contextual_parent_max_abs", "isolated_parent_max_abs",
                     "contextual_bin_fresh_mean_max_abs", "isolated_bin_fresh_mean_max_abs"):
            value = record["checks"][accession][name]
            require(np.isfinite(value) and 0 <= value <= TOLERANCE, "Invalid saved fresh reconstruction evidence")
        bins = distance_bins(length, index > 0, index < len(query["components"]) - 1)
        selected = [r for r in rows if r["component_index"] == index]
        require([r["bin_index"] for r in selected] == sorted(set(bins.tolist())), "Missing/repeated distance bin")
        for row in selected:
            mask = bins == row["bin_index"]
            expected_row_keys.append((index, row["bin_index"]))
            require(row["accession"] == accession and row["component_length"] == length and
                    row["residue_count"] == int(mask.sum()) and row["left_artificial_join"] == (index > 0) and
                    row["right_artificial_join"] == (index < len(query["components"]) - 1), "Wrong bin support/joins")
            label = BIN_NAMES[row["bin_index"]] if row["bin_index"] >= 0 else "clean_no_artificial_join"
            require(row["bin"] == label, "Wrong distance label")
            close_vectors(arrays["isolated_bin_means"][row["row_index"]],
                          mean_residues(isolated_states[accession][mask]), "Isolated bin mean")
            for metric in ("cosine", "relative_l2"):
                stats = row[metric]
                require(stats["count"] == row["residue_count"] and stats["quantile_probabilities"] == list(QUANTILES), "Scalar count/quantiles mismatch")
                values = [stats["sum"], stats["mean"], *stats["quantiles"]]
                require(np.isfinite(values).all() and len(stats["quantiles"]) == len(QUANTILES) and
                        stats["quantiles"] == sorted(stats["quantiles"]) and
                        abs(stats["sum"] / stats["count"] - stats["mean"]) <= 1e-12, "Invalid scalar summary")
                lo, hi = stats["quantiles"][0], stats["quantiles"][-1]
                require(lo <= stats["mean"] + 1e-12 and stats["mean"] <= hi + 1e-12 and
                        (lo >= -1 and hi <= 1 if metric == "cosine" else lo >= 0), "Scalar range mismatch")
        require(sum(r["residue_count"] for r in selected) == length, "Bins do not partition component")
        means = {}
        for kind in ("contextual", "isolated"):
            means[kind] = sum(r["residue_count"] * arrays[kind + "_bin_means"][r["row_index"]].astype(np.float64)
                              for r in selected) / length
            error = close_vectors(means[kind], parent_arrays[kind + "_" + accession], kind + " bin reconstruction versus parent")
            if kind == "isolated":
                close_vectors(means[kind], mean_residues(isolated_states[accession]), "Isolated bin reconstruction")
            check_value(record["checks"].setdefault(accession, {}), kind + "_bin_parent_max_abs", error)
        reconstructed_components.append((length, means["contextual"]))
    require([(r["component_index"], r["bin_index"]) for r in rows] == expected_row_keys, "Unexpected component rows")
    whole = sum(n * v for n, v in reconstructed_components) / query["length_aa"]
    check_value(record["checks"], "binned_whole_parent_max_abs",
                close_vectors(whole, parent_arrays["whole"], "All bins versus parent whole mean"))


def summarize_profiles(records):
    """Retain scalar observation summaries only; quantile averaging is forbidden."""
    import numpy as np
    strata, matched = defaultdict(list), defaultdict(list)
    for record in records:
        q = record["query"]
        by_component = defaultdict(dict)
        for row in record["bins"]:
            observation = {"target": row["accession"], "group": q["group_id"], "query": q["query_id"],
                           "count": row["residue_count"], "cosine": row["cosine"], "relative_l2": row["relative_l2"]}
            for suite in ("all", q["suite"]):
                strata[(record["scope"], suite, row["bin"])].append(observation)
            by_component[row["accession"]][row["bin_index"]] = row
        if record["scope"] == "combined":
            for target, bins in by_component.items():
                if 0 in bins and 4 in bins:
                    matched[(q["group_id"], target)].append({
                        "query": q["query_id"], "cosine_far_minus_near": bins[4]["cosine"]["mean"] - bins[0]["cosine"]["mean"],
                        "relative_l2_far_minus_near": bins[4]["relative_l2"]["mean"] - bins[0]["relative_l2"]["mean"]})
    summaries = []
    for (scope, suite, bin_name), observations in sorted(strata.items()):
        count = sum(o["count"] for o in observations)
        result = {"scope": scope, "suite": suite, "bin": bin_name, "residue_observations": count,
                  "component_bin_observations": len(observations), "queries": len({o["query"] for o in observations}),
                  "groups": len({o["group"] for o in observations}), "targets": sorted({o["target"] for o in observations}),
                  "target_count": len({o["target"] for o in observations}),
                  "group_target_pairs": len({(o["group"], o["target"]) for o in observations})}
        for metric in ("cosine", "relative_l2"):
            target_values, group_target_values = defaultdict(list), defaultdict(list)
            for o in observations:
                target_values[o["target"]].append(o[metric]["mean"])
                group_target_values[(o["group"], o["target"])].append(o[metric]["mean"])
            result[metric] = {"residue_weighted_mean": sum(o[metric]["sum"] for o in observations) / count,
                "component_bin_mean_distribution": scalar_summary([o[metric]["mean"] for o in observations]),
                "target_macro_mean": float(np.mean([np.mean(v) for v in target_values.values()])),
                "group_target_macro_mean": float(np.mean([np.mean(v) for v in group_target_values.values()]))}
        summaries.append(result)
    matched_rows = []
    for (group, target), values in sorted(matched.items()):
        matched_rows.append({"group_id": group, "accession": target, "matched_query_count": len(values),
                             "query_ids": sorted(v["query"] for v in values),
                             **{metric: float(np.mean([v[metric] for v in values])) for metric in
                                ("cosine_far_minus_near", "relative_l2_far_minus_near")}})
    return {"by_scope_suite_bin": summaries, "matched_far128plus_vs_near0to7_by_group_target": matched_rows,
            "quantile_scope": "Per-unit quantiles are over residues; aggregate quantiles are over component-bin means, never pooled residue quantiles.",
            "interpretation": "Exploratory descriptive profile of supplied component boundaries; repeated orders are correlated; contrasts are nonadditive and not causal rank-loss decompositions."}


def select_queries(suite, parent, design, subset):
    require(design["version"] == "context-profile-v1" and
            design["bins_distance_to_nearest_artificial_join_residues"] == [list(x) for x in BINS], "Frozen profile design changed")
    require(parent.identity["mode"] == "model", "Parent must be a completed model core job")
    if subset == "full":
        require(parent.root.name == design["model_first"], "Full profile is frozen to the first8M parent; use separately labelled repeat for another model")
        wanted = {q["query_id"] for q in suite["queries"]}
    else:
        require(subset == "repeat", "Unknown profile subset")
        wanted = set(suite["manifest"]["repeat_query_ids"])
    require(wanted <= set(parent.identity["query_ids"]), "Parent lacks required profile queries")
    return sorted([q for q in suite["queries"] if q["query_id"] in wanted],
                  key=lambda q: ({"clean": 0, "pairs": 1, "triples": 2, "quadruples": 3}[q["suite"]], q["query_id"]))


def validate_runtime_identity(encoder, parent):
    require(model_identity(encoder) == parent.identity["method"], "Loaded model/environment/representation differs from explicit parent")
    for name, version in parent.identity["environment"].items():
        require(importlib.metadata.version(name) == version, "Runtime package differs from parent: " + name)
    for name in ("src/iri_constructs/session_models.py", "src/iri_constructs/session_benchmark.py", "src/iri_constructs/session_store.py"):
        require(parent.identity["code_hashes"][name] == file_hash(ROOT / name), "Frozen parent dependency changed: " + name)


def validate_precision_record(record, length, dimension, device):
    require(record["dtype"] == "torch.float32" and record["device_type"] == device and
            record["shape"] == [1, length + 2, dimension] and record["requires_grad"] is False and
            record["model_training"] is False and record["autocast_enabled"] is False and
            record["adapter_return_dtype"] == "float32", "Missing or incompatible raw precision evidence")


def run(args):
    import numpy as np
    from iri_constructs.session_models import load_encoder
    session = args.session.resolve()
    check_deadline(session)
    require(file_hash(__file__) == args.expected_code_sha256, "Profile code differs from explicit precompute freeze")
    initial_code = code_identity()
    design_path = args.design or session / "inputs/context_profile_design_v1.json"
    design = json.loads(design_path.read_text())
    suite = load_suite(ROOT, session / "inputs")
    parent = open_store(args.parent_job, session)
    validate_source_suite(parent.identity, suite, file_hash(session / "inputs/challenge_manifest.json"))
    queries = select_queries(suite, parent, design, args.subset)
    targets = sorted(suite["manifest"]["targets"])
    require(len(targets) == 16, "Frozen standalone panel must contain16targets")
    encoder = precision = store = None
    start = time.monotonic()
    try:
        encoder = load_encoder(parent.identity["config"])
        validate_runtime_identity(encoder, parent)
        precision = PrecisionAudit(encoder)
        identity = {"mode": "context_profile_v1", "design_sha256": file_hash(design_path), "design": design,
            "parent_job": str(parent.root), "parent_identity_hash": parent.identity_hash,
            "parent_completion_sha256": parent.completion_sha256, "config": parent.identity["config"],
            "method": model_identity(encoder), "precision_audit": precision.metadata, "code_hashes": initial_code,
            "challenge_sha256": parent.identity["challenge_sha256"], "subset": args.subset,
            "query_ids": [q["query_id"] for q in queries], "targets": targets,
            "sequence_sha256": {a: digest(suite["sequences"][a].encode()) for a in targets},
            "quantiles": list(QUANTILES), "tolerance": TOLERANCE,
            "python_executable": sys.executable, "python_version": sys.version,
            "numerical_policy": "per-residue comparisons and reductions float64; stored means float32",
            "storage_policy": "standalone residue states only; joined matrices and per-residue deltas never persisted"}
        job_id = args.name + "-" + digest(canonical(identity))[:12]
        store = UnitStore(session / "jobs" / job_id, identity, session / "job_ledger.jsonl")
        store.event(event="job_started", state="running", job_id=job_id, pid=os.getpid(), query_count=len(queries))
        atomic_json(store.root / ("invocation-" + str(os.getpid()) + ".json"),
                    {"argv": sys.argv, "started_utc": utcnow(), "pid": os.getpid()})
        isolated, expected = {}, []
        for accession in targets:
            check_deadline(session)
            key = "standalone-" + accession
            expected.append(key)
            sequence = suite["sequences"][accession]
            source_record, source_arrays = completed_read(parent, "catalogue-" + accession, arrays=True)
            require(source_record["sequence_sha256"] == identity["sequence_sha256"][accession] and
                    source_record["length"] == len(sequence), "Parent standalone source identity mismatch")
            saved = store.read(key, arrays=True)
            if saved is None:
                states, raw = precision.encode(sequence)
                error = close_vectors(mean_residues(states), source_arrays["mean"], "Standalone versus parent catalogue")
                record = {"accession": accession, "sequence_sha256": identity["sequence_sha256"][accession],
                          "length": len(sequence), "precision": raw, "parent_mean_max_abs": error}
                store.publish(key, record, {"residues": states, "mean": mean_residues(states)})
            else:
                record, saved_arrays = saved
                require(set(saved_arrays) == {"residues", "mean"}, "Unexpected standalone arrays")
                states = saved_arrays["residues"]
                close_vectors(mean_residues(states), saved_arrays["mean"], "Standalone cache mean")
            require(record["accession"] == accession and record["sequence_sha256"] == identity["sequence_sha256"][accession]
                    and record["length"] == len(sequence) and states.ndim == 2 and states.shape[0] == len(sequence)
                    and states.dtype == np.float32 and np.isfinite(states).all(), "Invalid standalone cache")
            close_vectors(mean_residues(states), source_arrays["mean"], "Standalone cached versus parent catalogue")
            validate_precision_record(record["precision"], len(sequence), states.shape[1], encoder.config["device"])
            isolated[accession] = states
        for index, query in enumerate(queries):
            check_deadline(session)
            require(code_identity() == initial_code, "Profile source changed during run")
            key = "query-" + query["query_id"]
            expected.append(key)
            parent_record, parent_arrays = completed_read(parent, key, arrays=True)
            require(parent_record["query"] == public_query(query), "Parent query metadata mismatch")
            saved = store.read(key, arrays=True)
            if saved is None:
                tick = time.monotonic()
                states, raw = precision.encode(query["sequence"])
                record, arrays = profile_query(query, states, isolated, parent_arrays)
                del states
                record.update(precision=raw, seconds=time.monotonic() - tick, peak_process_rss_mib=peak_mib())
                store.publish(key, record, arrays)
            else:
                record, arrays = saved
                validate_profile_unit(record, arrays, query, isolated, parent_arrays)
            validate_precision_record(record["precision"], query["length_aa"], arrays["contextual_bin_means"].shape[1], encoder.config["device"])
            if index % 25 == 0 or index + 1 == len(queries):
                print(json.dumps({"job_id": job_id, "completed_queries": index + 1,
                    "expected_queries": len(queries), "elapsed_seconds": time.monotonic() - start,
                    "peak_rss_mib": peak_mib()}), flush=True)
        summary = summarize_profiles(store.read("query-" + q["query_id"]) for q in queries)
        expected.append("summary")
        old = store.read("summary")
        if old is None:
            store.publish("summary", {"summary": summary, "query_count": len(queries),
                "standalone_target_count": len(targets), "elapsed_seconds": time.monotonic() - start,
                "peak_process_rss_mib": peak_mib()})
        else:
            require(old["summary"] == summary, "Recomputed profile summary differs from immutable saved summary")
        require(code_identity() == initial_code, "Profile source changed before completion")
        done = store.complete(expected, {"query_count": len(queries), "parent_identity_hash": parent.identity_hash,
                                        "peak_process_rss_mib": peak_mib()})
        print(json.dumps({"status": "complete", "job_id": job_id, "units": done["expected_count"]}), flush=True)
    except BaseException as exc:
        if store:
            store.event(event="job_failed", state="failed", exception=f"{type(exc).__name__}: {exc}", traceback=traceback.format_exc())
        raise
    finally:
        if precision:
            precision.close()
        if encoder:
            encoder.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--session", type=Path, required=True)
    parser.add_argument("--parent-job", type=Path, required=True)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--subset", choices=("full", "repeat"), default="full")
    parser.add_argument("--name", required=True)
    parser.add_argument("--expected-code-sha256", required=True)
    run(parser.parse_args())
