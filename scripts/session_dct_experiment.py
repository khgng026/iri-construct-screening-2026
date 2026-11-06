#!/usr/bin/env python3
"""Bounded native DCTdomain extraction, interoperability probes and retrieval.

All execution must be scheduled inside the lead's resource guard. This script
never compiles code, downloads weights, installs packages, or repairs upstream
algorithm failures. Encoding and ranking are separate immutable jobs.
"""
from __future__ import annotations
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
for thread_variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[thread_variable] = "4"
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash, atomic_json, utcnow
from iri_constructs.session_benchmark import load_suite, rank_scores, evaluate_query, summarize
from iri_constructs.session_dct_native import (NativeDCTModel, ObservedModel, verify_sources, load_module,
    window_lengths, fingerprint_embedding, parse_domains, boundary_metrics, exhaustive_distances, require)
from iri_constructs.session_dctdomain_adapter import DCTdomainHFModel, MODEL_ID, REVISION, WEIGHT_SHA256
from session_experiment import check_deadline, public_query, peak_mib


def code_hashes():
    names = ["scripts/session_dct_experiment.py", "src/iri_constructs/session_dct_native.py",
             "src/iri_constructs/session_dctdomain_adapter.py", "src/iri_constructs/session_models.py",
             "src/iri_constructs/session_benchmark.py", "src/iri_constructs/session_store.py",
             "scripts/session_experiment.py"]
    return {name: file_hash(ROOT / name) for name in names}


def open_completed(path, session):
    path = Path(path).resolve()
    store = UnitStore(path, json.loads((path / "identity.json").read_text()), session / "job_ledger.jsonl")
    candidates = [(json.loads(p.read_text())["completed_utc"], p) for p in path.glob("completion-*.json")]
    require(bool(candidates), "Source DCT job is incomplete")
    completion_path = max(candidates)[1]
    completion = json.loads(completion_path.read_text())
    expected = store.identity["unit_keys"]
    require(len(set(expected)) == len(expected) and completion["identity_hash"] == store.identity_hash and
            set(completion["units"]) == set(expected) and completion["expected_count"] == len(expected), "DCT completion coverage differs")
    store.completed_paths = {}
    for key, entry in completion["units"].items():
        unit = (path / entry["path"]).resolve()
        require(unit.is_relative_to(path) and file_hash(unit / "complete.json") == entry["descriptor_sha256"], "DCT descriptor changed")
        require(store.verify(unit)["key"] == key, "DCT unit key changed")
        store.completed_paths[key] = unit
    store.completion_sha256 = file_hash(completion_path)
    return store


def completed_unit(store, key, arrays=False):
    import numpy as np
    path = store.completed_paths[key]
    record = json.loads((path / "record.json").read_text())
    if not arrays:
        return record
    require("arrays.npz" in store.verify(path)["files"], "DCT arrays not checksum covered")
    with np.load(path / "arrays.npz", allow_pickle=False) as data:
        return record, {name: data[name].copy() for name in data.files}


def build_gate(path, roots, repo):
    evidence = json.loads(path.read_text())
    require(evidence["status"] == "passed", "RecCut build/fixture gate has not passed")
    for fixture in ("continuous", "discontinuous", "constant_quantization", "span_validation", "ce_integerization"):
        require(evidence["fixtures"].get(fixture) in (True, "passed"), "Unpassed RecCut fixture: " + fixture)
    binary = roots["dctdomain"] / "src/RecCut"
    require(os.access(binary, os.X_OK) and file_hash(binary) == evidence["binary_sha256"], "RecCut binary differs from build gate")
    require(file_hash(roots["dctdomain"] / "src/RecCut.cpp") == evidence["source_sha256"], "RecCut source differs from build gate")
    require(file_hash(roots["dctdomain"].parent / "compat/bits/stdc++.h") == evidence["compatibility_header_sha256"], "RecCut compatibility include differs")
    require(bool(evidence.get("compiler")), "Build gate lacks compiler evidence")
    return {"evidence_sha256": file_hash(path), "binary_sha256": evidence["binary_sha256"],
            "compiler": evidence["compiler"], "fixtures": evidence["fixtures"]}


def setup(args):
    session = args.session.resolve()
    check_deadline(session)
    require(file_hash(__file__) == args.expected_code_sha256, "DCT runner differs from explicit code freeze")
    design_path = args.design or session / "inputs/dctdomain_design_v1.json"
    require(file_hash(design_path) == args.expected_design_sha256, "DCT design differs from explicit freeze")
    design = json.loads(design_path.read_text())
    require(design["version"] == "dctdomain-v1" and design["model_id"] == MODEL_ID and
            design["model_revision"] == REVISION and design["checkpoint_sha256"] == WEIGHT_SHA256 and
            design["layers"] == [15, 21] and design["maxlen"] == 500 and design["overlap"] == 200 and
            design["contact_multiplier"] == 2.6 and design["quantization_dimensions"] == [3, 80, 3, 80], "DCT method settings changed")
    require(design["device"] == "cpu" and design["dtype"] == "float32" and design["batch_size"] == 1
            and type(design["threads"]) is int and 1 <= design["threads"] <= 4, "DCT resources/configuration changed")
    suite = load_suite(ROOT, session / "inputs")
    require(file_hash(session / "inputs/challenge_manifest.json") == design["source_challenge_sha256"], "DCT challenge changed")
    return session, design, suite, design_path


def select_probe(suite, design):
    by_id = {q["query_id"]: q for q in suite["queries"]}
    queries = []
    for item in design["probe_inputs"]:
        query = by_id[item["query_id"]]
        require(query["length_aa"] == item["length"] and query["joined_sequence_sha256"] == item["sequence_sha256"], "Frozen probe input differs")
        queries.append(query)
    require([q["length_aa"] for q in queries] == [67, 248, 642, 1011], "Bounded probe panel differs")
    return queries


def validate_encoded(record, arrays, length, sequence_hash):
    import numpy as np
    require(record["status"] == "success" and record["length"] == length and record["sequence_sha256"] == sequence_hash, "DCT cached input mismatch")
    require(record["domains"] == parse_domains([d["domain"] for d in record["domains"]], length), "DCT domain interpretation changed")
    fingerprints = arrays["fingerprints"]
    require(fingerprints.shape == (len(record["domains"]), 480) and np.issubdtype(fingerprints.dtype, np.integer)
            and ((0 <= fingerprints) & (fingerprints <= 127)).all(), "DCT fingerprint cache invalid")
    require([w["residue_count"] for w in record["windows"]] == window_lengths(length), "DCT cached windows differ")
    require(all(w["raw_dtype"] == "torch.float32" and w["device"] == "cpu" and w["autocast"] is False for w in record["windows"]), "DCT raw precision evidence differs")


def encode_one(sequence, backend, embedding_class, fingerprint_class, temp_root, capture):
    import numpy as np
    sha = digest(sequence.encode())
    observed = ObservedModel(backend, capture=capture)
    embedding = embedding_class(pid="dct_" + sha[:24], seq=sequence)
    embedding.embed_seq(observed, "cpu", [15, 21], 500)
    require([w["residue_count"] for w in observed.windows] == window_lengths(len(sequence)), "Original window count differs")
    for layer in (15, 21):
        matrix = embedding.embed[layer]
        require(matrix.shape == (len(sequence), 640) and matrix.dtype == np.float32 and np.isfinite(matrix).all(), "Invalid stitched native state")
    require(embedding.contacts.shape == (len(sequence), len(sequence)) and np.isfinite(embedding.contacts).all(), "Invalid stitched contact map")
    record, arrays = fingerprint_embedding(embedding, fingerprint_class, temp_root, capture=capture)
    record.update(status="success", length=len(sequence), sequence_sha256=sha, windows=observed.windows,
                  parameter_dtype_evidence=observed.parameter_evidence)
    if capture:
        arrays.update(observed.arrays)
        arrays.update(stitched_layer15=embedding.embed[15], stitched_layer21=embedding.embed[21], stitched_contacts=embedding.contacts)
    validate_encoded(record, arrays, len(sequence), sha)
    return record, arrays


def hf_public_probe(backend, sequence, recorded_arrays):
    import numpy as np
    torch = backend.torch
    _, _, tokens = backend.esm_tokenizer([("public_api_probe", sequence)])
    mask = tokens.ne(backend.tokenizer.pad_token_id).long()
    result = {}
    with torch.no_grad():
        output = backend.model(input_ids=tokens, attention_mask=mask, output_hidden_states=True, output_attentions=False, return_dict=True)
        for layer in (15, 21):
            observed = output.hidden_states[layer].detach().numpy()
            expected = recorded_arrays[f"window0_layer{layer}"]
            delta = np.abs(observed.astype(np.float64) - expected.astype(np.float64))
            require(np.all(delta <= 1e-5 + 1e-5 * np.abs(expected)), "HF hook/public hidden state differs")
            result[f"layer{layer}_max_abs"] = float(delta.max())
        del output
        contacts = backend.model.predict_contacts(tokens, mask).detach().numpy()
        delta = np.abs(contacts.astype(np.float64) - recorded_arrays["window0_contacts"].astype(np.float64))
        require(np.all(delta <= 1e-5 + 1e-5 * np.abs(recorded_arrays["window0_contacts"])), "HF shim/public contacts differ")
        result["contacts_max_abs"] = float(delta.max())
    return {"status": "passed", "atol": 1e-5, "rtol": 1e-5, **result}


def check_admission(path, session, design_hash, sources, build):
    require(path is not None, "Native full encoding requires explicit completed probe/resource admission")
    gate = json.loads(path.read_text())
    require(gate["status"] == "passed" and gate["resource_gate"] == "passed", "Native resource gate not passed")
    probe = open_completed(Path(gate["native_probe_job"]), session)
    require(probe.completion_sha256 == gate["native_probe_completion_sha256"] and
            probe.identity["mode"] == "dct_probe_v1" and probe.identity["backend"] == "native" and
            probe.identity["design_sha256"] == design_hash and probe.identity["source_hashes"] == sources and
            probe.identity["build"] == build and probe.identity["code_hashes"] == code_hashes(), "Native admission does not match frozen source/probe")
    require(completed_unit(probe, "summary")["failed_units"] == [], "Native probe contains failures")
    return {"gate_sha256": file_hash(path), "native_probe_identity_hash": probe.identity_hash,
            "native_probe_completion_sha256": probe.completion_sha256, "native_probe_method": probe.identity["method"]}


def encode(args, session, design, suite, design_path):
    probe = args.mode.startswith("probe-")
    hf = args.mode == "probe-hf"
    roots, sources = verify_sources(ROOT, session)
    require(args.build_evidence is not None, "DCT extraction requires native build evidence")
    build = build_gate(args.build_evidence, roots, ROOT)
    admission = None if probe else check_admission(args.admission, session, file_hash(design_path), sources, build)
    queries = select_probe(suite, design) if probe else sorted(suite["queries"], key=lambda q: ({"clean": 0, "pairs": 1, "triples": 2, "quadruples": 3}[q["suite"]], q["query_id"]))
    targets = [] if probe else sorted(suite["sequences"])
    planned = [("catalogue-" + a, suite["sequences"][a], {"accession": a}) for a in targets]
    planned += [("query-" + q["query_id"], q["sequence"], {"query": public_query(q)}) for q in queries]
    backend = None
    started = time.monotonic()
    initial_code = code_hashes()
    try:
        import torch
        torch.set_num_interop_threads(1)
        backend = (DCTdomainHFModel({"threads": design["threads"]}) if hf
                   else NativeDCTModel(roots["fair_esm"], threads=design["threads"]))
        if hf:
            require(importlib.metadata.version("transformers") == "4.57.6", "Only independently reviewed HF4.57.6 probe is admitted")
        embedding_module = load_module(roots["dctdomain"] / "src/embedding.py", "_iri_pinned_dct_embedding")
        fingerprint_module = load_module(roots["dctdomain"] / "src/fingerprint.py", "_iri_pinned_dct_fingerprint")
        method = {k: v for k, v in backend.metadata.items() if k != "load_seconds"}
        method["packages"]["scipy"] = importlib.metadata.version("scipy")
        method["actual_thread_configuration"] = {"torch_intraop": torch.get_num_threads(),
                                                   "torch_interop": torch.get_num_interop_threads()}
        require(1 <= torch.get_num_threads() <= 4 and torch.get_num_interop_threads() == 1, "Actual DCT thread bound differs")
        if hf:
            method["token_vocabulary"] = [token for token, index in sorted(backend.tokenizer.get_vocab().items(), key=lambda x: x[1])]
        if admission:
            require(method == admission["native_probe_method"], "Native runtime differs from admitted probe")
        identity = {"mode": "dct_probe_v1" if probe else "dct_encode_v1", "backend": "hf" if hf else "native",
            "design_sha256": file_hash(design_path), "design": design, "method": method, "source_hashes": sources,
            "build": build, "admission": admission, "code_hashes": initial_code,
            "challenge_sha256": design["source_challenge_sha256"], "query_ids": [q["query_id"] for q in queries],
            "catalogue_accessions": targets, "catalogue_sequence_sha256": {a: digest(suite["sequences"][a].encode()) for a in targets},
            "unit_keys": [key for key, _, _ in planned] + ["summary"], "python_executable": sys.executable,
            "python_version": sys.version, "temporary_raw_sequence_policy": "Owned ignored temporary CE only, deleted after each input; never logs"}
        store = UnitStore(session / "jobs" / (args.name + "-" + digest(canonical(identity))[:12]), identity, session / "job_ledger.jsonl")
        store.event(event="job_started", state="running", pid=os.getpid(), mode=identity["mode"])
        atomic_json(store.root / ("invocation-" + str(os.getpid()) + ".json"), {"argv": sys.argv, "started_utc": utcnow()})
        for index, (key, sequence, public) in enumerate(planned):
            check_deadline(session)
            require(code_hashes() == initial_code, "DCT source changed during extraction")
            saved = store.read(key)
            if saved is not None:
                require(saved["status"] in ("success", "failed") and all(saved[k] == v for k, v in public.items()), "DCT cache public identity differs")
                require(saved["sequence_sha256"] == digest(sequence.encode()) and saved["length"] == len(sequence), "DCT cached failed/successful input differs")
                if saved["status"] == "success":
                    record, arrays = store.read(key, arrays=True)
                    validate_encoded(record, arrays, len(sequence), digest(sequence.encode()))
                continue
            tick = time.monotonic()
            try:
                record, arrays = encode_one(sequence, backend, embedding_module.Embedding, fingerprint_module.Fingerprint, session / "tmp/dctdomain", capture=probe)
                if hf and index == 0:
                    record["hf_public_interface"] = hf_public_probe(backend, sequence, arrays)
                record.update(public, seconds=time.monotonic() - tick, peak_process_rss_mib=peak_mib())
                store.publish(key, record, arrays)
            except Exception as exc:
                # No exception text/traceback: upstream diagnostics can contain inputs.
                store.publish(key, {**public, "status": "failed", "exception_type": type(exc).__name__,
                    "failure_stage": getattr(exc, "iri_dct_stage", "representation_or_output_validation"),
                    "failure_policy": "No silent algorithm repair or input omission; raw exception text intentionally excluded",
                    "length": len(sequence), "sequence_sha256": digest(sequence.encode()),
                    "seconds": time.monotonic() - tick, "peak_process_rss_mib": peak_mib()})
                store.event(event="dct_input_failed", key=key, exception_type=type(exc).__name__)
                if isinstance(exc, MemoryError) or "out of memory" in str(exc).lower():
                    raise RuntimeError("DCT resource failure; completed input checkpoints retained") from None
            if index % 25 == 0 or index + 1 == len(planned):
                print(json.dumps({"job": str(store.root), "processed": index + 1, "planned": len(planned), "peak_rss_mib": peak_mib()}), flush=True)
        failed = [key for key, _, _ in planned if store.read(key)["status"] != "success"]
        summary = {"planned_inputs": len(planned), "successful_inputs": len(planned) - len(failed), "failed_units": failed,
                   "query_count": len(queries), "catalogue_count": len(targets), "scientific_coverage_complete": not failed}
        old = store.read("summary")
        if old is None:
            store.publish("summary", {**summary, "elapsed_seconds": time.monotonic() - started, "peak_process_rss_mib": peak_mib()})
        else:
            require(all(old[k] == v for k, v in summary.items()), "DCT extraction summary changed")
        store.complete(identity["unit_keys"], {"failed_unit_count": len(failed), "scientific_coverage_complete": not failed})
        print(json.dumps({"status": "execution_complete", "job": str(store.root), **summary}), flush=True)
        return store.root
    finally:
        if backend is not None:
            backend.close()


def whole_fingerprint(record, arrays):
    matches = [i for i, d in enumerate(record["domains"]) if d["domain"] == f"1-{record['length']}"]
    return arrays["fingerprints"][matches] if len(matches) == 1 else None


def ranking_from_distances(minima):
    accessions = sorted(minima)
    ranking = rank_scores([-minima[a] for a in accessions], accessions, score_name="negative_raw_l1")
    for row in ranking:
        row["raw_l1_distance"] = minima[row["accession"]]
        row["unclipped_display_similarity"] = 1 - minima[row["accession"]] / 17000
    return ranking


def rank(args, session, design, suite, design_path, encoding_path):
    source = open_completed(encoding_path, session)
    require(source.identity["mode"] == "dct_encode_v1" and source.identity["backend"] == "native" and
            source.identity["design_sha256"] == file_hash(design_path) and
            source.identity["catalogue_accessions"] == sorted(suite["sequences"]) and
            source.identity["catalogue_sequence_sha256"] == {a: digest(s.encode()) for a, s in suite["sequences"].items()}, "DCT encoding source differs from full frozen suite")
    queries = {q["query_id"]: q for q in suite["queries"]}
    require(set(source.identity["query_ids"]) == set(queries), "DCT query coverage differs")
    reference_arrays, global_arrays, failed_references, missing_global = {}, {}, [], []
    reference_domain_records = {}
    for accession in sorted(suite["sequences"]):
        key = "catalogue-" + accession
        record = completed_unit(source, key)
        if record["status"] != "success":
            failed_references.append(accession)
            continue
        record, arrays = completed_unit(source, key, arrays=True)
        validate_encoded(record, arrays, len(suite["sequences"][accession]), digest(suite["sequences"][accession].encode()))
        reference_arrays[accession] = arrays["fingerprints"]
        reference_domain_records[accession] = record["domains"]
        whole = whole_fingerprint(record, arrays)
        if whole is None:
            missing_global.append(accession)
        else:
            global_arrays[accession] = whole
    identity = {"mode": "dct_ranking_v1", "source_job": str(source.root), "source_identity_hash": source.identity_hash,
        "source_completion_sha256": source.completion_sha256, "design_sha256": file_hash(design_path), "design": design,
        "code_hashes": code_hashes(), "catalogue_accessions": sorted(suite["sequences"]),
        "query_ids": source.identity["query_ids"], "unit_keys": ["query-" + qid for qid in source.identity["query_ids"]] + ["summary"],
        "reference_domains": reference_domain_records, "failed_reference_accessions": failed_references,
        "missing_global_reference_accessions": missing_global,
        "budget": "Fixed total unique-accessionK per original query; adapted from native domain-pair output budget"}
    store = UnitStore(session / "jobs" / (args.name + "-ranking-" + digest(canonical(identity))[:12]), identity, session / "job_ledger.jsonl")
    for query_id in identity["query_ids"]:
        check_deadline(session)
        key = "query-" + query_id
        if store.find(key):
            continue
        query = queries[query_id]
        encoded = completed_unit(source, key)
        record = {"query": public_query(query), "rankings": {}, "evaluations": {}, "unavailable_arms": {}}
        arrays = {}
        if encoded["status"] != "success" or failed_references:
            reason = "query_encoding_failed" if encoded["status"] != "success" else "incomplete_frozen_reference_catalogue"
            record.update(status="unscorable", reason=reason)
            record["unavailable_arms"] = {"domain_pairs": reason, "global_only": reason}
        else:
            encoded, fingerprints = completed_unit(source, key, arrays=True)
            validate_encoded(encoded, fingerprints, query["length_aa"], query["joined_sequence_sha256"])
            minima, pairs = exhaustive_distances(fingerprints["fingerprints"], reference_arrays)
            record["rankings"]["domain_pairs"] = ranking_from_distances(minima)
            arrays["domain_pair_distances"] = pairs
            record.update(status="success", fingerprint_count=encoded["fingerprint_count"],
                query_domains=encoded["domains"], domain_pair_comparisons=len(pairs),
                distance_array_columns=["reference_accession_index", "query_fingerprint_index", "reference_fingerprint_index", "raw_l1_distance"],
                distance_row_order="raw_distance,accession,query_fingerprint_index,reference_fingerprint_index",
                boundary_evaluation=boundary_metrics(encoded["domains"], query["spans"]))
            whole = whole_fingerprint(encoded, fingerprints)
            if whole is None or missing_global:
                record["unavailable_arms"]["global_only"] = "missing_existing_whole_fingerprint"
            else:
                global_minima, global_pairs = exhaustive_distances(whole, global_arrays)
                record["rankings"]["global_only"] = ranking_from_distances(global_minima)
                arrays["global_pair_distances"] = global_pairs
            for arm, ranking in record["rankings"].items():
                record["evaluations"][arm] = evaluate_query(ranking, query, suite["equivalence"])
        store.publish(key, record, arrays if arrays else None)
    records = [store.read("query-" + qid) for qid in identity["query_ids"]]
    summary = {"planned_queries": len(records), "planned_component_observations": sum(len(r["query"]["components"]) for r in records),
               "failed_reference_accessions": failed_references, "missing_global_reference_accessions": missing_global, "arms": {}}
    for arm in ("domain_pairs", "global_only"):
        values = [r["evaluations"][arm] for r in records if arm in r["evaluations"]]
        excluded = [r["query"]["query_id"] for r in records if arm not in r["evaluations"]]
        summary["arms"][arm] = {"scorable_query_count": len(values), "unscorable_query_ids": excluded,
            "unscorable_component_observations": sum(len(r["query"]["components"]) for r in records if arm not in r["evaluations"]),
            "by_suite": {s: summarize([v for v in values if v["suite"] == s]) for s in ("clean", "pairs", "triples", "quadruples")}}
    old = store.read("summary")
    if old is None:
        store.publish("summary", summary)
    else:
        require(old == summary, "DCT ranking summary changed")
    store.complete(identity["unit_keys"], {"planned_queries": len(records), "full_scientific_coverage": all(not a["unscorable_query_ids"] for a in summary["arms"].values())})
    print(json.dumps({"status": "ranking_execution_complete", "job": str(store.root), "planned_queries": len(records)}), flush=True)
    return store.root


def compare_probes(args, session, design, design_path):
    """Sequential saved-output comparison; never loads either model."""
    import numpy as np
    require(args.native_probe_job is not None and args.hf_probe_job is not None, "Both explicit probe jobs are required")
    native = open_completed(args.native_probe_job, session)
    hf = open_completed(args.hf_probe_job, session)
    for source, backend in ((native, "native"), (hf, "hf")):
        require(source.identity["mode"] == "dct_probe_v1" and source.identity["backend"] == backend and
                source.identity["design_sha256"] == file_hash(design_path), "Probe source/design differs")
    require(native.identity["query_ids"] == hf.identity["query_ids"] and native.identity["source_hashes"] == hf.identity["source_hashes"]
            and native.identity["build"] == hf.identity["build"] and
            native.identity["method"]["token_vocabulary"] == hf.identity["method"]["token_vocabulary"], "Probe vocabulary/source/build differs")
    for package in ("torch", "numpy", "scipy", "safetensors"):
        require(native.identity["method"]["packages"][package] == hf.identity["method"]["packages"][package],
                "Interoperability probes use different runtime packages: " + package)
    identity = {"mode": "dct_interoperability_v1", "native_identity_hash": native.identity_hash,
        "native_completion_sha256": native.completion_sha256, "hf_identity_hash": hf.identity_hash,
        "hf_completion_sha256": hf.completion_sha256, "design_sha256": file_hash(design_path),
        "code_hashes": code_hashes(), "unit_keys": native.identity["unit_keys"], "query_ids": native.identity["query_ids"]}
    store = UnitStore(session / "jobs" / (args.name + "-" + digest(canonical(identity))[:12]), identity, session / "job_ledger.jsonl")
    atol, rtol = design["equivalence_gate"]["atol"], design["equivalence_gate"]["rtol"]
    for query_id in identity["query_ids"]:
        check_deadline(session)
        key = "query-" + query_id
        if store.find(key):
            continue
        nr, hr = completed_unit(native, key), completed_unit(hf, key)
        if nr["status"] != "success" or hr["status"] != "success":
            store.publish(key, {"query_id": query_id, "status": "unavailable", "native_status": nr["status"], "hf_status": hr["status"]})
            continue
        nr, na = completed_unit(native, key, arrays=True)
        hr, ha = completed_unit(hf, key, arrays=True)
        require(nr["sequence_sha256"] == hr["sequence_sha256"] and set(na) == set(ha), "Probe input/array contract differs")
        discrete = {"domains_equal": nr["domains"] == hr["domains"], "tokens_equal": nr["windows"] == hr["windows"]}
        errors = {}
        for name in sorted(na):
            a, b = ha[name], na[name]
            require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), "Probe tensor shape/finiteness differs")
            if name.startswith("stitched_") or name.startswith("window"):
                require(a.dtype == b.dtype == np.float32, "Probe raw saved state dtype differs")
                delta = np.abs(a.astype(np.float64) - b.astype(np.float64))
                bound = atol + rtol * np.abs(b.astype(np.float64))
                errors[name] = {"max_abs": float(delta.max()), "p99_abs": float(np.quantile(delta, .99)),
                    "max_tolerance_ratio": float((delta / bound).max()), "failing_elements": int(np.sum(delta > bound)),
                    "element_count": delta.size, "passed": bool(np.all(delta <= bound))}
                if name.startswith("window") and "layer" in name:
                    actual = delta[:, 1:-1]
                    errors[name]["actual_residue_max_abs"] = float(actual.max())
                    errors[name]["actual_residue_failing_elements"] = int(np.sum(actual > bound[:, 1:-1]))
            else:
                discrete[name + "_equal"] = bool(np.array_equal(a, b))
        passed = all(discrete.values()) and all(row["passed"] for row in errors.values())
        store.publish(key, {"query_id": query_id, "status": "equivalent" if passed else "not_equivalent",
            "float_errors": errors, "discrete_comparisons": discrete, "atol": atol, "rtol": rtol,
            "reference": "native", "native_admission_independent": True})
    records = [store.read("query-" + qid) for qid in identity["query_ids"]]
    summary = {"query_count": len(records), "equivalent": all(r["status"] == "equivalent" for r in records),
               "status_by_query": {r["query_id"]: r["status"] for r in records},
               "interpretation": "HF interoperability only; failure does not automatically invalidate independently admitted native extraction."}
    if store.find("summary") is None:
        store.publish("summary", summary)
    else:
        require(store.read("summary") == summary, "Interoperability summary changed")
    store.complete(identity["unit_keys"], {"equivalent": summary["equivalent"]})
    print(json.dumps({"job": str(store.root), **summary}), flush=True)


def main(args):
    session, design, suite, design_path = setup(args)
    if args.mode == "compare-probes":
        compare_probes(args, session, design, design_path)
    elif args.mode == "rank":
        require(args.encoding_job is not None, "Ranking requires an explicit completed encoding job")
        rank(args, session, design, suite, design_path, args.encoding_job)
    else:
        job = encode(args, session, design, suite, design_path)
        if args.mode == "full":
            rank(args, session, design, suite, design_path, job)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--session", type=Path, required=True)
    parser.add_argument("--mode", choices=("probe-native", "probe-hf", "compare-probes", "encode-native", "rank", "full"), required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--expected-code-sha256", required=True)
    parser.add_argument("--expected-design-sha256", required=True)
    parser.add_argument("--build-evidence", type=Path)
    parser.add_argument("--admission", type=Path)
    parser.add_argument("--encoding-job", type=Path)
    parser.add_argument("--native-probe-job", type=Path)
    parser.add_argument("--hf-probe-job", type=Path)
    main(parser.parse_args())
