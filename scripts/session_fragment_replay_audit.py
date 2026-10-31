#!/usr/bin/env python3
"""Independent, guarded representative fresh-fragment alignment replay.

``freeze`` is a stdlib-only inventory of existing local evidence. ``run`` needs
explicit lead admission under the prepared guard. No model or production
alignment/selection implementation is imported. Each fixed case is attempted at
most once; complete accounting is distinct from an all-passed audit outcome.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import session_segmentation_audit as saved

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.session_store import UnitStore, atomic_json, canonical, digest, safe_name, utcnow

base = saved.base
read, sha, require = base.read, base.sha, base.require
VERSION = "independent-fresh-fragment-replay-v1"
DESIGN_NAME = "fragment_replay_audit_design_v1.json"
PREDICTION_ID = "d4c478a398b2a2dfa0c3614c0e8fc120a4d6fbe9a4523edb248a8563e55d2259"
PREDICTION_COMPLETION = "aa75152146f00c6ae3bf4798be5238c4020a0e9f9fa901e29e61f912c88158f2"
STORED_AUDIT_SHA = "92b80d886168a077f52d6029b71adbb3f91a3db5ef8f2b14aee1fb67e50b4394"
PROPOSAL_SHA = "43b793a272b399b0fe77a5524e683fc601c3e42aa4a1b93698763e7611a54609"
GUARD_SHA = "a5435e14ce302c8d7a78844baed1e015ce526931965a33f960552716d439673d"
CONFIG = "alignment_cap4_max_score"
POLICIES = ("exact_removed", "near90_removed")
SUITES = ("pairs", "triples", "quadruples")
LIMITS = {"maximum_cases": 7, "case_seconds": 60, "overall_seconds": 600,
          "maximum_combined_parent_child_rss_bytes": 512 * 1024**2,
          "maximum_threads": 4, "rss_poll_seconds": .1}
CODE_FILES = ("scripts/session_fragment_replay_audit.py", "scripts/session_segmentation_audit.py",
              "scripts/session_extension_audit.py", "scripts/session_independent_audit.py",
              "src/iri_constructs/session_store.py", "tests/test_session_fragment_replay_audit.py")
LIMITATION = ("Representative replay only. Unselected fresh-fragment attempts retain the v3 stored-audit limitation. "
             "Old full fresh-fragment rankings were not recorded, so their row-for-row agreement cannot be asserted. "
             "This is neither an audit of all 1,933 attempts nor a new retrieval-performance estimate.")


def exact(actual, expected, where="value"):
    """Finite numeric equality, including exact half-integer scores and ties."""
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and set(actual) == set(expected), where + ": fields differ")
        for key in expected:
            exact(actual[key], expected[key], where + "." + key)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), where + ": length differs")
        for i, value in enumerate(expected):
            exact(actual[i], value, where + "." + str(i))
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        require(isinstance(actual, (int, float)) and not isinstance(actual, bool) and
                math.isfinite(actual) and math.isfinite(expected) and actual == expected, where + ": number differs")
    else:
        require(type(actual) is type(expected) and actual == expected, where + ": value differs")


def peak_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


class Budget(base.Deadline):
    def check(self):
        super().check()
        require(peak_bytes() <= LIMITS["maximum_combined_parent_child_rss_bytes"], "Parent memory bound exceeded")


def selection_hash(case):
    value = "iri-fresh-fragment-replay-v1|{}|{}|{}|{}|{}".format(
        case["policy"], case["suite"], case["query_id"], case["start0"], case["end0"])
    return digest(value.encode())


def select_cases(pool):
    """Selection depends only on frozen identifiers/lengths, never outcomes."""
    unique = {}
    for case in pool:
        key = (case["fragment_sequence_sha256"], case["included_accessions_sha256"])
        if key not in unique or selection_hash(case) < selection_hash(unique[key]):
            unique[key] = case
    representatives = sorted(unique.values(), key=selection_hash)
    chosen, strata = [], []
    for policy in POLICIES:
        for suite in SUITES:
            candidates = [c for c in representatives if c["policy"] == policy and c["suite"] == suite]
            choice = candidates[0] if candidates else None
            strata.append({"policy": policy, "suite": suite, "available_distinct_cases": len(candidates),
                           "selected_hash": selection_hash(choice) if choice else None,
                           "status": "selected" if choice else "unavailable_no_replacement"})
            if choice is not None:
                chosen.append(choice)
    longest = min(representatives, key=lambda c: (-c["length_aa"], selection_hash(c))) if representatives else None
    stress_added = longest is not None and longest not in chosen
    if stress_added:
        chosen.append(longest)
    require(len(chosen) <= 7, "Too many fixed cases")
    return [{**c, "case_id": "case-" + selection_hash(c)[:16], "selection_sha256": selection_hash(c),
             "selection_roles": (["stratum"] if any(s["selected_hash"] == selection_hash(c) for s in strata) else []) +
                                (["longest_fragment"] if c is longest else [])} for c in chosen], {
        "fresh_attempts": len(pool), "distinct_fragment_mask_cases": len(representatives),
        "deduplication": "global (fragment sequence SHA256, included-accession SHA256); minimum selection hash is canonical representative",
        "strata": strata, "longest_selection_sha256": selection_hash(longest) if longest else None,
        "longest_already_selected": longest is not None and not stress_added,
        "selected_cases": len(chosen), "no_outcome_dependent_selection": True}


def environment_fingerprint():
    """Read installed distribution metadata/files without importing Bio/numpy."""
    files = {}
    versions = {}
    for package in ("biopython", "numpy"):
        dist = importlib.metadata.distribution(package)
        versions[package] = dist.version
        for relative in dist.files or []:
            name = str(relative)
            relevant = (name.endswith((".dist-info/METADATA", ".dist-info/RECORD")) or
                        (package == "biopython" and (name in ("Bio/__init__.py", "Bio/Align/__init__.py",
                            "Bio/Align/substitution_matrices/__init__.py", "Bio/Align/substitution_matrices/data/BLOSUM62") or
                            (name.startswith("Bio/Align/") and name.endswith(".so")))) or
                        (package == "numpy" and name.startswith("numpy/_core/") and name.endswith(".so")))
            if relevant:
                path = Path(dist.locate_file(relative)).resolve()
                files[str(path)] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    require(versions == {"biopython": "1.88", "numpy": "2.5.4"}, "Original alignment package versions required")
    require(any("_pairwisealigner" in p for p in files) and any("_alignmentcounts" in p for p in files), "Missing official alignment binaries")
    executable = Path(sys.executable).resolve()
    return {"versions": versions, "python_version": sys.version, "python_executable": str(executable),
            "python_executable_sha256": sha(executable), "prefix": str(Path(sys.prefix).resolve()), "files": files}


def all_sequences(ctx):
    result = {}
    for entry in read(ROOT / "data/manifests/t05_expanded_catalogue.json")["records"]:
        path = base.safe_file(ROOT, entry["raw_path"])
        require(sha(path) == entry["checksum_sha256"], "Approved raw source changed")
        document = read(path)
        if "results" in document:
            require(len(document["results"]) == 1, "Ambiguous source record")
            document = document["results"][0]
        value = document["sequence"].get("value", document["sequence"].get("sequence"))
        a = entry["accession"]
        require(isinstance(value, str) and set(value) <= set("ACDEFGHIKLMNPQRSTVWY") and
                digest(value.encode()) == ctx["hashes"][a] and len(value) == ctx["lengths"][a], "Approved sequence identity differs")
        result[a] = value
    return result


def inventory(session, deadline):
    """Verify complete source snapshots and rebuild the score-blind sample."""
    old_path = session / "audit/independent_segmentation_audit_v3.json"
    require(sha(old_path) == STORED_AUDIT_SHA and read(old_path)["status"] == "passed", "Passed v3 audit changed")
    proposal = session / "audit/fresh_fragment_replay_proposal_v1.json"
    require(sha(proposal) == PROPOSAL_SHA, "Approved replay proposal changed")
    ctx = base.context(session, deadline)
    design, transfer, selection = saved.frozen_designs(session, ctx)
    sources = saved.Sources(session, deadline, ctx)
    prediction = sources.open(PREDICTION_ID, PREDICTION_COMPLETION)
    raw = sources.open(prediction.identity["raw_alignment_identity"], prediction.identity["raw_alignment_completion_sha256"])
    saved.verify_production_source(prediction.identity)
    saved.verify_production_source(raw.identity)
    require(raw.identity["mode"] == "alignment" and set(raw.identity["query_ids"]) == set(ctx["queries"]), "Incomplete original query lookup")
    exact(raw.identity["environment"], {"biopython": "1.88", "numpy": "2.5.4"}, "Original environment")
    config = next(c for c in prediction.identity["configs"] if c["configuration_id"] == CONFIG)
    require(config in saved.selected_configs(design, selection) and len(design["test_query_ids"]) == 950,
            "Frozen development winner or held-out coverage differs")
    require(prediction.identity["primary_design_sha256"] == saved.DESIGN_SHA and
            prediction.identity["transfer_design_sha256"] == saved.TRANSFER_SHA and
            prediction.identity["selection_sha256"] == saved.SELECTION_SHA, "Source transfer pins differ")
    known = {q["joined_sequence_sha256"] for q in ctx["queries"].values()}
    sequences = all_sequences(ctx)
    pool, full_panel_fresh, excluded_clean_fresh = [], 0, 0
    for qid in design["test_query_ids"] + design["clean_sanity_query_ids"]:
        deadline.check()
        query = ctx["queries"][qid]
        sequence = "".join(sequences[a] for a in query["components"])
        for policy in POLICIES:
            key = qid + "--" + policy + "--" + CONFIG
            record = prediction.record(key)
            exact(record["query"], query, "Source query")
            mask = saved.common_mask(ctx, query, policy)
            exact(record["common_mask"], mask, "Rebuilt all-component mask")
            exact(record["prediction_catalogue_accessions"], mask["included_accessions"], "Prediction reference set")
            exact(record["prediction_config"], saved.geometry(config), "Prediction configuration")
            require(record["boundary_labels_used_for_prediction"] is False, "Source prediction used boundary labels")
            segments = record["prediction"]["segments"]
            for attempt in record["prediction"]["interval_attempts"]:
                a, b = attempt["start0"], attempt["end0"]
                require(type(a) is int and type(b) is int and 0 <= a < b <= len(sequence), "Invalid recorded fragment")
                fragment_sha = digest(sequence[a:b].encode())
                if fragment_sha in known:
                    continue
                full_panel_fresh += 1
                if query["suite"] == "clean":
                    excluded_clean_fresh += 1
                    continue
                accepted = [s for s in segments if s.get("seed_hit") is not None and s["source_interval"] == [a, b]]
                require(len(accepted) <= 1, "Fragment accepted more than once")
                unit = prediction.units[key]
                pool.append({"query_id": qid, "suite": query["suite"], "policy": policy,
                    "query_sequence_sha256": query["joined_sequence_sha256"], "query_components": query["components"],
                    "start0": a, "end0": b, "length_aa": b - a, "fragment_sequence_sha256": fragment_sha,
                    "included_accessions_sha256": mask["included_accessions_sha256"], "common_mask": mask,
                    "expected_attempt": attempt, "expected_accepted_segment": accepted[0] if accepted else None,
                    "source_unit_key": key, "source_unit_path": str(unit.relative_to(session)),
                    "source_descriptor_sha256": sha(unit / "complete.json"), "source_record_sha256": sha(unit / "record.json"),
                    "reference_lengths": {a: ctx["lengths"][a] for a in mask["included_accessions"]},
                    "reference_sequence_sha256": {a: ctx["hashes"][a] for a in mask["included_accessions"]},
                    "estimated_dp_cells": (b - a) * sum(ctx["lengths"][a] for a in mask["included_accessions"])})
    chosen, selection_evidence = select_cases(pool)
    # The prior v3 total includes the 16 clean sanity queries; the approved
    # representative sample is restricted to the 950 combined test queries.
    prior_family = next(j for j in read(old_path)["jobs"] if Path(j["job"]).name == "family-test-8m-infer-v1-d0424de04de8")
    require(full_panel_fresh == prior_family["trace_evidence"]["fresh_fragment_attempts_without_saved_full_ranking"] == 1933 and
            full_panel_fresh == len(pool) + excluded_clean_fresh and chosen, "Expected full-panel fresh-fragment accounting differs")
    selection_evidence.update(original_full_panel_fresh_attempts=full_panel_fresh,
                              excluded_clean_sanity_fresh_attempts=excluded_clean_fresh,
                              eligible_combined_test_queries=950, excluded_clean_sanity_queries=16)
    sources.unchanged()
    source_evidence = [snapshot.evidence() for snapshot in sources.snapshots.values()]
    return {"cases": chosen, "selection": selection_evidence, "parameters": config["parameters"],
            "original_alignment_method": raw.identity["method"], "source_snapshots": source_evidence,
            "frozen_input_evidence": ctx["evidence"], "original_stored_audit_sha256": STORED_AUDIT_SHA,
            "proposal_sha256": PROPOSAL_SHA}, sequences


def freeze(args):
    session = args.session.resolve()
    destination = session / "inputs" / DESIGN_NAME
    require(not destination.exists(), "Frozen design already exists; no replacement permitted")
    started = time.monotonic()
    data, _ = inventory(session, Budget(600, session))
    value = {"version": VERSION, "status": "frozen_before_any_fresh_alignment", "frozen_utc": utcnow(),
             "inventory_seconds": time.monotonic() - started, "limits": LIMITS, "guard_source_sha256": GUARD_SHA,
             "environment": environment_fingerprint(), "code_sha256": {p: sha(ROOT / p) for p in CODE_FILES},
             "data": data, "information_scope": "Source intervals and all-component masks replayed; no model or new evaluation labels",
             "limitation": LIMITATION, "no_replacements_or_retuning": True}
    atomic_json(destination, value)
    print(json.dumps({"status": value["status"], "design": str(destination), "sha256": sha(destination),
                      "cases": [{k: c[k] for k in ("case_id", "query_id", "suite", "policy", "start0", "end0", "length_aa")}
                                for c in data["cases"]], "selection": data["selection"]}), flush=True)
    return 0


def guard_admission():
    path, token = os.environ.get("IRI_GUARD_OWNERSHIP_PATH"), os.environ.get("IRI_GUARD_OWNERSHIP_TOKEN")
    require(path and token, "Explicit root-owned prepared guard is required")
    until = time.monotonic() + 5
    while True:
        try:
            record = read(path)
        except FileNotFoundError:
            record = None
        if record is not None:
            require(record["token"] == token and record["guard_source_sha256"] == GUARD_SHA, "Prepared guard identity differs")
            if any(e["pgid"] == os.getpgrp() and e["pid"] == e["pgid"] and e["start_new_session"]
                   for e in record["child_groups"]):
                return {"ownership_path": str(Path(path).resolve()), "guard_source_sha256": GUARD_SHA,
                        "owned_pgid": os.getpgrp(), "wrapper_pid": record["wrapper_pid"]}
        require(time.monotonic() < until, "Process group is not recorded by prepared guard")
        time.sleep(.02)


def independent_pair(aligner, query, reference):
    """First official optimum; metrics from Alignment.counts, not our old loop."""
    alignments = aligner.align(query, reference)
    score = float(alignments.score)
    require(math.isfinite(score) and score >= 0, "Invalid official alignment score")
    result = {"score": score, "query_length": len(query), "reference_length": len(reference),
              "query_span": None, "reference_span": None, "coordinates": [[], []], "aligned_length": 0,
              "aligned_residue_pairs": 0, "identities": 0, "identity_fraction": 0.0, "query_coverage": 0.0,
              "reference_coverage": 0.0, "query_paired_coverage": 0.0, "reference_paired_coverage": 0.0,
              "gap_columns": 0, "gap_runs": 0}
    if score == 0:
        return result
    alignment = next(iter(alignments))
    counts = alignment.counts(aligner)
    coords = alignment.coordinates.tolist()
    require(len(coords) == 2 and len(coords[0]) == len(coords[1]) >= 2, "Invalid official coordinates")
    for row, length in zip(coords, (len(query), len(reference))):
        require(all(type(v) is int for v in row) and 0 <= row[0] < row[-1] <= length and
                all(a <= b for a, b in zip(row, row[1:])), "Nonmonotone/out-of-range official coordinates")
    pairs, gaps, runs, identities = int(counts.aligned), int(counts.gaps), int(counts.open_gaps), int(counts.identities)
    columns = pairs + gaps
    require(pairs > 0 and 0 <= identities <= pairs and columns == alignment.shape[1], "Official count/shape mismatch")
    reconstructed = float(counts.substitution_score) - 10.0 * runs - .5 * (gaps - runs)
    require(float(counts.score) == score == reconstructed, "Official count/score reconstruction differs")
    qa, qb = coords[0][0], coords[0][-1]
    ra, rb = coords[1][0], coords[1][-1]
    result.update(query_span=[qa, qb], reference_span=[ra, rb], coordinates=coords,
                  aligned_length=columns, aligned_residue_pairs=pairs, identities=identities,
                  identity_fraction=identities / columns, query_coverage=(qb - qa) / len(query),
                  reference_coverage=(rb - ra) / len(reference), query_paired_coverage=pairs / len(query),
                  reference_paired_coverage=pairs / len(reference), gap_columns=gaps, gap_runs=runs,
                  reconstructed_score=reconstructed)
    return result


def rank_rows(rows):
    require(len({r["accession"] for r in rows}) == len(rows) and all(math.isfinite(r["score"]) for r in rows), "Invalid ranking candidates")
    result = sorted((dict(r) for r in rows), key=lambda r: (-r["score"], r["accession"]))
    scores = [r["score"] for r in result]
    for index, row in enumerate(result, 1):
        positions = [i + 1 for i, value in enumerate(scores) if value == row["score"]]
        row.update(rank=index, exact_score_tie_count=len(positions), tie_best_rank=min(positions),
                   tie_worst_rank=max(positions), near_tie_count_1e6=sum(abs(v - row["score"]) <= 1e-6 for v in scores))
    return result


def independent_candidates(rows, interval, parameters):
    result = []
    for row in rows:
        span = row["query_span"]
        if span is None:
            continue
        a, b = span
        require(0 <= a < b <= interval[1] - interval[0], "Candidate fragment coordinates invalid")
        if (b - a >= parameters["minimum_span"] and row["score"] >= parameters["minimum_score"] and
                row["reference_coverage"] >= parameters["minimum_reference_coverage"] and
                row["identity_fraction"] >= parameters["minimum_identity_fraction"]):
            result.append({"start0": interval[0] + a, "end0": interval[0] + b,
                           "source_interval": list(interval), "seed_hit": row})
    return sorted(result, key=lambda r: (-r["seed_hit"]["score"], r["seed_hit"]["accession"], r["start0"], r["end0"]))


def compare_case(case, rows, parameters):
    exact(sorted(r["accession"] for r in rows), case["common_mask"]["included_accessions"], "Full masked candidate set")
    exact(rows, rank_rows(rows), "Independent rank and tie metadata")
    candidates = independent_candidates(rows, (case["start0"], case["end0"]), parameters)
    best = candidates[0] if candidates else None
    observed = {"start0": case["start0"], "end0": case["end0"], "catalogue_count": len(rows),
                "eligible_reference_count": len(candidates), "best_eligible_accession": best["seed_hit"]["accession"] if best else None,
                "best_eligible_score": best["seed_hit"]["score"] if best else None}
    errors = []
    try:
        exact(observed, case["expected_attempt"], "Saved attempt versus fresh official replay")
    except ValueError as exc:
        errors.append(str(exc))
    accepted = case["expected_accepted_segment"]
    if accepted is not None:
        try:
            require(best is not None, "Saved accepted seed has no fresh eligible candidate")
            for key in ("start0", "end0", "source_interval", "seed_hit"):
                exact(best[key], accepted[key], "Accepted seed " + key)
        except ValueError as exc:
            errors.append(str(exc))
    return {"status": "failed" if errors else "passed", "errors": errors, "observed_attempt": observed,
            "best_eligible_candidate": best, "accepted_seed_geometry_checked": accepted is not None,
            "all_remaining_references_aligned": len(rows), "old_full_ranking_available": False,
            "fresh_complete_ranking_sha256": digest(canonical(rows))}


def official_aligner(method):
    from Bio import Align, __version__
    from Bio.Align import substitution_matrices
    require(__version__ == method["implementation_version"] == "1.88", "Biopython version differs")
    require(method["matrix"] == "BLOSUM62" and method["mode"] == "local" and method["gap_open"] == 10.0 and
            method["gap_extend"] == .5 and method["epsilon"] == 1e-6, "Original alignment settings differ")
    aligner = Align.PairwiseAligner(mode="local")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -.5
    aligner.epsilon = 1e-6
    require(aligner.algorithm == method["algorithm"] and
            digest(aligner.substitution_matrix.tobytes()) == method["matrix_sha256"] and
            aligner.substitution_matrix.alphabet == method["matrix_alphabet"], "Official algorithm/matrix identity differs")
    return aligner


def worker(request_path):
    """Private child: no subprocesses; inherited guard-owned process group."""
    guard_admission()
    request = read(request_path)
    require(os.getppid() == request["parent_pid"] == int(os.environ.get("IRI_FRAGMENT_REPLAY_PARENT_PID", "-1")), "Worker has no owning parent")
    require(os.getpgrp() == request["owned_pgid"] and sha(__file__) == request["code_sha256"], "Worker identity/group differs")
    design = read(request["design_path"])
    require(sha(request["design_path"]) == request["design_sha256"], "Worker design changed")
    exact(environment_fingerprint(), design["environment"], "Worker package/source fingerprint")
    case = next(c for c in design["data"]["cases"] if c["case_id"] == request["case_id"])
    session = Path(request["session"]).resolve()
    deadline = Budget(60, session)
    ctx = base.context(session, deadline)
    sequences = all_sequences(ctx)
    full = "".join(sequences[a] for a in case["query_components"])
    fragment = full[case["start0"]:case["end0"]]
    require(digest(full.encode()) == case["query_sequence_sha256"] and digest(fragment.encode()) == case["fragment_sequence_sha256"], "Worker input changed")
    aligner = official_aligner(design["data"]["original_alignment_method"])
    rows = []
    started = time.monotonic()
    for accession in case["common_mask"]["included_accessions"]:
        deadline.check()
        require(digest(sequences[accession].encode()) == case["reference_sequence_sha256"][accession] and
                len(sequences[accession]) == case["reference_lengths"][accession], "Worker reference changed")
        rows.append({"accession": accession, **independent_pair(aligner, fragment, sequences[accession])})
    rows = rank_rows(rows)
    comparison = compare_case(case, rows, design["data"]["parameters"])
    deadline.check()
    atomic_json(Path(request["output_path"]), {"case_id": case["case_id"], "comparison": comparison, "rankings": rows,
                "alignment_seconds": time.monotonic() - started, "peak_process_rss_bytes": peak_bytes(),
                "official_count_method": "Alignment.counts(aligner)", "code_sha256": sha(__file__)})
    return 0


def rss_sum(pids):
    output = subprocess.run(["/bin/ps", "-o", "pid=,rss=", "-p", ",".join(str(p) for p in pids)],
                            capture_output=True, text=True, timeout=1, check=True).stdout
    rows = {int(parts[0]): int(parts[1]) * 1024 for line in output.splitlines() if (parts := line.split())}
    require(os.getpid() in rows and set(rows) <= set(pids), "Invalid owned-process RSS observation")
    return sum(rows.values())


def execute_case(case, design_path, store, guard, deadline):
    tmp = store.root / "tmp" / case["case_id"]
    tmp.mkdir(parents=True, exist_ok=False)
    output = tmp / "observed.json"
    request = {"parent_pid": os.getpid(), "owned_pgid": guard["owned_pgid"], "code_sha256": sha(__file__),
               "design_path": str(design_path), "design_sha256": sha(design_path), "case_id": case["case_id"],
               "session": str(design_path.parent.parent), "output_path": str(output)}
    request_path = tmp / "request.json"
    atomic_json(request_path, request)  # identifiers/hashes only, no raw sequence
    env = dict(os.environ, IRI_FRAGMENT_REPLAY_PARENT_PID=str(os.getpid()))
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = "1"
    command = [sys.executable, str(Path(__file__).resolve()), "worker", "--request", str(request_path)]
    started, observed_rss, reason, child = time.monotonic(), peak_bytes(), None, None
    try:
        # No new session: the unchanged root guard owns this entire group.
        child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        require(os.getpgid(child.pid) == os.getpgrp(), "Child escaped owned process group")
        while child.poll() is None:
            deadline.check()
            if time.monotonic() - started >= LIMITS["case_seconds"]:
                reason = "case_wall_time_exceeded"
                break
            observed_rss = max(observed_rss, rss_sum([os.getpid(), child.pid]))
            if observed_rss > LIMITS["maximum_combined_parent_child_rss_bytes"]:
                reason = "combined_memory_bound_exceeded"
                break
            time.sleep(LIMITS["rss_poll_seconds"])
        if reason is None and child.poll() != 0:
            reason = "owned_child_failed"
    except (TimeoutError, ValueError, OSError, subprocess.SubprocessError) as exc:
        reason = "supervisor_" + type(exc).__name__
    finally:
        if child is not None and child.poll() is None:
            child.kill()
            child.wait(timeout=5)
    result = {"case_id": case["case_id"], "status": "failed", "reason": reason,
              "wall_seconds": time.monotonic() - started, "maximum_observed_combined_rss_bytes": observed_rss,
              "child_pid": child.pid if child else None, "child_exit_code": child.returncode if child else None,
              "child_reaped": child is not None and child.returncode is not None, "owned_pgid": guard["owned_pgid"],
              "request_sha256": sha(request_path), "command_sha256": digest(canonical(command)),
              "case_manifest_sha256": digest(canonical(case))}
    if reason is None and output.exists():
        observed = read(output)
        if "worker_error_type" in observed:
            result.update(reason="worker_" + observed["worker_error_type"], observed_file_sha256=sha(output))
            return result
        require(observed["case_id"] == case["case_id"] and observed["code_sha256"] == sha(__file__), "Worker output identity differs")
        comparison = compare_case(case, observed["rankings"], read(design_path)["data"]["parameters"])
        exact(comparison, observed["comparison"], "Parent independent comparison")
        result.update(observed, status=comparison["status"], observed_file_sha256=sha(output))
        # Sum of separate peaks is conservative and supplements sampled RSS.
        result["conservative_parent_child_peak_sum_bytes"] = peak_bytes() + observed["peak_process_rss_bytes"]
        if result["conservative_parent_child_peak_sum_bytes"] > LIMITS["maximum_combined_parent_child_rss_bytes"]:
            result.update(status="failed", reason="conservative_memory_bound_exceeded")
        if result["wall_seconds"] > LIMITS["case_seconds"]:
            result.update(status="failed", reason="case_wall_time_exceeded")
    elif reason is None:
        result["reason"] = "missing_owned_child_output"
    return result


def accounting(cases, records, pool):
    require(cases, "No fixed cases to account for")
    require([r["case_id"] for r in records] == [c["case_id"] for c in cases], "Missing/reordered case accounting")
    counts = Counter(r["status"] for r in records)
    require(set(counts) <= {"passed", "failed", "unattempted"}, "Unknown case status")
    return {"status": "passed" if counts["passed"] == len(cases) else "incomplete_or_failed",
            "expected_cases": len(cases), "attempted_cases": counts["passed"] + counts["failed"],
            "passed_cases": counts["passed"], "failed_cases": counts["failed"], "unattempted_cases": counts["unattempted"],
            "original_fresh_attempts": pool["original_full_panel_fresh_attempts"],
            "eligible_combined_fresh_attempts": pool["fresh_attempts"],
            "excluded_clean_sanity_fresh_attempts": pool["excluded_clean_sanity_fresh_attempts"],
            "distinct_fragment_mask_cases": pool["distinct_fragment_mask_cases"],
            "representative_cases_passed": counts["passed"], "unreplayed_distinct_cases": pool["distinct_fragment_mask_cases"] - counts["passed"],
            "all_1933_attempts_validated": False, "retrieval_accuracy_estimated": False, "limitation": LIMITATION}


def run(args):
    require(args.execute_fresh_replay, "Explicit fresh-alignment execution flag required")
    guard = guard_admission()
    started = time.monotonic()
    session = args.session.resolve()
    design_path = session / "inputs" / DESIGN_NAME
    require(sha(__file__) == args.expected_code_sha256 and sha(design_path) == args.expected_design_sha256, "Explicit replay freeze differs")
    design = read(design_path)
    require(design["version"] == VERSION and design["status"] == "frozen_before_any_fresh_alignment" and
            design["no_replacements_or_retuning"] is True, "Unfrozen replay design")
    exact(design["limits"], LIMITS, "Frozen resource limits")
    exact({p: sha(ROOT / p) for p in CODE_FILES}, design["code_sha256"], "Frozen auditor code")
    exact(environment_fingerprint(), design["environment"], "Pinned runtime")
    deadline = Budget(600, session)
    deadline.end = min(deadline.end, started + LIMITS["overall_seconds"])
    rebuilt, _ = inventory(session, deadline)
    exact(rebuilt, design["data"], "Frozen cases, masks and source completion")
    cases = design["data"]["cases"]
    identity = {"mode": VERSION, "design_sha256": sha(design_path), "code_sha256": design["code_sha256"],
                "source_snapshots": rebuilt["source_snapshots"], "environment_sha256": digest(canonical(design["environment"])),
                "case_ids": [c["case_id"] for c in cases], "limits": LIMITS, "information_scope": design["information_scope"]}
    job = session / "jobs" / (safe_name(args.name) + "-" + digest(canonical(identity))[:12])
    require(not job.exists(), "A replay attempt already exists; retries/replacements are not permitted")
    store = UnitStore(job, identity)
    store.event(event="guarded_replay_started", guard=guard, case_ids=identity["case_ids"])
    records, stopped = [], None
    for case in cases:
        if stopped is not None:
            record = {"case_id": case["case_id"], "status": "unattempted", "reason": stopped,
                      "case_manifest_sha256": digest(canonical(case))}
        else:
            try:
                deadline.check()
                require(time.monotonic() - started < 600, "Overall replay budget exhausted")
            except (ValueError, TimeoutError, OSError) as exc:
                record = {"case_id": case["case_id"], "status": "unattempted", "reason": "prelaunch_" + type(exc).__name__,
                          "case_manifest_sha256": digest(canonical(case))}
            else:
                try:
                    record = execute_case(case, design_path, store, guard, deadline)
                except (ValueError, TimeoutError, OSError, KeyError) as exc:
                    record = {"case_id": case["case_id"], "status": "failed", "reason": "case_" + type(exc).__name__,
                              "case_manifest_sha256": digest(canonical(case))}
            # No outcome-dependent backfill or retuning; fail closed on the first
            # discrepancy/resource failure and account for every remaining case.
            if record["status"] != "passed":
                stopped = "stopped_after_case_failure_or_budget"
        store.publish(case["case_id"], record)
        records.append(record)
        print(json.dumps({"case_id": case["case_id"], "status": record["status"], "reason": record.get("reason")}), flush=True)
    summary = accounting(cases, records, rebuilt["selection"])
    summary.update(elapsed_seconds=time.monotonic() - started, parent_peak_rss_bytes=peak_bytes(), guard=guard,
                   source_validation="All source units and frozen masks revalidated before any alignment",
                   memory_enforcement="100 ms owned-process RSS polling plus conservative sum of parent/worker peak RSS; no OS hard RSS reservation",
                   case_units=[{"case_id": r["case_id"], "status": r["status"]} for r in records])
    store.publish("summary", summary)
    store.complete(identity["case_ids"] + ["summary"], {"audit_status": summary["status"],
                              "accounting_complete": True, "representative_only": True})
    completion = next(store.root.glob("completion-*.json"))
    print(json.dumps({"job": str(store.root), "completion": str(completion), "completion_sha256": sha(completion), "status": summary["status"],
                      "passed_cases": summary["passed_cases"], "expected_cases": len(cases)}), flush=True)
    return 0 if summary["status"] == "passed" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("freeze", "run", "worker"))
    parser.add_argument("--session", type=Path)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--name", default="fragment-replay-audit-v1")
    parser.add_argument("--execute-fresh-replay", action="store_true")
    parser.add_argument("--expected-code-sha256")
    parser.add_argument("--expected-design-sha256")
    args = parser.parse_args()
    if args.mode == "worker":
        require(args.request is not None, "Missing private worker request")
        try:
            return worker(args.request)
        except Exception as exc:
            # Persist only the exception class: upstream exceptions need not be
            # safe to quote when input sequences were transiently in memory.
            request = read(args.request)
            atomic_json(Path(request["output_path"]), {"case_id": request["case_id"], "worker_error_type": type(exc).__name__})
            return 1
    require(args.session is not None, "Missing session directory")
    return freeze(args) if args.mode == "freeze" else run(args)


if __name__ == "__main__":
    sys.exit(main())
