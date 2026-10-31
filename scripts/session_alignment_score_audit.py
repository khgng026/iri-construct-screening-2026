#!/usr/bin/env python3
"""Independently audit full saved raw-alignment scores; run only in a lead slot.

All 3,536 constructed targets are checked against an attainable matrix-diagonal
upper bound. Eight deterministically selected real pairs are checked with the
existing independent Python affine DP, never PairwiseAligner scoring. Selection
and source evidence are saved before DP; budget-limited work cannot pass.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import datetime as dt
from fractions import Fraction
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import time
import uuid


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DP_SOURCE_SHA256 = "07e7dbc6c82be152b469d9e2d07c3c4132420439b56c4e3c980254a91b33f4d3"
EXPECTED_QUERY_COUNT = 1146
EXPECTED_TARGET_COUNT = 3536
SAMPLE_SIZE = 8
RECIPE = {
    "version": "independent-real-alignment-score-audit-v1",
    "frozen_before_any_dp_execution": True,
    "sample_size": SAMPLE_SIZE,
    "sampling_scope": "Complete raw105-reference job; all1146 queries. Audit-only stratification, not a performance estimate.",
    "off_target": "Exclude every reference whose sequence hash equals any component of this query, including duplicate accessions.",
    "unique_pair": "Deduplicate by (joined query sequence SHA256, reference sequence SHA256), not accession alone.",
    "selection": "In listed category order choose the first unused unique pair. Retain the best8 distinct sequence pairs per category while streaming; no fallback or resampling after scores.",
    "categories": [
        {"name": "shortest_target_query", "sort": ["query_length ascending", "reference_length ascending", "query_id", "accession"]},
        {"name": "longest_target_query", "sort": ["query_length descending", "reference_length descending", "query_id", "accession"]},
        {"name": "shortest_target_reference", "sort": ["reference_length ascending", "query_length descending", "query_id", "accession"]},
        {"name": "longest_target_reference", "sort": ["reference_length descending", "query_length descending", "query_id", "accession"]},
        {"name": "most_unequal_target", "sort": ["exact query/reference length ratio descending", "query_length descending", "query_id", "accession"]},
        {"name": "highest_score_off_target", "sort": ["saved raw score descending", "query_length descending", "query_id", "accession"]},
        {"name": "highest_score_gapped_off_target", "sort": ["saved raw score descending", "gap_columns descending", "query_id", "accession"]},
        {"name": "most_gapped_off_target", "sort": ["gap_columns descending", "gap_runs descending", "saved raw score descending", "query_id", "accession"]}
    ],
    "gapped_eligibility": "saved gap_columns >0",
    "dp_order": "First the upper median by cell count among frozen8 (tie category/query/accession); then remaining category order.",
    "timing_rule": "Forecast remaining DP cost from median-pair seconds/cell, multiplied by2 plus5seconds. Stop if forecast exceeds remaining DP budget; no smaller replacement cases.",
    "numerics": "BLOSUM62 matrix lookup only; local affine gap open10 extend0.5, first gap open; exact score equality because inputs and penalties are half-integer representable.",
    "default_dp_budget_seconds": 120,
    "default_per_pair_timeout_seconds": 30,
    "bound_proof": "For each reference symbol b, require positive B[b,b] and B[a,b]<=B[b,b] for all canonical a. A local alignment uses each reference position at most once and gaps incur negative costs. Therefore sum_b B[b,b] bounds every local alignment. An exact full reference at the constructed true span attains that bound.",
    "scope_limit": "Supplementary numerical verification, not a new biological result or random audit of every off-target optimum."
}


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def canonical_hash(value):
    return digest_bytes(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode())


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    with temporary.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def utcnow():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def category_key(category, row):
    q, r, score = row["query_length"], row["reference_length"], row["saved_score"]
    tail = (row["query_id"], row["accession"])
    if category == "shortest_target_query": return (q, r, *tail)
    if category == "longest_target_query": return (-q, -r, *tail)
    if category == "shortest_target_reference": return (r, -q, *tail)
    if category == "longest_target_reference": return (-r, -q, *tail)
    if category == "most_unequal_target": return (-Fraction(q, r), -q, *tail)
    if category == "highest_score_off_target": return (-score, -q, *tail)
    if category == "highest_score_gapped_off_target": return (-score, -row["gap_columns"], *tail)
    if category == "most_gapped_off_target": return (-row["gap_columns"], -row["gap_runs"], -score, *tail)
    raise ValueError("Unknown frozen selection category")


def retain_candidate(pool, category, candidate):
    pair = (candidate["query_sequence_sha256"], candidate["reference_sequence_sha256"])
    old = pool.get(pair)
    if old is None or category_key(category, candidate) < category_key(category, old):
        pool[pair] = candidate
    if len(pool) > SAMPLE_SIZE:
        worst = max(pool, key=lambda key: category_key(category, pool[key]))
        del pool[worst]


def inspect_target(query, reference, span, row, selfscore, source):
    n, start, end = len(reference), span["start0"], span["end0"]
    checks = {
        "constructed_span_is_reference": query["sequence"][start:end] == reference and end - start == n,
        "raw_score_equals_selfscore": row["score"] == selfscore,
        "saved_query_span_is_true_span": row["query_span"] == [start, end],
        "saved_reference_span_is_full": row["reference_span"] == [0, n],
        "lengths": row["query_length"] == len(query["sequence"]) and row["reference_length"] == n,
        "identities": row["identities"] == n and row["identity_fraction"] == 1.0,
        "aligned_lengths": row["aligned_length"] == row["aligned_residue_pairs"] == n,
        "no_gaps": row["gap_columns"] == row["gap_runs"] == 0,
        "reference_coverage": row["reference_coverage"] == row["reference_paired_coverage"] == 1.0,
        "query_coverage": abs(row["query_coverage"] - n / len(query["sequence"])) <= 1e-12 and
                          abs(row["query_paired_coverage"] - n / len(query["sequence"])) <= 1e-12,
    }
    coordinates = row.get("coordinates")
    try:
        require(isinstance(coordinates, list) and len(coordinates) == 2, "Bad coordinate rows")
        qcoords, rcoords = coordinates
        require(len(qcoords) == len(rcoords) >= 2, "Bad coordinate width")
        require(all(type(v) is int for values in coordinates for v in values), "Noninteger coordinates")
        require(qcoords[0] == start and qcoords[-1] == end and rcoords[0] == 0 and rcoords[-1] == n,
                "Coordinate endpoints differ from exact true span")
        require(all(qp - start == rp for qp, rp in zip(qcoords, rcoords)), "Coordinates leave exact diagonal")
        require(all(a < b for a, b in zip(qcoords, qcoords[1:])), "Nonpositive coordinate step")
        checks["exact_coordinate_path"] = True
    except (TypeError, ValueError) as error:
        checks["exact_coordinate_path"] = False
        checks["coordinate_error"] = str(error)
    return {"query_id": query["query_id"], "suite": query["suite"], "accession": span["accession"],
            "source": source, "expected_span": [start, end], "reference_selfscore": selfscore,
            "saved_raw_score": row["score"], "score_difference": row["score"] - selfscore,
            "saved_alignment_row": row, "checks": checks,
            "passed": all(value is True for name, value in checks.items() if name != "coordinate_error")}


@contextmanager
def timed_dp(seconds):
    """Bound this pure-Python DP without an additional concurrent process."""
    require(seconds > 0 and signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0), "Invalid/nested DP timer")
    previous = signal.getsignal(signal.SIGALRM)
    def timeout(signum, frame):
        raise TimeoutError("Independent DP exceeded its recorded wall-time cap")
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def run(args):
    session, raw_job = args.session.resolve(), args.raw_job.resolve()
    session.relative_to((ROOT / "results").resolve())
    raw_job.relative_to(session / "jobs")
    require(1 <= args.dp_budget_seconds <= 300 and 1 <= args.per_pair_seconds <= 60, "Requested audit time cap is outside allowed bounds")
    output = session / "audit" / ("independent_alignment_score_" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:10])
    output.mkdir(parents=True, exist_ok=False)
    (output / "source_records").mkdir()
    report = {"schema_version": 1, "status": "running", "started_utc": utcnow(), "output": str(output),
              "recipe": RECIPE, "recipe_sha256": canonical_hash(RECIPE), "dp_results": [],
              "target_observations_checked": 0, "target_failures": 0,
              "target_scores_equal_selfscore": 0, "target_exact_true_span_count": 0,
              "target_counts_by_suite": {},
              "dp_budget_seconds": args.dp_budget_seconds, "per_pair_seconds": args.per_pair_seconds,
              "script_sha256": sha(Path(__file__)), "scope": RECIPE["scope_limit"]}
    atomic_json(output / "recipe.json", RECIPE)
    atomic_json(output / "report.json", report)
    exit_code = 1
    try:
        deadline = dt.datetime.fromisoformat(read(session / "launch.json")["hard_stop_utc"])
        require(deadline.tzinfo is not None and dt.datetime.now(dt.timezone.utc) < deadline, "Session hard stop reached")
        from iri_constructs.session_benchmark import load_suite
        from session_independent_audit import verify_unit, validate_ranking
        from iri_constructs.session_alignment import independent_affine_score
        from Bio.Align import substitution_matrices
        import Bio.Align as bio_align
        dp_path = Path(inspect.getsourcefile(independent_affine_score)).resolve()
        require(sha(dp_path) == DP_SOURCE_SHA256, "Independent DP implementation changed from frozen full-raw job source")
        report["dp_source_sha256"] = sha(dp_path)
        report["dp_function_sha256"] = digest_bytes(inspect.getsource(independent_affine_score).encode())
        report["validation_code_sha256"] = {str(path): sha(path) for path in
            (ROOT / "scripts/session_independent_audit.py", ROOT / "src/iri_constructs/session_benchmark.py")}
        suite = load_suite(ROOT, session / "inputs")
        queries = {q["query_id"]: q for q in suite["queries"]}
        references = suite["sequences"]
        require(len(queries) == EXPECTED_QUERY_COUNT and len(references) == 105 and
                sum(len(q["spans"]) for q in queries.values()) == EXPECTED_TARGET_COUNT, "Full frozen suite denominator changed")
        identity_path = raw_job / "identity.json"
        identity = read(identity_path)
        identity_hash = canonical_hash(identity)
        require(identity["mode"] == "alignment" and set(identity["query_ids"]) == set(queries)
                and len(identity["query_ids"]) == len(queries), "Expected complete full raw-alignment job")
        require(identity["catalogue_accessions"] == sorted(references), "Source catalogue changed")
        reference_hashes = {a: digest_bytes(s.encode()) for a, s in references.items()}
        require(identity["catalogue_sequence_sha256"] == reference_hashes, "Reference sequence hashes differ")
        require(identity["challenge_sha256"] == sha(session / "inputs/challenge_manifest.json"), "Source challenge differs")
        require(identity["code_hashes"]["src/iri_constructs/session_alignment.py"] == DP_SOURCE_SHA256,
                "Saved alignment implementation differs from frozen source")
        method = identity["method"]
        require(method["mode"] == "local" and method["matrix"] == "BLOSUM62" and
                method["gap_open"] == 10 and method["gap_extend"] == 0.5 and
                method["gap_cost"] == "gap_open + (gap_length - 1) * gap_extend" and
                method["normalization"] == "none; raw scores favour longer high-scoring matches",
                "Saved scoring method differs from audit recurrence")
        completions = [(read(path)["completed_utc"], str(path), path) for path in raw_job.glob("completion-*.json")]
        require(completions, "Source has no completion descriptor")
        completion_path = max(completions)[2]
        completion = read(completion_path)
        expected_units = {"query-" + q for q in queries} | {"summary"}
        require(completion["identity_hash"] == identity_hash and set(completion["units"]) == expected_units
                and completion["expected_count"] == len(expected_units), "Source completion coverage differs")
        matrix = substitution_matrices.load("BLOSUM62")
        require(digest_bytes(matrix.tobytes()) == method["matrix_sha256"], "Matrix bytes differ from scored job")
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        require(all(float(matrix[b, b]) > 0 and all(float(matrix[a, b]) <= float(matrix[b, b]) for a in amino_acids)
                    for b in amino_acids), "Matrix column-diagonal upper-bound premise fails")
        require(all(float(matrix[a, b]).is_integer() for a in amino_acids for b in amino_acids), "Matrix no longer integer valued")
        selfscores = {a: math.fsum(float(matrix[x, x]) for x in sequence) for a, sequence in references.items()}
        report["upper_bound"] = {"columnwise_diagonal_maximum_verified": True, "positive_diagonals_verified": True,
                                 "negative_gap_costs": True, "matrix_sha256": method["matrix_sha256"],
                                 "canonical_alphabet_checked": amino_acids,
                                 "proof": RECIPE["bound_proof"], "reference_selfscores": selfscores}
        source_manifest = {"job": str(raw_job), "identity_sha256": sha(identity_path), "identity_hash": identity_hash,
                           "completion_path": str(completion_path), "completion_sha256": sha(completion_path),
                           "challenge_sha256": identity["challenge_sha256"], "data_source_files_sha256": suite["manifest"]["source_files_sha256"],
                           "catalogue_sequence_sha256": reference_hashes, "units": {}}
        shutil.copyfile(identity_path, output / "source-identity.json")
        shutil.copyfile(completion_path, output / "source-completion.json")
        category_names = [category["name"] for category in RECIPE["categories"]]
        pools = {name: {} for name in category_names}
        total_bounds_checked = 0
        target_path = output / "all-target-checks.jsonl"
        with target_path.open("x") as target_stream:
            for qid in sorted(queries):
                require(dt.datetime.now(dt.timezone.utc) < deadline, "Session hard stop reached during source verification")
                query = queries[qid]
                key = "query-" + qid
                path = verify_unit(raw_job, key, completion["units"][key], identity_hash)
                descriptor = read(path / "complete.json")
                source = {"unit": str(path), "descriptor_sha256": completion["units"][key]["descriptor_sha256"],
                          "record_sha256": descriptor["files"]["record.json"]["sha256"],
                          "record_bytes": descriptor["files"]["record.json"]["bytes"]}
                source_manifest["units"][key] = source
                record = read(path / "record.json")
                require(record["query"] == {k: v for k, v in query.items() if k != "sequence"}, "Saved query differs from reviewed construction")
                require(set(record["rankings"]) == {"local_alignment"}, "Unexpected source score arm")
                rows = record["rankings"]["local_alignment"]
                validate_ranking(rows, sorted(references))
                by_accession = {row["accession"]: row for row in rows}
                target_hashes = {reference_hashes[a] for a in query["components"]}
                for span in query["spans"]:
                    accession = span["accession"]
                    check = inspect_target(query, references[accession], span, by_accession[accession], selfscores[accession], source)
                    target_stream.write(json.dumps(check, sort_keys=True, allow_nan=False) + "\n")
                    report["target_observations_checked"] += 1
                    report["target_failures"] += int(not check["passed"])
                    report["target_scores_equal_selfscore"] += int(check["checks"]["raw_score_equals_selfscore"])
                    report["target_exact_true_span_count"] += int(all(check["checks"][name] for name in
                        ("constructed_span_is_reference", "saved_query_span_is_true_span", "saved_reference_span_is_full",
                         "identities", "aligned_lengths", "no_gaps", "exact_coordinate_path")))
                    report["target_counts_by_suite"][query["suite"]] = report["target_counts_by_suite"].get(query["suite"], 0) + 1
                for row in rows:
                    accession = row["accession"]
                    require(row["query_length"] == len(query["sequence"]) and
                            row["reference_length"] == len(references[accession]), "Saved alignment sequence lengths differ")
                    require(type(row["gap_columns"]) is int and type(row["gap_runs"]) is int and
                            0 <= row["gap_runs"] <= row["gap_columns"], "Invalid saved gap counts")
                    require(0 <= row["score"] <= selfscores[accession] and float(2 * row["score"]).is_integer(),
                            "Saved raw score violates nonnegative half-integer/diagonal bound")
                    total_bounds_checked += 1
                    candidate = {"query_id": qid, "accession": accession, "query_length": len(query["sequence"]),
                                 "reference_length": len(references[accession]), "cells": len(query["sequence"]) * len(references[accession]),
                                 "query_sequence_sha256": query["joined_sequence_sha256"], "reference_sequence_sha256": reference_hashes[accession],
                                 "saved_score": row["score"], "gap_columns": row["gap_columns"], "gap_runs": row["gap_runs"],
                                 "saved_rank": row["rank"], "source": source, "saved_alignment_row": row}
                    if accession in query["components"]:
                        eligible = category_names[:5]
                    elif reference_hashes[accession] not in target_hashes:
                        eligible = category_names[5:] if row["gap_columns"] > 0 else category_names[5:6]
                    else:
                        eligible = []
                    for category in eligible:
                        retain_candidate(pools[category], category, candidate)
                if len(source_manifest["units"]) % 200 == 0:
                    target_stream.flush()
                    atomic_json(output / "sources.json", source_manifest)
                    atomic_json(output / "report.json", report)
            target_stream.flush()
            os.fsync(target_stream.fileno())
        summary_path = verify_unit(raw_job, "summary", completion["units"]["summary"], identity_hash)
        summary_descriptor = read(summary_path / "complete.json")
        source_manifest["units"]["summary"] = {"unit": str(summary_path),
            "descriptor_sha256": completion["units"]["summary"]["descriptor_sha256"],
            "record_sha256": summary_descriptor["files"]["record.json"]["sha256"]}
        report["all_target_checks_sha256"] = sha(target_path)
        report["all_raw_rows_bound_checked"] = total_bounds_checked
        require(report["target_observations_checked"] == EXPECTED_TARGET_COUNT, "Target denominator mismatch")
        require(total_bounds_checked == EXPECTED_QUERY_COUNT * 105, "Raw-row bound-check denominator mismatch")
        selected, used = [], set()
        for category in category_names:
            candidates = sorted(pools[category].items(), key=lambda item: category_key(category, item[1]))
            choice = next(((key, row) for key, row in candidates if key not in used), None)
            require(choice is not None, "Frozen category has no unused eligible pair: " + category)
            key, candidate = choice
            used.add(key)
            selected.append({"category": category, **candidate})
        require(len(selected) == SAMPLE_SIZE, "Numerical sample size differs")
        for candidate in selected:
            source = candidate["source"]
            for filename, checksum in (("record.json", source["record_sha256"]), ("complete.json", source["descriptor_sha256"])):
                destination = output / "source_records" / (checksum + "-" + filename)
                if not destination.exists():
                    shutil.copyfile(Path(source["unit"]) / filename, destination)
                require(sha(destination) == checksum, "Copied selected source evidence differs")
        atomic_json(output / "sources.json", source_manifest)
        representative = sorted(selected, key=lambda row: (row["cells"], row["category"], row["query_id"], row["accession"]))[SAMPLE_SIZE // 2]
        order = [representative] + [row for row in selected if row is not representative]
        selection = {"selected_utc": utcnow(), "recipe_sha256": report["recipe_sha256"], "sources_sha256": sha(output / "sources.json"),
                     "selection_frozen_before_dp": True, "pairs": selected, "dp_order_categories": [r["category"] for r in order]}
        atomic_json(output / "selection.json", selection)
        report["selection_sha256"] = sha(output / "selection.json")
        report["dp_results"] = [{"category": row["category"], "query_id": row["query_id"], "accession": row["accession"],
                                  "cells": row["cells"], "saved_score": row["saved_score"], "status": "not_run"} for row in order]
        atomic_json(output / "report.json", report)
        original_aligner = bio_align.PairwiseAligner
        def forbidden_aligner(*unused_args, **unused_kwargs):
            raise RuntimeError("PairwiseAligner is forbidden in this independent score audit")
        bio_align.PairwiseAligner = forbidden_aligner
        report["production_scoring_guard"] = "PairwiseAligner constructor replaced with a raising stub during every independent DP call."
        dp_started = time.monotonic()
        try:
            for position, row in enumerate(order):
                remaining = min(args.dp_budget_seconds - (time.monotonic() - dp_started),
                                (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds())
                if remaining <= 0:
                    report["budget_stop_reason"] = "No DP wall-time budget remains"
                    break
                if position == 1:
                    probe = report["dp_results"][0]
                    require(probe["status"] == "passed", "Representative DP did not pass; remaining DP not admitted")
                    forecast = 2 * probe["elapsed_seconds"] * sum(r["cells"] for r in order[1:]) / probe["cells"] + 5
                    report["representative_forecast_remaining_seconds"] = forecast
                    if forecast > remaining:
                        report["budget_stop_reason"] = "Frozen representative timing forecast exceeds remaining DP budget"
                        break
                result = report["dp_results"][position]
                result["status"] = "running"
                result["timeout_seconds"] = min(args.per_pair_seconds, remaining)
                atomic_json(output / "report.json", report)
                start = time.monotonic()
                try:
                    with timed_dp(result["timeout_seconds"]):
                        value = independent_affine_score(queries[row["query_id"]]["sequence"], references[row["accession"]])
                    result["independent_score"] = value
                    result["difference"] = value - row["saved_score"]
                    result["status"] = "passed" if math.isfinite(value) and value == row["saved_score"] else "failed"
                except TimeoutError as error:
                    result["status"] = "timeout"
                    result["error"] = str(error)
                finally:
                    result["elapsed_seconds"] = time.monotonic() - start
                    atomic_json(output / "report.json", report)
                print(json.dumps({"stage": "independent_dp", "pair_completed": position + 1,
                                  "pairs_frozen": SAMPLE_SIZE, "status": result["status"],
                                  "seconds": result["elapsed_seconds"]}), flush=True)
                if result["status"] == "timeout":
                    report["budget_stop_reason"] = "A frozen pair exceeded its recorded cap; no replacement selected"
                    break
        finally:
            bio_align.PairwiseAligner = original_aligner
            report["dp_elapsed_seconds"] = time.monotonic() - dp_started
        require(sha(identity_path) == source_manifest["identity_sha256"] and sha(completion_path) == source_manifest["completion_sha256"],
                "Source identity/completion changed during audit")
        complete_dp = all(result["status"] == "passed" for result in report["dp_results"])
        numerical_failure = report["target_failures"] > 0 or any(result["status"] == "failed" for result in report["dp_results"])
        report["status"] = "failed" if numerical_failure else ("passed" if complete_dp else
                           ("budget_limited" if "budget_stop_reason" in report else "failed"))
        exit_code = 0 if report["status"] == "passed" else 1
    except BaseException as error:
        report["status"] = "interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        exit_code = int(error.code) if isinstance(error, SystemExit) and isinstance(error.code, int) else 1
    finally:
        report["finished_utc"] = utcnow()
        atomic_json(output / "report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output / "report.json")}))
    return exit_code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", type=Path, required=True)
    parser.add_argument("--raw-job", type=Path, required=True)
    parser.add_argument("--dp-budget-seconds", type=float, default=120)
    parser.add_argument("--per-pair-seconds", type=float, default=30)
    args = parser.parse_args()
    def interrupted(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
