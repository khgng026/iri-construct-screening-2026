#!/usr/bin/env python3
"""Exploratory reference-self-score normalization of immutable raw alignments.

Primary raw scores remain untouched. Prediction uses no target boundaries/labels.
Run only after the raw worker completes and the lead allocates an analysis slot.
"""
from __future__ import annotations
import argparse
from collections import Counter
import datetime as dt
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash, utcnow


def read(path):
    return json.loads(Path(path).read_text())


def require(value, message):
    if not value: raise ValueError(message)


def latest_completion(job):
    candidates = [(read(p)["completed_utc"], p) for p in job.glob("completion-*.json")]
    require(candidates, "Source alignment job has not completed")
    return max(candidates)[1]


def reference_selfscores(sequences, *, check_examples=True):
    from Bio.Align import substitution_matrices
    from iri_constructs.session_alignment import make_aligner
    matrix = substitution_matrices.load("BLOSUM62")
    aa = "ACDEFGHIKLMNPQRSTVWY"
    require(all(matrix[a, a] > 0 and matrix[a, b] <= matrix[a, a] for a in aa for b in aa),
            "Canonical diagonal no longer bounds substitution scores")
    scores = {a: math.fsum(float(matrix[x, x]) for x in s) for a, s in sequences.items()}
    checks = []
    if check_examples:
        by_length = sorted(sequences, key=lambda a: (len(sequences[a]), a))
        chosen = list(dict.fromkeys([by_length[0], by_length[len(by_length) // 2], by_length[-1]]))
        aligner = make_aligner()
        for accession in chosen:
            actual = float(aligner.score(sequences[accession], sequences[accession]))
            require(actual == scores[accession], "Direct self-alignment disagrees with diagonal sum")
            checks.append({"accession": accession, "length": len(sequences[accession]),
                           "diagonal_sum": scores[accession], "independent_aligner_selfscore": actual})
    return scores, checks, digest(matrix.tobytes())


def normalized_ranking(raw_ranking, selfscores):
    from iri_constructs.session_benchmark import rank_scores, filter_ranking
    raw = {r["accession"]: r for r in raw_ranking}
    require(len(raw) == len(raw_ranking) and set(raw) == set(selfscores), "Raw catalogue coverage mismatch")
    accessions = sorted(raw)
    values = [float(raw[a]["score"]) / selfscores[a] for a in accessions]
    require(all(math.isfinite(v) and -1e-12 <= v <= 1 + 1e-12 for v in values),
            "Normalized score outside diagonal upper bound; investigate without clipping")
    rows = rank_scores(values, accessions, "normalized_alignment_score")
    for row in rows:
        a = row["accession"]
        row.update(raw_alignment_score=float(raw[a]["score"]), reference_selfscore=selfscores[a])
    # Recompute exact and near-tie metadata after normalization, never carry raw ties.
    return filter_ranking(rows, accessions)


def run(args):
    from iri_constructs.session_benchmark import (load_suite, catalogue_masks, filter_ranking,
        evaluate_query, summarize, common_reference_removal_mask)
    from session_independent_audit import (verify_unit, validate_ranking, validate_evaluation,
                                           validate_summary, canonical_hash)
    started = time.monotonic(); session, rawjob = args.session.resolve(), args.raw_job.resolve()
    launch = read(session / "launch.json")
    deadline = dt.datetime.fromisoformat(launch["hard_stop_utc"])
    def within_deadline():
        require(dt.datetime.now(dt.timezone.utc) < deadline, "Session hard stop reached")
    within_deadline()
    design = read(args.design)
    require(design["version"] == "alignment-reference-selfscore-normalization-v1"
            and design["status"] == "frozen_before_normalized_scoring", "Unrecognized/unfrozen normalization design")
    require(design["secondary_score"] == "S(query, reference) / S(reference, reference)", "Changed normalization formula")
    suite = load_suite(ROOT, session / "inputs")
    require(design["source_challenge_sha256"] == file_hash(session / "inputs/challenge_manifest.json"), "Normalization challenge changed")
    raw_identity = read(rawjob / "identity.json"); raw_hash = canonical_hash(raw_identity)
    require(raw_identity["mode"] == "alignment", "Source must be the primary raw alignment arm")
    require(raw_identity["challenge_sha256"] == design["source_challenge_sha256"], "Source challenge mismatch")
    accessions = sorted(suite["sequences"])
    require(raw_identity["catalogue_accessions"] == accessions, "Source catalogue changed")
    require(raw_identity["catalogue_sequence_sha256"] == {a: digest(suite["sequences"][a].encode()) for a in accessions},
            "Raw/source reference sequence hashes disagree")
    completion_path = latest_completion(rawjob); completion = read(completion_path)
    require(completion["identity_hash"] == raw_hash, "Raw completion identity mismatch")
    query_ids = raw_identity["query_ids"]
    require(len(query_ids) == len(set(query_ids)), "Duplicate raw query ID")
    if not args.allow_partial:
        require(set(query_ids) == {q["query_id"] for q in suite["queries"]}, "Full normalization requires completed full raw suite")
    expected_raw = {"query-" + q for q in query_ids} | {"summary"}
    require(set(completion["units"]) == expected_raw and completion["expected_count"] == len(expected_raw), "Raw completion has missing/extra units")
    paths = {key: verify_unit(rawjob, key, index, raw_hash) for key, index in completion["units"].items()}
    selfscores, selfchecks, matrix_hash = reference_selfscores(suite["sequences"])
    require(matrix_hash == design["matrix_sha256"], "BLOSUM62 matrix changed")
    masks = catalogue_masks(suite)
    identity = {"version": 1, "mode": "alignment_normalization", "design": design,
                "design_sha256": file_hash(args.design), "source_raw_job": str(rawjob),
                "source_identity_hash": raw_hash, "source_completion_sha256": file_hash(completion_path),
                "challenge_sha256": design["source_challenge_sha256"], "query_ids": query_ids,
                "catalogue_accessions": accessions, "catalogue_sequence_sha256": raw_identity["catalogue_sequence_sha256"],
                "reference_selfscores": selfscores, "matrix_sha256": matrix_hash,
                "catalogue_masks": masks, "partial_source": args.allow_partial,
                "code_hashes": {str(p.relative_to(ROOT)): file_hash(p) for p in
                    [Path(__file__), ROOT / "scripts/session_independent_audit.py", ROOT / "src/iri_constructs/session_alignment.py",
                     ROOT / "src/iri_constructs/session_benchmark.py", ROOT / "src/iri_constructs/session_store.py"]},
                "interpretation": "Outcome-informed exploratory ranking diagnostic; no calibrated significance or biological function"}
    jobid = args.name + "-" + digest(canonical(identity))[:12]
    store = UnitStore(session / "jobs" / jobid, identity, session / "job_ledger.jsonl")
    store.event(event="job_started", state="running", job_id=jobid, pid=os.getpid(), query_count=len(query_ids))
    expected = []
    try:
        key = "reference-selfscores"; expected.append(key)
        if store.find(key) is None:
            store.publish(key, {"scores": selfscores, "self_alignment_checks": selfchecks,
                               "matrix_sha256": matrix_hash, "no_target_labels_used": True})
        queries = {q["query_id"]: q for q in suite["queries"]}
        for i, qid in enumerate(query_ids):
            within_deadline(); key = "query-" + qid; expected.append(key)
            if store.find(key): continue
            q = queries[qid]; source_record = read(paths[key] / "record.json")
            require(source_record["query"] == {k: v for k, v in q.items() if k != "sequence"}, "Raw/source query metadata disagree")
            raw = source_record["rankings"]["local_alignment"]
            validate_ranking(raw, accessions)
            ranking = normalized_ranking(raw, selfscores)
            by_accession = {r["accession"]: r for r in raw}
            evaluations, rankings, independent = {}, {}, {}
            for mask_id, included in masks["masks"].items():
                masked = filter_ranking(ranking, included)
                evaluation = evaluate_query(masked, q, suite["equivalence"])
                validate_ranking(masked, included)
                eq = {a: [b for b in suite["equivalence"][a] if b in included] for a in included}
                independent[mask_id] = validate_evaluation(evaluation, q, masked, eq)
                rankings[mask_id] = masked; evaluations[mask_id] = evaluation
            targets = []
            for span in q["spans"]:
                row = by_accession[span["accession"]]
                exact = (row["identity_fraction"] == 1 and row["reference_coverage"] == 1
                         and row["query_span"] == [span["start0"], span["end0"]]
                         and row["reference_span"] == [0, span["end0"] - span["start0"]])
                targets.append({"accession": span["accession"], "fraction": span["fraction"], "raw_rank": row["rank"],
                                "raw_score": row["score"], "full_exact_match_at_true_span": exact,
                                "reference_coverage": row["reference_coverage"], "identity_fraction": row["identity_fraction"],
                                "query_span": row["query_span"], "expected_span": [span["start0"], span["end0"]]})
            removal = common_reference_removal_mask(q, suite["equivalence"], accessions)
            store.publish(key, {"query": {k: v for k, v in q.items() if k != "sequence"},
                "rankings": rankings, "evaluations": evaluations, "raw_target_diagnostics": targets,
                "source_raw_unit": str(paths[key]), "source_raw_unit_descriptor_sha256": file_hash(paths[key] / "complete.json"),
                "reference_removal_diagnostic": {"removed_accessions": sorted(set(accessions) - set(removal)),
                    "remaining_count": len(removal), "ranking": filter_ranking(ranking, removal),
                    "exact_accuracy_defined": False, "reason": "All exact-equivalent target references removed together; no independently defined surviving labels"}})
            if (i + 1) % 100 == 0 or i + 1 == len(query_ids):
                print(json.dumps({"job": jobid, "queries_completed": i + 1, "queries_expected": len(query_ids)}), flush=True)
        by_mask = {mask_id: [] for mask_id in masks["masks"]}; diagnostics = []
        for qid in query_ids:
            record = store.read("query-" + qid)
            for mask_id, evaluation in record["evaluations"].items():
                by_mask[mask_id].append(evaluation)
            diagnostics.extend({"query_id": qid, "suite": record["query"]["suite"], **t} for t in record["raw_target_diagnostics"])
        summaries = {m: {"all": summarize(evs), "by_suite": {s: summarize([e for e in evs if e["suite"] == s])
                       for s in sorted({e["suite"] for e in evs})}} for m, evs in by_mask.items()}
        # Independently recompute every summary count from saved component ranks.
        for mask_id, evs in by_mask.items():
            validate_summary(summaries[mask_id]["all"], evs)
            for suite_name, summary in summaries[mask_id]["by_suite"].items():
                validate_summary(summary, [e for e in evs if e["suite"] == suite_name])
        report = {"summary": summaries, "query_count": len(query_ids), "catalogue_masks": masks,
                  "raw_component_observations": len(diagnostics),
                  "raw_exact_true_span_count": sum(t["full_exact_match_at_true_span"] for t in diagnostics),
                  "raw_top5_misses_with_exact_true_span": [t for t in diagnostics if t["raw_rank"] > 5 and t["full_exact_match_at_true_span"]],
                  "raw_top5_miss_count": sum(t["raw_rank"] > 5 for t in diagnostics),
                  "all_denominators_independently_checked": True}
        expected.append("summary")
        old = store.read("summary")
        if old is None: store.publish("summary", report)
        else: require(old == report, "Existing normalized summary differs from recomputed saved results")
        complete = store.complete(expected, {"query_count": len(query_ids), "elapsed_seconds": time.monotonic() - started,
                                             "exploratory": True, "primary_raw_unchanged": True})
        print(json.dumps({"status": "complete", "job": str(store.root), "units": complete["expected_count"],
                          "elapsed_seconds": time.monotonic() - started}), flush=True)
    except BaseException as exc:
        store.event(event="job_failed", state="failed", job_id=jobid, exception=f"{type(exc).__name__}: {exc}")
        raise


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--session", type=Path, required=True); p.add_argument("--raw-job", type=Path, required=True)
    p.add_argument("--design", type=Path, required=True); p.add_argument("--name", default="alignment-normalized-exploratory")
    p.add_argument("--allow-partial", action="store_true", help="Explicitly label a completed trial subset instead of full suite")
    run(p.parse_args())
