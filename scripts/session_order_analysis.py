#!/usr/bin/env python3
"""Descriptive order sensitivity from complete immutable saved model outputs.

No encoder, alignment, production metric accumulator, parameter selection or
IID inference. Stream one model/group; keep exact observed mixture discrepancies
visible. A failed sanity check produces accounted-for outputs with failed status.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import itertools
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

import session_extension_audit as base

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.session_store import UnitStore, canonical, digest, safe_name

require, read, sha = base.require, base.read, base.sha
DESIGN_SHA = "86625d80fe1d49927e7571daf19d0515d947a0e281e200899141c14373527335"
ARMS = ("whole", "isolated_mixture", "contextual_component")
MODELS = ("esm2_8m", "esm2_150m", "esm2_650m", "esmc_300m", "esmc_600m")
SUITES = ("pairs", "triples", "quadruples")
QUANTILES = (0., .25, .5, .75, 1.)
VECTOR_TOLERANCE = 1e-6
SCORE_TOLERANCE = 1e-6
NEAR_TIE_GAP = 2e-6


def peak_mib():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (1024**2 if sys.platform == "darwin" else 1024)


class Deadline(base.Deadline):
    def check(self):
        super().check()
        require(peak_mib() <= 512, "Saved-order analysis exceeded 512 MiB")


def mean(values):
    values = list(values)
    return math.fsum(values) / len(values) if values else None


def distribution(values):
    values = sorted(values)
    require(all(math.isfinite(v) for v in values), "Nonfinite descriptive value")
    if not values:
        return {"count": 0, "mean": None, "quantiles": {str(q): None for q in QUANTILES}}
    result = {}
    for q in QUANTILES:
        position = (len(values) - 1) * q
        lo, hi = math.floor(position), math.ceil(position)
        result[str(q)] = values[lo] + (position - lo) * (values[hi] - values[lo])
    return {"count": len(values), "mean": mean(values), "quantiles": result}


def allowed_orders(components):
    values = tuple(sorted(components))
    require(len(values) == len(set(values)) and 1 <= len(values) <= 4, "Invalid component set")
    if len(values) <= 3:
        return set(itertools.permutations(values))
    rotations = [values[i:] + values[:i] for i in range(4)]
    return set(rotations + [order[::-1] for order in rotations])


def group_queries(queries):
    groups = defaultdict(list)
    require(len({q["query_id"] for q in queries}) == len(queries), "Duplicate query IDs")
    for query in queries:
        groups[query["group_id"]].append(query)
    for gid, rows in groups.items():
        rows.sort(key=lambda q: q["query_id"])
        first = rows[0]
        components = set(first["components"])
        lengths = {s["accession"]: s["end0"] - s["start0"] for s in first["spans"]}
        require(len(components) == first["component_count"] and set(lengths) == components, "Group components differ from spans")
        for query in rows:
            require(set(query["components"]) == components and len(query["components"]) == len(components) and
                    query["component_count"] == len(components) and query["length_aa"] == first["length_aa"] and
                    query["suite"] == first["suite"], "Group composition/length/suite changed across orders")
            offset = 0
            require([s["accession"] for s in query["spans"]] == query["components"], "Span order differs")
            for span in query["spans"]:
                a = span["accession"]
                require(span["start0"] == offset and span["end0"] == offset + lengths[a] and
                        span["fraction"] == lengths[a] / first["length_aa"], "Target fraction/geometry changes within group")
                offset = span["end0"]
            require(offset == query["length_aa"], "Group spans do not cover input")
        observed = [tuple(q["components"]) for q in rows]
        require(len(observed) == len(set(observed)) and set(observed) == allowed_orders(components), "Missing, duplicated or unplanned frozen order")
    return dict(sorted(groups.items()))


def coverage(groups):
    rows = [q for values in groups.values() for q in values]
    combined = {g: values for g, values in groups.items() if values[0]["component_count"] > 1}
    return {"queries": len(rows), "suite_queries": dict(sorted(Counter(q["suite"] for q in rows).items())),
            "combined_groups": dict(sorted(Counter(v[0]["suite"] for v in combined.values()).items())),
            "combined_group_target_cells": sum(v[0]["component_count"] for v in combined.values()),
            "combined_query_component_observations": sum(q["component_count"] for q in rows if q["component_count"] > 1)}


def require_full_panel(identity, expected_ids, accessions):
    ids = identity["query_ids"]
    require(identity["mode"] == "model" and len(ids) == len(set(ids)) == 1146 and set(ids) == set(expected_ids),
            "Every model must have the exact complete 1146-query panel; intersections are forbidden")
    require(identity["catalogue_accessions"] == sorted(accessions) and len(accessions) == 105, "Full catalogue differs")


def rank_map(rows, accessions):
    base.validate_ranking(rows, accessions)
    mapping = {row["accession"]: row for row in rows}
    require([r["accession"] for r in rows] == sorted(accessions, key=lambda a: (-base.score(mapping[a]), a)),
            "Source score/accession ranking differs")
    return mapping


def mixture_changes(anchor_rows, observed_rows, accessions):
    """Retain exact ranks; classify candidate-pair inversions without reranking."""
    anchor, observed = rank_map(anchor_rows, accessions), rank_map(observed_rows, accessions)
    maximum = max(abs(base.score(anchor[a]) - base.score(observed[a])) for a in accessions)
    changed = [{"accession": a, "anchor_rank": anchor[a]["rank"], "observed_rank": observed[a]["rank"],
                "anchor_score": base.score(anchor[a]), "observed_score": base.score(observed[a])}
               for a in accessions if anchor[a]["rank"] != observed[a]["rank"]]
    inversions = []
    for a, b in itertools.combinations(accessions, 2):
        if (anchor[a]["rank"] - anchor[b]["rank"]) * (observed[a]["rank"] - observed[b]["rank"]) < 0:
            left, right = abs(base.score(anchor[a]) - base.score(anchor[b])), abs(base.score(observed[a]) - base.score(observed[b]))
            inversions.append({"accessions": [a, b], "anchor_score_gap": left, "observed_score_gap": right,
                               "within_near_tie_gap": left <= NEAR_TIE_GAP and right <= NEAR_TIE_GAP})
    return {"maximum_absolute_score_difference": maximum, "scores_exactly_equal": maximum == 0,
            "changed_candidate_ranks": changed, "inverted_candidate_pairs": inversions,
            "all_rank_changes_within_near_ties": all(row["within_near_tie_gap"] for row in inversions),
            "scores_within_tolerance": maximum <= SCORE_TOLERANCE}


def position_means(observations, component_count):
    overall = {key: mean(row[key] for row in observations) for key in ("rank", "score", "hit_at_5")}
    result = []
    for position in range(1, component_count + 1):
        rows = [row for row in observations if row["position1"] == position]
        values = {key: mean(row[key] for row in rows) for key in overall}
        result.append({"position1": position, "observation_count": len(rows), "query_ids": [r["query_id"] for r in rows],
                       "means": values, "within_cell_centered_means": {key: values[key] - overall[key] if rows else None for key in overall}})
    return result


def describe_cell(group, accession, arm, observations):
    require(observations and len({r["query_id"] for r in observations}) == len(observations), "Missing/duplicate cell observations")
    first = group[0]
    require({r["query_id"] for r in observations} == {q["query_id"] for q in group}, "Cell lost an observed order")
    ranks, scores, hits = ([row[key] for row in observations] for key in ("rank", "score", "hit_at_5"))
    positions = position_means(observations, first["component_count"])
    span = next(s for s in first["spans"] if s["accession"] == accession)
    return {"group_id": first["group_id"], "accession": accession, "arm": arm, "suite": first["suite"],
            "component_count": first["component_count"], "component_fraction": span["fraction"],
            "target_length": span["end0"] - span["start0"], "group_length": first["length_aa"],
            "observed_orders": len(observations), "possible_permutations": math.factorial(first["component_count"]),
            "all_permutations_observed": len(observations) == math.factorial(first["component_count"]),
            "rank_min": min(ranks), "rank_max": max(ranks), "rank_range": max(ranks) - min(ranks), "rank_mean": mean(ranks),
            "rank_min_query_ids": [r["query_id"] for r in observations if r["rank"] == min(ranks)],
            "rank_max_query_ids": [r["query_id"] for r in observations if r["rank"] == max(ranks)],
            "score_min": min(scores), "score_max": max(scores), "score_range": max(scores) - min(scores), "score_mean": mean(scores),
            "r5_hits": sum(hits), "r5_hit_fraction": mean(hits), "r5_flip": min(hits) != max(hits),
            "single_order_control": len(observations) == 1,
            "position_support_balanced": len({row["observation_count"] for row in positions}) == 1,
            "position_means": positions, "observations": observations}


def aggregate_cells(cells):
    require(all(not c["single_order_control"] for c in cells), "Clean/single-order cells cannot enter permutation summaries")
    groups = defaultdict(list)
    for cell in cells:
        groups[cell["group_id"]].append(cell)
    metrics = ("rank_range", "score_range", "r5_flip")
    return {"group_target_cells": len(cells), "groups": len(groups), "targets": len({c["accession"] for c in cells}),
            "query_component_observations": sum(c["observed_orders"] for c in cells),
            "queries": len({r["query_id"] for c in cells for r in c["observations"]}),
            "cell_weighted": {key: distribution([c[key] for c in cells]) for key in metrics},
            "group_weighted_means": {key: mean(mean(c[key] for c in group) for group in groups.values()) for key in metrics},
            "r5_flip_cells": sum(c["r5_flip"] for c in cells),
            "incomplete_permutation_cells": sum(not c["all_permutations_observed"] for c in cells),
            "unbalanced_ordinal_position_cells": sum(not c["position_support_balanced"] for c in cells)}


def aggregate_positions(cells):
    result = []
    for count in (2, 3, 4):
        for position in range(1, count + 1):
            supported = [(c, next(p for p in c["position_means"] if p["position1"] == position))
                         for c in cells if c["component_count"] == count]
            supported = [(c, p) for c, p in supported if p["observation_count"]]
            result.append({"component_count": count, "position1": position, "group_target_cells": len(supported),
                           "groups": len({c["group_id"] for c, _ in supported}), "targets": len({c["accession"] for c, _ in supported}),
                           "order_observations": sum(p["observation_count"] for _, p in supported),
                           "queries": len({q for _, p in supported for q in p["query_ids"]}),
                           "equal_cell_position_means": {key: mean(p["means"][key] for _, p in supported) for key in ("rank", "score", "hit_at_5")},
                           "equal_cell_within_cell_centered_means": {key: mean(p["within_cell_centered_means"][key] for _, p in supported)
                                                                    for key in ("rank", "score", "hit_at_5")},
                           "causal_position_effect": False})
    return result


def read_arrays(snapshot, key):
    import numpy as np
    completion = read(snapshot.completion_path)
    path = base.verify_unit(snapshot.path, key, completion["units"][key], snapshot.identity_hash)
    with np.load(path / "arrays.npz", allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def references(snapshot, ctx, deadline):
    import numpy as np
    values = {}
    for accession in sorted(ctx["hashes"]):
        deadline.check()
        record = snapshot.record("catalogue-" + accession)
        require(record["accession"] == accession and record["sequence_sha256"] == ctx["hashes"][accession] and
                record["length"] == ctx["lengths"][accession], "Reference source metadata differs")
        arrays = read_arrays(snapshot, "catalogue-" + accession)
        require(set(arrays) == {"mean"}, "Unexpected catalogue vector schema")
        vector = arrays["mean"]
        require(vector.ndim == 1 and len(vector) > 0 and vector.dtype == np.float32 and np.isfinite(vector).all(), "Invalid catalogue vector")
        values[accession] = vector
    require(len({len(vector) for vector in values.values()}) == 1, "Catalogue vector widths differ")
    return values


def independent_cosine(left, right):
    numerator = math.fsum(float(a) * float(b) for a, b in zip(left, right))
    denominator = math.sqrt(math.fsum(float(a)**2 for a in left) * math.fsum(float(b)**2 for b in right))
    require(denominator > 0, "Undefined zero-norm mixture/reference")
    return numerator / denominator


def summarize_group(snapshot, group, refs, deadline):
    import numpy as np
    accessions = sorted(refs)
    components = sorted(group[0]["components"])
    lengths = {s["accession"]: s["end0"] - s["start0"] for s in group[0]["spans"]}
    width = len(refs[accessions[0]])
    reconstructed = [math.fsum(lengths[a] * float(refs[a][i]) for a in components) / group[0]["length_aa"] for i in range(width)]
    observations = {(arm, a): [] for arm in ARMS for a in components}
    anchor_vector = anchor_rows = None
    checks = []
    failures = []
    for query in group:
        deadline.check()
        record = snapshot.record("query-" + query["query_id"])
        require(record["query"] == query, "Saved order/query metadata differs")
        arrays = read_arrays(snapshot, "query-" + query["query_id"])
        expected = {"whole", "isolated_mixture"} | {arm + "_" + a for arm in ("isolated", "contextual") for a in components}
        require(set(arrays) == expected and all(v.shape == (width,) and v.dtype == np.float32 and np.isfinite(v).all()
                                              for v in arrays.values()), "Saved query vector schema/dtype/finite check failed")
        require(set(record["rankings"]) == {"whole", "isolated_mixture", "contextual_component", "isolated_component"}, "Unknown source ranking schema")
        require(set(record["rankings"]["contextual_component"]) == set(record["rankings"]["isolated_component"]) == set(components), "Missing oracle component ranking")
        for a in components:
            require(np.array_equal(arrays["isolated_" + a], refs[a]), "Saved isolated vector differs from catalogue")
            rank_map(record["rankings"]["isolated_component"][a], accessions)
        mappings = {arm: rank_map(record["rankings"][arm], accessions) for arm in ("whole", "isolated_mixture")}
        mappings["contextual_component"] = {a: rank_map(record["rankings"]["contextual_component"][a], accessions) for a in components}
        mixture = arrays["isolated_mixture"]
        reconstruction_error = max(abs(float(a) - b) for a, b in zip(mixture, reconstructed))
        if anchor_vector is None:
            anchor_vector, anchor_rows = mixture.copy(), record["rankings"]["isolated_mixture"]
            numeric_error = max(abs(independent_cosine(anchor_vector, refs[a]) - base.score(mappings["isolated_mixture"][a])) for a in accessions)
            if numeric_error > SCORE_TOLERANCE:
                failures.append({"query_id": query["query_id"], "check": "anchor_independent_cosines", "error": numeric_error})
        change = mixture_changes(anchor_rows, record["rankings"]["isolated_mixture"], accessions)
        vector_error = float(np.max(np.abs(anchor_vector.astype(np.float64) - mixture.astype(np.float64))))
        check = {"query_id": query["query_id"], "anchor_query_id": group[0]["query_id"], **change,
                 "vector_exactly_equal": bool(np.array_equal(anchor_vector, mixture)),
                 "maximum_absolute_vector_difference": vector_error, "reconstruction_maximum_absolute_error": reconstruction_error}
        check["passed"] = vector_error <= VECTOR_TOLERANCE and reconstruction_error <= VECTOR_TOLERANCE and change["scores_within_tolerance"] and change["all_rank_changes_within_near_ties"]
        if not check["passed"]:
            failures.append({"query_id": query["query_id"], "check": "mixture_invariance", "details": check})
        checks.append(check)
        spans = {s["accession"]: s for s in query["spans"]}
        for a in components:
            for arm in ARMS:
                row = mappings[arm][a][a] if arm == "contextual_component" else mappings[arm][a]
                observations[(arm, a)].append({"query_id": query["query_id"], "order": query["components"],
                    "position1": query["components"].index(a) + 1, "start0": spans[a]["start0"], "end0": spans[a]["end0"],
                    "rank": row["rank"], "score": base.score(row), "hit_at_5": int(row["rank"] <= 5)})
    cells = [describe_cell(group, a, arm, observations[(arm, a)]) for arm in ARMS for a in components]
    contrasts = []
    for a in components:
        whole = next(c for c in cells if c["arm"] == "whole" and c["accession"] == a)
        mixture = next(c for c in cells if c["arm"] == "isolated_mixture" and c["accession"] == a)
        contrasts.append({"group_id": whole["group_id"], "accession": a, "suite": whole["suite"],
                          **{key + "_whole_minus_mixture": whole[key] - mixture[key] for key in ("rank_range", "score_range", "r5_flip")}})
    return {"group_id": group[0]["group_id"], "suite": group[0]["suite"], "query_ids": [q["query_id"] for q in group],
            "cells": cells, "whole_minus_mixture": contrasts,
            "mixture_sanity": {"status": "failed" if failures else "passed", "anchor_independent_cosine_max_abs_error": numeric_error,
                               "checks": checks, "failures": failures}, "no_rank_or_vector_equality_forced": True}


def model_summary(records, design):
    combined = [r for r in records if r["suite"] != "clean"]
    clean = [r for r in records if r["suite"] == "clean"]
    cells = [c for r in combined for c in r["cells"]]
    result = {"status": "passed" if all(r["mixture_sanity"]["status"] == "passed" for r in records) else "failed_sanity",
              "combined_groups": len(combined), "clean_queries": len(clean), "arms": {},
              "failed_sanity_groups": [r["group_id"] for r in records if r["mixture_sanity"]["status"] != "passed"],
              "clean_single_order_sanity_not_permutation_evidence": True,
              "clean_summary": {arm: {"observations": len(clean),
                   "mean_rank": mean(c["rank_mean"] for r in clean for c in r["cells"] if c["arm"] == arm),
                   "hit_at_5_count": sum(c["r5_hits"] for r in clean for c in r["cells"] if c["arm"] == arm)} for arm in ARMS}}
    for arm in ARMS:
        selected = [c for c in cells if c["arm"] == arm]
        require(len(selected) == 680 and sum(c["observed_orders"] for c in selected) == 3520, "Missing combined group-target cells")
        result["arms"][arm] = {"information_and_budget": design["arms"][arm],
            "by_suite": {scope: aggregate_cells([c for c in selected if scope == "combined" or c["suite"] == scope]) for scope in (*SUITES, "combined")},
            "supported_positions": aggregate_positions(selected)}
    contrasts = [c for r in combined for c in r["whole_minus_mixture"]]
    result["whole_minus_mixture"] = {scope: {key: distribution([c[key] for c in contrasts if scope == "combined" or c["suite"] == scope])
        for key in ("rank_range_whole_minus_mixture", "score_range_whole_minus_mixture", "r5_flip_whole_minus_mixture")}
        for scope in (*SUITES, "combined")}
    result["sanity_totals"] = {"orders_checked": sum(len(r["mixture_sanity"]["checks"]) for r in records),
        "orders_with_changed_mixture_ranks": sum(bool(c["changed_candidate_ranks"]) for r in records for c in r["mixture_sanity"]["checks"]),
        "inverted_mixture_candidate_pairs": sum(len(c["inverted_candidate_pairs"]) for r in records for c in r["mixture_sanity"]["checks"]),
        "non_near_tie_inversions": sum(not pair["within_near_tie_gap"] for r in records for c in r["mixture_sanity"]["checks"] for pair in c["inverted_candidate_pairs"])}
    return result


def run(args):
    require(args.execute_saved_analysis, "Explicit lead-admitted saved-output execution flag required")
    session = args.session.resolve()
    require(session.is_relative_to(ROOT / "results"), "Session escaped results")
    design_path = session / "inputs/order_sensitivity_design_v1.json"
    admission_path = session / "inputs/order_sensitivity_model_admission_v1.json"
    require(sha(__file__) == args.expected_code_sha256 and sha(design_path) == DESIGN_SHA == args.expected_design_sha256 and
            sha(admission_path) == args.expected_admission_sha256, "Frozen code/design/roster invocation differs")
    design, admission = read(design_path), read(admission_path)
    require(design["version"] == "iri-saved-order-sensitivity-v1" and admission["design_sha256"] == DESIGN_SHA and
            admission["required_model_labels"] == list(MODELS) and [s["label"] for s in admission["sources"]] == list(MODELS) and
            admission["allow_missing_model_or_query_intersection"] is False, "Incomplete prospective model roster")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    started = time.monotonic()
    deadline = Deadline(600, session)
    ctx = base.context(session, deadline)
    groups = group_queries(list(ctx["queries"].values()))
    support = coverage(groups)
    require(support == {"queries": 1146, "suite_queries": design["population"]["suite_queries"],
                        "combined_groups": design["population"]["combined_groups"], "combined_group_target_cells": 680,
                        "combined_query_component_observations": 3520}, "Frozen order/group support differs")
    snapshots = {}
    for declared, pin in zip(design["models"], admission["sources"]):
        deadline.check()
        require(declared["label"] == pin["label"] and sha(session / declared["config_path"]) == declared["config_sha256"], "Pinned model configuration changed")
        snapshot = base.Snapshot(session / pin["job"], session, deadline, pin["completion_sha256"])
        identity = snapshot.identity
        require(snapshot.identity_hash == pin["identity_sha256"] and len(snapshot.units) == 1252, "Source identity/units differ")
        require_full_panel(identity, ctx["queries"], ctx["hashes"])
        require(identity["challenge_sha256"] == design["challenge_sha256"] and identity["catalogue_sequence_sha256"] == ctx["hashes"] and
                identity["config"] == read(session / declared["config_path"]) and identity["config"]["model_id"] == declared["model_id"], "Source representation/sequence identity differs")
        for relative, expected in identity["code_hashes"].items():
            require(sha(base.safe_file(ROOT, relative)) == expected, "Source implementation changed")
        snapshots[declared["label"]] = snapshot
    code = {p: sha(ROOT / p) for p in ("scripts/session_order_analysis.py", "tests/test_session_order_analysis.py",
        "scripts/session_extension_audit.py", "scripts/session_independent_audit.py", "src/iri_constructs/session_store.py")}
    identity = {"mode": "saved_order_sensitivity_v1", "design_sha256": DESIGN_SHA, "roster_sha256": sha(admission_path),
                "sources": {label: snapshot.evidence() for label, snapshot in snapshots.items()}, "code_sha256": code,
                "query_ids": sorted(ctx["queries"]), "support": support, "interpretation": design["interpretation"]}
    store = UnitStore(session / "analysis" / (safe_name(args.name) + "-" + digest(canonical(identity))[:12]), identity)
    expected, summaries = [], {}
    for label, snapshot in snapshots.items():
        refs = references(snapshot, ctx, deadline)
        records = []
        for gid, queries in groups.items():
            deadline.check()
            key = label + "--" + gid
            result = summarize_group(snapshot, queries, refs, deadline)
            result["model_label"] = label
            if store.read(key) is None:
                store.publish(key, result)
            else:
                require(store.read(key) == result, "Resumed order analysis differs")
            expected.append(key)
            records.append(result)
        summaries[label] = model_summary(records, design)
        deadline.check()
        summary_key = label + "--summary"
        if store.read(summary_key) is None:
            store.publish(summary_key, summaries[label])
        else:
            require(store.read(summary_key) == summaries[label], "Model summary changed")
        expected.append(summary_key)
        del refs, records
        print(json.dumps({"model": label, "status": summaries[label]["status"], "groups": len(groups)}), flush=True)
    for snapshot in snapshots.values():
        require(base.canonical_hash(read(snapshot.path / "identity.json")) == snapshot.identity_hash and
                sha(snapshot.completion_path) == snapshot.completion_sha256, "Source snapshot changed during analysis")
    summary = {"status": "passed" if all(v["status"] == "passed" for v in summaries.values()) else "failed_sanity",
               "models": summaries, "common_support": support, "model_count": 5,
               "no_cross_model_cosine_calibration": True, "pooled_iid_inference": False,
               "interpretation": design["interpretation"], "design": design, "source_roster": admission}
    deadline.check()
    if store.read("summary") is None:
        store.publish("summary", summary)
    else:
        require(store.read("summary") == summary, "Final order summary changed")
    expected.append("summary")
    store.complete(expected, {"analysis_status": summary["status"], "accounting_complete": True,
                             "elapsed_seconds": time.monotonic() - started, "peak_rss_mib": peak_mib()})
    print(json.dumps({"job": str(store.root), "status": summary["status"], "models": 5, "queries_per_model": 1146}), flush=True)
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", type=Path, required=True)
    parser.add_argument("--name", default="order-sensitivity-v1")
    parser.add_argument("--execute-saved-analysis", action="store_true")
    parser.add_argument("--expected-code-sha256", required=True)
    parser.add_argument("--expected-design-sha256", required=True)
    parser.add_argument("--expected-admission-sha256", required=True)
    sys.exit(run(parser.parse_args()))
