#!/usr/bin/env python3
"""Independent saved-result audit; deliberately does not import benchmark metrics.

Reads a completed job, its immutable descriptors and frozen query manifest. Checks
every result/denominator and a deterministic bounded sample of saved vector scores.
Run after the model worker releases memory, as required by the hardware contract.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
KS = (1, 5, 10)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def equal(observed, expected, context, tolerance=1e-10):
    if isinstance(expected, float):
        require(isinstance(observed, (int, float)) and math.isfinite(observed)
                and abs(observed - expected) <= tolerance, context)
    else:
        require(observed == expected, context)


def compare_fields(observed, expected, context):
    for field, value in expected.items():
        require(field in observed, context + ": missing " + field)
        equal(observed[field], value, context + ": " + field)


def score(row):
    value = float(row["score"])
    require(math.isfinite(value), "Nonfinite saved ranking score")
    if "cosine_similarity" in row:
        equal(row["cosine_similarity"], value, "Conflicting cosine/score aliases")
    return value


def validate_ranking(rows, accessions):
    require(len(rows) == len(accessions), "Incomplete ranking length")
    require([row["rank"] for row in rows] == list(range(1, len(rows) + 1)), "Ranking positions are not contiguous")
    names = [row["accession"] for row in rows]
    require(len(set(names)) == len(names) and set(names) == set(accessions), "Ranking catalogue mismatch")
    ordered = sorted(rows, key=lambda row: (-score(row), row["accession"]))
    require(rows == ordered, "Ranking order disagrees with score/accession rule")
    counts = Counter(score(row) for row in rows)
    for row in rows:
        equal(row["exact_score_tie_count"], counts[score(row)], "Incorrect exact tie count")


def expected_match(rows, acceptable):
    surviving = set(acceptable) & {r["accession"] for r in rows}
    require(surviving, "Primary benchmark component unexpectedly absent")
    match = next(r for r in rows if r["accession"] in surviving)
    rank, value = match["rank"], score(match)
    tied = [r["rank"] for r in rows if score(r) == value]
    near = [r["rank"] for r in rows if abs(score(r) - value) <= 1e-6]
    return {"scorable": True, "accepted_surviving": sorted(surviving),
            "best_accession": match["accession"], "rank": rank, "score": value,
            "reciprocal_rank": 1 / rank, "tie_best_rank": min(tied), "tie_worst_rank": max(tied),
            "near_tie_count": len(near), "near_tie_best_rank": min(near), "near_tie_worst_rank": max(near),
            **{f"hit_at_{k}": int(rank <= k) for k in KS}}


def bin_name(value):
    if value < .25: return "lt_0.25"
    if value < .4: return "0.25_to_lt_0.40"
    if value < .6: return "0.40_to_lt_0.60"
    if value < .75: return "0.60_to_lt_0.75"
    return "ge_0.75"


def validate_evaluation(saved, query, ranking, equivalence):
    components = query["components"]
    oracle = isinstance(ranking, dict)
    require(len(saved["components"]) == len(components), "Evaluation omitted a component")
    compare_fields(saved, {"query_id": query["query_id"], "suite": query["suite"],
                   "group_id": query["group_id"], "component_count": len(components),
                   "catalogue_count": len(equivalence),
                   "retrieval_budget": "k per known component" if oracle else "k total per original query"},
                   "Evaluation identity")
    independent = []
    for i, (accession, span) in enumerate(zip(components, query["spans"])):
        rows = ranking[accession] if oracle else ranking
        position = "single" if len(components) == 1 else ("N" if i == 0 else "C" if i == len(components) - 1 else "internal")
        expected = {"accession": accession, "component_index": i, "position": position,
                    "start0": span["start0"], "end0": span["end0"], "fraction": span["fraction"],
                    "fraction_bin": bin_name(span["fraction"])}
        compare_fields(saved["components"][i], expected, "Component identity")
        for kind, accepted in (("exact", [accession]), ("duplicate_equivalent", equivalence[accession])):
            expected[kind] = expected_match(rows, accepted)
            compare_fields(saved["components"][i][kind], expected[kind], "Component " + kind)
        independent.append(expected)
    for kind in ("exact", "duplicate_equivalent"):
        fields = {"component_count": len(components), "scorable_components": len(components),
                  "unscorable_components": 0, "all_components_scorable": True}
        for k in KS:
            hits = sum(c[kind][f"hit_at_{k}"] for c in independent)
            fields.update({f"component_hits_at_{k}": hits, f"component_recall_at_{k}": hits / len(components),
                           f"any_component_hit_at_{k}": int(hits > 0),
                           f"any_scorable_component_hit_at_{k}": int(hits > 0),
                           f"all_components_hit_at_{k}": int(hits == len(components))})
            equal(saved["exact_recall_ceiling"][str(k)], 1.0 if oracle else min(1, k / len(components)), "Recall ceiling")
            equal(saved["returned_hits_ceiling"][str(k)], min(k, len(equivalence)) * (len(components) if oracle else 1), "Review budget ceiling")
        compare_fields(saved["query_metrics"][kind], fields, "Query " + kind)
    return {"query_id": query["query_id"], "group_id": query["group_id"], "suite": query["suite"],
            "component_count": len(components), "retrieval_budget": saved["retrieval_budget"],
            "components": independent}


def row_summary(entries, kind):
    n = len(entries)
    targets = defaultdict(list)
    for query, component in entries:
        targets[component["accession"]].append(component[kind])
    result = {"original_components": n, "scorable_components": n, "unscorable_components": 0,
              "parent_queries": len({q["query_id"] for q, _ in entries}),
              "distinct_groups": len({q["group_id"] for q, _ in entries}), "scorable_targets": len(targets),
              "mean_rank": math.fsum(c[kind]["rank"] for _, c in entries) / n if n else None,
              "mrr": math.fsum(1 / c[kind]["rank"] for _, c in entries) / n if n else None}
    for k in KS:
        hits = sum(c[kind]["rank"] <= k for _, c in entries)
        result.update({f"hits_at_{k}": hits, f"recall_at_{k}": hits / n if n else None,
                       f"per_target_macro_recall_at_{k}": math.fsum(sum(c["rank"] <= k for c in cs) / len(cs)
                              for cs in targets.values()) / len(targets) if targets else None})
    return result


def validate_summary(saved, evaluations):
    if not evaluations:
        compare_fields(saved, {"query_count": 0, "component_count": 0, "kinds": {}}, "Empty summary")
        return
    compare_fields(saved, {"query_count": len(evaluations),
                   "distinct_group_count": len({e["group_id"] for e in evaluations}),
                   "component_count": sum(e["component_count"] for e in evaluations),
                   "retrieval_budget": evaluations[0]["retrieval_budget"]}, "Summary totals")
    entries = [(q, c) for q in evaluations for c in q["components"]]
    for kind in ("exact", "duplicate_equivalent"):
        observed = saved["kinds"][kind]
        compare_fields(observed, row_summary(entries, kind), "Summary components " + kind)
        query_values = {"original_queries": len(evaluations), "all_components_scorable_queries": len(evaluations),
                        "excluded_unscorable_queries": 0}
        for k in KS:
            any_count = sum(any(c[kind]["rank"] <= k for c in q["components"]) for q in evaluations)
            all_count = sum(all(c[kind]["rank"] <= k for c in q["components"]) for q in evaluations)
            query_values.update({f"any_component_hit_count_at_{k}": any_count,
                                 f"all_components_hit_count_at_{k}": all_count,
                                 f"any_component_hit_at_{k}": any_count / len(evaluations),
                                 f"all_components_hit_at_{k}": all_count / len(evaluations)})
        compare_fields(observed["query_metrics"], query_values, "Summary query denominators")
        groupers = {"by_component_count": lambda q, c: str(q["component_count"]),
                    "by_fraction": lambda q, c: c["fraction_bin"], "by_position": lambda q, c: c["position"],
                    "by_target": lambda q, c: c["accession"]}
        for name, get in groupers.items():
            groups = defaultdict(list)
            for q, c in entries:
                groups[get(q, c)].append((q, c))
            require(set(observed[name]) == set(groups), "Summary stratum keys mismatch")
            for value, group in groups.items():
                compare_fields(observed[name][value], row_summary(group, kind), "Summary " + name + "/" + value)
                # Later summary schema adds query-denominator fields within
                # component-count strata. Audit them when present while retaining
                # compatibility with already frozen earlier primary summaries.
                if name == "by_component_count" and "query_metrics" in observed[name][value]:
                    distinct = {q["query_id"]: q for q, _ in group}
                    n = len(distinct)
                    nested = {"original_queries": n, "all_components_scorable_queries": n,
                              "excluded_unscorable_queries": 0}
                    for k in KS:
                        for metric, predicate in (("any_component_hit", any), ("all_components_hit", all)):
                            hits = sum(predicate(c[kind]["rank"] <= k for c in q["components"]) for q in distinct.values())
                            nested[f"{metric}_count_at_{k}"] = hits
                            nested[f"{metric}_at_{k}"] = hits / n
                    compare_fields(observed[name][value]["query_metrics"], nested, "Count-stratum query denominators")


def verify_unit(job, key, index, identity_hash):
    path = (job / index["path"]).resolve()
    require(path.is_relative_to(job.resolve()), "Completion path escaped job root")
    require(sha(path / "complete.json") == index["descriptor_sha256"], "Unit descriptor hash mismatch")
    descriptor = load(path / "complete.json")
    require(descriptor["key"] == key and descriptor["identity_hash"] == identity_hash, "Unit descriptor identity mismatch")
    require("record.json" in descriptor["files"], "Unit missing saved record")
    for filename, expected in descriptor["files"].items():
        require(Path(filename).name == filename and filename not in (".", ".."), "Unsafe descriptor filename")
        target = path / filename
        require(not target.is_symlink(), "Unit contains unexpected symlink")
        require(target.stat().st_size == expected["bytes"] and sha(target) == expected["sha256"], "Unit bytes/hash mismatch")
    return path


def numeric_sample(query, record, path, reference_vectors, accessions, tolerance):
    import numpy as np
    with np.load(path / "arrays.npz", allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    components = query["components"]
    required = {"whole", "isolated_mixture"} | {f"{arm}_{a}" for arm in ("contextual", "isolated") for a in components}
    require(set(arrays) == required, "Query vector keys mismatch")
    width = len(reference_vectors[accessions[0]])
    require(all(v.shape == (width,) and v.dtype == np.float32 and np.isfinite(v).all() for v in arrays.values()), "Invalid saved vector shapes/dtype")
    total = query["length_aa"]
    lengths = {s["accession"]: s["end0"] - s["start0"] for s in query["spans"]}
    mixture = sum(lengths[a] * reference_vectors[a].astype(np.float64) for a in sorted(components)) / total
    reconstructed = sum(lengths[a] * arrays["contextual_" + a].astype(np.float64) for a in sorted(components)) / total
    require(np.max(np.abs(mixture - arrays["isolated_mixture"])) <= tolerance, "Saved mixture disagrees with isolated vectors")
    require(np.max(np.abs(reconstructed - arrays["whole"])) <= tolerance, "Weighted contextual identity failed independently")
    for a in components:
        require(np.array_equal(arrays["isolated_" + a], reference_vectors[a]), "Isolated vector differs from catalogue vector")
    vectors = [("whole", None, arrays["whole"]), ("isolated_mixture", None, arrays["isolated_mixture"])]
    for a in components:
        vectors += [("contextual_component", a, arrays["contextual_" + a]), ("isolated_component", a, arrays["isolated_" + a])]
    maximum, comparisons, disagreements, near_disagreements = 0.0, 0, [], []
    for arm, accession, vector in vectors:
        # Python fsum dot products independently check the NumPy ranking arithmetic.
        q = [float(x) for x in vector]
        qnorm = math.sqrt(math.fsum(x * x for x in q))
        recalculated = {}
        for a in accessions:
            r = [float(x) for x in reference_vectors[a]]
            recalculated[a] = math.fsum(x * y for x, y in zip(q, r)) / (qnorm * math.sqrt(math.fsum(y * y for y in r)))
        ranking = record["rankings"][arm][accession] if accession else record["rankings"][arm]
        for row in ranking:
            error = abs(recalculated[row["accession"]] - score(row))
            maximum = max(maximum, error); comparisons += 1
            require(error <= tolerance, "Independent cosine tolerance exceeded")
        order = sorted(accessions, key=lambda a: (-recalculated[a], a))
        saved_order = [r["accession"] for r in ranking]
        if order != saved_order:
            changed = [a for a, b in zip(order, saved_order) if a != b]
            numerical_tie_only = all(abs(recalculated[a] - recalculated[b]) <= tolerance
                                     for a, b in zip(order, saved_order) if a != b)
            item = {"arm": arm, "component": accession, "changed_accessions": changed,
                    "within_score_tolerance": numerical_tie_only}
            disagreements.append(item)
            if numerical_tie_only: near_disagreements.append(item)
            require(numerical_tie_only, "Independent ranking disagreement exceeds near-tie tolerance")
    return {"query_id": query["query_id"], "cosine_comparisons": comparisons,
            "max_abs_cosine_error": maximum, "ranking_disagreements": disagreements,
            "near_tie_only_disagreements": near_disagreements}


def audit(job, manifest_path=None, numerical_samples=8, tolerance=1e-6, completion_path=None):
    start = time.monotonic(); job = Path(job).resolve()
    identity = load(job / "identity.json"); identity_hash = canonical_hash(identity)
    query_ids, accessions = identity["query_ids"], identity["catalogue_accessions"]
    require(len(set(query_ids)) == len(query_ids) and query_ids, "Invalid query identity list")
    require(len(set(accessions)) == len(accessions) == 105, "Primary catalogue is not 105 unique accessions")
    manifest_path = Path(manifest_path) if manifest_path else job.parent.parent / "inputs/challenge_manifest.json"
    require(sha(manifest_path) == identity["challenge_sha256"], "Frozen query manifest mismatch")
    manifest = load(manifest_path); queries = {q["query_id"]: q for q in manifest["queries"]}
    require(set(query_ids) <= set(queries), "Job contains query outside frozen suite")
    if completion_path is None:
        candidates = [(load(p)["completed_utc"], p) for p in job.glob("completion-*.json")]
        require(candidates, "Job has no completion descriptor")
        completion_path = max(candidates)[1]
    completion = load(completion_path)
    require(completion["identity_hash"] == identity_hash, "Job completion identity mismatch")
    expected = {"query-" + q for q in query_ids} | {"summary"}
    if identity["mode"] == "model": expected |= {"catalogue-" + a for a in accessions}
    require(set(completion["units"]) == expected and completion["expected_count"] == len(expected), "Completion coverage mismatch")
    paths = {key: verify_unit(job, key, index, identity_hash) for key, index in completion["units"].items()}
    seq_hashes = identity["catalogue_sequence_sha256"]
    require(set(seq_hashes) == set(accessions), "Catalogue hash coverage mismatch")
    equivalence = {a: [b for b in accessions if seq_hashes[a] == seq_hashes[b]] for a in accessions}
    expected_arms = {"whole", "isolated_mixture", "contextual_component", "isolated_component"} if identity["mode"] == "model" else {"local_alignment"}
    by_arm = {arm: [] for arm in expected_arms}
    rankings_checked = component_metrics_checked = 0
    # Fixed query-ID hashes, with one per component-count suite before extras.
    selected = []
    for suite in sorted({queries[q]["suite"] for q in query_ids}):
        candidates = [q for q in query_ids if queries[q]["suite"] == suite]
        selected.append(min(candidates, key=lambda q: hashlib.sha256(q.encode()).hexdigest()))
    selected += [q for q in sorted(query_ids, key=lambda q: hashlib.sha256(q.encode()).hexdigest()) if q not in selected]
    selected = selected[:max(0, numerical_samples)] if identity["mode"] == "model" else []
    records = {}
    for qid in query_ids:
        query, record = queries[qid], load(paths["query-" + qid] / "record.json")
        for key, value in query.items(): equal(record["query"][key], value, "Saved query disagrees with frozen input: " + key)
        require("sequence" not in record["query"], "Raw sequence unexpectedly saved in public query metadata")
        require(set(record["rankings"]) == set(record["evaluations"]) == expected_arms, "Missing or unexpected arm")
        if identity["mode"] == "model":
            descriptor = load(paths["query-" + qid] / "complete.json")
            require("arrays.npz" in descriptor["files"], "Model query lacks checksum-covered arrays")
            encode = record["encode"]
            compare_fields(encode, {"residue_count": query["length_aa"], "token_count": query["length_aa"] + 2,
                           "truncated": False, "model_training": False, "inference_mode": True}, "Query representation metadata")
            require(encode["residue_shape"][0] == query["length_aa"] and encode["hidden_shape"][:2] == [1, query["length_aa"] + 2],
                    "Query token/residue shape mismatch")
            for field in ("whole_weighted_context_max_abs", "mixture_reverse_max_abs"):
                value = record["controls"][field]
                require(math.isfinite(value) and 0 <= value <= 1e-6, "Recorded control identity failed")
            if len(query["components"]) == 1:
                value = record["controls"]["independent_clean_reference_max_abs"]
                require(math.isfinite(value) and 0 <= value <= 1e-6, "Clean/reference control failed")
        for arm, ranking in record["rankings"].items():
            if isinstance(ranking, dict):
                require(set(ranking) == set(query["components"]), "Oracle component coverage mismatch")
                lists = list(ranking.values())
            else: lists = [ranking]
            for rows in lists: validate_ranking(rows, accessions); rankings_checked += 1
            by_arm[arm].append(validate_evaluation(record["evaluations"][arm], query, ranking, equivalence))
            component_metrics_checked += len(query["components"]) * 2
        if qid in selected:
            records[qid] = record
    summary_record = load(paths["summary"] / "record.json")
    require(summary_record["query_count"] == len(query_ids), "Summary query count mismatch")
    require(set(summary_record["summary"]) == expected_arms, "Summary arm coverage mismatch")
    for arm, evaluations in by_arm.items():
        summary = summary_record["summary"][arm]
        validate_summary(summary["all"], evaluations)
        suites = {e["suite"] for e in evaluations}
        require(set(summary["by_suite"]) == suites, "Summary suite coverage mismatch")
        for suite in suites:
            validate_summary(summary["by_suite"][suite], [e for e in evaluations if e["suite"] == suite])
    numerical = []
    if identity["mode"] == "model":
        import numpy as np
        references = {}
        for accession in accessions:
            path = paths["catalogue-" + accession]; record = load(path / "record.json")
            require("arrays.npz" in load(path / "complete.json")["files"], "Catalogue lacks checksum-covered arrays")
            require(record["accession"] == accession and record["sequence_sha256"] == seq_hashes[accession], "Cached reference identity mismatch")
            with np.load(path / "arrays.npz", allow_pickle=False) as z:
                require(z.files == ["mean"], "Unexpected catalogue array keys")
                vector = z["mean"].copy()
            shape = record["encode"]["residue_shape"]
            require(vector.shape == (shape[1],) and record["length"] == shape[0]
                    and vector.dtype == np.float32 and np.isfinite(vector).all(), "Invalid cached catalogue dimensions")
            references[accession] = vector
        require(len({v.shape for v in references.values()}) == 1, "Catalogue widths differ")
        for qid in selected:
            numerical.append(numeric_sample(queries[qid], records[qid], paths["query-" + qid], references, accessions, tolerance))
    result = {"status": "passed", "job": str(job), "identity_hash": identity_hash,
              "completion_descriptor": str(completion_path), "completion_sha256": sha(completion_path),
              "query_count": len(query_ids), "catalogue_count": len(accessions), "valid_units": len(paths),
              "rankings_checked": rankings_checked, "component_metric_records_checked": component_metrics_checked,
              "all_summary_denominators_independently_checked": True,
              "numerical_samples": numerical, "numeric_tolerance": tolerance,
              "numerical_scope": "Saved vector cosine, mixture and contextual weighted identities; no new inference or residue-state reconstruction",
              "elapsed_seconds": time.monotonic() - start,
              "auditor_code_sha256": sha(__file__), "audited_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("job", type=Path); p.add_argument("--manifest", type=Path)
    p.add_argument("--completion", type=Path); p.add_argument("--numerical-samples", type=int, default=8)
    p.add_argument("--tolerance", type=float, default=1e-6); p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    result = audit(a.job, a.manifest, a.numerical_samples, a.tolerance, a.completion)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=a.output.parent, delete=False) as f:
        json.dump(result, f, indent=2, allow_nan=False); f.write("\n"); f.flush(); os.fsync(f.fileno()); temp = Path(f.name)
    os.replace(temp, a.output)
    print(json.dumps({key: value for key, value in result.items() if key != "numerical_samples"}))


if __name__ == "__main__":
    main()
