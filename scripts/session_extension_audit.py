#!/usr/bin/env python3
"""Independent streaming audit of frozen catalogue/family and total-budget controls.

No production metric, mask, fusion or SummaryAccumulator imports. No inference,
alignment or NumPy import. Counts are rebuilt from outcome-frequency tables.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import datetime as dt
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import signal
import tempfile
import time

from session_independent_audit import verify_unit, validate_ranking, score, bin_name

ROOT = Path(__file__).resolve().parents[1]
CHALLENGE_SHA256 = "7e671e4d281f3f107be36e88bf1c0eeb7a487ac40525ffbabd026b1c7e9608b0"
FAMILY_SHA256 = "92a2bb852426948fb627f4029e46af7c97dcdbafbf7ad579f57fa5e7a5155d94"
ORACLE_DESIGN_SHA256 = "228a6a2d19b9637c04e0bb62a3aa5862833fd6bfc80e4abb03c2f8c5906f5316"
KS = (1, 5, 10)
MODES = {"catalogue_and_family_analysis", "oracle_total_budget_v1"}
DEFAULT_JOBS = (
    "catalogue-family-swnorm-v1-49730919f2c6", "catalogue-family-8m-v1-0c972ca22a95",
    "catalogue-family-swraw-v1-aaf772ce6584", "catalogue-family-esmc300m-v1-532f89eecdc8",
    "catalogue-family-650m-v1-983fdbdf3d01", "oracle-total-650m-v1-b39d17272eb6",
    "oracle-total-esmc300m-v1-ef761ba35367", "oracle-total-8m-v1-4ea2ff4ee781",
)


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def same(observed, expected, where="root"):
    """Strict structural comparison; only finite numerical roundoff is tolerated."""
    if isinstance(expected, dict):
        require(isinstance(observed, dict) and set(observed) == set(expected), where + ": keys differ")
        for key, value in expected.items():
            same(observed[key], value, where + "/" + str(key))
    elif isinstance(expected, list):
        require(isinstance(observed, list) and len(observed) == len(expected), where + ": length differs")
        for index, value in enumerate(expected):
            same(observed[index], value, where + "/" + str(index))
    elif isinstance(expected, float):
        require(type(observed) in (int, float) and math.isfinite(observed) and abs(observed - expected) <= 1e-10,
                where + ": finite numerical value differs")
    else:
        require(type(observed) is type(expected) and observed == expected, where + ": value/type differs")


def fields(observed, expected, where):
    for key, value in expected.items():
        require(key in observed, where + ": missing " + key)
        same(observed[key], value, where + "/" + key)


class Deadline:
    def __init__(self, seconds, session):
        require(0 < seconds <= 1800, "Audit timeout must be between 0 and 1800 seconds")
        stop = dt.datetime.fromisoformat(read(session / "launch.json")["hard_stop_utc"])
        require(stop.tzinfo is not None, "Session stop lacks timezone")
        self.end = time.monotonic() + min(seconds, (stop - dt.datetime.now(dt.timezone.utc)).total_seconds())
    def check(self):
        if time.monotonic() >= self.end:
            raise TimeoutError("Independent audit or session deadline reached")


def safe_file(root, relative):
    path = (root / relative).resolve()
    require(path.is_relative_to(root.resolve()), "Input path escaped its owning root")
    return path


def rebuild_catalogues(manifest, sequence_hashes):
    groups = defaultdict(list)
    for accession in sorted(sequence_hashes):
        groups[sequence_hashes[accession]].append(accession)
    targets = {sequence_hashes[a] for a in manifest["targets"]}
    decoys = set(groups) - targets
    expected = []
    for replicate in range(3):
        ordered = sorted(decoys, key=lambda value: hashlib.sha256(
            (manifest["selection_seed"] + "|decoys:" + str(replicate) + "|" + value).encode()).hexdigest())
        for size in (25, 50, len(groups)):
            chosen = sorted(targets) + ordered[:size - len(targets)]
            expected.append({"subset_id": f"unique{size}_rep{replicate}", "replicate": replicate,
                "unique_sequence_count": size,
                "entries": [{"sequence_sha256": h, "representative_accession": min(groups[h]),
                             "equivalent_accessions": groups[h]} for h in chosen]})
    same(manifest["catalogue_subsets"], expected, "Rebuilt catalogue subsets")
    masks, aliases, seen = {"accession105": sorted(sequence_hashes)}, {}, {}
    for subset in expected:
        accessions = sorted(row["representative_accession"] for row in subset["entries"])
        key = tuple(accessions)
        if key in seen:
            aliases[subset["subset_id"]] = seen[key]
        else:
            masks[subset["subset_id"]] = accessions
            seen[key] = subset["subset_id"]
    return {"masks": masks, "aliases": aliases}


def rebuild_removals(labels, sequence_hashes, lengths, pair_rows):
    names = sorted(sequence_hashes)
    targets = sorted(labels["acceptable_matches"])
    exact = {a: sorted(b for b in names if sequence_hashes[a] == sequence_hashes[b]) for a in targets}
    near = {a: {a} for a in targets}
    seen = set()
    for row in pair_rows:
        a, b = row["a"], row["b"]
        require(a in sequence_hashes and b in sequence_hashes and a < b and (a, b) not in seen, "Invalid duplicate pairwise row")
        seen.add((a, b))
        identical, columns = int(row["identical_residues"]), int(row["alignment_columns"])
        identity = float(row["global_identity"])
        require(0 <= identical <= min(lengths[a], lengths[b]) and max(lengths[a], lengths[b]) <= columns <= lengths[a] + lengths[b] and
                math.isfinite(identity) and abs(identity - identical / columns) <= 1e-12 and
                row["exact_duplicate"] == str(sequence_hashes[a] == sequence_hashes[b]), "Invalid frozen pairwise arithmetic")
        if identity >= 0.9:
            if a in near:
                near[a].add(b)
            if b in near:
                near[b].add(a)
    require(seen == set(itertools.combinations(names, 2)), "Frozen pairwise coverage incomplete")
    result = {"exact_removed": exact, "near90_removed": {a: sorted(values) for a, values in near.items()}}
    for policy, values in result.items():
        same(labels["removal_policies"][policy]["removed_accessions_by_target"], values, "Rebuilt " + policy)
    return result


def family_masks(query, accessions, removals):
    result = {}
    for policy, by_target in removals.items():
        removed = sorted({a for target in query["components"] for a in by_target[target]})
        result["family_" + policy] = {"included_accessions": sorted(set(accessions) - set(removed)),
                                    "removed_accessions": removed, "policy": policy}
    return result


def denominator_rows(queries, accessions, families, removals, policy):
    result = []
    for query in queries:
        included = set(family_masks(query, accessions, removals)["family_" + policy]["included_accessions"])
        counts = {a: len(included & set(families[a])) for a in query["components"]}
        n = sum(count > 0 for count in counts.values())
        result.append({"query_id": query["query_id"], "group_id": query["group_id"], "suite": query["suite"],
                       "component_count": len(counts), "acceptable_surviving_counts": counts,
                       "scorable_components": n, "unscorable_components": len(counts) - n,
                       "all_components_scorable": n == len(counts)})
    return result


def denominator_totals(rows):
    def totals(selected):
        n = len(selected)
        full = sum(row["all_components_scorable"] for row in selected)
        return {"original_queries": n, "original_component_observations": sum(r["component_count"] for r in selected),
                "scorable_component_observations": sum(r["scorable_components"] for r in selected),
                "unscorable_component_observations": sum(r["unscorable_components"] for r in selected),
                "all_components_scorable_queries": full, "excluded_from_all_components_metric": n - full}
    return {"all": totals(rows), "combined": totals([r for r in rows if r["component_count"] > 1]),
            "by_suite": {suite: totals([r for r in rows if r["suite"] == suite]) for suite in sorted({r["suite"] for r in rows})}}


def context(session, deadline):
    manifest_path = session / "inputs/challenge_manifest.json"
    label_path = session / "inputs/family_labels_v1.json"
    require(sha(manifest_path) == CHALLENGE_SHA256 and sha(label_path) == FAMILY_SHA256, "Frozen audit input hash differs")
    manifest, labels = read(manifest_path), read(label_path)
    for relative, expected in manifest["source_files_sha256"].items():
        deadline.check()
        require(sha(safe_file(ROOT, relative)) == expected, "Frozen source input changed")
    catalogue = read(ROOT / "data/manifests/t05_expanded_catalogue.json")
    sequences, hashes, lengths = {}, {}, {}
    for entry in catalogue["records"]:
        deadline.check()
        a = entry["accession"]
        require(a not in sequences and all(entry[key] == "approved" for key in ("review_status", "licensing_status", "safety_review")), "Unreviewed or duplicate source record")
        raw_path = safe_file(ROOT, entry["raw_path"])
        require(sha(raw_path) == entry["checksum_sha256"], "Original catalogue record changed")
        require(sha(safe_file(ROOT, entry["license_evidence_path"])) == entry["license_evidence_sha256"], "Licence evidence changed")
        document = read(raw_path)
        if "results" in document:
            require(len(document["results"]) == 1, "Ambiguous raw record")
            document = document["results"][0]
        seq = document["sequence"].get("value", document["sequence"].get("sequence"))
        require(isinstance(seq, str) and set(seq) <= set("ACDEFGHIKLMNPQRSTVWY") and
                hashlib.sha256(seq.encode()).hexdigest() == entry["sequence_sha256"] and len(seq) == entry["sequence_length_aa"], "Original sequence identity differs")
        sequences[a], hashes[a], lengths[a] = seq, entry["sequence_sha256"], len(seq)
        fields(labels["source_evidence"][a], {"raw_sha256": entry["checksum_sha256"], "sequence_sha256": hashes[a]}, "Family source evidence")
    require(len(hashes) == 105 and manifest["catalogue_accessions"] == sorted(hashes), "Original catalogue coverage differs")
    queries = {}
    for item in manifest["queries"]:
        sequence = "".join(sequences[a] for a in item["components"])
        require(item["query_id"] not in queries and len(sequence) == item["length_aa"] and
                len(set(item["components"])) == item["component_count"] == len(item["spans"]), "Frozen query invalid")
        offset = 0
        for a, span in zip(item["components"], item["spans"]):
            same(span, {"accession": a, "start0": offset, "end0": offset + lengths[a],
                        "fraction": lengths[a] / len(sequence), "sequence_sha256": hashes[a]}, "Frozen span")
            offset += lengths[a]
        queries[item["query_id"]] = {**item, "joined_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}
    del sequences
    require(len(queries) == 1146, "Frozen suite size differs")
    require(labels["version"] == "iri-family-labels-v1" and labels["status"] == "frozen_before_family_scoring" and
            labels["source_challenge_manifest_sha256"] == CHALLENGE_SHA256 and
            labels["source_catalogue_manifest_sha256"] == sha(ROOT / "data/manifests/t05_expanded_catalogue.json"), "Family freeze identity differs")
    families = labels["acceptable_matches"]
    require(set(families) == set(manifest["targets"]), "Family labels do not cover targets")
    for a, accepted in families.items():
        require(a in accepted and len(accepted) == len(set(accepted)) and set(accepted) <= set(hashes), "Invalid frozen acceptable family")
    near = labels["near90_policy"]
    require(near["threshold"] == 0.9 and near["comparison"] == ">=", "Near-duplicate policy differs")
    pair_path = safe_file(ROOT, near["audit_csv_path"])
    require(sha(pair_path) == near["audit_csv_sha256"] and sha(safe_file(ROOT, near["audit_summary_path"])) == near["audit_summary_sha256"], "Frozen identity audit changed")
    with pair_path.open() as stream:
        removals = rebuild_removals(labels, hashes, lengths, csv.DictReader(stream))
    review_path = session / "audit/independent_family_label_review.json"
    review = read(review_path)
    require(review["status"] == "passed" and review["family_label_sha256"] == FAMILY_SHA256, "Independent family review differs")
    for policy in removals:
        rows = denominator_rows(list(queries.values()), hashes, families, removals, policy)
        frozen_rows = labels["removal_policies"][policy]["query_scoring_status"]
        same({r["query_id"]: r for r in frozen_rows}, {r["query_id"]: r for r in rows}, "Frozen query denominators")
        totals = denominator_totals(rows)
        same(labels["removal_policies"][policy]["denominators"], totals, "Frozen family denominator totals")
        same(review["masks_and_denominators_recomputed"][policy], totals, "Independent family review denominators")
        for query in queries.values():
            same(labels["removal_policies"][policy]["common_removed_accessions_by_group"][query["group_id"]],
                 family_masks(query, hashes, removals)["family_" + policy]["removed_accessions"], "Common all-component removal")
    return {"manifest": manifest, "queries": queries, "hashes": hashes, "lengths": lengths, "labels": labels,
            "families": families, "removals": removals, "catalogues": rebuild_catalogues(manifest, hashes),
            "equivalence": {a: sorted(b for b in hashes if hashes[b] == hashes[a]) for a in hashes},
            "evidence": {"challenge_sha256": CHALLENGE_SHA256, "family_label_sha256": FAMILY_SHA256,
                         "independent_family_review_sha256": sha(review_path), "pairwise_csv_sha256": sha(pair_path)}}


class Snapshot:
    def __init__(self, path, session, deadline, completion_sha256=None, completion_path=None):
        self.path = Path(path).resolve()
        require(self.path.is_relative_to(Path(session).resolve() / "jobs"), "Completed source escaped session jobs")
        self.identity = read(self.path / "identity.json")
        self.identity_hash = canonical_hash(self.identity)
        require(self.path.name.endswith("-" + self.identity_hash[:12]), "Job directory differs from immutable identity")
        ids = self.identity["query_ids"]
        require(ids and len(ids) == len(set(ids)), "Invalid source query IDs")
        mode = self.identity["mode"]
        require(mode in MODES | {"model", "alignment", "alignment_normalization"}, "Unsupported audit source mode")
        keys = {"query-" + q for q in ids} | {"summary"}
        if mode == "model":
            keys |= {"catalogue-" + a for a in self.identity["catalogue_accessions"]}
        elif mode == "alignment_normalization":
            keys.add("reference-selfscores")
        if completion_path is not None:
            path = Path(completion_path).resolve()
            require(path.parent == self.path, "Completion descriptor escaped its source")
        else:
            choices = [p for p in self.path.glob("completion-*.json") if completion_sha256 is None or sha(p) == completion_sha256]
            require(choices, "No matching completed source snapshot")
            path = max(choices, key=lambda p: read(p)["completed_utc"])
        self.completion_path, self.completion_sha256 = path, sha(path)
        require(completion_sha256 is None or self.completion_sha256 == completion_sha256, "Pinned completion hash differs")
        descriptor = read(path)
        require(descriptor["identity_hash"] == self.identity_hash and descriptor["expected_count"] == len(keys) and
                set(descriptor["units"]) == keys, "Completed source unit coverage differs")
        self.units, self.record_hashes = {}, {}
        self.deadline = deadline
        for key, entry in descriptor["units"].items():
            deadline.check()
            unit = verify_unit(self.path, key, entry, self.identity_hash)
            self.units[key] = unit
            self.record_hashes[key] = read(unit / "complete.json")["files"]["record.json"]["sha256"]
    def record(self, key):
        self.deadline.check()
        path = self.units[key] / "record.json"
        require(sha(path) == self.record_hashes[key], "Record changed after snapshot verification")
        return read(path)
    def evidence(self):
        return {"job": str(self.path), "identity_sha256": self.identity_hash, "completion": str(self.completion_path),
                "completion_sha256": self.completion_sha256, "verified_units": len(self.units),
                "source_code_hashes": self.identity.get("code_hashes", self.identity.get("code", {})),
                "source_model_configuration": self.identity.get("config"), "source_method": self.identity.get("method")}


def filtered(rows, included):
    included = set(included)
    chosen = [r for r in rows if r["accession"] in included]
    require(len(chosen) == len(included), "Mask contains reference absent from original complete ranking")
    return [{"accession": r["accession"], "score": score(r), "rank": index} for index, r in enumerate(chosen, 1)]


def match(rows, accepted):
    surviving = sorted(set(accepted) & {row["accession"] for row in rows})
    if not surviving:
        return {"scorable": False, "accepted_surviving": [], "rank": None, "reciprocal_rank": None,
                **{f"hit_at_{k}": None for k in KS}}
    selected = next(row for row in rows if row["accession"] in surviving)
    rank, value = selected["rank"], score(selected)
    exact = [row["rank"] for row in rows if score(row) == value]
    near = [row["rank"] for row in rows if abs(score(row) - value) <= 1e-6]
    return {"scorable": True, "accepted_surviving": surviving, "rank": rank, "best_accession": selected["accession"],
            "score": value, "reciprocal_rank": 1 / rank, "tie_best_rank": min(exact), "tie_worst_rank": max(exact),
            "near_tie_count": len(near), "near_tie_best_rank": min(near), "near_tie_worst_rank": max(near),
            **{f"hit_at_{k}": int(rank <= k) for k in KS}}


def expected_evaluation(query, ranking, equivalence, families=None):
    oracle = isinstance(ranking, dict)
    if oracle:
        require(set(ranking) == set(query["components"]), "Oracle component lists differ")
        catalogues = [{r["accession"] for r in rows} for rows in ranking.values()]
        require(all(c == catalogues[0] for c in catalogues), "Oracle components use different masks")
        catalogue_count = len(catalogues[0])
    else:
        catalogue_count = len(ranking)
    components = []
    for index, (a, span) in enumerate(zip(query["components"], query["spans"])):
        rows = ranking[a] if oracle else ranking
        kinds = {"exact": [a], "duplicate_equivalent": equivalence[a]}
        if families is not None:
            kinds["labelled"] = families[a]
        position = "single" if len(query["components"]) == 1 else "N" if index == 0 else "C" if index == len(query["components"]) - 1 else "internal"
        components.append({"accession": a, "component_index": index, "position": position,
                           "start0": span["start0"], "end0": span["end0"], "fraction": span["fraction"],
                           "fraction_bin": bin_name(span["fraction"]), **{kind: match(rows, accepted) for kind, accepted in kinds.items()}})
    metrics = {}
    for kind in kinds:
        scorable = [c[kind] for c in components if c[kind]["scorable"]]
        n, total = len(scorable), len(components)
        values = {"component_count": total, "scorable_components": n, "unscorable_components": total - n,
                  "all_components_scorable": n == total}
        for k in KS:
            hits = sum(row[f"hit_at_{k}"] for row in scorable)
            values.update({f"component_hits_at_{k}": hits, f"component_recall_at_{k}": hits / n if n else None,
                           f"any_component_hit_at_{k}": int(hits > 0) if n == total else None,
                           f"any_scorable_component_hit_at_{k}": int(hits > 0) if n else None,
                           f"all_components_hit_at_{k}": int(hits == total) if n == total else None})
        metrics[kind] = values
    return {"query_id": query["query_id"], "suite": query["suite"], "group_id": query["group_id"],
            "component_count": len(components), "catalogue_count": catalogue_count, "components": components,
            "query_metrics": metrics, "exact_recall_ceiling": {str(k): 1.0 if oracle else min(1, k / len(components)) for k in KS},
            "retrieval_budget": "k per known component" if oracle else "k total per original query",
            "returned_hits_ceiling": {str(k): min(k, catalogue_count) * (len(components) if oracle else 1) for k in KS},
            "near_tie_tolerance": 1e-6}


class OutcomeCell:
    """A frequency table of target/rank outcomes, including explicit missing labels."""
    def __init__(self):
        self.outcomes, self.parents, self.groups = Counter(), set(), set()
    def add(self, query, component, kind):
        self.outcomes[(component["accession"], component[kind]["rank"])] += 1
        self.parents.add(query["query_id"])
        self.groups.add(query["group_id"])
    def result(self):
        original = sum(self.outcomes.values())
        present = {(a, r): n for (a, r), n in self.outcomes.items() if r is not None}
        n = sum(present.values())
        targets = {a for a, _ in present}
        result = {"original_components": original, "scorable_components": n, "unscorable_components": original - n,
                  "parent_queries": len(self.parents), "distinct_groups": len(self.groups), "scorable_targets": len(targets),
                  "mean_rank": sum(rank * count for (_, rank), count in present.items()) / n if n else None,
                  "mrr": math.fsum(count / rank for (_, rank), count in present.items()) / n if n else None}
        for k in KS:
            hits = sum(count for (_, rank), count in present.items() if rank <= k)
            macro = [sum(count for (a, rank), count in present.items() if a == target and rank <= k) /
                     sum(count for (a, _), count in present.items() if a == target) for target in targets]
            result.update({f"hits_at_{k}": hits, f"recall_at_{k}": hits / n if n else None,
                           f"per_target_macro_recall_at_{k}": math.fsum(macro) / len(macro) if macro else None})
        return result


def query_frequency_summary(frequency):
    original = sum(frequency.values())
    full = {key: count for key, count in frequency.items() if key[0] == key[1]}
    n = sum(full.values())
    result = {"original_queries": original, "all_components_scorable_queries": n, "excluded_unscorable_queries": original - n}
    for position, k in enumerate(KS, 2):
        for name, test in (("any_component_hit", lambda key: key[position] > 0), ("all_components_hit", lambda key: key[position] == key[0])):
            hits = sum(count for key, count in full.items() if test(key))
            result[f"{name}_count_at_{k}"] = hits
            result[f"{name}_at_{k}"] = hits / n if n else None
    return result


class Rollup:
    def __init__(self):
        self.ids, self.groups = set(), set()
        self.components, self.budget, self.data = 0, None, {}
    def add(self, evaluation):
        require(evaluation["query_id"] not in self.ids, "Repeated query in independent rollup")
        self.ids.add(evaluation["query_id"])
        self.groups.add(evaluation["group_id"])
        self.components += evaluation["component_count"]
        require(self.budget is None or self.budget == evaluation["retrieval_budget"], "Mixed retrieval budgets")
        self.budget = evaluation["retrieval_budget"]
        for kind, qm in evaluation["query_metrics"].items():
            data = self.data.setdefault(kind, {"cells": {}, "queries": Counter(), "by_count": defaultdict(Counter)})
            frequency_key = (evaluation["component_count"], qm["scorable_components"], *(qm[f"component_hits_at_{k}"] for k in KS))
            data["queries"][frequency_key] += 1
            data["by_count"][str(evaluation["component_count"])][frequency_key] += 1
            for component in evaluation["components"]:
                for cell in (("total", ""), ("by_component_count", str(evaluation["component_count"])),
                             ("by_fraction", component["fraction_bin"]), ("by_position", component["position"]),
                             ("by_target", component["accession"])):
                    data["cells"].setdefault(cell, OutcomeCell()).add(evaluation, component, kind)
    def result(self):
        if not self.ids:
            return {"query_count": 0, "component_count": 0, "kinds": {}}
        result = {"query_count": len(self.ids), "distinct_group_count": len(self.groups), "component_count": self.components,
                  "retrieval_budget": self.budget, "kinds": {}}
        for kind, data in self.data.items():
            values = data["cells"][("total", "")].result()
            values["query_metrics"] = query_frequency_summary(data["queries"])
            for field in ("by_component_count", "by_fraction", "by_position", "by_target"):
                values[field] = {key: cell.result() for (name, key), cell in data["cells"].items() if name == field}
            for count, frequency in data["by_count"].items():
                values["by_component_count"][count]["query_metrics"] = query_frequency_summary(frequency)
            result["kinds"][kind] = values
        return result


class Sources:
    """Cache verified descriptors and hashes, never query records or arrays."""
    def __init__(self, session, deadline, ctx):
        self.session, self.deadline, self.ctx, self.cache = session.resolve(), deadline, ctx, {}
    def open(self, path, pinned=None, completion_path=None):
        path = Path(path).resolve()
        key = (str(path), pinned)
        if key not in self.cache:
            snapshot = Snapshot(path, self.session, self.deadline, pinned, completion_path)
            self.cache[key] = snapshot
            identity = snapshot.identity
            fields(identity, {"challenge_sha256": CHALLENGE_SHA256,
                              "catalogue_accessions": sorted(self.ctx["hashes"])}, "Source suite identity")
            require(set(identity["query_ids"]) == set(self.ctx["queries"]), "This audit requires complete 1146-query source coverage")
            if identity["mode"] not in MODES:
                same(identity["catalogue_sequence_sha256"], self.ctx["hashes"], "Source sequence hashes")
            if identity["mode"] == "alignment_normalization":
                raw = self.open(identity["source_raw_job"], identity["source_completion_sha256"])
                require(raw.identity["mode"] == "alignment" and raw.identity_hash == identity["source_identity_hash"], "Normalized source chain differs")
                same(identity["query_ids"], raw.identity["query_ids"], "Normalized source query coverage")
                same(identity["catalogue_masks"], self.ctx["catalogues"], "Normalized source catalogue masks")
        snapshot = self.cache[key]
        if completion_path is not None:
            require(Path(completion_path).resolve() == snapshot.completion_path, "Conflicting pinned completion path")
        return snapshot
    def unchanged(self):
        for snapshot in self.cache.values():
            self.deadline.check()
            require(canonical_hash(read(snapshot.path / "identity.json")) == snapshot.identity_hash and
                    sha(snapshot.completion_path) == snapshot.completion_sha256, "Source identity/completion changed during audit")


def source_arms(snapshot, record, query, accessions):
    """Validate complete original rankings before any mask or fusion operation."""
    same(record["query"], query, "Source query metadata")
    mode = snapshot.identity["mode"]
    if mode == "model":
        expected = {"whole", "isolated_mixture", "contextual_component", "isolated_component"}
        require(set(record["rankings"]) == expected, "Original model arms differ")
        result = {a: (record["rankings"][a], "/rankings/" + a) for a in sorted(expected)}
    elif mode == "alignment":
        require(set(record["rankings"]) == {"local_alignment"}, "Original alignment arms differ")
        result = {"local_alignment": (record["rankings"]["local_alignment"], "/rankings/local_alignment")}
    elif mode == "alignment_normalization":
        result = {"normalized_alignment": (record["rankings"]["accession105"], "/rankings/accession105")}
    else:
        raise ValueError("Unsupported original ranking source")
    for name, (rows, _) in result.items():
        oracle = mode == "model" and name in ("contextual_component", "isolated_component")
        require(isinstance(rows, dict) == oracle, "Unexpected source ranking shape")
        if oracle:
            require(set(rows) == set(query["components"]), "Original oracle component coverage differs")
        for values in rows.values() if oracle else [rows]:
            validate_ranking(values, accessions)
    return result


def expected_recipes(ctx, query):
    masks = {name: {"included_accessions": accessions, "label_kind": "exact_and_duplicate",
                    "source": "frozen_challenge.catalogue_subsets"} for name, accessions in ctx["catalogues"]["masks"].items()}
    for name, mask in family_masks(query, ctx["hashes"], ctx["removals"]).items():
        masks[name] = {**mask, "label_kind": "labelled_family", "source": "family_labels.removal_policies." + mask["policy"]}
    return {name: {**mask, "included_accessions_sha256": canonical_hash(mask["included_accessions"]),
                  "ordered_ranking_recipe": "Read immutable source ranking pointer; retain included accessions in existing score/accession order; renumber from 1; recompute exact/near ties at tolerance 1e-6."}
            for name, mask in masks.items()}


def mask_rankings(raw, included):
    return {a: filtered(rows, included) for a, rows in raw.items()} if isinstance(raw, dict) else filtered(raw, included)


def removal_diagnostics(raw, masked, query):
    components = {}
    for a in query["components"]:
        original, surviving = (raw[a], masked[a]) if isinstance(raw, dict) else (raw, masked)
        exact = next(row for row in original if row["accession"] == a)
        original_ranks = {row["accession"]: row["rank"] for row in original}
        top = [{"accession": row["accession"], "masked_rank": row["rank"],
                "source_rank": original_ranks[row["accession"]], "score": score(row)} for row in surviving[:10]]
        components[a] = {"removed_exact_source_rank": exact["rank"], "removed_exact_score": score(exact),
                         "surviving_top10": top,
                         "best_survivor_minus_removed_exact_score": top[0]["score"] - score(exact) if top else None}
    return {"exact_accession_accuracy_defined": False,
            "interpretation": "Exact references are absent; surviving identities/scores are diagnostics, not functional correctness.",
            "components": components}


def check_family_query(record, raw_arms, query, ctx):
    same(record["query"], query, "Derived query metadata")
    recipes = expected_recipes(ctx, query)
    same(record["ranking_recipes"], recipes, "Independently rebuilt ranking masks/recipes")
    same(record["source_ranking_pointers"], {name: pointer for name, (_, pointer) in raw_arms.items()}, "Original ranking pointers")
    require(set(record["evaluations"]) == set(recipes), "Derived condition coverage differs")
    evaluations, diagnostics = {}, {}
    for condition, recipe in recipes.items():
        require(set(record["evaluations"][condition]) == set(raw_arms), "Derived arm coverage differs")
        evaluations[condition] = {}
        for name, (raw, _) in raw_arms.items():
            masked = mask_rankings(raw, recipe["included_accessions"])
            evaluation = expected_evaluation(query, masked, ctx["equivalence"],
                                             ctx["families"] if recipe["label_kind"] == "labelled_family" else None)
            same(record["evaluations"][condition][name], evaluation, condition + "/" + name + "/evaluation")
            evaluations[condition][name] = evaluation
            if condition == "family_exact_removed":
                diagnostics[name] = removal_diagnostics(raw, masked, query)
    same(record["exact_reference_removal_diagnostics"], diagnostics, "Removed-exact diagnostics")
    return evaluations


def audit_family(snapshot, sources):
    ctx, identity = sources.ctx, snapshot.identity
    source = sources.open(identity["source_job"], identity["source_completion_sha256"], identity["source_completion_path"])
    fields(identity, {"source_identity_hash": source.identity_hash, "source_mode": source.identity["mode"],
                      "source_model_configuration": source.identity.get("config"), "query_ids": source.identity["query_ids"],
                      "catalogue_sequence_sha256": ctx["hashes"], "catalogue_masks": ctx["catalogues"], "partial_source": False,
                      "family_labels_path": str(sources.session / "inputs/family_labels_v1.json"),
                      "family_label_validation": {"label_sha256": FAMILY_SHA256,
                          "independent_review_sha256": ctx["evidence"]["independent_family_review_sha256"],
                          "independent_review_path": str(sources.session / "audit/independent_family_label_review.json")}}, "Family analysis identity")
    bundles = {}
    for qid in identity["query_ids"]:
        sources.deadline.check()
        query, key = ctx["queries"][qid], "query-" + qid
        raw, record = source.record(key), snapshot.record(key)
        same(record["source_record"], {"unit_path": str(source.units[key]), "record_sha256": source.record_hashes[key],
             "descriptor_sha256": sha(source.units[key] / "complete.json"), "source_identity_hash": source.identity_hash,
             "source_completion_sha256": source.completion_sha256}, "Derived record source pins")
        arms = source_arms(source, raw, query, sorted(ctx["hashes"]))
        values = check_family_query(record, arms, query, ctx)
        for condition, arms in values.items():
            for arm, evaluation in arms.items():
                bundle = bundles.setdefault((condition, arm), {"all": Rollup(), "combined": Rollup(), "by_suite": {}})
                bundle["all"].add(evaluation)
                if query["component_count"] > 1:
                    bundle["combined"].add(evaluation)
                bundle["by_suite"].setdefault(query["suite"], Rollup()).add(evaluation)
    summary = {}
    for (condition, arm), bundle in bundles.items():
        sources.deadline.check()
        summary.setdefault(condition, {})[arm] = {"all": bundle["all"].result(), "combined": bundle["combined"].result(),
                                                  "by_suite": {name: value.result() for name, value in bundle["by_suite"].items()}}
    denominator_evidence = {policy: denominator_totals(denominator_rows(
        [ctx["queries"][qid] for qid in identity["query_ids"]], ctx["hashes"], ctx["families"], ctx["removals"], policy))
        for policy in ctx["removals"]}
    fields(snapshot.record("summary"), {"summary": summary, "query_count": len(identity["query_ids"]),
        "catalogue_masks": ctx["catalogues"], "family_frozen_denominators_for_source_subset": denominator_evidence,
        "all_family_denominators_match_freeze": True, "unique_catalogue_conditions": len(ctx["catalogues"]["masks"]),
        "prepared_subset_design_count": len(ctx["manifest"]["catalogue_subsets"]), "full101_shared_endpoint_scored_once": True}, "Independent family summaries")
    return {"status": "passed", **snapshot.evidence(), "original_source": source.evidence(),
            "queries_checked": len(identity["query_ids"]), "conditions_checked": sorted(summary),
            "arms_checked": sorted({arm for _, arm in bundles}), "family_denominators": denominator_evidence,
            "checks": ["independent_common_all_component_masks", "independent25_50_101_catalogue_masks",
                       "all_component_ranks_ties_hits", "absent_label_nulls", "all_summary_strata_and_support",
                       "all_scorable_query_denominators_including_zero", "removed_exact_diagnostics"]}


def fused_ranking(component_rankings, method, accessions):
    require(method in ("max_score", "rrf") and component_rankings, "Unsupported or empty fixed fusion")
    contributions = {a: [] for a in accessions}
    for rows in component_rankings:
        validate_ranking(rows, accessions)
        for row in rows:
            require("cosine_similarity" in row, "Oracle source lacks explicit same-model cosine")
            contributions[row["accession"]].append(1 / (60 + row["rank"]) if method == "rrf" else score(row))
    values = {a: math.fsum(v) if method == "rrf" else max(v) for a, v in contributions.items()}
    ties = Counter(values.values())
    return [{"accession": a, "rank": rank, "score": values[a], "fusion_score": values[a], "exact_score_tie_count": ties[values[a]]}
            for rank, a in enumerate(sorted(values, key=lambda a: (-values[a], a)), 1)]


def audit_oracle(snapshot, sources):
    ctx, identity = sources.ctx, snapshot.identity
    design_path = sources.session / "inputs/oracle_total_budget_design_v1.json"
    require(sha(design_path) == identity["design_sha256"] == ORACLE_DESIGN_SHA256, "Frozen oracle design changed")
    design = read(design_path)
    fields(design, {"version": "iri-oracle-total-budget-v1", "status": "frozen_before_oracle_fusion_scoring",
                   "challenge_sha256": CHALLENGE_SHA256, "arms": ["contextual_component", "isolated_component"],
                   "fusions": ["max_score", "rrf"], "rrf_constant": 60, "budgets": [1, 5, 10], "known_boundaries": True}, "Oracle design")
    source = sources.open(identity["source_path"], identity["source_completion_sha256"])
    require(source.identity["mode"] == "model", "Oracle source is not a model job")
    fields(identity, {"source_identity": source.identity_hash, "query_ids": source.identity["query_ids"],
                      "coverage": {"full_suite": True, "observed_queries": len(ctx["queries"]), "expected_suite_queries": len(ctx["queries"])},
                      "information_budget": "Known complete component boundaries; k total unique accessions per query"}, "Oracle identity")
    bundles = {arm + "__" + method: {"clean": Rollup(), "combined": Rollup()} for arm in design["arms"] for method in design["fusions"]}
    for qid in identity["query_ids"]:
        query, key = ctx["queries"][qid], "query-" + qid
        record = snapshot.record(key)
        original = source_arms(source, source.record(key), query, sorted(ctx["hashes"]))
        fields(record, {"query": query, "known_component_count": query["component_count"], "boundaries_supplied": True}, "Oracle query information")
        require(set(record["rankings"]) == set(record["evaluations"]) == set(bundles), "Oracle output arm coverage differs")
        for arm in design["arms"]:
            component_rows = [original[arm][0][a] for a in query["components"]]
            for method in design["fusions"]:
                name = arm + "__" + method
                expected = fused_ranking(component_rows, method, sorted(ctx["hashes"]))
                same(record["rankings"][name], expected, "Independent total-budget fusion")
                evaluation = expected_evaluation(query, expected, ctx["equivalence"])
                same(record["evaluations"][name], evaluation, "Independent total-budget evaluation")
                bundles[name]["clean" if query["component_count"] == 1 else "combined"].add(evaluation)
    same(snapshot.record("summary"), {"summaries": {name: {scope: value.result() for scope, value in scopes.items()}
         for name, scopes in bundles.items()}, "design": design}, "Independent oracle summary strata")
    return {"status": "passed", **snapshot.evidence(), "original_source": source.evidence(),
            "queries_checked": len(identity["query_ids"]), "arms_checked": sorted(bundles), "design_sha256": sha(design_path),
            "checks": ["known_boundary_information_label", "independent_maxcos_and_rrf60", "single_totalK_unique_accession_list",
                       "all_exact_duplicate_ranks_hits_ties", "all_clean_combined_summary_strata_and_support"]}


def atomic_report(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix="." + path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def run(args):
    session = args.session.resolve()
    output = (args.output or session / "audit/independent_extension_audit.json").resolve()
    require(output.is_relative_to(session / "audit"), "Audit output must remain in session audit directory")
    require(not output.exists(), "Preserve previous audit evidence; choose a new output filename")
    started = time.monotonic()
    report = {"status": "running", "version": "iri-independent-extension-audit-v1", "jobs": [],
              "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "timeout_seconds": args.timeout_seconds,
              "audit_code_sha256": sha(Path(__file__)), "independent_helper_sha256": sha(Path(__file__).with_name("session_independent_audit.py")),
              "scope": "Saved-result derivation audit; no inference, alignment, biological-label validation or representation recomputation"}
    sources = None
    def interrupted(signum, frame):
        raise InterruptedError("Audit interrupted by signal " + str(signum))
    previous = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        deadline = Deadline(args.timeout_seconds, session)
        deadline.check()
        ctx = context(session, deadline)
        report["input_evidence"] = ctx["evidence"]
        sources = Sources(session, deadline, ctx)
        jobs = args.job or [session / "jobs" / name for name in DEFAULT_JOBS]
        require(len(jobs) == args.expected_jobs and len({Path(p).resolve() for p in jobs}) == len(jobs), "Completed audit job coverage/count differs")
        report["requested_jobs"] = [str(Path(p).resolve()) for p in jobs]
        for path in jobs:
            deadline.check()
            report["active_job"] = str(Path(path).resolve())
            snapshot = sources.open(path)
            require(snapshot.identity["mode"] in MODES, "Requested job is not a supported extension")
            audit = audit_family if snapshot.identity["mode"] == "catalogue_and_family_analysis" else audit_oracle
            report["jobs"].append(audit(snapshot, sources))
            atomic_report(output, report)
            print(json.dumps({"status": "job_passed", "job": str(snapshot.path), "jobs_checked": len(report["jobs"])}), flush=True)
        sources.unchanged()
        report.pop("active_job", None)
        report["status"] = "passed"
    except (Exception, KeyboardInterrupt) as error:
        report["status"] = "failed"
        report["failure"] = {"type": type(error).__name__, "message": str(error)}
    finally:
        if sources is not None:
            report["verified_snapshots"] = [source.evidence() for source in sources.cache.values()]
        report["elapsed_seconds"] = time.monotonic() - started
        report["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        atomic_report(output, report)
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    print(json.dumps({"status": report["status"], "report": str(output), "jobs_checked": len(report["jobs"])}), flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", type=Path, required=True)
    parser.add_argument("--job", type=Path, action="append")
    parser.add_argument("--expected-jobs", type=int, default=8)
    parser.add_argument("--timeout-seconds", type=float, default=900)
    parser.add_argument("--output", type=Path)
    raise SystemExit(run(parser.parse_args()))
