#!/usr/bin/env python3
"""Independent held-out segmentation audit; no production metric/model imports.

Rebuilds development selection, frozen masks, fusion, interval diagnostics and
all metrics. Fresh fragment DP candidate lists were not persisted: their global
optimality cannot be re-proved here. Saved trace eligibility and available raw
source rankings are checked separately and the limitation is retained in reports.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter, defaultdict
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import signal
import struct
import time
import zipfile

import session_extension_audit as base

ROOT = Path(__file__).resolve().parents[1]
DESIGN_SHA = "a91092feb094dfe46bb9f157b90758933bce0ff7fe58d0f19634de5188e435fe"
TRANSFER_SHA = "f2e74790784074bf2326056787a12bd62674e5eebe73012f352cadb85ec7c209"
SELECTION_SHA = "833fe8774c38eddb0c476de54800ff3d6b39b6b118b699cfec370bf2f3e5d1e7"
DEFAULT_JOBS = ("segmentation-test-8m-infer-v1-42b54ea33217", "family-test-8m-infer-v1-d0424de04de8")
SEGMENT_MODES = {"automatic_spans_v1", "automatic_segment_retrieval_v1", "family_automatic_spans_v1", "family_automatic_retrieval_v1"}
require, same, fields = base.require, base.same, base.fields
read, sha, canonical_hash = base.read, base.sha, base.canonical_hash


def expected_keys(identity):
    ids = identity["query_ids"]
    require(ids and len(ids) == len(set(ids)), "Invalid segmentation query coverage")
    mode, configs = identity["mode"], identity["configs"]
    if mode == "automatic_spans_v1":
        require(isinstance(configs, dict) and configs, "Missing prediction configurations")
        return {qid + "--" + key for qid in ids for key in configs}
    names = [cfg["configuration_id"] for cfg in configs]
    require(names and len(names) == len(set(names)), "Invalid retrieval configurations")
    if mode == "automatic_segment_retrieval_v1":
        return {qid + "--" + key for qid in ids for key in names} | {"summary"}
    require(mode in ("family_automatic_spans_v1", "family_automatic_retrieval_v1"), "Unknown segmentation mode")
    require(len(identity["policies"]) == 2 and set(identity["policies"]) == {"exact_removed", "near90_removed"}, "Family policy coverage differs")
    keys = {qid + "--" + policy + "--" + key for qid in ids for policy in identity["policies"] for key in names}
    return keys | {"summary"} if mode == "family_automatic_retrieval_v1" else keys


class SegmentSnapshot(base.Snapshot):
    """Same independent checksum contract with segmentation-specific unit keys."""
    def __init__(self, path, session, deadline, pinned=None):
        self.path = Path(path).resolve()
        require(self.path.is_relative_to(session.resolve() / "jobs"), "Segmentation source escaped jobs")
        self.identity = read(self.path / "identity.json")
        self.identity_hash = canonical_hash(self.identity)
        require(self.path.name.endswith("-" + self.identity_hash[:12]), "Segmentation job identity differs from directory")
        keys = expected_keys(self.identity)
        choices = [p for p in self.path.glob("completion-*.json") if pinned is None or sha(p) == pinned]
        require(choices, "Missing pinned segmentation completion")
        self.completion_path = max(choices, key=lambda p: read(p)["completed_utc"])
        self.completion_sha256 = sha(self.completion_path)
        descriptor = read(self.completion_path)
        require(descriptor["identity_hash"] == self.identity_hash and descriptor["expected_count"] == len(keys) and
                set(descriptor["units"]) == keys, "Segmentation complete-unit coverage differs")
        self.deadline, self.units, self.record_hashes = deadline, {}, {}
        for key, index in descriptor["units"].items():
            deadline.check()
            path = base.verify_unit(self.path, key, index, self.identity_hash)
            self.units[key] = path
            self.record_hashes[key] = read(path / "complete.json")["files"]["record.json"]["sha256"]


class Sources:
    def __init__(self, session, deadline, ctx):
        self.session, self.deadline, self.ctx = session.resolve(), deadline, ctx
        self.index, self.snapshots = {}, {}
        for path in sorted((session / "jobs").glob("*/identity.json")):
            deadline.check()
            identity = read(path)
            digest = canonical_hash(identity)
            self.index.setdefault(digest, []).append(path.parent.resolve())
    def open(self, digest, pinned=None, path=None):
        require(digest in self.index, "Pinned source identity is missing")
        candidates = [value for value in self.index[digest] if path is None or value == Path(path).resolve()]
        if pinned is not None:
            candidates = [value for value in candidates if any(sha(p) == pinned for p in value.glob("completion-*.json"))]
        require(len(candidates) == 1, "Ambiguous or missing pinned source job")
        path = candidates[0]
        key = (digest, pinned, str(path))
        if key not in self.snapshots:
            mode = read(path / "identity.json")["mode"]
            cls = SegmentSnapshot if mode in SEGMENT_MODES else base.Snapshot
            snapshot = cls(path, self.session, self.deadline, pinned)
            require(snapshot.identity_hash == digest, "Source identity changed")
            fields(snapshot.identity, {"challenge_sha256": base.CHALLENGE_SHA256,
                   "catalogue_accessions": sorted(self.ctx["hashes"]), "catalogue_sequence_sha256": self.ctx["hashes"]}, "Segmentation source catalogue")
            require(set(snapshot.identity["query_ids"]) <= set(self.ctx["queries"]), "Unknown source query")
            self.snapshots[key] = snapshot
        return self.snapshots[key]
    def by_path(self, path):
        path = Path(path).resolve()
        digest = canonical_hash(read(path / "identity.json"))
        require(path in self.index.get(digest, []), "Source path/identity differs")
        return self.open(digest, path=path)
    def unchanged(self):
        for snapshot in self.snapshots.values():
            self.deadline.check()
            require(canonical_hash(read(snapshot.path / "identity.json")) == snapshot.identity_hash and
                    sha(snapshot.completion_path) == snapshot.completion_sha256, "Pinned source changed during audit")


def geometry(config):
    return {key: value for key, value in config.items() if key not in ("configuration_id", "fusion", "rrf_constant")}


def independently_ranked(rows, included, *, diagnostic_ties=False):
    allowed = set(included)
    chosen = [row for row in rows if row["accession"] in allowed]
    require(len(chosen) == len(allowed), "Missing candidate before masked reranking")
    ties = Counter(base.score(row) for row in chosen)
    result = [{**row, "rank": i, "exact_score_tie_count": ties[base.score(row)]} for i, row in enumerate(chosen, 1)]
    if diagnostic_ties:
        for row in result:
            exact = [r["rank"] for r in result if base.score(r) == base.score(row)]
            near = [r["rank"] for r in result if abs(base.score(r) - base.score(row)) <= 1e-6]
            row.update(tie_best_rank=min(exact), tie_worst_rank=max(exact), near_tie_count=len(near),
                       near_tie_best_rank=min(near), near_tie_worst_rank=max(near), near_tie_tolerance=1e-6)
            if "near_tie_count_1e6" in row:
                row["near_tie_count_1e6"] = len(near)
    return result


def interval_metrics(segments, targets, length):
    ranges = []
    for value in segments:
        a, b = value["start0"], value["end0"]
        require(type(a) is int and type(b) is int and 0 <= a < b <= length, "Invalid predicted interval")
        ranges.append((a, b))
    rows = []
    for target in targets:
        a, b = target["start0"], target["end0"]
        require(0 <= a < b <= length, "Invalid reference interval")
        # At most1011 positions: direct set union is independent of production's
        # interval sweep and avoids double-counting overlapping windows.
        covered = {position for x, y in ranges for position in range(max(a, x), min(b, y))}
        ious = [len(set(range(a, b)) & set(range(x, y))) / len(set(range(a, b)) | set(range(x, y))) for x, y in ranges]
        rows.append({"accession": target.get("accession"), "start0": a, "end0": b,
                     "best_interval_iou": max(ious, default=0), "covered_fraction_by_union": len(covered) / (b - a),
                     "exact_boundary_match": (a, b) in ranges})
    return {"components": rows, "predicted_segment_count": len(ranges), "target_component_count": len(rows), "evaluation_only": True}


def expected_windows(length, config):
    width, stride = config["width"], config["stride"]
    require(0 < stride <= width, "Invalid frozen window geometry")
    starts = sorted(set([0] if length <= width else list(range(0, length - width + 1, stride)) + [length - width]))
    segments = [{"start0": a, "end0": min(length, a + width)} for a in starts]
    return {"segments": segments, "segment_count": len(segments), "query_length": length,
            "total_queried_residues": sum(row["end0"] - row["start0"] for row in segments),
            "width": width, "stride": stride, "cap_bound": False, "method": "fixed_overlapping_windows", "boundaries_supplied": False}


def union_result(rankings, k, query, families=None, included=None):
    occurrences = defaultdict(list)
    for segment, rows in enumerate(rankings):
        for row in rows[:k]:
            occurrences[row["accession"]].append({"segment_index": segment, "rank": row["rank"]})
    names, n = sorted(occurrences), sum(len(values) for values in occurrences.values())
    result = {"accessions": names, "occurrences": dict(occurrences), "per_segment_k": k, "segment_count": len(rankings),
              "returned_hit_occurrences": n, "unique_returned_hits": len(names), "duplicate_returned_hits": n - len(names),
              "retrieval_budget": "k per predicted segment"}
    present = set(names)
    if families is None:
        hits = len(present & set(query["components"]))
        result["evaluation"] = {"components": len(query["components"]), "component_hits": hits,
            "component_recall": hits / len(query["components"]), "any_components": int(hits > 0),
            "all_components": int(hits == len(query["components"])), "extraneous_to_exact_component_labels": len(present - set(query["components"])),
            "extraneous_note": "Outside exact target labels, not necessarily a false functional match"}
    else:
        remaining, components = set(included), []
        for a in query["components"]:
            accepted = remaining & set(families[a])
            components.append({"accession": a, "scorable": bool(accepted), "accepted_surviving": sorted(accepted),
                "hit": int(bool(accepted & present)) if accepted else None, "exact_reference_present": a in remaining,
                "exact_hit": int(a in present) if a in remaining else None})
        scorable = [c for c in components if c["scorable"]]
        hits, count = sum(c["hit"] for c in scorable), len(scorable)
        all_scorable = count == len(components)
        result["evaluation"] = {"components": components, "original_components": len(components), "scorable_components": count,
            "unscorable_components": len(components) - count, "component_hits": hits, "component_recall": hits / count if count else None,
            "all_components_scorable": all_scorable, "any_component_hit": int(hits > 0) if all_scorable else None,
            "all_components_hit": int(hits == len(components)) if all_scorable else None,
            "unlabelled_returned_hits": len(present - {a for c in components for a in c["accepted_surviving"]}),
            "unlabelled_note": "Outside accepted source-family sets, not independently verified functional negatives."}
    return result


class UnionRollup:
    def __init__(self):
        self.ids, self.targets, self.queries, self.volumes = set(), Counter(), Counter(), Counter()
    def add(self, qid, value):
        require(qid not in self.ids, "Repeated union query")
        self.ids.add(qid)
        evaluation = value["evaluation"]
        for component in evaluation["components"]:
            self.targets[(component["accession"], component["hit"])] += 1
        self.queries[(evaluation["all_components_scorable"], evaluation["any_component_hit"], evaluation["all_components_hit"])] += 1
        self.volumes[(value["returned_hit_occurrences"], value["unique_returned_hits"], value["duplicate_returned_hits"])] += 1
    def result(self):
        original = sum(self.targets.values())
        scorable = {(a, h): n for (a, h), n in self.targets.items() if h is not None}
        count, hits = sum(scorable.values()), sum(h * n for (a, h), n in scorable.items())
        targets = {a for a, _ in scorable}
        full = sum(n for (eligible, _, _), n in self.queries.items() if eligible)
        volume = [sum(key[i] * n for key, n in self.volumes.items()) for i in range(3)]
        macro = [sum(h * n for (a, h), n in scorable.items() if a == target) / sum(n for (a, _), n in scorable.items() if a == target) for target in targets]
        return {"original_queries": len(self.ids), "original_components": original, "scorable_components": count,
            "unscorable_components": original - count, "component_hits": hits, "component_recall": hits / count if count else None,
            "per_target_macro_component_recall": math.fsum(macro) / len(macro) if macro else None,
            "all_components_scorable_queries": full, "excluded_all_component_queries": len(self.ids) - full,
            "any_component_hit": sum(n * anyhit for (eligible, anyhit, _), n in self.queries.items() if eligible) / full if full else None,
            "all_components_hit": sum(n * allhit for (eligible, _, allhit), n in self.queries.items() if eligible) / full if full else None,
            "returned_hit_occurrences": volume[0], "unique_hit_occurrences": volume[1], "duplicate_hit_occurrences": volume[2],
            "mean_returned_hit_occurrences": volume[0] / len(self.ids) if self.ids else None, "retrieval_budget": "k per predicted segment"}


def frozen_designs(session, ctx):
    primary_path = session / "inputs/segmentation_design_v1.json"
    transfer_path = session / "inputs/family_segmentation_design_v1.json"
    selection_path = session / "inputs/segmentation_selection_v1.json"
    require(sha(primary_path) == DESIGN_SHA and sha(transfer_path) == TRANSFER_SHA and sha(selection_path) == SELECTION_SHA,
            "Frozen segmentation design/selection changed")
    design, transfer, selection = read(primary_path), read(transfer_path), read(selection_path)
    fields(design, {"version": "iri-segmentation-design-v1", "status": "design_frozen_no_predictions_scored",
                   "source_challenge_sha256": base.CHALLENGE_SHA256}, "Primary design")
    fields(transfer, {"version": "family-segmentation-transfer-v1", "status": "prospective_frozen_before_transfer_scoring",
                     "primary_design_sha256": DESIGN_SHA, "family_labels_sha256": base.FAMILY_SHA256,
                     "required_combined_test_queries": 950, "retuning_allowed": False}, "Family transfer design")
    fields(selection, {"design_sha256": DESIGN_SHA}, "Frozen development selection")
    partitions, groups, grouped_ids = {}, {}, defaultdict(list)
    for qid, query in ctx["queries"].items():
        group = query["group_id"]
        digest = hashlib.sha256((design["partition"]["seed"] + "|" + group).encode()).hexdigest()
        part = "dev" if int(digest, 16) < (1 << 256) // 5 else "test"
        partitions[qid] = part
        item = {"components": sorted(query["components"]), "partition": part, "sha256": digest}
        require(group not in groups or groups[group] == item, "A group crosses component sets/partitions")
        groups[group] = item
        grouped_ids[group].append(qid)
    same(design["partition"]["queries"], partitions, "Independent grouped split")
    same(design["partition"]["groups"], groups, "Independent group split evidence")
    same({g: sorted(ids) for g, ids in design["group_queries"].items()}, {g: sorted(ids) for g, ids in grouped_ids.items()}, "Group membership")
    require(len(design["selection_query_ids"]) == 180 and len(design["test_query_ids"]) == 950 and
            len(design["clean_sanity_query_ids"]) == 16, "Frozen panel counts differ")
    require(set(design["selection_query_ids"]) == {q for q, p in partitions.items() if p == "dev" and ctx["queries"][q]["component_count"] > 1} and
            set(design["test_query_ids"]) == {q for q, p in partitions.items() if p == "test" and ctx["queries"][q]["component_count"] > 1} and
            set(design["clean_sanity_query_ids"]) == {q for q, value in ctx["queries"].items() if value["component_count"] == 1}, "Frozen panel membership differs")
    require(not ({ctx["queries"][q]["group_id"] for q in design["selection_query_ids"]} &
                 {ctx["queries"][q]["group_id"] for q in design["test_query_ids"]}), "Development/test group leakage")
    for relative, digest in design["implementation_sha256"].items():
        require(sha(base.safe_file(ROOT, relative)) == digest, "Frozen segmentation implementation changed")
    return design, transfer, selection


def selection_from_outcomes(candidates, design, ctx):
    require(set(candidates) == {c["configuration_id"] for c in design["candidate_configurations"]}, "Missing development candidate")
    scored = []
    for name, rows in candidates.items():
        require(set(rows) == set(design["selection_query_ids"]), "Selection includes missing/test queries")
        grouped = defaultdict(list)
        for qid, values in rows.items():
            require(design["partition"]["queries"][qid] == "dev", "Test query influenced development selection")
            require(len(values) == 3 and math.isfinite(values[0]) and 0 <= values[0] <= 1 and values[1] in (0, 1) and
                    math.isfinite(values[2]) and values[2] >= 0, "Invalid independent development outcome")
            grouped[ctx["queries"][qid]["group_id"]].append(values)
        group_metrics = [[math.fsum(row[i] for row in values) / len(values) for values in grouped.values()] for i in (0, 1)]
        scored.append({"configuration_id": name,
            "group_macro_component_recall_at_5": math.fsum(group_metrics[0]) / len(group_metrics[0]),
            "group_macro_all_components_hit_at_5": math.fsum(group_metrics[1]) / len(group_metrics[1]),
            "mean_total_queried_residues": math.fsum(row[2] for row in rows.values()) / len(rows),
            "dev_query_count": len(rows), "dev_group_count": len(grouped)})
    scored.sort(key=lambda row: (-row["group_macro_component_recall_at_5"], -row["group_macro_all_components_hit_at_5"],
                                 row["mean_total_queried_residues"], row["configuration_id"]))
    families = {c["configuration_id"]: c["segmenter"] for c in design["candidate_configurations"]}
    winners = {}
    for row in scored:
        winners.setdefault(families[row["configuration_id"]], row["configuration_id"])
    return {"selected_configuration_id": scored[0]["configuration_id"], "selected_by_segmenter": winners,
            "development_scores": scored, "selection_scope": "combined dev queries only; test remains untouched",
            "must_freeze_before_test_scoring": True}


def trusted_target_sequences(ctx):
    result = {}
    catalogue = read(ROOT / "data/manifests/t05_expanded_catalogue.json")
    for entry in catalogue["records"]:
        if entry["accession"] not in ctx["manifest"]["targets"]:
            continue
        path = base.safe_file(ROOT, entry["raw_path"])
        require(sha(path) == entry["checksum_sha256"], "Target source changed")
        document = read(path)
        if "results" in document:
            document = document["results"][0]
        sequence = document["sequence"].get("value", document["sequence"].get("sequence"))
        require(hashlib.sha256(sequence.encode()).hexdigest() == ctx["hashes"][entry["accession"]], "Target sequence digest differs")
        result[entry["accession"]] = sequence
    return result


def eligible_candidates(rows, interval, config):
    start, end = interval
    eligible = []
    for row in rows:
        span = row.get("query_span")
        if span is None:
            continue
        a, b = span
        require(type(a) is int and type(b) is int and 0 <= a < b <= end - start, "Invalid alignment seed coordinates")
        coverage, identity = row["reference_coverage"], row["identity_fraction"]
        require(math.isfinite(coverage) and math.isfinite(identity) and 0 <= coverage <= 1 and 0 <= identity <= 1, "Invalid seed quality metadata")
        if b - a >= config["minimum_span"] and base.score(row) >= config["minimum_score"] and coverage >= config["minimum_reference_coverage"] and identity >= config["minimum_identity_fraction"]:
            eligible.append({"start0": start + a, "end0": start + b, "source_interval": list(interval), "seed_hit": row})
    return sorted(eligible, key=lambda row: (-base.score(row["seed_hit"]), row["seed_hit"]["accession"], row["start0"], row["end0"]))


def audit_alignment_trace(prediction, sequence, config, included, raw_source, ctx, counters):
    params, length = config["parameters"], len(sequence)
    fields(prediction, {"method": "iterative_alignment_candidates", "config": params, "query_length": length,
                       "boundaries_supplied": False}, "Alignment prediction inputs")
    attempts = prediction["interval_attempts"]
    intervals = [(row["start0"], row["end0"]) for row in attempts]
    require(len(set(intervals)) == len(intervals) == prediction["intervals_aligned"], "Duplicate/missing attempted interval")
    require(all(0 <= a < b <= length and b - a >= params["minimum_span"] for a, b in intervals), "Invalid attempted interval")
    require(prediction["total_queried_residues"] == sum(b - a for a, b in intervals) and
            prediction["alignment_comparisons"] == len(included) * len(intervals), "Alignment scan accounting differs")
    known_queries = {q["joined_sequence_sha256"]: qid for qid, q in ctx["queries"].items()}
    known_candidates, accepted_by_interval = {}, {}
    segments = prediction["segments"]
    fallback = prediction["whole_query_fallback_used"]
    accepted = [] if fallback else sorted(segments, key=lambda row: row["selection_order"])
    require([s["selection_order"] for s in accepted] == list(range(len(accepted))) and
            prediction["accepted_alignment_segment_count"] == len(accepted) <= params["maximum_segments"], "Accepted seed order/cap differs")
    for segment in accepted:
        interval = tuple(segment["source_interval"])
        require(interval in intervals and interval not in accepted_by_interval, "Seed refers to an unknown/reused source interval")
        require(segment["seed_hit"]["accession"] in included, "Removed reference informed prediction")
        candidate = eligible_candidates([segment["seed_hit"]], interval, params)
        require(len(candidate) == 1, "Chosen seed violates frozen eligibility")
        fields(segment, candidate[0], "Seed coordinate translation")
        accepted_by_interval[interval] = segment
    for attempt, interval in zip(attempts, intervals):
        require(attempt["catalogue_count"] == len(included) and 0 <= attempt["eligible_reference_count"] <= len(included), "Attempt used a different catalogue")
        fragment_hash = hashlib.sha256(sequence[interval[0]:interval[1]].encode()).hexdigest()
        qid = known_queries.get(fragment_hash)
        if qid is not None:
            raw = raw_source.record("query-" + qid)
            same(raw["query"], ctx["queries"][qid], "Cached prediction source sequence metadata")
            rows = raw["rankings"]["local_alignment"]
            base.validate_ranking(rows, sorted(ctx["hashes"]))
            candidates = eligible_candidates(independently_ranked(rows, included, diagnostic_ties=set(included) != set(ctx["hashes"])), interval, params)
            fields(attempt, {"eligible_reference_count": len(candidates),
                "best_eligible_accession": candidates[0]["seed_hit"]["accession"] if candidates else None,
                "best_eligible_score": base.score(candidates[0]["seed_hit"]) if candidates else None}, "Independent cached fragment candidates")
            known_candidates[interval] = candidates[0] if candidates else None
            if interval in accepted_by_interval:
                require(candidates, "Selected source interval has no eligible candidate")
                fields(accepted_by_interval[interval], candidates[0], "Masked cached fragment selected seed")
            counters["attempts_replayed_from_complete_original_rankings"] += 1
        else:
            counters["fresh_fragment_attempts_without_saved_full_ranking"] += 1
        if attempt["eligible_reference_count"] == 0:
            require(attempt["best_eligible_accession"] is None and attempt["best_eligible_score"] is None, "Empty attempt has a best hit")
        else:
            require(attempt["best_eligible_accession"] in included and math.isfinite(attempt["best_eligible_score"]), "Invalid attempt best hit")
    # Replay active intervals and strongest saved best-hit scores without labels.
    by_interval = dict(zip(intervals, attempts))
    active, discarded, visited = [(0, length)], [], set()
    for segment in accepted:
        candidates = []
        for interval in list(active):
            a, b = interval
            if b - a < params["minimum_span"]:
                active.remove(interval); discarded.append({"start0": a, "end0": b, "reason": "below_minimum_span"}); continue
            require(interval in by_interval, "Active candidate interval was never aligned")
            visited.add(interval)
            attempt = by_interval[interval]
            if attempt["eligible_reference_count"] == 0:
                active.remove(interval); discarded.append({"start0": a, "end0": b, "reason": "no_eligible_reference"})
            else:
                candidates.append((attempt["best_eligible_score"], attempt["best_eligible_accession"], interval))
        chosen_interval = tuple(segment["source_interval"])
        require(chosen_interval in active and candidates, "Selection is outside active candidates")
        best_value = min((-value, accession) for value, accession, _ in candidates)
        require((-base.score(segment["seed_hit"]), segment["seed_hit"]["accession"]) == best_value, "Saved global candidate choice disagrees with score/accession order")
        active.remove(chosen_interval)
        for interval in ((chosen_interval[0], segment["start0"]), (segment["end0"], chosen_interval[1])):
            if interval[1] > interval[0]:
                if interval[1] - interval[0] >= params["minimum_span"]:
                    active.append(interval)
                else:
                    discarded.append({"start0": interval[0], "end0": interval[1], "reason": "below_minimum_span"})
    if len(accepted) < params["maximum_segments"]:
        for interval in list(active):
            a, b = interval
            if b - a < params["minimum_span"]:
                reason = "below_minimum_span"
            else:
                require(interval in by_interval and by_interval[interval]["eligible_reference_count"] == 0, "Unfinished eligible interval below cap")
                visited.add(interval); reason = "no_eligible_reference"
            active.remove(interval); discarded.append({"start0": a, "end0": b, "reason": reason})
    require(visited == set(intervals), "Trace contains unused/missing alignment attempts")
    same(prediction["discarded_intervals"], discarded, "Independent discarded interval trace")
    same(prediction["unprocessed_intervals"], [list(interval) for interval in sorted(active)], "Unprocessed intervals")
    require(prediction["cap_bound"] == (len(accepted) == params["maximum_segments"] and bool(active)), "Cap-bound status differs")
    if fallback:
        require(not accepted, "Fallback has accepted seeds")
        same(segments, [{"start0": 0, "end0": length, "source_interval": [0, length], "seed_hit": None,
                         "selection_order": None, "fallback": True}], "Whole-query fallback")
    else:
        require(segments == sorted(segments, key=lambda row: (row["start0"], row["end0"])) and
                all(a["end0"] <= b["start0"] for a, b in zip(segments, segments[1:])), "Accepted spans overlap or reorder")
    require(prediction["total_predicted_segment_residues"] == sum(s["end0"] - s["start0"] for s in segments), "Predicted residue accounting differs")


def audit_proposal(record, query, config, sequence, included, raw, ctx, counters, common_mask=None):
    same(record["query"], query, "Prediction query metadata")
    same(record["prediction_config"], geometry(config), "Prediction parameters")
    prediction = record["prediction"]
    require(prediction["boundaries_supplied"] is False and prediction["segment_count"] == len(prediction["segments"]) > 0, "Prediction information/coverage differs")
    if config["segmenter"] == "fixed_windows":
        same(prediction, expected_windows(len(sequence), config), "Independent fixed windows")
        scanned = 0
    else:
        audit_alignment_trace(prediction, sequence, config, included, raw, ctx, counters)
        scanned = prediction["total_queried_residues"]
    encoded = sum(s["end0"] - s["start0"] for s in prediction["segments"])
    fields(record, {"logical_candidate_scan_residues": scanned, "logical_encoded_segment_residues": encoded,
                   "logical_total_queried_residues": scanned + encoded,
                   "boundary_evaluation": interval_metrics(prediction["segments"], query["spans"], query["length_aa"])}, "Independent prediction work/boundaries")
    if common_mask is not None:
        fields(record, {"common_mask": common_mask, "prediction_catalogue_accessions": included,
                       "boundary_labels_used_for_prediction": False}, "Prediction common mask/information")


def float32_mean(path):
    """Read one small saved NumPy vector with stdlib; never load model libraries."""
    with zipfile.ZipFile(path) as archive:
        require(archive.namelist() == ["mean.npy"], "Vector archive has unexpected arrays")
        info = archive.getinfo("mean.npy")
        require(info.file_size <= 1024 * 1024, "Vector archive exceeds bounded size")
        stream = io.BytesIO(archive.read(info))
    require(stream.read(6) == b"\x93NUMPY", "Invalid NumPy vector header")
    version = tuple(stream.read(2))
    require(version in ((1, 0), (2, 0)), "Unsupported NumPy vector version")
    size_bytes, fmt = (2, "<H") if version == (1, 0) else (4, "<I")
    header_size = struct.unpack(fmt, stream.read(size_bytes))[0]
    require(header_size <= 16384, "Unbounded NumPy header")
    header = ast.literal_eval(stream.read(header_size).decode("latin1"))
    shape = header["shape"]
    require(header["descr"] == "<f4" and header["fortran_order"] is False and
            type(shape) is tuple and len(shape) == 1 and type(shape[0]) is int and 0 < shape[0] <= 4096, "Saved vector dtype/shape differs")
    data = stream.read()
    require(len(data) == 4 * shape[0], "Saved vector byte count differs")
    values = struct.unpack("<" + "f" * shape[0], data)
    require(all(math.isfinite(v) for v in values), "Nonfinite saved vector")
    return values


def cosine(vector, reference):
    require(len(vector) == len(reference), "Vector dimensions differ")
    denominator = math.sqrt(math.fsum(v * v for v in vector)) * math.sqrt(math.fsum(v * v for v in reference))
    require(denominator > 0, "Zero-norm saved vector")
    return math.fsum(a * b for a, b in zip(vector, reference)) / denominator


class CacheEvidence:
    def __init__(self, sources, model, parent, sample_limit=32):
        require(0 < sample_limit <= 32, "Numeric sample limit must be1 through32")
        self.sources, self.model, self.parent = sources, model, parent
        self.cache_index, self.verified, self.samples = {}, {}, {}
        self.sample_limit = sample_limit
        for path in (sources.session / "caches").glob("*/identity.json"):
            identity = read(path)
            digest = canonical_hash(identity)
            require(digest not in self.cache_index, "Ambiguous cache identity")
            self.cache_index[digest] = (path.parent.resolve(), identity)
        self.catalogue_by_hash = {}
        for a in sorted(sources.ctx["hashes"]):
            key = "catalogue-" + a
            metadata = parent.record(key)
            fields(metadata, {"accession": a, "sequence_sha256": sources.ctx["hashes"][a], "length": sources.ctx["lengths"][a]}, "Parent catalogue vector metadata")
            self.catalogue_by_hash.setdefault(sources.ctx["hashes"][a], parent.units[key])
    def unit(self, cache_hash, sequence_hash, unit_path=None, descriptor_hash=None):
        require(cache_hash in self.cache_index, "Segment cache identity is missing")
        root, identity = self.cache_index[cache_hash]
        require(root.is_relative_to(self.sources.session / "caches") and root.name.endswith(cache_hash[:16]), "Cache path/identity differs")
        same(identity["model"], self.model, "Exact segment representation identity")
        require(identity["pooling"] == "float64 residue mean stored float32; cosine float64", "Cache pooling differs")
        if "catalogue_sequence_sha256" in identity:
            same(identity["catalogue_sequence_sha256"], self.sources.ctx["hashes"], "Cache catalogue sequence hashes")
        if unit_path is None:
            candidates = [p.parent for p in (root / "units" / sequence_hash).glob("*/complete.json") if not p.parent.name.startswith(".")]
            require(len(candidates) == 1, "Unpinned primary segment cache has missing/ambiguous versions")
            path = candidates[0].resolve()
            descriptor_hash = sha(path / "complete.json")
        else:
            path = Path(unit_path).resolve()
        require(path.is_relative_to(root), "Segment evidence escaped cache root")
        key = (str(path), descriptor_hash)
        if key not in self.verified:
            self.sources.deadline.check()
            base.verify_unit(root, sequence_hash, {"path": str(path.relative_to(root)), "descriptor_sha256": descriptor_hash}, cache_hash)
            descriptor = read(path / "complete.json")
            require(set(descriptor["files"]) == {"record.json", "arrays.npz"}, "Segment cache lacks exactly record and vector")
            self.verified[key] = {"path": str(path), "descriptor_sha256": descriptor_hash, "record_sha256": descriptor["files"]["record.json"]["sha256"],
                                  "arrays_sha256": descriptor["files"]["arrays.npz"]["sha256"], "cache_identity_sha256": cache_hash}
        require(sha(path / "record.json") == self.verified[key]["record_sha256"], "Segment cache record changed")
        return path, read(path / "record.json")
    def offer_sample(self, sequence_hash, vector_unit, ranking_unit, pointer):
        priority = hashlib.sha256(("iri-independent-segment-vector-v1|" + sequence_hash).encode()).hexdigest()
        if sequence_hash in self.samples:
            return
        self.samples[sequence_hash] = {"priority": priority, "vector_unit": str(vector_unit), "ranking_unit": str(ranking_unit), "pointer": pointer,
            "vector_descriptor_sha256": sha(vector_unit / "complete.json"), "vector_array_sha256": sha(vector_unit / "arrays.npz"),
            "ranking_record_sha256": sha(ranking_unit / "record.json")}
        if len(self.samples) > self.sample_limit:
            last = max(self.samples, key=lambda key: (self.samples[key]["priority"], key))
            del self.samples[last]
    def family_rows(self, record, query, sequence, output_unit):
        references, segments = record["segment_ranking_references"], record["prediction"]["segments"]
        require(len(references) == len(segments), "Segment cache reference coverage differs")
        result = []
        for evidence, span in zip(references, segments):
            fragment = sequence[span["start0"]:span["end0"]]
            digest = hashlib.sha256(fragment.encode()).hexdigest()
            fields(evidence, {"sequence_sha256": digest, "length": len(fragment), "ranking_pointer": "/full_catalogue_ranking"}, "Fragment cache sequence/pointer")
            require(evidence["source"] in ("parent_catalogue", "new_encode", "segment_cache"), "Unknown cache work source")
            path, metadata = self.unit(evidence["cache_identity_hash"], digest, evidence["unit_path"], evidence["descriptor_sha256"])
            fields(metadata, {"sequence_sha256": digest, "length": len(fragment), "representation_identity_hash": canonical_hash(self.model)}, "Family segment cache identity")
            rows = metadata["full_catalogue_ranking"]
            base.validate_ranking(rows, sorted(self.sources.ctx["hashes"]))
            result.append(rows)
            self.offer_sample(digest, path, path, ["full_catalogue_ranking"])
        same(record["cache_counts"], dict(Counter(e["source"] for e in references)), "Per-query cache work accounting")
        return result
    def primary_rows(self, record, sequence, identity, output_unit):
        rows, segments = record["segment_rankings"], record["prediction"]["segments"]
        require(len(rows) == len(segments), "Primary segment ranking coverage differs")
        cache_identity = {"model": self.model, "pooling": identity["pooling"], "code_sha256": identity["code_sha256"], "pooling_code_sha256": identity["pooling_code_sha256"]}
        cache_hash = canonical_hash(cache_identity)
        for index, (ranking, span) in enumerate(zip(rows, segments)):
            base.validate_ranking(ranking, sorted(self.sources.ctx["hashes"]))
            fragment = sequence[span["start0"]:span["end0"]]
            digest = hashlib.sha256(fragment.encode()).hexdigest()
            if digest in self.catalogue_by_hash:
                vector_unit = self.catalogue_by_hash[digest]
            else:
                vector_unit, metadata = self.unit(cache_hash, digest)
                fields(metadata, {"sequence_sha256": digest, "length": len(fragment)}, "Primary segment cache identity")
            self.offer_sample(digest, vector_unit, output_unit, ["segment_rankings", index])
        return rows
    def numerical_report(self):
        completion = read(self.parent.completion_path)
        for a in sorted(self.sources.ctx["hashes"]):
            key = "catalogue-" + a
            base.verify_unit(self.parent.path, key, completion["units"][key], self.parent.identity_hash)
        references = {a: float32_mean(self.parent.units["catalogue-" + a] / "arrays.npz") for a in sorted(self.sources.ctx["hashes"])}
        checked = []
        for digest, sample in sorted(self.samples.items(), key=lambda item: (item[1]["priority"], item[0])):
            self.sources.deadline.check()
            vector_unit, ranking_unit = Path(sample["vector_unit"]), Path(sample["ranking_unit"])
            require(sha(vector_unit / "complete.json") == sample["vector_descriptor_sha256"] and
                    sha(vector_unit / "arrays.npz") == sample["vector_array_sha256"] and
                    sha(ranking_unit / "record.json") == sample["ranking_record_sha256"], "Numeric sample source changed")
            vector = float32_mean(vector_unit / "arrays.npz")
            rows = read(ranking_unit / "record.json")
            for key in sample["pointer"]:
                rows = rows[key]
            errors = [abs(cosine(vector, references[row["accession"]]) - base.score(row)) for row in rows]
            require(len(errors) == len(references) and max(errors) <= 1e-6, "Independent saved-vector cosine differs")
            checked.append({"sequence_sha256": digest, "maximum_absolute_score_error": max(errors),
                "vector_sha256": sha(vector_unit / "arrays.npz"), "ranking_record_sha256": sha(ranking_unit / "record.json"), "ranking_pointer": sample["pointer"]})
        for value in self.verified.values():
            self.sources.deadline.check()
            path = Path(value["path"])
            require(sha(path / "complete.json") == value["descriptor_sha256"] and sha(path / "record.json") == value["record_sha256"] and
                    sha(path / "arrays.npz") == value["arrays_sha256"], "Verified segment cache changed during audit")
        return {"sample_policy": "Up to32 unique segment hashes with lowest SHA256(iri-independent-segment-vector-v1|sequence_hash), independent of scores",
                "samples_checked": len(checked), "score_tolerance": 1e-6, "evidence": checked,
                "cache_units_verified": len(self.verified), "verified_cache_evidence_sha256": canonical_hash(sorted(self.verified.values(), key=lambda row: row["path"]))}


def verify_production_source(identity):
    pinned = dict(identity.get("code", identity.get("code_hashes", {})))
    for field, path in (("code_sha256", "scripts/session_segment_experiment.py"),
                        ("segmentation_code_sha256", "src/iri_constructs/session_segmentation.py"),
                        ("pooling_code_sha256", "scripts/session_experiment.py"),
                        ("metrics_code_sha256", "src/iri_constructs/session_benchmark.py")):
        if field in identity:
            pinned[path] = identity[field]
    require(pinned, "Segmentation source hashes are absent")
    for relative, expected in pinned.items():
        require(sha(base.safe_file(ROOT, relative)) == expected, "Production source differs from recorded job identity: " + relative)
    return pinned


def static_input_contract():
    """Check reviewed source signatures/calls without importing their modules."""
    path = ROOT / "src/iri_constructs/session_segmentation.py"
    tree = ast.parse(path.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "alignment_guided_segments")
    require([arg.arg for arg in function.args.args] == ["sequence", "catalogue", "config"] and
            [arg.arg for arg in function.args.kwonlyargs] == ["ranker"], "Prediction API accepts unexpected information")
    require(not ({node.id for node in ast.walk(function) if isinstance(node, ast.Name)} &
                 {"targets", "target_spans", "acceptable_matches", "expected_components"}), "Label variables entered predictor")
    evidence = {str(path.relative_to(ROOT)): sha(path)}
    for relative, family in (("scripts/session_segment_experiment.py", False), ("scripts/session_family_segmentation.py", True)):
        file = ROOT / relative
        calls = [node for node in ast.walk(ast.parse(file.read_text())) if isinstance(node, ast.Call) and
                 isinstance(node.func, ast.Name) and node.func.id == "alignment_guided_segments"]
        require(len(calls) == 1, "Unexpected prediction call sites")
        args = calls[0].args
        require(len(args) == 3 and isinstance(args[0], ast.Subscript) and
                isinstance(args[0].value, ast.Name) and args[0].value.id == "q" and
                isinstance(args[0].slice, ast.Constant) and args[0].slice.value == "sequence", "Prediction did not receive only query sequence")
        if family:
            require(isinstance(args[1], ast.Name) and args[1].id == "remaining", "Family prediction did not receive masked catalogue")
        else:
            require(isinstance(args[1], ast.Subscript) and isinstance(args[1].slice, ast.Constant) and args[1].slice.value == "sequences", "Primary prediction catalogue expression changed")
        evidence[relative] = sha(file)
    return {"status": "passed_static_signature_and_callsite_review", "source_sha256": evidence,
            "limitation": "Static source contract and saved traces; not a new execution trace of every fresh alignment candidate."}


def time_order(source, successor, keys=None):
    completed = dt.datetime.fromisoformat(read(source.completion_path)["completed_utc"])
    require(completed.tzinfo is not None, "Source completion time lacks timezone")
    for key in keys or successor.units:
        timestamp = dt.datetime.fromisoformat(read(successor.units[key] / "complete.json")["completed_utc"])
        require(timestamp.tzinfo is not None and timestamp >= completed, "Source stage was not complete before successor output")


def parent_and_predictions(snapshot, sources, family=False):
    identity = snapshot.identity
    parent = sources.open(identity["parent_model_identity"], identity["parent_completion_sha256"])
    require(parent.identity["mode"] == "model" and set(parent.identity["query_ids"]) == set(sources.ctx["queries"]), "Incomplete parent representation job")
    same(identity["model"], parent.identity["method"], "Exact parent representation/backend/library identity")
    predictions = sources.open(identity["prediction_identity" if family else "predictions_identity"],
                               identity["prediction_completion_sha256" if family else "predictions_completion_sha256"])
    expected_mode = "family_automatic_spans_v1" if family else "automatic_spans_v1"
    require(predictions.identity["mode"] == expected_mode, "Wrong prediction stage")
    same(predictions.identity["query_ids"], identity["query_ids"], "Prediction/retrieval query coverage")
    raw = sources.open(predictions.identity["raw_alignment_identity" if family else "alignment_identity"],
                       predictions.identity["raw_alignment_completion_sha256" if family else "alignment_completion_sha256"])
    require(raw.identity["mode"] == "alignment" and set(raw.identity["query_ids"]) == set(sources.ctx["queries"]), "Incomplete original alignment source")
    verify_production_source(identity)
    verify_production_source(predictions.identity)
    time_order(predictions, snapshot)
    time_order(parent, snapshot)
    time_order(raw, predictions)
    return parent, predictions, raw


def selected_configs(design, selection):
    winners = selection["selection"]["selected_by_segmenter"]
    require(set(winners) == {c["segmenter"] for c in design["candidate_configurations"]}, "Missing selected segmenter family")
    by_id = {c["configuration_id"]: c for c in design["candidate_configurations"]}
    require(all(name in by_id and by_id[name]["segmenter"] == family for family, name in winners.items()), "Wrong-family selected configuration")
    return [c for c in design["candidate_configurations"] if c["configuration_id"] in set(winners.values())]


def audit_primary(snapshot, sources, design, selection, sequences, *, development=False):
    identity, ctx = snapshot.identity, sources.ctx
    ids = design["selection_query_ids"] if development else design["test_query_ids"] + design["clean_sanity_query_ids"]
    configs = design["candidate_configurations"] if development else selected_configs(design, selection)
    fields(identity, {"mode": "automatic_segment_retrieval_v1", "design_sha256": DESIGN_SHA,
                      "partition": "dev" if development else "test", "query_ids": ids, "configs": configs,
                      "selection_sha256": None if development else SELECTION_SHA}, "Primary segmentation identity")
    parent, predictions, raw = parent_and_predictions(snapshot, sources)
    unique = {canonical_hash(geometry(config))[:12]: geometry(config) for config in configs}
    fields(predictions.identity, {"design_sha256": DESIGN_SHA, "partition": identity["partition"],
                                  "selection_sha256": identity["selection_sha256"], "configs": unique}, "Primary prediction identity")
    method = identity["model"]
    expected_method = design["selection_representation"]
    fields(method, {key: expected_method[key] for key in ("model_id", "revision", "device", "dtype")}, "Frozen development representation")
    require(type(method["layer"]) is int and method["layer"] == method["num_layers"] > 0 and method["truncated"] is False, "Final-layer/nontruncating representation changed")
    cache = CacheEvidence(sources, method, parent)
    bundles = {c["configuration_id"]: {"all": base.Rollup(), "by_suite": {}} for c in configs}
    candidates = {c["configuration_id"]: {} for c in configs}
    counters, proposal_count = Counter(), 0
    for qid in ids:
        sources.deadline.check()
        query = ctx["queries"][qid]
        sequence = "".join(sequences[a] for a in query["components"])
        checked_proposals = set()
        for config in configs:
            name, geometric_hash = config["configuration_id"], canonical_hash(geometry(config))[:12]
            key, prediction_key = qid + "--" + name, qid + "--" + geometric_hash
            proposal = predictions.record(prediction_key)
            if prediction_key not in checked_proposals:
                audit_proposal(proposal, query, config, sequence, sorted(ctx["hashes"]), raw, ctx, counters)
                checked_proposals.add(prediction_key); proposal_count += 1
            record = snapshot.record(key)
            fields(record, {"query": query, "configuration": config, "prediction_unit": prediction_key,
                            "prediction": proposal["prediction"], "boundary_evaluation": proposal["boundary_evaluation"],
                            "logical_total_queried_residues": proposal["logical_total_queried_residues"]}, "Saved retrieval/proposal identity")
            rankings = cache.primary_rows(record, sequence, identity, snapshot.units[key])
            fused = base.fused_ranking(rankings, config["fusion"], sorted(ctx["hashes"]))
            same(record["ranking"], fused, "Independent primary total-budget fusion")
            evaluation = base.expected_evaluation(query, fused, ctx["equivalence"])
            same(record["evaluation"], evaluation, "Independent primary component evaluation")
            same(record["larger_budget"], {str(k): union_result(rankings, k, query) for k in base.KS}, "Independent primary larger segment budget")
            bundles[name]["all"].add(evaluation)
            bundles[name]["by_suite"].setdefault(query["suite"], base.Rollup()).add(evaluation)
            if development:
                metrics = evaluation["query_metrics"]["exact"]
                candidates[name][qid] = (metrics["component_recall_at_5"], metrics["all_components_hit_at_5"], proposal["logical_total_queried_residues"])
    expected_selection = selection_from_outcomes(candidates, design, ctx) if development else None
    scientific = {"summaries": {name: value["all"].result() for name, value in bundles.items()},
                  "summaries_by_suite": {name: {suite: values.result() for suite, values in value["by_suite"].items()} for name, value in bundles.items()},
                  "selection": expected_selection}
    fields(snapshot.record("summary"), scientific, "Independent primary summary strata/selection")
    if development:
        same(selection["selection"], expected_selection, "Independently recomputed development winner")
        require(selection["development_identity"] == snapshot.identity_hash and
                Path(selection["development_job"]).resolve() == snapshot.path and
                sha(snapshot.units["summary"] / "complete.json") == selection["development_summary_descriptor_sha256"], "Selection lacks exact completed development evidence")
    return {"status": "passed", **snapshot.evidence(), "queries_checked": len(ids), "configurations_checked": len(configs),
            "prediction_units_checked": proposal_count, "trace_evidence": dict(counters), "summary_scientific_sha256": canonical_hash(scientific),
            "saved_vector_cosines": cache.numerical_report(), "independent_selection": expected_selection,
            "checks": ["exact_complete_unit_coverage", "frozen_group_partition", "source_stage_completion_order", "fixed_or_transferred_parameters",
                       "independent_window_geometry_and_alignment_trace", "prediction_input_contract", "all_boundary_interval_metrics",
                       "segment_cache_sequence_and_representation", "single_totalK_fusion", "larger_union_budget_and_review_volume", "all_summary_strata"]}


def common_mask(ctx, query, policy):
    mask = base.family_masks(query, ctx["hashes"], ctx["removals"])["family_" + policy]
    return {**mask, "included_accessions_sha256": canonical_hash(mask["included_accessions"]), "remaining_catalogue_count": len(mask["included_accessions"])}


def audit_family(snapshot, sources, design, transfer, selection, sequences, development):
    identity, ctx = snapshot.identity, sources.ctx
    require(type(identity["include_clean_sanity"]) is bool, "Ambiguous clean sanity flag")
    ids = design["test_query_ids"] + (design["clean_sanity_query_ids"] if identity["include_clean_sanity"] else [])
    configs = selected_configs(design, selection)
    shared = {"primary_design_sha256": DESIGN_SHA, "transfer_design_sha256": TRANSFER_SHA,
              "selection_sha256": SELECTION_SHA, "family_labels_sha256": base.FAMILY_SHA256,
              "query_ids": ids, "configs": configs, "policies": transfer["policies"],
              "parameter_transfer": "unchanged winners from completed original development only",
              "prediction_mask_policy": "same all-components removal mask is applied during prediction and retrieval; no target-boundary reuse"}
    fields(identity, {**shared, "mode": "family_automatic_retrieval_v1",
                      "fusion_order": "filter and rerank each segment using common all-components mask BEFORE fusion"}, "Family transfer identity")
    parent, predictions, raw = parent_and_predictions(snapshot, sources, family=True)
    fields(predictions.identity, shared, "Family prediction transferred settings")
    same(identity["model"], development.identity["model"], "Exact selected development representation")
    time_order(development, predictions)
    cache = CacheEvidence(sources, identity["model"], parent)
    bundles, counters = {}, Counter()
    for policy in transfer["policies"]:
        for config in configs:
            bundles[(policy, config["configuration_id"])] = {
                "scopes": {scope: base.Rollup() for scope in ("all", "combined", "clean")}, "by_suite": {},
                "unions": {scope: {str(k): UnionRollup() for k in base.KS} for scope in ("all", "combined", "clean")},
                "resources": {"cache_counts": Counter(), "logical_candidate_scan_residues": 0, "logical_encoded_segment_residues": 0,
                              "cap_bound_queries": 0, "whole_query_fallback_queries": 0}}
    for qid in ids:
        sources.deadline.check()
        query = ctx["queries"][qid]
        sequence = "".join(sequences[a] for a in query["components"])
        for policy in transfer["policies"]:
            mask = common_mask(ctx, query, policy)
            included = mask["included_accessions"]
            for config in configs:
                name = config["configuration_id"]
                key = qid + "--" + policy + "--" + name
                proposal = predictions.record(key)
                audit_proposal(proposal, query, config, sequence, included, raw, ctx, counters, mask)
                record = snapshot.record(key)
                fields(record, {"query": query, "configuration": config, "common_mask": mask,
                    "prediction": proposal["prediction"], "boundary_evaluation": proposal["boundary_evaluation"], "prediction_unit": key,
                    "prediction_unit_descriptor_sha256": sha(predictions.units[key] / "complete.json"),
                    **{field: proposal[field] for field in ("logical_candidate_scan_residues", "logical_encoded_segment_residues", "logical_total_queried_residues")}}, "Family proposal/retrieval evidence")
                full = cache.family_rows(record, query, sequence, snapshot.units[key])
                masked = [independently_ranked(rows, included) for rows in full]
                fused = base.fused_ranking(masked, config["fusion"], included)
                same(record["ranking"], fused, "Independent mask-before-fusion ranking")
                evaluation = base.expected_evaluation(query, fused, ctx["equivalence"], ctx["families"])
                same(record["evaluation"], evaluation, "Independent family transfer component metrics")
                require(evaluation["query_metrics"]["exact"]["scorable_components"] == 0 and
                        evaluation["query_metrics"]["duplicate_equivalent"]["scorable_components"] == 0, "Removed exact label leaked into retrieval")
                unions = {str(k): union_result(masked, k, query, ctx["families"], included) for k in base.KS}
                same(record["larger_budget"], unions, "Independent family larger-budget metrics/volume")
                bundle, scope = bundles[(policy, name)], "clean" if query["component_count"] == 1 else "combined"
                for target in ("all", scope):
                    bundle["scopes"][target].add(evaluation)
                    for k, value in unions.items():
                        bundle["unions"][target][k].add(qid, value)
                bundle["by_suite"].setdefault(query["suite"], base.Rollup()).add(evaluation)
                resources = bundle["resources"]
                resources["cache_counts"].update(record["cache_counts"])
                for field in ("logical_candidate_scan_residues", "logical_encoded_segment_residues"):
                    resources[field] += proposal[field]
                resources["cap_bound_queries"] += int(proposal["prediction"]["cap_bound"])
                resources["whole_query_fallback_queries"] += int(proposal["prediction"].get("whole_query_fallback_used", False))
    summaries = {}
    for (policy, name), bundle in bundles.items():
        summaries.setdefault(policy, {})[name] = {**{scope: value.result() for scope, value in bundle["scopes"].items()},
            "by_suite": {suite: value.result() for suite, value in bundle["by_suite"].items()},
            "larger_budget_summaries": {scope: {k: value.result() for k, value in values.items()} for scope, values in bundle["unions"].items()},
            "resources": {**bundle["resources"], "cache_counts": dict(bundle["resources"]["cache_counts"])}}
    denominator_evidence = {policy: base.denominator_totals(base.denominator_rows([ctx["queries"][qid] for qid in ids],
        ctx["hashes"], ctx["families"], ctx["removals"], policy)) for policy in transfer["policies"]}
    scientific = {"summaries": summaries, "frozen_family_denominators_for_panel": denominator_evidence, "query_count": len(ids),
                  "transfer_only_no_retuning": True, "common_mask_applied_during_prediction_and_retrieval": True}
    fields(snapshot.record("summary"), scientific, "Independent family transfer summary strata and denominators")
    return {"status": "passed", **snapshot.evidence(), "queries_checked": len(ids), "combined_test_queries": len(design["test_query_ids"]),
            "configurations_checked": len(configs), "policies_checked": transfer["policies"], "trace_evidence": dict(counters),
            "summary_scientific_sha256": canonical_hash(scientific), "family_denominators": denominator_evidence,
            "saved_vector_cosines": cache.numerical_report(),
            "checks": ["completed_development_winner_unchanged", "exact_development_model_identity", "heldout_group_partition",
                       "source_stage_completion_order", "all_component_common_mask_before_prediction", "all_seed_eligibility_and_available_source_rankings",
                       "checksum_pinned_segment_cache", "each_segment_masked_and_reranked_before_fusion", "totalK_and_larger_union_budgets",
                       "all_exact_duplicate_family_metrics", "all_summary_strata_and_null_denominators", "boundary_and_resource_accounting"]}


def run(args):
    session = args.session.resolve()
    output = (args.output or session / "audit/independent_segmentation_audit_v1.json").resolve()
    require(output.is_relative_to(session / "audit") and not output.exists(), "Use a new owned audit report path")
    started = time.monotonic()
    report = {"status": "running", "version": "iri-independent-segmentation-audit-v1", "jobs": [],
              "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "timeout_seconds": args.timeout_seconds,
              "audit_code_sha256": sha(Path(__file__)), "independent_extension_helper_sha256": sha(Path(base.__file__)),
              "independent_unit_helper_sha256": sha(Path(__file__).with_name("session_independent_audit.py")),
              "scope": "Independent saved-result metrics, selection, masks, traces, cache identity and bounded saved-vector arithmetic; no model inference or fresh alignment",
              "limitations": ["Full fresh-fragment alignment candidate lists were not retained; their numerical DP optimality and every coordinate tie cannot be independently replayed from saved artifacts.",
                              "The source and prediction interface review plus saved trace checks do not constitute a new runtime information-flow proof.",
                              "Primary exact segmentation persisted full segment rankings but no per-output cache descriptor references; matching vector cache versions must be unique and are checked as current supporting evidence.",
                              "Boundary metrics describe artificial component intervals, not source-family domain/function correctness.",
                              "Transferred family larger-budget summaries have all/combined/clean scopes; production did not save union-by-suite summaries."]}
    sources = None
    def interrupted(signum, frame):
        raise InterruptedError("Segmentation audit interrupted by signal " + str(signum))
    previous = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        deadline = base.Deadline(args.timeout_seconds, session)
        deadline.check()
        ctx = base.context(session, deadline)
        sources = Sources(session, deadline, ctx)
        design, transfer, selection = frozen_designs(session, ctx)
        report["frozen_inputs"] = {**ctx["evidence"], "segmentation_design_sha256": DESIGN_SHA,
                                   "family_transfer_design_sha256": TRANSFER_SHA, "development_selection_sha256": SELECTION_SHA}
        report["prediction_static_source_review"] = static_input_contract()
        sequences = trusted_target_sequences(ctx)
        development = sources.by_path(selection["development_job"])
        report["active_job"] = str(development.path)
        report["development_audit"] = audit_primary(development, sources, design, selection, sequences, development=True)
        base.atomic_report(output, report)
        print(json.dumps({"status": "development_audit_passed", "job": str(development.path)}), flush=True)
        jobs = args.job or [session / "jobs" / name for name in DEFAULT_JOBS]
        require(len(jobs) == args.expected_jobs and len({Path(p).resolve() for p in jobs}) == len(jobs), "Heldout audit job coverage differs")
        report["requested_jobs"] = [str(Path(path).resolve()) for path in jobs]
        for path in jobs:
            deadline.check()
            snapshot = sources.by_path(path)
            report["active_job"] = str(snapshot.path)
            mode = snapshot.identity["mode"]
            same(snapshot.identity["model"], development.identity["model"], "Heldout representation differs from selected development")
            if mode == "automatic_segment_retrieval_v1":
                result = audit_primary(snapshot, sources, design, selection, sequences)
                predictions = sources.open(snapshot.identity["predictions_identity"], snapshot.identity["predictions_completion_sha256"])
                time_order(development, predictions)
            elif mode == "family_automatic_retrieval_v1":
                result = audit_family(snapshot, sources, design, transfer, selection, sequences, development)
            else:
                raise ValueError("Requested job is not heldout exact/family segment retrieval")
            report["jobs"].append(result)
            base.atomic_report(output, report)
            print(json.dumps({"status": "job_passed", "job": str(snapshot.path), "jobs_checked": len(report["jobs"])}), flush=True)
        sources.unchanged()
        report["status"] = "passed"
        report.pop("active_job", None)
    except (Exception, KeyboardInterrupt) as error:
        report["status"] = "failed"
        report["failure"] = {"type": type(error).__name__, "message": str(error)}
    finally:
        if sources is not None:
            report["verified_snapshots"] = [snapshot.evidence() for snapshot in sources.snapshots.values()]
        report["elapsed_seconds"] = time.monotonic() - started
        report["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        base.atomic_report(output, report)
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    print(json.dumps({"status": report["status"], "report": str(output), "jobs_checked": len(report["jobs"])}), flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", type=Path, required=True)
    parser.add_argument("--job", type=Path, action="append")
    parser.add_argument("--expected-jobs", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=float, default=900)
    parser.add_argument("--output", type=Path)
    raise SystemExit(run(parser.parse_args()))
