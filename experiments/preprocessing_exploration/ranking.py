"""Replay geometry-only ranking policies on archived, sequence-free candidate pools."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
from itertools import islice
import json
from pathlib import Path
from statistics import mean
import sys
from time import perf_counter

from experiments.paired_preprocessing.verify import _covered, _location, _score_construct


ARMS = ("minimal", "start_ranked", "gene_caller")
SCORERS = ("esm2", "blosum62")
QUERY_CAPS = (1, 2, 3, 4, 6, 8, 12)
RESIDUE_CAPS = (512, 1024, 2048, 4096)
REPLAY_CAPS = (1, 2, 3, 4)
MODEL_LENGTH_CEILING = 1022
LENGTH_CENTERS = (250, 450, 650)
POLICIES = (
    "recorded", "defer_ceiling",
    *(f"length_center_{center}" for center in LENGTH_CENTERS),
    "novel_span_any", "novel_span_same",
    "soft_diversity_any", "soft_diversity_same", "balanced_diversity",
)
COHORT = Path("results/paired-cohort-v2")


@dataclass(frozen=True)
class Candidate:
    strand: int
    spans: tuple[tuple[int, int], ...]
    residues: int

    @property
    def locus(self):
        return self.strand, self.spans


def candidate(row):
    residues = row["residues"] if "residues" in row else row["input_residues"]
    return Candidate(row["strand"], tuple(sorted(map(tuple, row["spans"]))), residues)


def merged_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(end, merged[-1][1])
        else:
            merged.append((start, end))
    return merged


def length_quality(residues, center):
    return min(residues / center, center / residues)


def _overlap_fraction(unit, covered):
    size = sum(end - start for start, end in unit.spans)
    overlap = sum(
        max(0, min(end, existing_end) - max(start, existing_start))
        for start, end in unit.spans
        for existing_start, existing_end in covered
    )
    return overlap / size


def _ranked_indices(pool, policy):
    if policy == "keep_first_then_250":
        if pool:
            yield 0
            yield from sorted(range(1, len(pool)), key=lambda index: -length_quality(pool[index].residues, 250))
        return
    if policy not in POLICIES:
        raise ValueError(f"Unrecognized ranking policy: {policy}")
    indices = range(len(pool))
    if policy == "recorded":
        yield from indices
        return
    if policy == "defer_ceiling":
        yield from sorted(indices, key=lambda index: pool[index].residues >= MODEL_LENGTH_CEILING)
        return
    if policy.startswith("length_center_"):
        center = int(policy.removeprefix("length_center_"))
        yield from sorted(indices, key=lambda index: -length_quality(pool[index].residues, center))
        return

    same_strand = policy.endswith("_same")
    covered = defaultdict(list)
    remaining = list(indices)
    while remaining:
        def score(index):
            unit = pool[index]
            strand = unit.strand if same_strand else 0
            overlap = _overlap_fraction(unit, covered[strand])
            if policy.startswith("novel_span_"):
                return unit.residues * (1 - overlap)
            quality = (
                length_quality(unit.residues, 450)
                if policy == "balanced_diversity"
                else 1 / (index + 1) ** 0.5
            )
            return quality * (1 - overlap / 2)

        picked = max(remaining, key=score)
        yield picked
        remaining.remove(picked)
        picked_unit = pool[picked]
        strand = picked_unit.strand if same_strand else 0
        covered[strand] = merged_intervals((*covered[strand], *picked_unit.spans))


def rank_indices(pool, policy, *, limit=None):
    """Order input indices without any evaluator-side annotation or scorer inputs."""
    return list(islice(_ranked_indices(pool, policy), limit))


def _budget_curves(pool, policy):
    by_count = []
    by_residues = {cap: [] for cap in RESIDUE_CAPS}
    remaining = dict.fromkeys(RESIDUE_CAPS, 0)
    remaining.update({cap: cap for cap in RESIDUE_CAPS})
    shortest = min((unit.residues for unit in pool), default=1)
    for index in _ranked_indices(pool, policy):
        if len(by_count) < QUERY_CAPS[-1]:
            by_count.append(index)
        cost = pool[index].residues
        for cap, allowance in remaining.items():
            if cost <= allowance:
                by_residues[cap].append(index)
                remaining[cap] -= cost
        if len(by_count) == QUERY_CAPS[-1] and max(remaining.values()) < shortest:
            break
    return {cap: by_count[:cap] for cap in QUERY_CAPS}, by_residues


def _frame(locus):
    strand, spans = locus
    if len(spans) != 1:
        return None
    first, last = spans[0]
    return (first if strand == 1 else last) % 3


def _reach_for_pool(pool, loci, overlap):
    all_hits, frame_hits = [], []
    frames = [_frame(locus) for locus in loci]
    for unit in pool:
        unit_frame = _frame(unit.locus)
        covered = sum(1 << index for index, locus in enumerate(loci) if _covered(locus, unit.locus) >= overlap)
        same_frame = sum(
            1 << index
            for index, frame in enumerate(frames)
            if covered & (1 << index) and frame is not None and frame == unit_frame
        )
        all_hits.append(covered)
        frame_hits.append(same_frame)
    return all_hits, frame_hits


def _hits(selected, candidate_hits):
    covered = 0
    for index in selected:
        covered |= candidate_hits[index]
    return covered


def summarize(rows):
    lineage = defaultdict(lambda: [0, 0])
    strata = defaultdict(lambda: [0, 0])
    reached = frame_denominator = frame_reached = 0
    for truth, outcome in rows:
        component_count = len(truth["components"])
        group_total = lineage[truth["lineage"]]
        group_total[0] += component_count
        group_total[1] += outcome["hits"].bit_count()
        reached += outcome["hits"].bit_count()
        frame_reached += outcome["frame_hits"].bit_count()
        frame_denominator += truth["frame_denominator"]
        for index, tier in enumerate(truth["tiers"]):
            hit = bool(outcome["hits"] & 1 << index)
            for field, label in (("layout", truth["layout"]), ("similarity", tier)):
                totals = strata[field, label]
                totals[0] += 1
                totals[1] += hit
    return {
        "selected_candidates": sum(outcome["candidates"] for _, outcome in rows),
        "input_residues": sum(outcome["residues"] for _, outcome in rows),
        "components": sum(n for n, _ in lineage.values()),
        "reached": reached,
        "lineage_groups": sum(n > 0 for n, _ in lineage.values()),
        "equal_lineage_recall": mean(hits / n for n, hits in lineage.values() if n),
        "single_span_components": frame_denominator,
        "single_span_same_frame_reached": frame_reached,
        "strata": {
            f"{field}:{label}": {"components": n, "reached": hits}
            for (field, label), (n, hits) in sorted(strata.items())
        },
    }


def _replayed_predictions(queries, rankings, limit):
    ranks = max((len(rankings[query["id"]]) for query in queries), default=0)
    return list(islice((
        {"query_id": query["id"], "family": rankings[query["id"]][rank]}
        for rank in range(ranks)
        for query in queries
        if rank < len(rankings[query["id"]])
    ), limit))


def _replay(construct, max_tier, selected, pool_size, cap, protocol):
    queries = [max_tier["queries"][index] for index in selected]
    predictions = _replayed_predictions(
        queries, max_tier["query_family_rankings"], protocol["prediction_limit"]
    )
    trace = {
        "status": "ok", "queries": queries, "predictions": predictions,
        "candidate_count": pool_size, "preprocess_seconds": 0, "characterize_seconds": 0,
    }
    result = _score_construct(
        construct, trace, cap, protocol["prediction_limit"], protocol["minimum_component_overlap"]
    )
    recovered = set(result["recovered_component_ids"])
    return sum(1 << index for index, component in enumerate(construct["components"]) if component["id"] in recovered)


def _truth(construct, metadata):
    loci = [_location(row, construct["length"], row["id"]) for row in construct["components"]]
    tiers = {row["id"]: row["similarity"]["tier"] for row in metadata["components"]}
    return {
        "components": loci, "layout": construct["archetype"], "lineage": metadata["lineage_group"],
        "frame_denominator": sum(_frame(locus) is not None for locus in loci),
        "tiers": [tiers[row["id"]] for row in construct["components"]],
    }


def evaluate_arm(traces, protocol, constructs, metadata, policies=POLICIES):
    """Evaluate selections from each trace and replay labels only for already scored candidates."""
    primary = traces[SCORERS[0]]
    pools = primary["candidate_pools"]
    if pools.keys() != constructs.keys() or any(
        trace["candidate_pools"] != pools
        or trace["selected_inputs_sha256"] != primary["selected_inputs_sha256"]
        or trace["input_sha256"] != primary["input_sha256"]
        or trace["reference_panel_sha256"] != primary["reference_panel_sha256"]
        for trace in traces.values()
    ):
        raise ValueError("Trace sources disagree on the selected inputs, pools, or references")
    tiers = {
        name: {
            row["query_cap"]: {item["id"]: item for item in row["constructs"]}
            for row in trace["tiers"]
        }
        for name, trace in traces.items()
    }

    reach, replay, full_pool = defaultdict(list), defaultdict(list), []
    ranking_seconds = 0.0
    for key, construct in constructs.items():
        truth = _truth(construct, metadata[key])
        pool = [candidate(row) for row in pools[key]]
        all_hits, frame_hits = _reach_for_pool(pool, truth["components"], protocol["minimum_component_overlap"])
        full_pool.append((truth, {
            "hits": _hits(range(len(pool)), all_hits),
            "frame_hits": _hits(range(len(pool)), frame_hits),
            "candidates": len(pool), "residues": sum(unit.residues for unit in pool),
        }))
        scored = [candidate(row) for row in tiers[SCORERS[0]][4][key]["queries"]]
        if scored != pool[:4] or any(
            [candidate(row) for row in tiers[name][4][key]["queries"]] != scored for name in SCORERS
        ):
            raise ValueError(f"Scored candidates do not match the first four pool entries for {key}")
        for policy in policies:
            started = perf_counter()
            candidate_caps, residue_caps = _budget_curves(pool, policy)
            scored_order = rank_indices(scored, policy)
            ranking_seconds += perf_counter() - started

            for budget, runs in (("candidate_cap", candidate_caps), ("residue_cap", residue_caps)):
                for cap, indices in runs.items():
                    reach[policy, budget, cap].append((truth, {
                        "hits": _hits(indices, all_hits), "frame_hits": _hits(indices, frame_hits),
                        "candidates": len(indices), "residues": sum(pool[index].residues for index in indices),
                    }))
            for cap in REPLAY_CAPS:
                indices = scored_order[:cap]
                for scorer in SCORERS:
                    max_tier = tiers[scorer][4][key]
                    if policy == "recorded" and cap in tiers[scorer]:
                        original = tiers[scorer][cap][key]
                        queries = [max_tier["queries"][index] for index in indices]
                        if original["queries"] != queries or original["predictions"] != _replayed_predictions(
                            queries, max_tier["query_family_rankings"], protocol["prediction_limit"]
                        ):
                            raise ValueError(f"Replay disagrees with recorded predictions for {key} at cap {cap}")
                    replay[policy, scorer, cap].append((truth, {
                        "hits": _replay(construct, max_tier, indices, len(pool), cap, protocol),
                        "frame_hits": 0, "candidates": len(indices),
                        "residues": sum(scored[index].residues for index in indices),
                    }))

    def scored_summary(rows):
        result = summarize(rows)
        result.pop("single_span_components")
        result.pop("single_span_same_frame_reached")
        result["identified"] = result.pop("reached")
        for group in result["strata"].values():
            group["identified"] = group.pop("reached")
        return result

    return {
        "pool_coverage": summarize(full_pool),
        "reach_curves": {
            policy: {
                budget: {str(cap): summarize(reach[policy, budget, cap]) for cap in caps}
                for budget, caps in (("candidate_cap", QUERY_CAPS), ("residue_cap", RESIDUE_CAPS))
            }
            for policy in policies
        },
        "identification_replay_limited_to_recorded_top_four": {
            policy: {
                scorer: {str(cap): scored_summary(replay[policy, scorer, cap]) for cap in REPLAY_CAPS}
                for scorer in SCORERS
            }
            for policy in policies
        },
    }, ranking_seconds


def run(root, splits, policies=POLICIES):
    inputs = {}

    def read(path):
        contents = (root / path).read_bytes()
        inputs[str(path)] = sha256(contents).hexdigest()
        return json.loads(contents)

    protocol = read(Path("experiments/paired_preprocessing/cohort_protocol.json"))
    results = {}
    for split in splits:
        labels = COHORT / "prepared" / split / "evaluator"
        constructs = {
            row["id"]: row for row in read(labels / "answer-key.corroborated.json")["constructs"]
        }
        metadata = {row["id"]: row for row in read(labels / "construct-map.private.json")}
        trace_run = COHORT / ("dev-run3" if split == "development" else split)
        arms = {}
        for arm in ARMS:
            traces = {scorer: read(trace_run / scorer / f"{arm}.json") for scorer in SCORERS}
            arms[arm], elapsed = evaluate_arm(traces, protocol, constructs, metadata, policies)
            print(f"{split} / {arm}: {elapsed:.2f} s for ranking and budget selection", file=sys.stderr)
        results[split] = arms
    return {"schema": "iri-preprocessing-ranking-exploration-1", "input_sha256": inputs, "splits": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--split", choices=("development", "holdout"), action="append")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--policy", choices=(*POLICIES, "keep_first_then_250"), action="append")
    args = parser.parse_args()
    result = run(args.root, args.split or ("development", "holdout"), args.policy or POLICIES)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
