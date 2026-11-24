"""Retrospective cross-construct candidate budgeting on sequence-free traces."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
import heapq
import json
import math
from pathlib import Path

from experiments.paired_preprocessing.verify import _covered, _location, _score_construct


SCORERS = ("esm2", "blosum62")
METRICS = ("coverage", *SCORERS)
METHODS = ("minimal", "gene_caller")
CAPS = (1, 2, 3, 4)
SAVED_CAPS = (1, 2, 4)
PRIORITIES = (
    "cheapest",
    "longest",
    "pool_demand",
    "pool_demand_per_residue",
    "new_region",
    "balanced_cheapest",
    "balanced_demand",
)


@dataclass(frozen=True)
class Candidate:
    residues: int
    locus: tuple

    @classmethod
    def from_metadata(cls, row):
        return cls(row["residues"], (row["strand"], tuple(map(tuple, row["spans"]))))


@dataclass(frozen=True)
class Budget:
    candidates: int
    residues: int

    def __add__(self, other):
        return Budget(self.candidates + other.candidates, self.residues + other.residues)

    def fits(self, limit):
        return self.candidates <= limit.candidates and self.residues <= limit.residues


def read_json(path):
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def allocation_cost(pools, prefix_lengths):
    return Budget(
        sum(prefix_lengths),
        sum(row.residues for pool, count in zip(pools, prefix_lengths, strict=True) for row in pool[:count]),
    )


def priority_key(policy, pool, index):
    candidate = pool[index]
    demand = (len(pool) - index) / (index + 1)
    if policy == "new_region":
        overlap = max((_covered(candidate.locus, previous.locus) for previous in pool[:index]), default=0.0)
        return (overlap, candidate.residues)
    priorities = {
        "cheapest": (candidate.residues,),
        "longest": (-candidate.residues,),
        "pool_demand": (-demand, candidate.residues),
        "pool_demand_per_residue": (-demand / candidate.residues,),
        "balanced_cheapest": (index, candidate.residues),
        "balanced_demand": (index, -demand, candidate.residues),
    }
    return priorities[policy]


def allocate(pools, construct_ids, limits, policy, max_cap=4):
    """Keep each pool's prefix order while assigning effort across the batch."""
    prefix = [min(1, len(pool)) for pool in pools]
    spent = allocation_cost(pools, prefix)
    if not spent.fits(limits):
        raise ValueError("Budget is below the required first candidate per nonempty construct")
    tie = [hashlib.sha256(cid.encode("ascii")).digest() for cid in construct_ids]
    frontier = []

    def offer(i):
        if prefix[i] < min(max_cap, len(pools[i])):
            heapq.heappush(frontier, (priority_key(policy, pools[i], prefix[i]), tie[i], i))

    for i in range(len(pools)):
        offer(i)
    while frontier and spent.candidates < limits.candidates:
        _, _, i = heapq.heappop(frontier)
        cost = Budget(1, pools[i][prefix[i]].residues)
        if (spent + cost).fits(limits):
            prefix[i] += 1
            spent += cost
            offer(i)
    return prefix


def predictions(query_ids, rankings, limit):
    """Use the recorded round-robin family assignment at the original output cap."""
    labels = [rankings[qid] for qid in query_ids]
    calls = []
    for rank in range(max(map(len, labels), default=0)):
        for qid, families in zip(query_ids, labels, strict=True):
            if rank < len(families):
                calls.append({"query_id": qid, "family": families[rank]})
                if len(calls) == limit:
                    return calls
    return calls


def coverage(pools, truth, overlap, limit):
    """Return sets of documented component IDs at a uniform locus-coverage limit."""
    found = []
    for pool, construct in zip(pools, truth, strict=True):
        loci = [
            (component["id"], _location(component, construct["length"], component["id"]))
            for component in construct["components"]
        ]
        found.append({
            cid for cid, locus in loci
            if any(_covered(locus, candidate.locus) >= overlap for candidate in pool[:limit])
        })
    return found


def parse_scored_run(run, metadata, truth, output_limit, overlap):
    """Replay all prefixes and check them against saved ranks, loci, and calls."""
    tiers = {
        tier["query_cap"]: {row["id"]: row for row in tier["constructs"]}
        for tier in run["tiers"]
    }
    if set(tiers) != set(SAVED_CAPS):
        raise ValueError("Expected saved tiers at caps 1, 2, and 4")
    ids = {construct["id"] for construct in truth}
    if any(ids != set(tier) for tier in tiers.values()):
        raise ValueError("Scored construct IDs differ from the answer key")
    recovered = {cap: [] for cap in CAPS}
    for construct, pool in zip(truth, metadata, strict=True):
        cid = construct["id"]
        full = tiers[4][cid]
        for cap in CAPS:
            queries = full["queries"][:cap]
            actual = [
                Candidate.from_metadata({**row, "residues": row["input_residues"]})
                for row in queries
            ]
            if actual != list(pool[:cap]):
                raise ValueError("Scored queries differ from the recorded prefix")
            calls = predictions([row["id"] for row in queries], full["query_family_rankings"], output_limit)
            if cap in tiers:
                saved = tiers[cap][cid]
                if saved["queries"] != queries or saved["predictions"] != calls:
                    raise ValueError(f"Saved ranking replay differs at cap {cap}")
            trace = {**full, "queries": queries, "predictions": calls}
            report = _score_construct(construct, trace, cap, output_limit, overlap)
            recovered[cap].append(set(report["recovered_component_ids"]))
    return recovered


def choose_prefix_outcomes(prefix_outcomes, allocation):
    return [
        prefix_outcomes[prefix][i] if prefix else set()
        for i, prefix in enumerate(allocation)
    ]


def outcome_summary(outcomes, truth, group_ids):
    """Calculate pooled and equal-lineage recall from complete per-record outcomes."""
    denominators, hits = Counter(), Counter()
    for tracked, construct, group in zip(outcomes, truth, group_ids, strict=True):
        denominators[group] += len(construct["components"])
        hits[group] += len(tracked)
    groups = sorted(group for group, count in denominators.items() if count)
    return {
        "recovered": hits.total(),
        "eligible": denominators.total(),
        "pooled": hits.total() / denominators.total(),
        "equal_lineage": sum(hits[group] / denominators[group] for group in groups) / len(groups),
    }


def paired_delta(first, second, truth, group_ids):
    """Replay the project's paired lineage bootstrap for one policy comparison."""
    from experiments.paired_preprocessing.cohort_analysis import paired_effect

    rows = [
        {
            "lineage_group": group,
            "outcome": {
                "policy": {"recovered": component["id"] in policy},
                "uniform": {"recovered": component["id"] in baseline},
            },
        }
        for policy, baseline, construct, group in zip(first, second, truth, group_ids, strict=True)
        for component in construct["components"]
    ]
    effect = paired_effect(rows, "policy", "uniform", family_bootstrap=False)
    return {key: effect[key] for key in (
        "pooled_gain", "equal_group_gain", "components_only_a", "components_only_b", "group_bootstrap",
    )}


def prefix_oracle(pools, outcomes, limit, time_limit=10):
    """Find a label-informed ceiling at both hard budgets for one objective."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix

    options, rows, columns, values, scores = [], [], [], [], []
    n = len(pools)
    for i, pool in enumerate(pools):
        for prefix in range(min(1, len(pool)), min(4, len(pool)) + 1):
            column = len(options)
            options.append((i, prefix))
            rows.extend((i, n, n + 1))
            columns.extend((column,) * 3)
            values.extend((1, prefix, sum(row.residues for row in pool[:prefix])))
            scores.append(len(outcomes[prefix][i]) if prefix else 0)
    matrix = coo_matrix((values, (rows, columns)), shape=(n + 2, len(scores))).tocsc()
    solver = milp(
        -np.asarray(scores, dtype=float),
        integrality=np.ones(len(scores)),
        bounds=Bounds(0, 1),
        constraints=LinearConstraint(
            matrix, [1] * n + [0, 0], [1] * n + [limit.candidates, limit.residues],
        ),
        options={"time_limit": time_limit, "mip_rel_gap": 0.0},
    )
    result = {
        "optimal": solver.status == 0,
        "objective_upper_bound": (
            -float(solver.mip_dual_bound) if getattr(solver, "mip_dual_bound", None) is not None else None
        ),
    }
    if solver.x is None:
        return result, None
    allocation = {i: prefix for (i, prefix), value in zip(options, solver.x, strict=True) if value > 0.5}
    if len(allocation) != n:
        raise ValueError("Oracle did not choose an admissible prefix for every construct")
    prefixes = [allocation[i] for i in range(n)]
    if not allocation_cost(pools, prefixes).fits(limit):
        raise ValueError("Oracle exceeds the original budget")
    return result, prefixes


def allocation_summary(pools, allocation):
    return {
        **vars(allocation_cost(pools, allocation)),
        "prefix_counts": dict(sorted(Counter(allocation).items())),
        "sum_squared_residues": sum(
            row.residues ** 2 for pool, count in zip(pools, allocation, strict=True) for row in pool[:count]
        ),
    }


def method_report(pools, truth, groups, prefix_outcomes, oracle_seconds):
    result = {"caps": {}}
    ids = [construct["id"] for construct in truth]
    for cap in CAPS:
        uniform = [min(cap, len(pool)) for pool in pools]
        limit = allocation_cost(pools, uniform)
        baseline = {metric: prefix_outcomes[metric][cap] for metric in METRICS}
        result["caps"][cap] = {
            "uniform": {
                "cost": allocation_summary(pools, uniform),
                "outcomes": {metric: outcome_summary(rows, truth, groups) for metric, rows in baseline.items()},
            },
        }
        if cap not in (2, 3):
            continue
        policies = result["caps"][cap]["policies"] = {}
        for name in PRIORITIES:
            allocation = allocate(pools, ids, limit, name)
            selected = {
                metric: choose_prefix_outcomes(prefix_outcomes[metric], allocation)
                for metric in METRICS
            }
            policies[name] = {
                "cost": allocation_summary(pools, allocation),
                "outcomes": {metric: outcome_summary(rows, truth, groups) for metric, rows in selected.items()},
                "paired_gain": {
                    metric: paired_delta(rows, baseline[metric], truth, groups)
                    for metric, rows in selected.items()
                },
            }
        if oracle_seconds:
            oracle = result["caps"][cap]["oracle"] = {}
            for metric in METRICS:
                bound, allocation = prefix_oracle(pools, prefix_outcomes[metric], limit, oracle_seconds)
                if allocation is not None:
                    rows = choose_prefix_outcomes(prefix_outcomes[metric], allocation)
                    bound.update(cost=allocation_summary(pools, allocation), outcome=outcome_summary(rows, truth, groups))
                oracle[metric] = bound
    return result


def survey(root, oracle_seconds):
    source_hashes = {}

    def load(relative):
        payload, source_hashes[relative] = read_json(root / relative)
        return payload

    base = "results/paired-cohort-v2"
    evaluator = f"{base}/prepared/holdout/evaluator"
    truth = load(f"{evaluator}/answer-key.corroborated.json")["constructs"]
    group_rows = load(f"{evaluator}/groups.lineage.json")["constructs"]
    lineage = {row["id"]: row["group"] for row in group_rows}
    groups = [lineage[construct["id"]] for construct in truth]
    protocols = {scorer: load(f"{base}/holdout/{scorer}/protocol.json") for scorer in SCORERS}
    geometry = {(row["minimum_component_overlap"], row["prediction_limit"]) for row in protocols.values()}
    if len(geometry) != 1:
        raise ValueError("Expected identical locus and assignment rules for both scorers")
    overlap, output_limit = geometry.pop()
    suite = {}
    for method in METHODS:
        runs = {scorer: load(f"{base}/holdout/{scorer}/{method}.json") for scorer in SCORERS}
        frozen = runs[SCORERS[0]]
        if any(
            run[field] != frozen[field] for run in runs.values()
            for field in ("candidate_pools", "reference_panel_sha256", "selected_inputs_sha256")
        ):
            raise ValueError("The scorer traces use different candidate orders or reference libraries")
        pools = [
            tuple(Candidate.from_metadata(row) for row in frozen["candidate_pools"][construct["id"]])
            for construct in truth
        ]
        prefix_outcomes = {"coverage": {cap: coverage(pools, truth, overlap, cap) for cap in CAPS}}
        for scorer, run in runs.items():
            prefix_outcomes[scorer] = parse_scored_run(run, pools, truth, output_limit, overlap)
        report = method_report(pools, truth, groups, prefix_outcomes, oracle_seconds)
        report["coverage_above_saved_scoring_cap"] = {
            cap if cap else "whole_pool": {
                "cost": allocation_summary(pools, [min(cap, len(pool)) if cap else len(pool) for pool in pools]),
                "outcome": outcome_summary(coverage(pools, truth, overlap, cap), truth, groups),
            }
            for cap in (8, 16, None)
        }
        suite[method] = report
    return {
        "schema": "preprocessing-budget-allocation-v1",
        "constructs": len(truth),
        "lineages": len(set(groups)),
        "source_sha256": source_hashes,
        "within_construct_output_limit": output_limit,
        "minimum_component_overlap": overlap,
        "methods": suite,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Repository checkout containing the traces")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--oracle-seconds", type=float, default=10, help="Time limit per oracle, 0 skips oracle")
    args = parser.parse_args()
    if not math.isfinite(args.oracle_seconds) or args.oracle_seconds < 0:
        parser.error("--oracle-seconds must be a finite nonnegative number")
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(survey(args.root, args.oracle_seconds), handle, indent=2, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
