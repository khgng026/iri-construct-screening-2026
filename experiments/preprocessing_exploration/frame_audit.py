"""Audit candidate coordinates in the archived sequence-free cohort traces."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from experiments.paired_preprocessing.verify import _covered, _location

ARMS = ("minimal", "start_ranked", "gene_caller")
COVERAGE = 0.8
METRICS = ("overlap", "in_frame_overlap", "contains_coding_span", "exact_coding_span")
SCHEMA = "preprocessing-frame-audit-1"


def boundary_metrics(component, candidates):
    """Compare coordinates with and without the documented terminal stop."""
    strand, ((start, end),) = component
    if end - start < 6 or (end - start) % 3:
        raise ValueError("Expected a complete single-span CDS including its terminal stop")
    translation_start = start if strand == 1 else end
    coding_start, coding_end = (start, end - 3) if strand == 1 else (start + 3, end)
    metrics = dict.fromkeys(METRICS, False)
    for candidate in candidates:
        candidate_strand, ((candidate_start, candidate_end),) = candidate
        if strand != candidate_strand:
            continue
        overlaps = _covered(component, candidate) >= COVERAGE
        metrics["overlap"] |= overlaps
        candidate_translation_start = candidate_start if strand == 1 else candidate_end
        if (translation_start - candidate_translation_start) % 3:
            continue
        metrics["in_frame_overlap"] |= overlaps
        if candidate_start <= coding_start and candidate_end >= coding_end:
            metrics["contains_coding_span"] = True
        if (candidate_start, candidate_end) == (coding_start, coding_end):
            metrics["exact_coding_span"] = True
    return metrics


def _read(path, hashes):
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    return json.loads(data)


def _candidate_positions(candidates, length):
    positions = []
    for candidate in candidates:
        position = _location(candidate, length, "candidate")
        _, spans = position
        residues = candidate.get("input_residues", candidate.get("residues"))
        if len(spans) != 1 or type(residues) is not int or residues * 3 != spans[0][1] - spans[0][0]:
            raise ValueError("Expected a contiguous translated candidate excluding the stop")
        positions.append(position)
    return positions


def _components_by_construct(key):
    eligible = {}
    excluded = Counter()
    for construct in key["constructs"]:
        cid = construct["id"]
        length = construct["length"]
        archetype = construct["archetype"]
        if cid in eligible or archetype == "overall":
            raise ValueError("Duplicate construct or reserved archetype")
        components = []
        for component in construct["components"]:
            location = _location(component, length, "documented component")
            _, spans = location
            if len(spans) != 1:
                excluded["joined"] += 1
                continue
            start, end = spans[0]
            if end - start < 6 or (end - start) % 3:
                excluded["unsupported_length"] += 1
                continue
            components.append(location)
        eligible[cid] = (length, archetype, components)
    return eligible, excluded


def _summarize(rows, eligible):
    if rows.keys() != eligible.keys():
        raise ValueError("Trace and answer key construct sets differ")
    counts = defaultdict(Counter)
    for cid, (length, archetype, components) in eligible.items():
        positions = _candidate_positions(rows[cid], length)
        for component in components:
            result = boundary_metrics(component, positions)
            for group in (archetype, "overall"):
                counts[group]["components"] += 1
                counts[group].update(result)
    fields = ("components", *METRICS)
    return {group: {field: tally[field] for field in fields} for group, tally in sorted(counts.items())}


def audit(cohort):
    hashes = {}
    key = _read(cohort / "prepared/holdout/evaluator/answer-key.corroborated.json", hashes)
    eligible, excluded = _components_by_construct(key)
    arms = {}
    for arm in ARMS:
        trace = _read(cohort / "holdout/esm2" / f"{arm}.json", hashes)
        if trace["method"] != arm:
            raise ValueError(f"Unexpected trace method for {arm}")
        stages = {"pool": _summarize(trace["candidate_pools"], eligible)}
        for tier in trace["tiers"]:
            cap = str(tier["query_cap"])
            if cap in stages:
                raise ValueError(f"Duplicate budget for {arm}: {cap}")
            rows = {row["id"]: row["queries"] for row in tier["constructs"]}
            if len(rows) != len(tier["constructs"]):
                raise ValueError(f"Duplicate constructs for {arm} at cap {cap}")
            if any(len(row) > tier["query_cap"] for row in rows.values()):
                raise ValueError(f"Candidate limit exceeded for {arm} at cap {cap}")
            stages[cap] = _summarize(rows, eligible)
        if stages.keys() != {"1", "2", "4", "pool"}:
            raise ValueError(f"Expected the archived caps of one, two and four for {arm}")
        arms[arm] = stages
    return {
        "schema": SCHEMA,
        "minimum_overlap": COVERAGE,
        "construct_count": len(eligible),
        "excluded_components": {reason: excluded[reason] for reason in ("joined", "unsupported_length")},
        "input_sha256": hashes,
        "arms": arms,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", default=Path("results/paired-cohort-v2"), type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = json.dumps(audit(args.cohort), indent=2, sort_keys=True, allow_nan=False)
    if args.out is None:
        print(result)
    else:
        with args.out.open("x", encoding="utf-8") as stream:
            stream.write(result + "\n")


if __name__ == "__main__":
    main()
