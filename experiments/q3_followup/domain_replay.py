"""Compare paired structural-domain outcomes across two recorded reference libraries."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from statistics import mean

from .evidence_audit import paired_counts, read


MODEL = "esm2_150m_local"
CONVENTIONAL = ("raw_alignment", "normalized_alignment", "fine_profiles", "superfamily_profiles")
METHODS = (*CONVENTIONAL, MODEL)
UNION = "conventional_oracle"
IDENTITY_LIMIT = 0.25


def rows_by_query(report):
    rows = {}
    for row in report["cases"]:
        if row["query"] in rows:
            raise ValueError("Duplicate query in domain report")
        if any(type(row["correct"].get(method)) is not bool for method in METHODS):
            raise ValueError("Missing correct/incorrect decision for a required method")
        rows[row["query"]] = row
    return rows


def matched_rows(original, wider):
    first, second = rows_by_query(original), rows_by_query(wider)
    if first.keys() != second.keys():
        raise ValueError("The two reference policies must cover the same domains")
    eligible = []
    for query, a in first.items():
        b = second[query]
        if any(a[key] != b[key] for key in ("fold", "fine_family", "label")):
            raise ValueError(f"Label or group changed for query {query}")
        if a["max_global_identity"] < IDENTITY_LIMIT:
            if b["max_global_identity"] >= IDENTITY_LIMIT:
                raise ValueError(f"Original low-identity query {query} no longer qualifies with wider references")
            eligible.append((a, b))
    return eligible


def decision(row):
    scores = row["correct"]
    return {**scores, UNION: any(scores[method] for method in CONVENTIONAL)}


def fold_family_accuracy(rows, method):
    families = defaultdict(lambda: defaultdict(list))
    for row in rows:
        families[row["label"]][row["fine_family"]].append(decision(row)[method])
    return mean(mean(mean(values) for values in label.values()) for label in families.values())


def summarize(rows):
    decisions = [decision(row) for row in rows]
    return {
        "domains": len(rows),
        "correct": {method: sum(row[method] for row in decisions) for method in (*METHODS, UNION)},
        "paired_150m_local": {
            method: paired_counts(decisions, MODEL, method) for method in (*CONVENTIONAL, UNION)
        },
    }


def policy_summary(rows):
    folds = defaultdict(list)
    for row in rows:
        folds[row["fold"]].append(row)
    total = summarize(rows)
    total["by_fold"] = {name: summarize(fold) for name, fold in folds.items()}
    if not folds:
        return total
    methods = (*METHODS, UNION)
    accuracies = {
        method: [fold_family_accuracy(fold, method) for fold in folds.values()] for method in methods
    }
    total["equal_fold_accuracy"] = {method: mean(values) for method, values in accuracies.items()}
    best = max(total["equal_fold_accuracy"][method] for method in CONVENTIONAL)
    total["best_conventional_equal_fold"] = [
        method for method in CONVENTIONAL if total["equal_fold_accuracy"][method] == best
    ]
    if len(folds) > 1:
        total["leave_one_fold_out_gain"] = {}
        for comparator in (*CONVENTIONAL, UNION):
            diffs = [a - b for a, b in zip(accuracies[MODEL], accuracies[comparator])]
            omitted = [(sum(diffs) - diff) / (len(diffs) - 1) for diff in diffs]
            total["leave_one_fold_out_gain"][comparator] = [min(omitted), max(omitted)]
    return total


def relative_result(row):
    flags = decision(row)
    if flags[MODEL] == flags[UNION]:
        return "both_correct" if flags[MODEL] else "both_missed"
    return "model_only" if flags[MODEL] else "conventional_only"


def replay(original, wider):
    pairs = matched_rows(original, wider)
    bands = {
        "under_25_both": pairs,
        "under_25_and_no_covered_same_label_hit_both": [
            (a, b) for a, b in pairs
            if a["covered_positive_local_identity"] is None and b["covered_positive_local_identity"] is None
        ],
    }
    results = {}
    for name, band in bands.items():
        transitions = defaultdict(lambda: defaultdict(int))
        for a, b in band:
            transitions[relative_result(a)][relative_result(b)] += 1
        results[name] = {
            "original": policy_summary([a for a, _ in band]),
            "wider": policy_summary([b for _, b in band]),
            "original_to_wider_vs_conventional_union": {
                old: dict(new_counts) for old, new_counts in transitions.items()
            },
        }
    return {"schema": "iri-q3-domain-robustness-v1", "bands": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--wider", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = replay(read(args.original), read(args.wider))
    result["source_sha256"] = {
        policy: hashlib.sha256(path.read_bytes()).hexdigest()
        for policy, path in (("original", args.original), ("wider", args.wider))
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
