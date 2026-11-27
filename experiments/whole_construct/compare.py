"""Calibrate score thresholds separately, then compare paired held-out recovery."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from math import inf, isfinite, nextafter
from pathlib import Path

from experiments.paired_preprocessing.verify import _index
from experiments.preprocessing_exploration.hybrid import paired_gain

from .evaluate import match_calls, summarize


ARMS = ("blastx", "minimal", "gene_caller")


def cutoff_grid(calls):
    """Use at most 101 empirical quantiles plus a cutoff that emits no calls."""
    scores = sorted({row["score"] for group in calls.values() for row in group})
    if any(not isfinite(score) for score in scores):
        raise ValueError("Nonfinite scores are unsupported")
    if not scores:
        return [0.0]
    return sorted({scores[i * (len(scores) - 1) // 100] for i in range(101)}) + [nextafter(scores[-1], inf)]


def calibrate(development, calls, burdens, overlap, **locus_rules):
    if not development["constructs"]:
        raise ValueError("Calibration requires annotated development constructs")
    report = []
    for cutoff in cutoff_grid(calls):
        values = summarize(development, calls, cutoff, overlap, **locus_rules)
        burden = (values["unmatched_complete"] + values["unverified_partial"]) / values["constructs"]
        report.append({"cutoff": cutoff, "unsupported_per_construct": burden, **values})
    # At matched burden, keep the most permissive point supported by development.
    return {burden: next(row for row in report if row["unsupported_per_construct"] <= burden) for burden in burdens}


def _load_stage(paths):
    arms = {}
    for path in paths:
        blob = path.read_bytes()
        row = json.loads(blob)
        arm = row["arm"]
        if row["schema"] != "iri-whole-construct-calls-v1" or arm in arms:
            raise ValueError("Unexpected or repeated scoring arm")
        arms[arm] = row | {"trace_sha256": sha256(blob).hexdigest()}
    if set(arms) != set(ARMS):
        raise ValueError("Each stage must contain all three scoring arms")
    for field in ("input_sha256", "references_sha256"):
        if len({row[field] for row in arms.values()}) != 1:
            raise ValueError(f"Scoring arms disagree on {field}")
    if arms["minimal"]["work"]["model"] != arms["gene_caller"]["work"]["model"]:
        raise ValueError("ESM-2 arms disagree on the model or model configuration")
    return arms


def _family_leave_out(examples, family_names, baseline):
    effects = []
    for family in family_names:
        filtered = []
        for row in examples:
            ids = {cid for cid in row["components"] if row["families"][cid] != family}
            if ids:
                filtered.append(row | {"components": ids, "runs": {
                    arm: ((), row["runs"][arm][1] & ids, ()) for arm in ARMS
                }})
        effect = paired_gain(filtered, "gene_caller", baseline, bootstrap=False)
        if effect is not None:
            effects.append(effect["lineage_gain"])
    return [min(effects), max(effects)] if effects else None


def _identity_summary(truth, metadata, hits):
    rows = defaultdict(list)
    families = defaultdict(Counter)
    for name, construct in truth.items():
        evaluator = metadata[name]
        tiers = {unit["id"]: unit["similarity"]["tier"] for unit in evaluator["components"]}
        components = defaultdict(set)
        family_by_id = {unit["id"]: unit["family"] for unit in construct["components"]}
        for unit in construct["components"]:
            tier = tiers[unit["id"]]
            groups = ("all", tier, "mid_low") if tier in ("low", "intermediate") else ("all", tier)
            for group in groups:
                components[group].add(unit["id"])
                families[group][unit["family"]] += 1
        for group, ids in components.items():
            rows[group].append({"lineage": evaluator["lineage_group"], "components": ids,
                                "families": family_by_id,
                                "runs": {arm: ((), hits[arm][name] & ids, ()) for arm in ARMS}})
    results = {}
    for group, examples in rows.items():
        n = sum(len(row["components"]) for row in examples)
        results[group] = {
            "documented": n, "families": len(families[group]),
            "constructs": len(examples), "lineages": len({row["lineage"] for row in examples}),
            "recovered": {arm: sum(len(row["runs"][arm][1]) for row in examples) for arm in ARMS},
            "gene_caller_versus": {arm: paired_gain(examples, "gene_caller", arm, bootstrap=True)
                                   for arm in ("minimal", "blastx")},
            "gene_caller_versus_leave_one_family_out_lineage_gain_range": {
                arm: _family_leave_out(examples, families[group], arm) for arm in ("minimal", "blastx")
            },
        }
    return results


def compare(dev, holdout, metadata, dev_arms, holdout_arms, burdens, overlap, **locus_rules):
    truth = _index(holdout["constructs"], "id", "held-out answer key")
    dev_ids = _index(dev["constructs"], "id", "development answer key")
    if truth.keys() & dev_ids.keys() or dev_arms["blastx"]["input_sha256"] == holdout_arms["blastx"]["input_sha256"]:
        raise ValueError("Development and holdout overlap")
    if holdout_arms["blastx"]["references_sha256"] != dev_arms["blastx"]["references_sha256"]:
        raise ValueError("Development and holdout must use one reference library")
    if set(metadata) != truth.keys():
        raise ValueError("Held-out identity metadata differs from the answer key")
    calibrated = {arm: calibrate(dev, dev_arms[arm]["calls"], burdens, overlap, **locus_rules) for arm in ARMS}
    output = {}
    for burden in burdens:
        selected = {arm: calibrated[arm][burden] for arm in ARMS}
        metrics = {arm: summarize(holdout, holdout_arms[arm]["calls"], selected[arm]["cutoff"], overlap, **locus_rules)
                   for arm in ARMS}
        hits = {arm: {name: match_calls(construct, holdout_arms[arm]["calls"][name],
                                      selected[arm]["cutoff"], overlap, **locus_rules)[0]
                      for name, construct in truth.items()}
                for arm in ARMS}
        output[burden] = {"development": selected, "holdout": metrics,
                          "held_out_by_identity": _identity_summary(truth, metadata, hits)}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dev-key", type=Path, required=True)
    parser.add_argument("--dev-arm", type=Path, action="append", required=True)
    parser.add_argument("--holdout-key", type=Path, required=True)
    parser.add_argument("--holdout-arm", type=Path, action="append", required=True)
    parser.add_argument("--holdout-metadata", type=Path, required=True)
    parser.add_argument("--max-unsupported-per-construct", type=float, action="append", required=True)
    parser.add_argument("--overlap", type=float, default=0.8)
    parser.add_argument("--join-linked-windows", action="store_true")
    parser.add_argument("--require-frame", action="store_true", help="Check reading frame for single spans")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    burdens = sorted(set(args.max_unsupported_per_construct))
    if not 0 < args.overlap <= 1 or any(not isfinite(x) or x < 0 for x in burdens):
        parser.error("Provide valid call-burden levels and locus overlap")
    dev, holdout = _load_stage(args.dev_arm), _load_stage(args.holdout_arm)
    if dev["minimal"]["work"]["model"] != holdout["minimal"]["work"]["model"]:
        raise ValueError("Development and holdout use different ESM-2 model configurations")
    keys = {name: path.read_bytes() for name, path in (("development", args.dev_key), ("holdout", args.holdout_key),
                                                       ("metadata", args.holdout_metadata))}
    metadata = _index(json.loads(keys["metadata"]), "id", "held-out identity metadata")
    result = compare(json.loads(keys["development"]), json.loads(keys["holdout"]),
                     metadata, dev, holdout, burdens, args.overlap, join_linked=args.join_linked_windows,
                     require_frame=args.require_frame)
    output = {
        "schema": "iri-whole-construct-comparison-v1", "locus_overlap": args.overlap,
        "join_linked_windows": args.join_linked_windows,
        "require_single_span_frame": args.require_frame,
        "evaluator_sha256": {name: sha256(data).hexdigest() for name, data in keys.items()},
        "traces_sha256": {name: {arm: row["trace_sha256"] for arm, row in stage.items()}
                          for name, stage in (("development", dev), ("holdout", holdout))},
        "comparisons": result,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
