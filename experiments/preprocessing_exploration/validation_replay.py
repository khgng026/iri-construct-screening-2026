"""Replay preprocessing policies with saved ESM-2 rankings and a fixed output budget."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from .hybrid import (
    REPLAY_CAPS, distinct_supplements, family_recovery, paired_gain, prepare, read,
    report_metrics, select_order, verify_family_replay,
)
from .ranking import rank_indices


SCORER = "esm2"
METHODS = ("minimal", "gene_caller", "mixed_ranking", "gene_then_fallback", "mixed_then_fallback")
MODES = ("candidate_cap", "gene_residue_ceiling")


def mixed_order(genes):
    return [genes[index] for index in rank_indices(genes, "keep_first_then_250")]


def evaluate(cohort):
    records, protocol, files = prepare(cohort, "holdout", scorers=(SCORER,))
    metadata = cohort / "prepared/holdout/evaluator/construct-map.private.json"
    files[str(metadata.relative_to(cohort))] = sha256(metadata.read_bytes()).hexdigest()
    tiers = {
        row["id"]: {unit["id"]: unit["similarity"]["tier"] for unit in row["components"]}
        for row in read(metadata)
    }
    overlap, output_limit = protocol["minimum_component_overlap"], protocol["prediction_limit"]
    capped_records = Counter()
    primary_fully_scored = []
    for record in records:
        verify_family_replay(record, overlap, output_limit)
        pools = record["pools"]
        full_genes = pools["gene_caller"]
        full_mixed = mixed_order(full_genes)
        genes = full_genes[:max(REPLAY_CAPS)]
        fallback = list(distinct_supplements(full_genes, pools["start_ranked"][:max(REPLAY_CAPS)], overlap))
        mixed = mixed_order(genes)
        orders = {
            "minimal": pools["minimal"][:max(REPLAY_CAPS)],
            "gene_caller": genes,
            "mixed_ranking": mixed,
            "gene_then_fallback": [*genes, *fallback],
            "mixed_then_fallback": [*mixed, *fallback],
        }
        record["similarity_tiers"] = {cid: tiers[record["truth"]["id"]][cid] for cid in record["components"]}
        for cap in REPLAY_CAPS:
            ceilings = (None, sum(candidate.residues for candidate in genes[:cap]))
            for mode, ceiling in zip(MODES, ceilings):
                full_selection = select_order(full_mixed, cap, ceiling, unique_loci=False)
                if any(candidate not in genes for candidate in full_selection):
                    capped_records[mode, cap] += 1
                elif mode == "candidate_cap" and cap == protocol["primary_query_cap"]:
                    primary_fully_scored.append(record)
                for method, order in orders.items():
                    selected = select_order(order, cap, ceiling, unique_loci=False)
                    hits, covered = family_recovery(record, selected, SCORER, cap, overlap, output_limit)
                    record["runs"][mode, cap, method] = selected, hits, covered

    denominators = Counter(tier for row in records for tier in row["similarity_tiers"].values())

    def summarize(mode, cap, method):
        key = mode, cap, method
        by_tier = Counter(
            record["similarity_tiers"][cid] for record in records for cid in record["runs"][key][1]
        )
        return {
            **report_metrics(records, key, include_coverage=True),
            "recovered_by_similarity_tier": {tier: by_tier[tier] for tier in denominators},
            "versus_minimal": paired_gain(records, key, (mode, cap, "minimal"), bootstrap=cap in (2, 4)),
            "versus_gene_caller": paired_gain(records, key, (mode, cap, "gene_caller"), bootstrap=cap in (2, 4)),
        }

    primary_key = "candidate_cap", protocol["primary_query_cap"]
    primary_mixed, primary_gene = (*primary_key, "mixed_ranking"), (*primary_key, "gene_caller")
    gained_family, lost_family = Counter(), Counter()
    for row in records:
        gained, lost = row["runs"][primary_mixed][1], row["runs"][primary_gene][1]
        families = {unit["id"]: unit["family"] for unit in row["truth"]["components"]}
        gained_family.update(families[cid] for cid in gained - lost)
        lost_family.update(families[cid] for cid in lost - gained)

    return {
        "schema": "iri-preprocessing-validation-esm2-replay-v1",
        "input_sha256": files,
        "scorer": SCORER,
        "output_limit": output_limit,
        "similarity_tier_components": dict(denominators),
        "full_mixed_would_select_unscored_records": {
            mode: {cap: capped_records[mode, cap] for cap in REPLAY_CAPS} for mode in MODES
        },
        "primary_budget_mixed_with_full_selection_scored": {
            "gene_caller": report_metrics(primary_fully_scored, primary_gene, include_coverage=True),
            "mixed_ranking": report_metrics(primary_fully_scored, primary_mixed, include_coverage=True),
        },
        "primary_budget_mixed_changes_by_family": {
            family: {"gained": gained_family[family], "lost": lost_family[family]}
            for family in gained_family.keys() | lost_family.keys()
        },
        "runs": {
            mode: {
                cap: {method: summarize(mode, cap, method) for method in METHODS}
                for cap in REPLAY_CAPS
            }
            for mode in MODES
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, default=Path("results/paired-cohort-v2"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args.cohort)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
