"""Audit full-pool candidate availability and cached ranking coverage without raw sequences."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path
from statistics import mean

from experiments.preprocessing_exploration.hybrid import prepare, read, verify_family_replay
from experiments.paired_preprocessing.verify import _covered


METHODS = ("minimal", "gene_caller")
SCORER = "esm2"


def component_evidence(locus, family, pool, rankings, overlap):
    """Return geometric reach and the observed top family from the cached subset."""
    compatible = [unit for unit in pool if _covered(locus, unit.locus) >= overlap]
    reciprocal = [unit for unit in compatible if _covered(unit.locus, locus) >= overlap]
    strand, spans = locus
    contiguous = len(spans) == 1
    frame = spans[0][0 if strand == 1 else 1] % 3 if contiguous else None
    scored = [unit for unit in compatible if unit in rankings]
    supported = any(rankings[unit][:1] == [family] for unit in scored)
    unknown = any(unit not in rankings for unit in compatible)
    return {
        "covering_candidate": bool(compatible),
        "covering_candidate_same_frame": contiguous and any(unit.frame == frame for unit in compatible),
        "reciprocal_covering_candidate_same_frame": contiguous and any(unit.frame == frame for unit in reciprocal),
        "single_span_component": contiguous,
        "cached_covering_candidate": bool(scored),
        "top_family_in_cached_covering_candidate": supported,
        "has_unscored_covering_candidate": unknown,
        "possible_top_family_in_full_pool": supported or unknown,
    }


def minimal_score_index(pools, rankings):
    """Use cached start-ranked scores for the same untrimmed translation candidates."""
    def key(unit):
        return unit.signature, unit.residues

    if Counter(map(key, pools["minimal"])) != Counter(map(key, pools["start_ranked"])):
        raise ValueError("The archived ranking arms do not contain the same candidate pool")
    by_locus = {}
    for unit, labels in rankings.items():
        if unit.source not in ("minimal", "start_ranked"):
            continue
        previous = by_locus.setdefault(key(unit), labels)
        if previous != labels:
            raise ValueError("Cached ranks disagree on identical archived candidates")
    return {unit: by_locus[key(unit)] for unit in pools["minimal"] if key(unit) in by_locus}


def audit(cohort):
    records, protocol, files = prepare(cohort, "holdout", scorers=(SCORER,))
    overlap = protocol["minimum_component_overlap"]
    metadata_path = cohort / "prepared/holdout/evaluator/construct-map.private.json"
    files[str(metadata_path.relative_to(cohort))] = sha256(metadata_path.read_bytes()).hexdigest()
    metadata = {row["id"]: row for row in read(metadata_path)}
    counts = defaultdict(Counter)
    work = defaultdict(Counter)
    strata = defaultdict(lambda: defaultdict(set))
    family_counts = defaultdict(Counter)
    reach_comparison = defaultdict(Counter)
    conservative_counts = defaultdict(Counter)
    bounds = defaultdict(Counter)
    group_bounds = defaultdict(lambda: defaultdict(Counter))
    for record in records:
        verify_family_replay(record, overlap, protocol["prediction_limit"])
        saved = record["ranks"][SCORER]
        rankings_by_method = {"minimal": minimal_score_index(record["pools"], saved), "gene_caller": saved}
        truth = record["truth"]
        record_metadata = metadata[truth["id"]]
        components = {unit["id"]: unit for unit in record_metadata["components"]}
        for method in METHODS:
            pool, rankings = record["pools"][method], rankings_by_method[method]
            work[method].update({
                "constructs": 1,
                "constructs_with_all_rankings": all(unit in rankings for unit in pool),
                "candidate_loci": len(pool),
                "candidate_residues": sum(unit.residues for unit in pool),
                "cached_candidate_loci": sum(unit in rankings for unit in pool),
                "cached_candidate_residues": sum(unit.residues for unit in pool if unit in rankings),
            })

        for component in truth["components"]:
            cid, family = component["id"], component["family"]
            meta = components[cid]
            tier = meta["similarity"]["tier"]
            evidence_by_method = {
                method: component_evidence(record["components"][cid], family, record["pools"][method],
                                           rankings_by_method[method], overlap)
                for method in METHODS
            }
            for group in ("all", tier, "mid_low") if tier in ("low", "intermediate") else ("all", tier):
                strata[group]["constructs"].add(truth["id"])
                strata[group]["submitted_lineages"].add(record["lineage"])
                strata[group]["conservative_groups"].add(record_metadata["conservative_group"])
                strata[group]["unique_protein_hashes"].add(meta["protein_sha256"])
                family_counts[group][family] += 1
                conservative_counts[group][record_metadata["conservative_group"]] += 1
                reached = tuple(evidence_by_method[method]["covering_candidate"] for method in METHODS)
                outcome = {(True, True): "both", (True, False): "minimal_only",
                           (False, True): "gene_only", (False, False): "neither"}[reached]
                reach_comparison[group][outcome] += 1
                for method, evidence in evidence_by_method.items():
                    counts[group, method].update(components=1, **evidence)
                minimal, genes = (evidence_by_method[method] for method in METHODS)
                limits = {
                    "lower": (int(genes["top_family_in_cached_covering_candidate"])
                              - int(minimal["possible_top_family_in_full_pool"])),
                    "upper": (int(genes["possible_top_family_in_full_pool"])
                              - int(minimal["top_family_in_cached_covering_candidate"])),
                }
                bounds[group].update(limits)
                for key, item_id in (("protein", meta["protein_sha256"]), ("family", family)):
                    group_bounds[group, key][item_id].update(components=1, **limits)

    return {
        "schema": "iri-whole-construct-availability-v1",
        "minimum_same_strand_overlap": overlap,
        "input_sha256": files,
        "candidate_work_before_cross_construct_caching": work,
        "strata": {
            group: {
                "evidence": {name: len(values) for name, values in populations.items()},
                "families": len(family_counts[group]),
                "largest_family_occurrences": max(family_counts[group].values()),
                "top_two_family_occurrences": sum(n for _, n in family_counts[group].most_common(2)),
                "largest_conservative_group_occurrences": max(conservative_counts[group].values()),
                "paired_geometric_reach": reach_comparison[group],
                "unfiltered_per_component_top_family_count_difference_bounds_gene_minus_minimal": bounds[group],
                "unfiltered_top_family_reweighted_fraction_difference_bounds": {
                    key: {side: mean(item[side] / item["components"] for item in group_bounds[group, key].values())
                          for side in ("lower", "upper")}
                    for key in ("protein", "family")
                },
                "methods": {method: counts[group, method] for method in METHODS},
            }
            for group, populations in strata.items()
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, default=Path("results/paired-cohort-v2"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.cohort)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(result, output, indent=2, sort_keys=True)
        output.write("\n")


if __name__ == "__main__":
    main()
