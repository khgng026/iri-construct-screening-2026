"""Replay conventional and combined family ranks on the same historical candidates."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from experiments.paired_preprocessing.verify import _covered, _location, _score_construct

from .evidence_audit import COHORT, read


SCORERS = ("blosum62", "esm2")
CAP = 2


def combine_ranks(conventional, model):
    model_positions = {label: position for position, label in enumerate(model)}
    if (
        not conventional
        or len(set(conventional)) != len(conventional)
        or len(model_positions) != len(model)
        or set(conventional) != model_positions.keys()
    ):
        raise ValueError("Expected complete ranked orders of the same families")
    return min(enumerate(conventional), key=lambda row: (row[0] + model_positions[row[1]], row[0]))[1]


def _choices(conventional, model):
    if conventional["queries"] != model["queries"]:
        raise ValueError("The two scorers must see identical candidates")
    queries = [query["id"] for query in conventional["queries"]]
    rankings = conventional["query_family_rankings"], model["query_family_rankings"]
    if any(set(queries) != ranks.keys() for ranks in rankings):
        raise ValueError("A candidate has no ranking or an unexpected ranking")
    baseline, model_ranks = rankings
    return {
        "conventional": [
            {"query_id": query, "family": baseline[query][0]} for query in queries
        ],
        "combined": [
            {"query_id": query, "family": combine_ranks(baseline[query], model_ranks[query])}
            for query in queries
        ],
    }


def _known_locus_calls(construct, trace, overlap):
    length = construct["length"]
    loci = {
        row["id"]: (_location(row, length, row["id"]), row["family"])
        for row in construct["components"]
    }
    qloc = {query["id"]: _location(query, length, query["id"]) for query in trace["queries"]}
    counts = Counter(correct=0, incompatible=0, unknown=0)
    for prediction in trace["predictions"]:
        predicted_locus = qloc[prediction["query_id"]]
        covered = [label for locus, label in loci.values() if _covered(locus, predicted_locus) >= overlap]
        if len(covered) == 1:
            category = "correct" if covered[0] == prediction["family"] else "incompatible"
        else:
            category = "unknown"
        counts[category] += 1
    return counts


def replay(root):
    from experiments.paired_preprocessing.cohort_analysis import check_matched_inputs, paired_effect

    run_dir = Path("results/paired-cohort-v2/holdout")
    protocol_paths = {scorer: run_dir / scorer / "protocol.json" for scorer in SCORERS}
    trace_paths = {scorer: run_dir / scorer / "gene_caller.json" for scorer in SCORERS}
    protocols = {scorer: read(root / path) for scorer, path in protocol_paths.items()}
    traces = {scorer: {"gene_caller": read(root / path)} for scorer, path in trace_paths.items()}
    matched = check_matched_inputs(protocols, traces, ("gene_caller",))
    transcript = {}
    for scorer in SCORERS:
        tier, = [run for run in traces[scorer]["gene_caller"]["tiers"] if run["query_cap"] == CAP]
        transcript[scorer] = {row["id"]: row for row in tier["constructs"]}
    evaluator = COHORT / "holdout" / "evaluator"
    key_path = evaluator / "answer-key.corroborated.json"
    metadata_path = evaluator / "construct-map.private.json"
    key = read(root / key_path)["constructs"]
    metadata = {row["id"]: row for row in read(root / metadata_path)}
    required_ids = {row["id"] for row in key}
    if any(rows.keys() != required_ids for rows in transcript.values()):
        raise ValueError("Construct set differs from the evaluator key")
    overlap = protocols["blosum62"]["minimum_component_overlap"]
    rows = []
    calls = {method: Counter() for method in ("conventional", "combined")}
    for construct in key:
        cid = construct["id"]
        conventional, model = (transcript[scorer][cid] for scorer in SCORERS)
        meta = metadata[cid]
        outcomes = {}
        for method, predictions in _choices(conventional, model).items():
            trace = {**conventional, "predictions": predictions}
            outcome = _score_construct(construct, trace, CAP, CAP, overlap)
            outcomes[method] = set(outcome["recovered_component_ids"])
            calls[method].update(_known_locus_calls(construct, trace, overlap))
        components = {row["id"]: row for row in meta["components"]}
        for component in construct["components"]:
            component_id = component["id"]
            rows.append({
                "family": component["family"],
                "lineage_group": meta["lineage_group"],
                "tier": components[component_id]["similarity"]["tier"],
                "outcome": {
                    method: {"recovered": component_id in recovered}
                    for method, recovered in outcomes.items()
                },
            })

    sections = {
        "all": rows,
        **{name: [row for row in rows if row["tier"] == name] for name in ("low", "insufficient_coverage")},
    }
    return {
        "schema": "iri-q3-historical-rank-replay-v1",
        "candidate_cap": CAP,
        "selected_input_sha256": matched,
        "effects_combined_minus_conventional": {
            label: paired_effect(section, "combined", "conventional", family_bootstrap=False)
            for label, section in sections.items()
        },
        "known_locus_calls": {method: dict(counts) for method, counts in calls.items()},
        "source_sha256": {
            path.as_posix(): hashlib.sha256((root / path).read_bytes()).hexdigest()
            for path in (*protocol_paths.values(), *trace_paths.values(), key_path, metadata_path)
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(replay(args.root), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
