"""Replay the Q3 feasibility audit from the committed sequence-free exports."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


COHORT = Path("results/paired-cohort-v2/prepared")
PILOTS = Path("experiments/paired_preprocessing")
CONVENTIONAL = ("alignment", "profile")
MODELS = ("esm2_8m_local", "esm2_150m_local")
METHODS = (*CONVENTIONAL, *MODELS)
PRIMARY_PANELS = ("alternate", "renilla", "gus")


def read(path):
    return json.loads(path.read_bytes())


def component_rows(root, split):
    directory = root / COHORT / split / "evaluator"
    documented = read(directory / "answer-key.corroborated.json")["constructs"]
    metadata = {row["id"]: row for row in read(directory / "construct-map.private.json")}
    rows = []
    for construct in documented:
        parts = {part["id"]: part for part in metadata[construct["id"]]["components"]}
        for component in construct["components"]:
            part = parts[component["id"]]
            if component["family"] != part["family"]:
                raise ValueError("Answer key and metadata disagree on a documented family")
            rows.append({
                **part,
                "construct": construct["id"],
                "archetype": construct["archetype"],
                "lineage_group": metadata[construct["id"]]["lineage_group"],
                "conservative_group": metadata[construct["id"]]["conservative_group"],
            })
    return rows


def group_summary(counts):
    total = counts.total()
    return {
        "groups": len(counts),
        "largest_group_size": max(counts.values(), default=0),
        "occurrence_weighted_effective_groups": (
            total ** 2 / sum(count ** 2 for count in counts.values()) if total else None
        ),
    }


def cohort_summary(rows, development_hashes):
    hashes = Counter(row["protein_sha256"] for row in rows)
    return {
        "occurrences": len(rows),
        "constructs": len({row["construct"] for row in rows}),
        "distinct_proteins": len(hashes),
        "proteins_also_in_development": len(hashes.keys() & development_hashes),
        "occurrences_also_in_development": sum(hashes[protein] for protein in hashes.keys() & development_hashes),
        "lineage": group_summary(Counter(row["lineage_group"] for row in rows)),
        "conservative": group_summary(Counter(row["conservative_group"] for row in rows)),
        "families": dict(Counter(row["family"] for row in rows)),
    }


def cohort_report(root):
    development_hashes = {row["protein_sha256"] for row in component_rows(root, "development")}
    holdout = component_rows(root, "holdout")
    low = [row for row in holdout if row["similarity"]["tier"] == "low"]
    bands = {
        "all": holdout,
        "original_low": low,
        "original_low_under_25": [row for row in low if row["similarity"]["best_same_family_identity"] < 0.25],
        "insufficient_coverage": [
            row for row in holdout if row["similarity"]["tier"] == "insufficient_coverage"
        ],
    }
    archetypes = sorted({row["archetype"] for row in holdout})
    return {
        name: {
            "all_archetypes": cohort_summary(rows, development_hashes),
            "by_archetype": {
                archetype: cohort_summary(
                    [row for row in rows if row["archetype"] == archetype], development_hashes
                )
                for archetype in archetypes
            },
        }
        for name, rows in bands.items()
    }


def paired_counts(decisions, a, b):
    table = Counter((row[a], row[b]) for row in decisions)
    return {
        "both_correct": table[True, True],
        "only_model": table[True, False],
        "only_comparator": table[False, True],
        "neither_correct": table[False, False],
    }


def decision_summary(decisions):
    if any(any(type(row.get(method)) is not bool for method in METHODS) for row in decisions):
        raise ValueError("Expected an explicit correct/incorrect decision for every scorer")
    return {
        "correct": {method: sum(row[method] for row in decisions) for method in METHODS},
        "both_conventional_miss": sum(not any(row[method] for method in CONVENTIONAL) for row in decisions),
        "model_correct_on_conventional_misses": {
            model: sum(row[model] and not any(row[method] for method in CONVENTIONAL) for row in decisions)
            for model in MODELS
        },
        "paired": {
            f"{model}_vs_{comparator}": paired_counts(decisions, model, comparator)
            for model in MODELS
            for comparator in CONVENTIONAL
        },
    }


def pilot_report(report):
    targets = report["targets"]
    ids = sorted(
        key for key, target in targets.items()
        if target["max_same_family_global_identity"] < 0.25
        and target["max_covered_same_family_local_identity"] is None
    )
    cap = str(report["primary_cap"])
    result = {
        "targets": len(ids),
        "ancestries": sorted({targets[key]["ancestry"] for key in ids}),
        "known_boundaries": decision_summary([report["known_boundaries_correct"][key] for key in ids]),
        "whole_dna": {},
    }
    for arm, tiers in report["construct_results"].items():
        rows = [tiers[cap]["cases"][key] for key in ids]
        if any(row["selected_covering"] and not row["pool_covering"] for row in rows):
            raise ValueError("A selected covering candidate is missing from the candidate pool")
        result["whole_dna"][arm] = {
            "pool_covering": sum(row["pool_covering"] for row in rows),
            "selected_covering": sum(row["selected_covering"] for row in rows),
            "selected_exact": sum(row["selected_exact_protein"] for row in rows),
            **decision_summary([row["correct"] for row in rows]),
        }
    return result


def zero_event_upper_bound(n, alpha=0.05):
    """One-sided exact upper limit with zero successes in n independent Bernoulli trials."""
    return 1 - alpha ** (1 / n)


def audit(root):
    pilot_path = PILOTS / "distant_reporter_results.json"
    gus_path = PILOTS / "gus_context_results.json"
    reports = read(root / pilot_path)["reports"]
    reports["gus"] = read(root / gus_path)
    pilots = {name: pilot_report(report) for name, report in reports.items()}
    ancestry = sorted({name for panel in PRIMARY_PANELS for name in pilots[panel]["ancestries"]})
    input_paths = [
        COHORT / split / "evaluator" / name
        for split in ("development", "holdout")
        for name in ("answer-key.corroborated.json", "construct-map.private.json")
    ] + [pilot_path, gus_path]
    return {
        "schema": "iri-q3-evidence-audit-v1",
        "sources": {
            path.as_posix(): hashlib.sha256((root / path).read_bytes()).hexdigest() for path in input_paths
        },
        "cohort": cohort_report(root),
        "pilots": pilots,
        "pilot_primary_unique_ancestries": ancestry,
        "illustrative_zero_event_upper_bound": {
            str(n): zero_event_upper_bound(n) for n in (3, 14, 30, 59, 60, 100)
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(audit(args.root), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
