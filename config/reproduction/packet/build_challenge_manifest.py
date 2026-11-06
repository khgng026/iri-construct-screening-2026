"""Prepare deterministic metadata-only challenges. No model inference or downloads."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent / "iri-construct-screening-2026"
SEED = "iri-overnight-2026-10-31-v1"
CAP = 1022


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_key(values, salt=""):
    return hashlib.sha256((SEED + "|" + salt + "|" + "|".join(values)).encode()).hexdigest()


def read_json(relative):
    return json.loads((REPO / relative).read_text())


def write_json(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    before = __import__("subprocess").check_output(
        ["git", "status", "--porcelain=v1"], cwd=REPO, text=True
    )
    previous_freeze = read_json("config/t05_expanded_freeze.json")
    for relative, expected in previous_freeze["file_sha256"].items():
        assert digest(REPO / relative) == expected, f"Original freeze mismatch: {relative}"
    catalogue = read_json("data/manifests/t05_expanded_catalogue.json")
    original = read_json("config/t05_expanded_queries.json")
    records = {r["accession"]: r for r in catalogue["records"]}
    targets = sorted(original["targets"])
    lengths, hashes = {}, {}
    for accession, r in records.items():
        assert all(r[k] == "approved" for k in
                   ("review_status", "licensing_status", "safety_review"))
        reference = (REPO / r["sequence_ref"]["path"]).resolve()
        assert reference.is_relative_to(REPO.resolve())
        lines = reference.read_text().splitlines()
        assert sum(line.startswith(">") for line in lines) == 1
        sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
        assert sequence and set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY")
        hashes[accession] = hashlib.sha256(sequence.encode()).hexdigest()
        lengths[accession] = len(sequence)
        assert hashes[accession] == r["sequence_sha256"]
        assert lengths[accession] == r["sequence_length_aa"]
    assert len({hashes[t] for t in targets}) == len(targets)

    queries = []

    def add_query(components, query_id, suite, group_id):
        start = 0
        spans = []
        total = sum(lengths[a] for a in components)
        assert total <= CAP and len(set(components)) == len(components)
        for accession in components:
            end = start + lengths[accession]
            spans.append({"accession": accession, "start0": start, "end0": end,
                          "fraction": lengths[accession] / total,
                          "sequence_sha256": hashes[accession]})
            start = end
        queries.append({"query_id": query_id, "suite": suite, "group_id": group_id,
                        "components": list(components), "length_aa": total,
                        "component_count": len(components), "spans": spans,
                        "artificial": len(components) > 1})

    for q in original["queries"]:
        add_query(q["components"], q["query_id"],
                  "clean" if len(q["components"]) == 1 else "pairs",
                  q.get("pair_id") or q["query_id"])

    def fraction_bin(value):
        if value < .10:
            return "lt_0.10"
        if value < .20:
            return "0.10_to_lt_0.20"
        return "ge_0.20"

    selection = {}
    for size, budget in ((3, 90), (4, 45)):
        bins = {name: [] for name in ("lt_0.10", "0.10_to_lt_0.20", "ge_0.20")}
        excluded = []
        for components in itertools.combinations(targets, size):
            total = sum(lengths[a] for a in components)
            group = {"group_id": "__".join(components), "components": list(components),
                     "length_aa": total,
                     "minimum_fraction": min(lengths[a] for a in components) / total}
            if total > CAP:
                excluded.append({**group, "reason": "shared_length_cap"})
            else:
                bins[fraction_bin(group["minimum_fraction"])].append(group)
        for name, groups in bins.items():
            groups.sort(key=lambda g: stable_key(g["components"], f"count{size}:{name}"))
        # Round-robin keeps available bins balanced; empty bins do not block selection.
        selected = []
        remaining = {name: list(groups) for name, groups in bins.items()}
        while len(selected) < budget and any(remaining.values()):
            for name in sorted(remaining):
                if remaining[name] and len(selected) < budget:
                    selected.append(remaining[name].pop(0))
        for group in selected:
            components = tuple(group["components"])
            if size == 3:
                orders = sorted(set(itertools.permutations(components)))
            else:
                forward = [components[i:] + components[:i] for i in range(size)]
                orders = sorted(set(forward + [tuple(reversed(order)) for order in forward]))
            for order in orders:
                add_query(order, f"join{size}_" + "__".join(order),
                          "triples" if size == 3 else "quadruples", group["group_id"])
        selection[str(size)] = {
            "group_budget": budget, "selected_groups": selected,
            "eligible_bin_counts": {k: len(v) for k, v in bins.items()},
            "selected_bin_counts": {k: sum(fraction_bin(g["minimum_fraction"]) == k
                                           for g in selected) for k in bins},
            "eligible_omitted_group_ids": [g["group_id"] for v in remaining.values() for g in v],
            "excluded_groups": excluded,
            "order_rule": "all_permutations" if size == 3 else "cyclic_rotations_and_reversals",
        }

    assert len({q["query_id"] for q in queries}) == len(queries)
    # Sequence-equivalent groups move together in the catalogue sensitivity analysis.
    equivalence = {}
    for accession in sorted(records):
        equivalence.setdefault(hashes[accession], []).append(accession)
    target_hashes = {hashes[t] for t in targets}
    decoys = sorted(set(equivalence) - target_hashes)
    subsets = []
    for replicate in range(3):
        ordered = sorted(decoys, key=lambda h: stable_key([h], f"decoys:{replicate}"))
        for size in (25, 50, len(equivalence)):
            group_hashes = sorted(target_hashes) + ordered[:size - len(target_hashes)]
            assert len(set(group_hashes)) == size
            subsets.append({"subset_id": f"unique{size}_rep{replicate}",
                            "replicate": replicate, "unique_sequence_count": size,
                            "entries": [{"sequence_sha256": h,
                                         "representative_accession": min(equivalence[h]),
                                         "equivalent_accessions": equivalence[h]}
                                        for h in group_hashes]})

    repeat_ids = list(original["repeat_query_ids"])
    for suite in ("triples", "quadruples"):
        eligible = [q for q in queries if q["suite"] == suite]
        repeat_ids.extend(q["query_id"] for q in sorted(
            eligible, key=lambda q: stable_key([q["query_id"]], "repeat"))[:4])
    source_files = ["config/t05_expanded_queries.json", "config/t05_expanded_freeze.json",
                    "config/t05_expanded_selection.json", "config/t05_model.json",
                    "data/manifests/t05_expanded_catalogue.json",
                    "docs/week1/T05_expanded_report.md"]
    counts = {suite: sum(q["suite"] == suite for q in queries)
              for suite in ("clean", "pairs", "triples", "quadruples")}
    manifest = {
        "version": SEED, "status": "prepared_not_run", "shared_max_residues": CAP,
        "selection_seed": SEED, "source_files_sha256": {p: digest(REPO / p) for p in source_files},
        "note": "Planned after observing the prior pilot; not an external holdout. No inference in preparation.",
        "catalogue_accessions": sorted(records), "targets": targets,
        "queries": queries, "group_selection": selection,
        "original_pair_exclusions": original["excluded_pairs"],
        "catalogue_subsets": subsets, "repeat_query_ids": repeat_ids,
        "suite_counts": counts,
    }
    write_json("challenge_manifest.json", manifest)
    after = __import__("subprocess").check_output(
        ["git", "status", "--porcelain=v1"], cwd=REPO, text=True
    )
    assert before == after, "Preparation unexpectedly changed experiment repository"
    checks = {
        "status": "passed", "no_model_runs": True, "no_network_requests": True,
        "experiment_worktree_status_unchanged": True,
        "existing_freeze_hashes_verified": len(previous_freeze["file_sha256"]),
        "reviewed_records_checked": len(records), "unique_sequences": len(equivalence),
        "target_count": len(targets), "suite_counts": counts, "query_count": len(queries),
        "combined_component_observations": sum(q["component_count"] for q in queries
                                               if q["component_count"] > 1),
        "max_query_length": max(q["length_aa"] for q in queries),
        "catalogue_subset_count": len(subsets), "repeat_query_count": len(repeat_ids),
        "challenge_manifest_sha256": digest(ROOT / "challenge_manifest.json"),
        "source_files_sha256": manifest["source_files_sha256"],
        "prepared_document_sha256": {name: digest(ROOT / name)
                                    for name in ("README.md", "PROTOCOL.md", "OPERATIONS.md", "LAUNCH_BRIEF.md", "EXTENSION_QUEUE.md")},
        "builder_sha256": digest(Path(__file__)),
    }
    write_json("preparation_checks.json", checks)
    print(json.dumps({k: v for k, v in checks.items() if not isinstance(v, dict)}, indent=2))
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
