"""Recount aggregate outcomes from pinned repository records with no model or network."""
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
REV = "43c5c1f4a0599700305de981bb78df09c656f754"


def archived(path):
    return json.loads(subprocess.check_output(["git", "show", f"{REV}:{path}"], cwd=ROOT))


def audit():
    cohort = "results/paired-cohort-v2/"
    scores = archived(cohort + "holdout/analysis.corroborated.verifier.json")
    mapping = archived(cohort + "prepared/holdout/evaluator/construct-map.private.json")
    source = archived("results/paired-corpus-batch-v1/audit.private.json")
    cds = {r["accession"]: {c["feature_index"]: c for c in r["cds"]} for r in source["constructs"]}
    runs = {a: {r["construct_id"]: r for r in v["2"]["per_construct"]}
            for a, v in scores["blosum62"]["runs"].items()}
    counts = {"evaluated": 0, "full_keys": 0, "full_key_gene_complete": 0,
              "full_key_start_complete": 0, "gene_only": 0, "start_only": 0,
              "gene_only_both_two_and_gene_shorter": 0}
    totals = {a: 0 for a in runs}
    ids = {r["id"] for r in mapping}
    assert len(ids) == len(mapping) == 777 and all(set(v) == ids for v in runs.values())
    for r in mapping:
        features = cds[r["accession"]]
        components = [c for c in r["components"] if c["evidence_tier"] == "corroborated"]
        if not components:
            continue
        used = set()
        for c in components:
            index = int(c["id"].split("-", 1)[0])
            assert features[index]["protein_sha256"] == c["protein_sha256"]
            used.add(index)
        rows = {a: v[r["id"]] for a, v in runs.items()}
        complete = {a: v["recovered_components"] == len(components) for a, v in rows.items()}
        for a in totals:
            totals[a] += complete[a]
        g, s = complete["gene_caller"], complete["start_ranked"]
        full = used == features.keys()
        counts["evaluated"] += 1
        counts["full_keys"] += full
        counts["full_key_gene_complete"] += full and g
        counts["full_key_start_complete"] += full and s
        counts["gene_only"] += g and not s
        counts["start_only"] += s and not g
        gene, start = rows["gene_caller"], rows["start_ranked"]
        counts["gene_only_both_two_and_gene_shorter"] += (
            g and not s and gene["queries"] == start["queries"] == 2
            and gene["input_residues"] < start["input_residues"]
        )
    return {"complete_inventories": totals, **counts}


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
