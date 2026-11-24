"""Audit the anonymous, already frozen predictions with the Python standard library."""

import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "sources"


def credited(calls, sites, total):
    """Maximum one-to-one matching of archived call positions to documented occurrences."""
    allowed = {}
    for site, features in enumerate(sites):
        for label, item in features:
            assert 0 <= item < total
            allowed.setdefault((site, label), []).append(item)
    assigned = {}

    def place(call, visited):
        for item in allowed.get(call, []):
            if item in visited:
                continue
            visited.add(item)
            if item not in assigned or place(assigned[item], visited):
                assigned[item] = call
                return True
        return False

    for call in calls:
        place(tuple(call), set())
    return len(assigned)


def recount(name, seal):
    payload = gzip.decompress((SOURCE / f"{name}-sealed-calls.json.gz").read_bytes())
    digest = hashlib.sha256(payload).hexdigest()
    assert digest == seal["prediction_sha256"]
    predictions = json.loads(payload)
    with gzip.open(SOURCE / f"{name}-keys.json.gz", "rt") as stream:
        keys = json.load(stream)
    counts = [row["n"] for row in keys["rows"]]
    assert len(counts) == seal["records"] and len(predictions) == seal["arms"]
    assert set(predictions) == set(keys["arms"])
    summary = {}
    earlier_counts = 0
    for arm, methods in predictions.items():
        assert set(methods) == set(seal["methods"])
        key = keys["arms"][arm]
        summary[arm] = {}
        for method, rows in methods.items():
            assert len(rows) == len(counts) == len(key["edges"])
            for calls in rows:
                assert len(calls) == len({tuple(c) for c in calls}) <= seal["largest_return_budget"]
                assert all(len(c) == 2 and all(type(i) is int and i >= 0 for i in c) for c in calls)
            by_budget = {}
            for budget in (5, 10):
                observed = [
                    credited(calls[:budget], sites, n)
                    for calls, sites, n in zip(rows, key["edges"], counts, strict=True)
                ]
                if budget == 5 and method in key["archived"]:
                    assert observed == key["archived"][method]
                    earlier_counts += len(observed)
                by_budget[str(budget)] = {
                    "credits": sum(observed),
                    "complete": sum(n > 0 and h == n for n, h in zip(counts, observed, strict=True)),
                }
            summary[arm][method] = by_budget
    return {
        "prediction_sha256": digest,
        "records": len(counts),
        "evaluable_records": sum(n > 0 for n in counts),
        "documented_occurrences": sum(counts),
        "ranked_lists_from_sealed_bytes": len(counts) * len(predictions) * len(seal["methods"]),
        "historical_per_record_counts_matched": earlier_counts,
        "outcomes": summary,
    }


if __name__ == "__main__":
    seals = json.loads((SOURCE / "prediction-seals.json").read_text())["partitions"]
    results = {k: recount(k, v) for k, v in seals.items()}
    archive = gzip.decompress((SOURCE / "committed-aggregate-results.json.gz").read_bytes())
    receipt = json.loads((SOURCE / "original-checkout-replay.json").read_text())
    assert hashlib.sha256(archive).hexdigest() == receipt["numerical_result_sha256"]
    original = json.loads(archive)
    for partition, result in results.items():
        matched = set()
        for budget, table in original[partition]["budgets"].items():
            for row in table["methods"]:
                arm, method = row["arm"], row["method"]
                actual = result["outcomes"][arm][method][budget]
                assert actual == {"credits": row["credited"], "complete": row["complete"]}
                assert row["plasmids"] == result["evaluable_records"]
                assert row["components"] == result["documented_occurrences"]
                matched.add((arm, method, budget))
        assert len(matched) == len(seals[partition]["methods"]) * seals[partition]["arms"] * 2
        result["original_aggregate_rows_matched"] = len(matched)
    (ROOT / "sealed-outcome-results.json").write_text(json.dumps(results, indent=2) + "\n")
    for partition, result in results.items():
        print(partition, result["ranked_lists_from_sealed_bytes"], result["historical_per_record_counts_matched"])
