"""Offline numerical replay of frozen conventional optical and yeast controls."""
import argparse
from bisect import bisect_right
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path


FOLDER = Path(__file__).resolve().parent
MANIFEST_SHA256 = "1247ffc46470d47db5d8d8c7ffcd7304a4cadba75351f19304bc1c027b6dfda3"
ARMS = ("blastx", "minimal", "gene_caller")
ALLOWANCES = ("0", "1", "2")


def verified_file(name, expected):
    data = (FOLDER / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Digest mismatch for {name}")
    return data


class CallPrefix:
    def __init__(self, nodes, n_components, named):
        ranks = [node[0] for node in nodes]
        assert ranks == sorted(ranks, reverse=True)
        retained = set()
        assignments = {0}
        self.negative_ranks = []
        self.results = [(0, 0, False)]
        for i, (rank, eligible, duplicates) in enumerate(nodes):
            assert 0 <= eligible < 1 << n_components
            assert all(0 <= earlier < i for earlier in duplicates)
            if not retained.isdisjoint(duplicates):
                continue
            retained.add(i)
            possible = {1 << index for index in range(n_components) if eligible & (1 << index)}
            # Distinct retained calls may each cover one distinct documented component.
            assignments |= {taken | feature for taken in assignments for feature in possible
                            if not taken & feature}
            best = max(mask.bit_count() for mask in assignments)
            got_named = named is not None and any(mask & (1 << named) for mask in assignments)
            self.negative_ranks.append(-rank)
            self.results.append((len(retained), best, got_named))

    def at(self, min_rank):
        return self.results[bisect_right(self.negative_ranks, -min_rank)]


def select_operating_points(rows, counts, training):
    ranks = sorted({node[0] for i in training for node in rows[i]})
    assert training and ranks
    quantiles = sorted({ranks[q * (len(ranks) - 1) // 100] for q in range(101)})
    grid = [(f"ge{r}", r) for r in quantiles]
    grid.append((f"gt{ranks[-1]}", ranks[-1] + 1))
    chosen = {}
    for token, min_rank in grid:
        observed = [counts[i].at(min_rank) for i in training]
        uncredited = sum(emitted - covered for emitted, covered, _ in observed)
        for n in ALLOWANCES:
            if n not in chosen and uncredited <= int(n) * len(training):
                chosen[n] = (token, min_rank, uncredited)
        if len(chosen) == len(ALLOWANCES):
            break
    assert len(chosen) == len(ALLOWANCES)
    return chosen


def summarize(records, decisions, indices):
    calls = [decisions[i][0] for i in indices]
    hits = [decisions[i][1] for i in indices]
    documented = [records[i]["components"] for i in indices]
    return {
        "stocks": len(indices), "documented": sum(documented), "recovered": sum(hits),
        "completed": sum(key == hit for key, hit in zip(documented, hits)),
        "exact": sum(key == hit == ret for key, hit, ret in zip(documented, hits, calls)),
        "returned": sum(calls), "uncredited": sum(calls) - sum(hits), "call_range": [min(calls), max(calls)],
    }


def rebuild(data):
    assert data["schema"] == 1 and set(data["arms"]) == set(ARMS)
    records = data["records"]
    by_role = {role: [i for i, r in enumerate(records) if r["role"] == role]
               for role in ("closed", "partial", "colony", "marker")}
    closed, colonies, yeast = by_role["closed"], by_role["colony"], by_role["marker"]
    named = [i for i, record in enumerate(records) if record["named"] is not None]
    assert [len(closed), len(by_role["partial"]), len(colonies), len(yeast), len(named)] == [11, 2, 4, 4, 7]
    assert sum(records[i]["components"] for i in closed) == 26
    assert all(records[i]["components"] == 1 for i in yeast)
    groups = {records[i]["fold"] for i in closed}
    assert len(groups) == 4
    training = {g: [i for i in closed if records[i]["fold"] != g] for g in groups}
    training["ALL"] = closed
    assert all(i not in training[r["fold"]] for i, r in enumerate(records))
    cutoffs = {n: {a: {} for a in ARMS} for n in ALLOWANCES}
    results = {n: {k: {} for k in ("primary", "blue_colonies", "named_reporter", "external_markers", "paired_primary")}
               for n in ALLOWANCES}
    complete = {n: {} for n in ALLOWANCES}
    for arm in ARMS:
        info = data["arms"][arm]
        rows = info["calls"]
        assert len(rows) == len(records)
        assert {node[0] for rr in rows for node in rr} == set(range(info["rank_count"]))
        profiles = [CallPrefix(nodes, record["components"], record["named"]) for nodes, record in zip(rows, records)]
        points = {g: select_operating_points(rows, profiles, ids) for g, ids in training.items()}
        for n in ALLOWANCES:
            for fold, ids in training.items():
                token, _, extra = points[fold][n]
                # Magnitudes enter only after the operating rule has selected an ordinal.
                score = info["reported_threshold_values"][token]
                cutoffs[n][arm][fold] = {"cutoff": score, "train_n": len(ids), "training_burden": extra / len(ids)}
            decisions = [profile.at(points[r["fold"]][n][1]) for r, profile in zip(records, profiles)]
            full = {i for i in closed if decisions[i][1] == records[i]["components"]}
            complete[n][arm] = full
            matched_series = []
            for g in groups:
                members = {i for i in closed if records[i]["fold"] == g}
                matched_series.append(Fraction(len(full & members), len(members)))
            row = results[n]
            row["primary"][arm] = summarize(records, decisions, closed)
            row["primary"][arm]["series_equal_completion"] = float(sum(matched_series) / len(matched_series))
            row["blue_colonies"][arm] = summarize(records, decisions, colonies)
            row["named_reporter"][arm] = {"stocks": len(named), "recovered": sum(decisions[i][2] for i in named),
                                          "returned_all_families": sum(decisions[i][0] for i in named)}
            marker = summarize(records, decisions, yeast)
            row["external_markers"][arm] = {"positive_loci_tested": len(yeast), "positive_markers": marker["recovered"],
                                            "all_calls": marker["returned"], "other_calls_unassessed": marker["uncredited"]}
    for n in ALLOWANCES:
        for a, b in combinations(ARMS, 2):
            first, second = complete[n][a], complete[n][b]
            results[n]["paired_primary"][a + "_vs_" + b] = {
                "both": len(first & second), "first_only": len(first - second),
                "second_only": len(second - first), "neither": len(set(closed) - first - second),
            }
    return {"cutoffs": cutoffs, "allowances": results}


def run():
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled.")
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--json", action="store_true", help="Print every verified aggregate and calibration.")
    config = args.parse_args()
    sources = json.loads(verified_file("provenance.json", MANIFEST_SHA256))
    data = json.loads(gzip.decompress(verified_file("evidence.json.gz", sources["evidence_sha256"])))
    result = rebuild(data)
    expected = json.loads(verified_file("expected.json", sources["expected_sha256"]))
    if result != expected:
        raise ValueError("Recalculated results differ from the recorded source tables.")
    if config.json:
        print(json.dumps(result, sort_keys=True, indent=2))
        return
    print("Verified all 45 training cutoffs and all published aggregate cells.")
    print("Allowance   Method        All keys /11   All & no extras /11   Reporter /7   Colonies /4   Yeast /4")
    for n, cells in result["allowances"].items():
        for arm in ARMS:
            primary = cells["primary"][arm]
            print(f"{n:9}   {arm:11} {primary['completed']:12} {primary['exact']:20} "
                  f"{cells['named_reporter'][arm]['recovered']:14} {cells['blue_colonies'][arm]['completed']:13} "
                  f"{cells['external_markers'][arm]['positive_markers']:9}")


if __name__ == "__main__":
    run()
