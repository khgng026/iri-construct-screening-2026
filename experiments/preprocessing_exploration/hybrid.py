"""Compare mixed candidate pools on the sequence-free benign cohort v2 traces."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
from itertools import zip_longest
import json
from pathlib import Path
import random
from statistics import mean
from typing import Iterable

from experiments.paired_preprocessing.verify import _covered, _location, _score_construct


POOLS = ("gene_caller", "minimal", "start_ranked")
METHODS = (
    *POOLS,
    "gene_first_fill",
    "gene_first_distinct",
    "gene_first_one_distinct",
    "alternate",
    "alternate_distinct",
    "reserve_one_distinct",
)
CAPS = (1, 2, 4, 8, 16, 32)
RESIDUE_BUDGETS = (384, 768, 1536, 3072, 6144)
REPLAY_CAPS = (1, 2, 4)
SCORERS = ("esm2", "blosum62")


@dataclass(frozen=True)
class Candidate:
    source: str
    index: int
    strand: int
    start: int
    end: int
    residues: int

    @property
    def signature(self):
        return self.strand, self.start, self.end

    @property
    def frame(self):
        return (self.start if self.strand == 1 else self.end) % 3

    @property
    def locus(self):
        return self.strand, ((self.start, self.end),)

    @property
    def query_id(self):
        return f"{self.source}:{self.index}"

    def query_row(self):
        return {
            "id": self.query_id,
            "strand": self.strand,
            "spans": [[self.start, self.end]],
            "input_residues": self.residues,
        }


def parse_pool(rows, source, construct_size):
    candidates = []
    for i, row in enumerate(rows):
        strand, spans = _location(row, construct_size, f"{source} candidate {i}")
        if len(spans) != 1:
            raise ValueError("Expected a contiguous sequence-only candidate")
        (start, end), = spans
        residues = row["residues"]
        if type(residues) is not int or residues <= 0 or residues * 3 != end - start:
            raise ValueError("Candidate length does not match its reported span")
        candidates.append(Candidate(source, i, strand, start, end, residues))
    return candidates


def same_frame_overlap(a, b):
    """Fraction of the shorter span shared by two same-frame candidates."""
    if a.strand != b.strand or a.frame != b.frame:
        return 0.0
    bases = max(0, min(a.end, b.end) - max(a.start, b.start))
    return bases / min(a.end - a.start, b.end - b.start)


def interleave(a, b):
    return (candidate for pair in zip_longest(a, b) for candidate in pair if candidate is not None)


def select_order(
    order: Iterable[Candidate], calls: int | None, residues: int | None, *,
    unique_loci: bool,
) -> tuple[Candidate, ...]:
    """Keep the order supplied by the policy, skipping candidates that do not fit."""
    selected = []
    seen = set()
    used = 0
    for candidate in order:
        if calls is not None and len(selected) == calls:
            break
        if unique_loci and candidate.signature in seen:
            continue
        if residues is not None and used + candidate.residues > residues:
            continue
        selected.append(candidate)
        if unique_loci:
            seen.add(candidate.signature)
        used += candidate.residues
    return tuple(selected)


def selection_orders(pools, overlap, reserve_at):
    genes, fallback = pools["gene_caller"], pools["start_ranked"]
    distinct = [
        unit for unit in fallback
        if all(same_frame_overlap(unit, gene) < overlap for gene in genes)
    ]
    return {
        **pools,
        "gene_first_fill": [*genes, *fallback],
        "gene_first_distinct": [*genes, *distinct],
        "gene_first_one_distinct": [*genes, *distinct[:1]],
        "alternate": list(interleave(genes, fallback)),
        "alternate_distinct": list(interleave(genes, distinct)),
        "reserve_one_distinct": [*genes[:reserve_at], *distinct[:1],
                                 *genes[reserve_at:], *distinct[1:]],
    }


def covered_ids(components, selected, overlap):
    return {
        name for name, locus in components.items()
        if any(_covered(locus, candidate.locus) >= overlap for candidate in selected)
    }


def summarize(rows, key, *, include_coverage=True):
    strata = defaultdict(list)
    for record in rows:
        strata[record["archetype"]].append(record)
    values = report_metrics(rows, key, include_coverage)
    values["by_archetype"] = {
        label: report_metrics(group, key, include_coverage) for label, group in sorted(strata.items())
    }
    return values


def report_metrics(rows, key, include_coverage):
    counts = defaultdict(lambda: [0, 0, 0])
    calls = residues = 0
    for record in rows:
        chosen, hits, coverage = record["runs"][key]
        joined = counts[record["lineage"]]
        joined[0] += len(record["components"])
        joined[1] += len(hits)
        joined[2] += len(coverage)
        calls += len(chosen)
        residues += sum(candidate.residues for candidate in chosen)
    components = sum(row[0] for row in counts.values())
    hits = sum(row[1] for row in counts.values())
    result = {
        "constructs": len(rows),
        "components": components,
        "groups": len(counts),
        "scorer_inputs": calls,
        "input_residues": residues,
        "recovered": hits,
        "pooled_recall": hits / components if components else None,
        "equal_lineage_recall": mean(hit / n for n, hit, _ in counts.values() if n) if components else None,
    }
    if include_coverage:
        result["covered"] = sum(row[2] for row in counts.values())
    return result


def paired_gain(rows, key, baseline, *, bootstrap):
    groups = defaultdict(lambda: [0, 0])
    gained = lost = 0
    for row in rows:
        _, alternative_hits, _ = row["runs"][key]
        _, baseline_hits, _ = row["runs"][baseline]
        gained += len(alternative_hits - baseline_hits)
        lost += len(baseline_hits - alternative_hits)
        count = groups[row["lineage"]]
        count[0] += len(row["components"])
        count[1] += len(alternative_hits) - len(baseline_hits)
    effects = [delta / n for n, delta in groups.values() if n]
    if not effects:
        return None
    result = {
        "gained": gained,
        "lost": lost,
        "lineage_gain": mean(effects),
        "groups_improved": sum(x > 0 for x in effects),
        "groups_declined": sum(x < 0 for x in effects),
    }
    if bootstrap:
        rng = random.Random(42)
        trials = sorted(mean(rng.choices(effects, k=len(effects))) for _ in range(2000))
        result["bootstrap_95"] = [trials[49], trials[1949]]
    return result


def read(path):
    return json.loads(path.read_bytes())


def prepare(cohort, partition):
    fold = "dev-run3" if partition == "development" else "holdout"
    trace_paths = {scorer: {method: cohort / fold / scorer / f"{method}.json" for method in POOLS}
                   for scorer in SCORERS}
    traces = {scorer: {method: read(path) for method, path in paths.items()}
              for scorer, paths in trace_paths.items()}
    protocol_paths = [cohort / fold / scorer / "protocol.json" for scorer in SCORERS]
    protocols = [read(path) for path in protocol_paths]
    echovalues = ("query_caps", "prediction_limit", "minimum_component_overlap",
                  "sequence_input_sha256", "reference_panel_sha256")
    if any(protocols[0][name] != protocols[1][name] for name in echovalues):
        raise ValueError("Historical scorers used different protocols")
    protocol = protocols[0]
    evaluator = cohort / "prepared" / partition / "evaluator"
    key_path = evaluator / "answer-key.corroborated.json"
    group_path = evaluator / "groups.lineage.json"
    truth = read(key_path)["constructs"]
    lineage = {row["id"]: row["group"] for row in read(group_path)["constructs"]}
    ids = {row["id"] for row in truth}
    replay_tiers = {}
    for scorer in SCORERS:
        for method in POOLS:
            trace = traces[scorer][method]
            if trace["input_sha256"] != protocol["sequence_input_sha256"]:
                raise ValueError("Historical trace input fingerprint differs from protocol")
            if trace["reference_panel_sha256"] != protocol["reference_panel_sha256"]:
                raise ValueError("Historical reference fingerprint differs from protocol")
            tier, = [item["constructs"] for item in trace["tiers"] if item["query_cap"] == max(REPLAY_CAPS)]
            replay_tiers[scorer, method] = {row["id"]: row for row in tier}
            if set(trace["candidate_pools"]) != ids or set(replay_tiers[scorer, method]) != ids:
                raise ValueError("Trace records differ from the answer-key records")

    records = []
    for construct in truth:
        cid, size = construct["id"], construct["length"]
        pools = {}
        ranks = {scorer: {} for scorer in SCORERS}
        historical = {scorer: {} for scorer in SCORERS}
        for method in POOLS:
            first, second = (traces[scorer][method]["candidate_pools"][cid] for scorer in SCORERS)
            if first != second:
                raise ValueError("Candidate pools differ across the saved scorers")
            pools[method] = parse_pool(first, method, size)
            for scorer in SCORERS:
                tier = replay_tiers[scorer, method][cid]
                queries = tier["queries"]
                if len(queries) != min(max(REPLAY_CAPS), len(pools[method])):
                    raise ValueError("Stored scorer ranks differ from the expected pool budget")
                for candidate, query in zip(pools[method], queries):
                    if query["spans"] != [[candidate.start, candidate.end]] \
                            or query["strand"] != candidate.strand or query["input_residues"] != candidate.residues:
                        raise ValueError("Stored scorer ranks and candidate order disagree")
                    ranks[scorer][candidate] = tier["query_family_rankings"][query["id"]]
                historical[scorer][method] = tier
        components = {row["id"]: _location(row, size, row["id"]) for row in construct["components"]}
        records.append({
            "truth": construct,
            "archetype": construct["archetype"],
            "lineage": lineage[cid],
            "pools": pools,
            "components": components,
            "ranks": ranks,
            "historical": historical,
            "runs": {},
        })
    paths = [key_path, group_path, *protocol_paths, *(path for rows in trace_paths.values() for path in rows.values())]
    files = {str(path.relative_to(cohort)): sha256(path.read_bytes()).hexdigest() for path in paths}
    return records, protocol, files


def family_predictions(candidates, ranks, cap):
    rows = [(candidate.query_id, ranks[candidate]) for candidate in candidates]
    predictions = []
    for rank in range(max((len(labels) for _, labels in rows), default=0)):
        for qid, labels in rows:
            if rank < len(labels):
                predictions.append({"query_id": qid, "family": labels[rank]})
                if len(predictions) == cap:
                    return predictions
    return predictions


def family_recovery(record, candidates, scorer, cap, overlap, call_limit):
    predictions = family_predictions(candidates, record["ranks"][scorer], call_limit)
    trace = {
        "queries": [candidate.query_row() for candidate in candidates],
        "predictions": predictions,
        "candidate_count": len(candidates),
        "status": "ok",
        "preprocess_seconds": 0,
        "characterize_seconds": 0,
    }
    outcome = _score_construct(record["truth"], trace, cap, call_limit, overlap)
    return set(outcome["recovered_component_ids"]), set(outcome["covered_component_ids"])


def verify_family_replay(record, overlap, call_limit):
    """Verify that the unchanged top-four arms exactly reproduce their saved predictions."""
    for scorer in SCORERS:
        for method in POOLS:
            tier = record["historical"][scorer][method]
            candidates = record["pools"][method][:max(REPLAY_CAPS)]
            qid_map = {query["id"]: candidate.query_id for query, candidate in zip(tier["queries"], candidates)}
            old = [{**row, "query_id": qid_map[row["query_id"]]} for row in tier["predictions"]]
            new = family_predictions(candidates, record["ranks"][scorer], call_limit)
            if old != new:
                raise ValueError("Original prediction scheduling did not replay exactly")


def add_coverage(record, overlap):
    source = record["pools"]
    loci = record["components"]
    runs = record["runs"]
    for cap in CAPS:
        orders = selection_orders(source, overlap, max(1, cap - 1))
        gene_residues = sum(candidate.residues for candidate in source["gene_caller"][:cap])
        for method, order in orders.items():
            for mode, residues in (("count", None), ("gene_residues", gene_residues)):
                selected = select_order(order, cap, residues, unique_loci=method not in POOLS)
                hits = covered_ids(loci, selected, overlap)
                runs[mode, cap, method] = selected, hits, hits
    orders = selection_orders(source, overlap, 1)
    for cap in RESIDUE_BUDGETS:
        for method, order in orders.items():
            selected = select_order(order, None, cap, unique_loci=method not in POOLS)
            hits = covered_ids(loci, selected, overlap)
            runs["residues", cap, method] = selected, hits, hits
    for method, order in orders.items():
        selected = select_order(order, None, None, unique_loci=method not in POOLS)
        hits = covered_ids(loci, selected, overlap)
        runs["pool", 0, method] = selected, hits, hits


def add_family_replay(record, overlap, call_limit):
    verify_family_replay(record, overlap, call_limit)
    source = {method: candidates[:max(REPLAY_CAPS)] for method, candidates in record["pools"].items()}
    runs = record["runs"]
    for cap in REPLAY_CAPS:
        orders = selection_orders(source, overlap, max(1, cap - 1))
        gene_residues = sum(candidate.residues for candidate in source["gene_caller"][:cap])
        for mode, residues in (("count", None), ("gene_residues", gene_residues)):
            for method, order in orders.items():
                selected = select_order(order, cap, residues, unique_loci=method not in POOLS)
                for scorer in SCORERS:
                    hits, coverage = family_recovery(record, selected, scorer, cap, overlap, call_limit)
                    runs[scorer, mode, cap, method] = selected, hits, coverage


def summarize_tier(records, prefix, *, bootstrap=False):
    baseline = (*prefix, "gene_caller")
    return {
        method: {
            **summarize(records, (*prefix, method), include_coverage=len(prefix) == 3),
            "versus_gene_caller": paired_gain(records, (*prefix, method), baseline, bootstrap=bootstrap),
        }
        for method in METHODS
    }


def evaluate(cohort, partition):
    records, protocol, files = prepare(cohort, partition)
    overlap, call_limit = protocol["minimum_component_overlap"], protocol["prediction_limit"]
    for record in records:
        add_coverage(record, overlap)
        add_family_replay(record, overlap, call_limit)
    results = {
        "input_sha256": files,
        "partition": partition,
        "minimum_strand_coverage": overlap,
        "replay_output_call_cap": call_limit,
        "coverage_full_pool": summarize_tier(records, ("pool", 0)),
        "coverage_fixed_candidate_cap": {
            cap: summarize_tier(records, ("count", cap), bootstrap=cap in (2, 4))
            for cap in CAPS
        },
        "coverage_fixed_residue_cap": {
            cap: summarize_tier(records, ("residues", cap)) for cap in RESIDUE_BUDGETS
        },
        "coverage_candidate_and_gene_residue_caps": {
            cap: summarize_tier(records, ("gene_residues", cap), bootstrap=cap in (2, 4))
            for cap in CAPS
        },
        "historical_top_four_family_replay": {
            scorer: {
                mode: {
                    cap: summarize_tier(records, (scorer, mode, cap), bootstrap=cap in (2, 4))
                    for cap in REPLAY_CAPS
                }
                for mode in ("count", "gene_residues")
            }
            for scorer in SCORERS
        },
    }
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, default=Path("results/paired-cohort-v2"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--partition", choices=("holdout", "development"), default="holdout")
    args = parser.parse_args()
    result = evaluate(args.cohort, args.partition)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
