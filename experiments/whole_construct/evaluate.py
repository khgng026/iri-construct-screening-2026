"""Evaluate uncapped family calls against held-out annotations at recorded cutoffs."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path

from experiments.paired_preprocessing.verify import _covered, _index, _location, _match
from experiments.preprocessing_exploration.ranking import merged_intervals


def _join_windows(ranked):
    parents = {}
    independent = []
    for score, family, (strand, spans), parent in ranked:
        if parent is None:
            independent.append((score, family, (strand, spans), parent))
            continue
        key = family, strand, parent
        prior, parts = parents.setdefault(key, (score, []))
        parts.extend(spans)
        parents[key] = max(prior, score), parts
    return independent + [(score, family, (strand, merged_intervals(parts)), parent)
                          for (family, strand, parent), (score, parts) in parents.items()]


def _same_apparent_frame(left, right):
    strand, spans = left
    other_strand, other_spans = right
    if strand != other_strand:
        return False
    if len(spans) != 1 or len(other_spans) != 1:
        return True  # A single genomic frame is undefined for a compound locus.
    boundary = 0 if strand == 1 else 1
    return (spans[0][boundary] - other_spans[0][boundary]) % 3 == 0


def match_calls(construct, calls, cutoff, overlap, *, join_linked=False, require_frame=False):
    """Score one call per documented component and track remaining call burden."""
    size, name = construct["length"], construct["id"]
    expected = [(unit, _location(unit, size, unit["id"])) for unit in construct["components"]]
    ranked = []
    for call in calls:
        score, family = call["score"], call["family"]
        if not isfinite(score) or not family:
            raise ValueError(f"{name}: invalid family score")
        if score >= cutoff:
            ranked.append((score, family, _location(call, size, name), call.get("parent_id")))
    if join_linked:
        ranked = _join_windows(ranked)
    # Consolidate repeated calls for the same family and nearly the same span.
    # A local domain and a much larger candidate can remain separate evidence.
    seen = {}
    edges = []
    for _, family, locus, _ in sorted(ranked, key=lambda row: (-row[0], row[1], row[2])):
        previous = seen.setdefault(family, [])
        if any(min(_covered(locus, prior), _covered(prior, locus)) >= 0.8
               and (not require_frame or _same_apparent_frame(locus, prior)) for prior in previous):
            continue
        previous.append(locus)
        edges.append([
            unit["id"] for unit, documented_locus in expected
            if unit["family"] == family and _covered(documented_locus, locus) >= overlap
            and (not require_frame or _same_apparent_frame(documented_locus, locus))
        ])
    return set(_match(edges)), len(edges), len(expected)


def summarize(answer, per_construct, cutoff, overlap, *, join_linked=False, require_frame=False):
    constructs = _index(answer["constructs"], "id", "answer key")
    if constructs.keys() != per_construct.keys():
        raise ValueError("Call records differ from the answer key")
    counts = Counter({"constructs": len(constructs), "documented": 0, "returned": 0, "recovered": 0,
                      "unmatched_complete": 0, "unverified_partial": 0})
    for name, construct in constructs.items():
        matched, emitted, documented = match_calls(construct, per_construct[name], cutoff, overlap,
                                                   join_linked=join_linked, require_frame=require_frame)
        unresolved = "unmatched_complete" if construct.get("annotation_complete") else "unverified_partial"
        counts.update(documented=documented, returned=emitted, recovered=len(matched))
        counts[unresolved] += emitted - len(matched)
    return dict(counts) | {
        "recall": counts["recovered"] / counts["documented"] if counts["documented"] else None,
        "observed_supported_fraction": counts["recovered"] / counts["returned"] if counts["returned"] else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--answer-key", type=Path, required=True)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--cutoff", type=float, action="append", required=True)
    parser.add_argument("--overlap", type=float, default=0.8)
    parser.add_argument("--join-linked-windows", action="store_true")
    parser.add_argument("--require-frame", action="store_true", help="Check reading frame for single spans")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.overlap <= 1 or not all(map(isfinite, args.cutoff)):
        parser.error("Provide valid score cutoffs and locus overlap")
    key_bytes, call_bytes = args.answer_key.read_bytes(), args.calls.read_bytes()
    answer, trace = json.loads(key_bytes), json.loads(call_bytes)
    if trace.get("schema") != "iri-whole-construct-calls-v1":
        raise ValueError("Unexpected call schema")
    response = {
        "schema": "iri-whole-construct-evaluation-v1", "arm": trace["arm"],
        "answer_key_sha256": sha256(key_bytes).hexdigest(), "prediction_sha256": sha256(call_bytes).hexdigest(),
        "input_sha256": trace["input_sha256"], "references_sha256": trace["references_sha256"],
        "minimum_documented_locus_overlap": args.overlap,
        "join_linked_windows": args.join_linked_windows,
        "require_single_span_frame": args.require_frame,
        "curve": {cutoff: summarize(answer, trace["calls"], cutoff, args.overlap, join_linked=args.join_linked_windows,
                                    require_frame=args.require_frame)
                  for cutoff in sorted(set(args.cutoff))},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(response, output, indent=2, sort_keys=True)
        output.write("\n")


if __name__ == "__main__":
    main()
