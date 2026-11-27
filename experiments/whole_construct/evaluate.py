"""Evaluate uncapped family calls against held-out annotations at recorded cutoffs."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path

from experiments.paired_preprocessing.verify import _covered, _index, _location, _match


def match_calls(construct, calls, cutoff, overlap):
    """Score one call per documented component and track remaining call burden."""
    size, name = construct["length"], construct["id"]
    expected = [(unit, _location(unit, size, unit["id"])) for unit in construct["components"]]
    ranked = []
    for call in calls:
        score, family = call["score"], call["family"]
        if not isfinite(score) or not family:
            raise ValueError(f"{name}: invalid family score")
        if score >= cutoff:
            ranked.append((score, family, _location(call, size, name)))
    # Consolidate repeated calls for the same family and nearly the same span.
    # A local domain and a much larger candidate can remain separate evidence.
    seen = {}
    edges = []
    for _, family, locus in sorted(ranked, key=lambda row: (-row[0], row[1], row[2])):
        previous = seen.setdefault(family, [])
        if any(min(_covered(locus, prior), _covered(prior, locus)) >= 0.8 for prior in previous):
            continue
        previous.append(locus)
        edges.append([
            unit["id"] for unit, documented_locus in expected
            if unit["family"] == family and _covered(documented_locus, locus) >= overlap
        ])
    return set(_match(edges)), len(edges), len(expected)


def summarize(answer, per_construct, cutoff, overlap):
    constructs = _index(answer["constructs"], "id", "answer key")
    if constructs.keys() != per_construct.keys():
        raise ValueError("Call records differ from the answer key")
    counts = Counter({"constructs": len(constructs), "documented": 0, "returned": 0, "recovered": 0,
                      "unmatched_complete": 0, "unverified_partial": 0})
    for name, construct in constructs.items():
        matched, emitted, documented = match_calls(construct, per_construct[name], cutoff, overlap)
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
        "curve": {cutoff: summarize(answer, trace["calls"], cutoff, args.overlap) for cutoff in sorted(set(args.cutoff))},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(response, output, indent=2, sort_keys=True)
        output.write("\n")


if __name__ == "__main__":
    main()
