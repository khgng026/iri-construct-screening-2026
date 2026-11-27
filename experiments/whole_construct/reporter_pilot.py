"""Exercise the uncapped raw-DNA baseline and profile candidate preparation on a verified archive."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter
from types import SimpleNamespace

from experiments.paired_preprocessing.cohort_study import _aligner, candidate_units, similarity_tier
from experiments.paired_preprocessing.reporter_pilot import PYRODIGAL_VERSION
from experiments.preprocessing_exploration.archive_coverage import _peptide_at_expected_position
from experiments.preprocessing_exploration.translated_search import E_VALUE, Region, prepare, search, write_fasta

from .evaluate import match_calls, summarize
from .scoring import blastx_calls


FAMILY = "corroborated_gfp_like"
EVALUES = (1e-20, 1e-10, 1e-5)
OVERLAPS = (0.5, 0.8)


def alignment_control(answer, panel, candidate_pools, blastp):
    """Hold family labels fixed and compare grading local hits and whole candidates."""
    result = {}
    with TemporaryDirectory() as scratch:
        work = Path(scratch)
        refs = work / "refs.fa"
        write_fasta(refs, panel)
        for method, records in candidate_pools.items():
            locations, sequences = {}, {}
            for record, pool in records.items():
                for index, unit in enumerate(pool):
                    (left, right), = unit["spans"]
                    query_id = f"q{len(sequences)}"
                    locations[query_id] = record, index, Region(left, right, unit["strand"])
                    sequences[query_id] = unit["protein"]
            inputs = work / f"{method}.fa"
            write_fasta(inputs, sequences)
            hits, elapsed = search(SimpleNamespace(blast_bin=blastp.parent), "blast", "blastp",
                                   inputs, refs, None, locations)
            calls = {endpoint: {record: [] for record in records} for endpoint in ("local_hit", "candidate_span")}
            for record, values in hits.items():
                for hit in values:
                    candidate = records[record][hit.candidate_index]
                    for endpoint, spans in (("local_hit", [[hit.region.left, hit.region.right]]),
                                            ("candidate_span", candidate["spans"])):
                        calls[endpoint][record].append({"family": FAMILY, "score": hit.score,
                                                       "strand": hit.region.strand, "spans": spans})
            result[method] = {
                "reporting_evalue_ceiling": float(E_VALUE), "wall_seconds": elapsed,
                "endpoint": {endpoint: {str(rule): summarize(answer, rows, 0, rule) for rule in OVERLAPS}
                             for endpoint, rows in calls.items()},
            }
    return result


def pilot(args):
    from Bio.Seq import Seq
    import pyrodigal

    if version("pyrodigal") != PYRODIGAL_VERSION:
        raise ValueError(f"Install pyrodigal=={PYRODIGAL_VERSION}")
    records, targets, panel, audit_sha = prepare(args)
    references = [{"accession": name, "family": FAMILY, "sequence": seq} for name, seq in panel.items()]
    aligner = _aligner()
    finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    answer = {"constructs": []}
    tiers, n_targets = {}, Counter()
    background = defaultdict(Counter)
    elapsed = defaultdict(float)
    distinct_inputs = defaultdict(set)
    candidate_pools = defaultdict(dict)
    translations = set()
    for name, sequence in records.items():
        label, region = targets[name]
        strand = region.strand
        spans = [[region.left, region.right]]
        component = {"id": name, "family": FAMILY, "strand": strand, "spans": spans}
        verified = label == FAMILY
        answer["constructs"].append({"id": name, "length": len(sequence),
                                     "components": [component] if verified else []})
        if verified:
            coding = Seq(sequence[region.left:region.right])
            protein = str((coding if strand == 1 else coding.reverse_complement()).translate())
            translations.add(protein)
            expected = {"protein": protein, "strand": strand, "spans": spans, "codon_offset": 0}
            tier = similarity_tier(aligner, protein, FAMILY, references)["tier"]
            tiers[name] = tier
            n_targets[tier] += 1
        for method in ("minimal", "gene_caller"):
            started = perf_counter()
            pool = candidate_units(method, sequence, finder)
            elapsed[method] += perf_counter() - started
            candidate_pools[method][name] = pool
            distinct_inputs[method].update(unit["protein"] for unit in pool)
            background[method, "all"].update(constructs=1, candidates=len(pool),
                                             residues=sum(len(unit["protein"]) for unit in pool))
            if verified:
                reached = any(_peptide_at_expected_position(expected, unit) for unit in pool)
                for group in ("verified", tier):
                    background[method, group].update(documented=1, exact_peptides_in_pool=reached)

    raw, work = blastx_calls(records, references, args.blastx, evalue=max(EVALUES))
    curve = {}
    for cutoff in EVALUES:
        passing = {name: [call for call in calls if call["evalue"] <= cutoff] for name, calls in raw.items()}
        for overlap in OVERLAPS:
            row = summarize(answer, passing, 0, overlap)
            strata = Counter()
            for construct in answer["constructs"]:
                name = construct["id"]
                if name in tiers:
                    matched, _, _ = match_calls(construct, passing[name], 0, overlap)
                    strata[tiers[name]] += len(matched)
            row["recovered_by_tier"] = dict(strata)
            curve.setdefault(str(cutoff), {})[str(overlap)] = row

    result = {
        "schema": "iri-whole-construct-open-reporters-pilot-v1",
        "evidence": {"archive_audit_sha256": audit_sha, "verified_target_occurrences": len(tiers),
                     "verified_unique_translations": len(translations), "other_catalogue_controls": len(records) - len(tiers),
                     "documented_by_similarity_tier": n_targets,
                     "references_sha256": sha256(json.dumps(references, sort_keys=True).encode()).hexdigest()},
        "uncapped_raw_dna": {"curve": curve, "work": work, "frame_aware_at_reporting_ceiling": {
            str(overlap): summarize(answer, raw, 0, overlap, require_frame=True) for overlap in OVERLAPS
        }},
        "candidate_preparation": {
            method: {"by_scope": {group: counts for (arm, group), counts in background.items() if arm == method},
                     "distinct_peptide_inputs": len(distinct_inputs[method]),
                     "distinct_input_residues": sum(map(len, distinct_inputs[method])),
                     "wall_seconds": elapsed[method]}
            for method in elapsed
        },
    }
    if args.blastp is not None:
        result["local_alignment_control"] = alignment_control(answer, panel, candidate_pools, args.blastp)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--construct-repo", type=Path, required=True)
    parser.add_argument("--protein-repo", type=Path, required=True)
    parser.add_argument("--blastx", type=Path, required=True)
    parser.add_argument("--blastp", type=Path, help="Run the optional coordinate-attribution control")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.blastp is not None and args.blastp.name != "blastp":
        parser.error("--blastp must name the blastp executable")
    data = pilot(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(data, output, indent=2, sort_keys=True)
        output.write("\n")


if __name__ == "__main__":
    main()
