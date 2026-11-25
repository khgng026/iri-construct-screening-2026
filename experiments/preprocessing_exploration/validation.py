"""Test preprocessing policies on archived benign constructs without a scorer comparison."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from .archive_coverage import (
    PROTOCOL, _peptide_at_expected_position, _source_record, _targets,
    exact_window_coverage, resolved_translation_segments,
)
from .hybrid import distinct_supplements, parse_pool
from .ranking import Candidate, rank_indices
from experiments.paired_preprocessing.reporter_pilot import (
    AMINO_ACIDS, MAX_RESIDUES, PYRODIGAL_VERSION, _append_chunks,
    _candidates, _gene_caller_candidates, _length_order,
)
from experiments.paired_preprocessing.verify import _covered


PANEL = Path(__file__).with_name("validation_panel.json")
CAPS = (1, 2, 4, 8, 16, "pool")
FALLBACK_OVERLAP = 0.8


def _geometry(pool, name, length):
    return parse_pool([
        {**unit, "residues": len(unit["protein"])} for unit in pool
    ], name, length)


def distinct_fallback(genes, fallback, length):
    """Discard fallback candidates that substantially repeat an existing same-frame locus."""
    regions = distinct_supplements(_geometry(genes, "genes", length),
                                   _geometry(fallback, "fallback", length), FALLBACK_OVERLAP)
    return [fallback[region.index] for region in regions]


def called_genes(dna, finder):
    parents = []
    for gene in finder.find_genes(dna.encode("ascii")):
        protein = gene.translate(include_stop=False)
        if not protein or set(protein) - AMINO_ACIDS:
            continue
        start, end = gene.begin - 1, gene.end
        if not (0 <= start < end <= len(dna) and gene.strand in (-1, 1)
                and 0 <= end - start - len(protein) * 3 <= 3):
            raise ValueError("Gene coordinates disagree with translation")
        anchor = start if gene.strand == 1 else end
        parents.append((protein, anchor, gene.strand, start))
    # Use the same parent order in the two window arms. Order the parent by
    # its longest standard window, as in the original shortlist.
    parents.sort(key=lambda row: (-min(len(row[0]), MAX_RESIDUES), row[3], -row[2]))
    return parents


def grouped_windows(parents, *, balanced):
    groups = []
    for protein, anchor, strand, _ in parents:
        output = []
        if balanced and len(protein) > MAX_RESIDUES:
            count = (len(protein) + MAX_RESIDUES - 1) // MAX_RESIDUES
            for index in range(count):
                lo = index * len(protein) // count
                hi = (index + 1) * len(protein) // count
                _append_chunks(output, protein[lo:hi], anchor + strand * 3 * lo, strand)
        else:
            _append_chunks(output, protein, anchor, strand)
        groups.append(output)
    return groups


def preserve_first_parent(mixed, groups):
    """Keep other windows of the first called gene directly after its lead window."""
    if not mixed:
        return mixed
    first = mixed[0]

    def locus(unit):
        return unit["strand"], tuple(map(tuple, unit["spans"]))

    lead = locus(first)
    siblings = next({locus(unit) for unit in group} for group in groups if any(locus(unit) == lead for unit in group))
    return [first, *(unit for unit in mixed[1:] if locus(unit) in siblings),
            *(unit for unit in mixed[1:] if locus(unit) not in siblings)]


def policies(dna, finder):
    minimal = _candidates(dna, False)
    genes = _gene_caller_candidates(dna, finder)
    resolved = resolved_translation_segments(dna)
    geometry = [Candidate(unit["strand"], tuple(map(tuple, unit["spans"])), len(unit["protein"]))
                for unit in genes]
    mixed = [genes[index] for index in rank_indices(geometry, "keep_first_then_250")]
    fallback = distinct_fallback(genes, resolved, len(dna))
    parents = called_genes(dna, finder)
    groups = grouped_windows(parents, balanced=False)
    balanced = grouped_windows(parents, balanced=True)
    archive_uncertainty = bool(set(dna) - set("ACGT"))
    gated = [*mixed[:2], *fallback[:1], *mixed[2:], *fallback[1:]] if archive_uncertainty else [*mixed, *fallback]
    return {
        "minimal": minimal,
        "resolved_translation": resolved,
        "gene_caller": genes,
        "mixed_ranking": mixed,
        "gene_then_fallback": [*genes, *fallback],
        "mixed_then_fallback": [*mixed, *fallback],
        "uncertainty_fallback": gated,
        "linked_gene_windows": [unit for group in groups for unit in group],
        "balanced_gene_windows": [unit for group in balanced for unit in group],
        "guarded_mixed_ranking": preserve_first_parent(mixed, groups),
    }


def select(pool, cap, residue_budget):
    chosen = []
    residues = 0
    for unit in pool:
        if cap != "pool" and len(chosen) == cap:
            break
        length = len(unit["protein"])
        if residue_budget is not None and residues + length > residue_budget:
            continue
        chosen.append(unit)
        residues += length
    return chosen, residues


def measure(chosen, target):
    locus = target["strand"], target["spans"]
    overlap = any(_covered(locus, (row["strand"], row["spans"])) >= 0.8 for row in chosen)
    contained = any(_peptide_at_expected_position(target, row) for row in chosen)
    long = len(target["protein"]) > MAX_RESIDUES
    window_complete = long and exact_window_coverage(target, chosen) == len(target["protein"])
    return Counter({
        "eligible": 1,
        "locus_reached": overlap,
        "peptide_in_one_candidate": contained,
        "long_eligible": long,
        "long_complete_across_windows": window_complete,
        "peptide_available": contained or window_complete,
    })


def run(new_manifest, reclone, lampra):
    from Bio import SeqIO
    import pyrodigal

    if version("pyrodigal") != PYRODIGAL_VERSION:
        raise ValueError(f"Install pyrodigal=={PYRODIGAL_VERSION}")
    new_bytes = new_manifest.read_bytes()
    new = json.loads(new_bytes)
    archived = json.loads(PROTOCOL.read_bytes())
    sources = [
        {"id": "open_enzyme", "revision": new["revision"], "records": new["records"]},
        *archived["sources"],
    ]
    totals = defaultdict(Counter)
    by_group = defaultdict(Counter)
    audit = {}
    finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    with TemporaryDirectory() as tmp:
        for source in sources:
            sid = source["id"]
            repo = lampra if sid == "lampra" else reclone
            disposition = Counter()
            hashes = set()
            group_counts = Counter()
            exclusions = []
            for entry in source["records"]:
                contents = _source_record(repo, source["revision"], entry["path"])
                if hashlib.sha256(contents).hexdigest() != entry["sha256"]:
                    raise ValueError(f"{entry['path']}: source differs from manifest")
                path = Path(tmp) / Path(entry["path"]).name
                path.write_bytes(contents)
                record = SeqIO.read(path, "genbank")
                dna = str(record.seq).upper()
                arms = policies(dna, finder)
                # The caller and ranker receive only DNA. Labels enter here.
                targets, excluded = _targets(record, path, set(entry["target_labels"]))
                disposition["records"] += 1
                disposition["records_with_uncertainty"] += bool(set(dna) - set("ACGT"))
                disposition["accepted_targets"] += len(targets)
                for target in targets:
                    hashes.add(hashlib.sha256(target["protein"].encode("ascii")).hexdigest())
                for row in excluded:
                    disposition[f"excluded_{row['reason']}"] += 1
                    exclusions.append({"path": entry["path"], **row})
                if not targets:
                    continue
                group = entry.get("group", sid)
                group_counts[group] += len(targets)
                for cap in CAPS:
                    _, baseline_residues = select(arms["gene_caller"], cap, None)
                    for mode, ceiling in (("candidate_cap", None), ("gene_residue_ceiling", baseline_residues)):
                        for method, pool in arms.items():
                            chosen, residues = select(pool, cap, ceiling)
                            key = sid, mode, str(cap), method
                            totals[key].update({"scorer_inputs": len(chosen), "residues": residues})
                            for target in targets:
                                metrics = measure(chosen, target)
                                totals[key].update(metrics)
                                if cap == 2 and mode == "candidate_cap":
                                    by_group[sid, group, method].update(metrics)
            audit[sid] = {
                **disposition, "unique_translations": len(hashes), "reference_groups": dict(group_counts),
                "excluded_targets": exclusions,
            }
    runs = {sid: {} for sid in audit}
    for (sid, mode, cap, method), metrics in totals.items():
        runs[sid].setdefault(mode, {}).setdefault(cap, {})[method] = dict(metrics)
    groups = {sid: {} for sid in audit}
    for (sid, group, method), metrics in by_group.items():
        groups[sid].setdefault(group, {})[method] = dict(metrics)
    return {
        "schema": "iri-preprocessing-validation-v1",
        "new_panel_sha256": hashlib.sha256(new_bytes).hexdigest(),
        "archive_panel_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        "environment": {pkg: version(pkg) for pkg in ("biopython", "pyrodigal")},
        "audit": audit,
        "runs": runs,
        "cap2_catalogue_group_outcomes": groups,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reclone", type=Path, required=True)
    parser.add_argument("--lampra", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run(PANEL, args.reclone, args.lampra)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
