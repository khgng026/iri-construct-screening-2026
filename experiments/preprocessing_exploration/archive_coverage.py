"""Evaluate archived yeast cargo and reporter candidate reach without classification."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory

from Bio import SeqIO
from Bio.Seq import Seq
import pyrodigal

from experiments.paired_preprocessing.reporter_pilot import (
    AMINO_ACIDS, MAX_RESIDUES, PYRODIGAL_VERSION, _append_chunks, _candidates,
    _gene_caller_candidates, _length_order, _rank_by_start,
)
from experiments.paired_preprocessing.verify import _covered
from iri_constructs import process_inputs


PROTOCOL = Path(__file__).with_name("archive_coverage_panel.json")
RESOLVED_PROTEIN = re.compile(f"[{''.join(sorted(AMINO_ACIDS))}]+")


def resolved_translation_segments(dna):
    """Retain individually resolved stretches in otherwise incomplete archive records."""
    candidates = []
    for strand, sequence in ((1, Seq(dna)), (-1, Seq(dna).reverse_complement())):
        for offset in range(3):
            end = len(sequence) - (len(sequence) - offset) % 3
            translated = str(sequence[offset:end].translate())
            for segment in RESOLVED_PROTEIN.finditer(translated):
                anchor = offset + 3 * segment.start()
                _append_chunks(candidates, segment.group(), anchor if strand == 1 else len(dna) - anchor, strand)
    return sorted(candidates, key=_length_order)


def _source_record(repo: Path, revision: str, path: str) -> bytes:
    return subprocess.check_output(
        ["git", "-C", str(repo), "show", f"{revision}:{path}"],
    )


def _targets(record, path: Path, labels: set[str]) -> tuple[list[dict], list[dict]]:
    # Reuse the toolkit's complete-CDS validation and sequence/coordinate checks
    # on the evaluator side. Candidate methods receive only the DNA string.
    result = process_inputs([path])
    if len(result.ok_records) != 1:
        raise ValueError(f"{path.name}: expected one valid annotated DNA record")
    parsed = result.ok_records[0]
    units = {unit.feature_id: unit for unit in parsed.units}
    targets = []
    excluded = []
    seen = Counter()
    for feature in parsed.features:
        label = feature.names["label"]
        if feature.source_feature_type != "CDS" or label not in labels:
            continue
        seen[label] += 1
        unit = units[feature.feature_id]
        if unit.extraction_state != "exported":
            excluded.append({"label": label, "reason": unit.reason_code})
            continue
        cds = record.features[feature.source_feature_index]
        targets.append({
            "label": label,
            "strand": {"+": 1, "-": -1}[feature.strand],
            "spans": [(part["start0"], part["end0"]) for part in feature.segments],
            "codon_offset": int(cds.qualifiers.get("codon_start", ["1"])[0]) - 1,
            "protein": unit.sequence,
        })
    if set(seen) != labels or any(count != 1 for count in seen.values()):
        raise ValueError(f"{path.name}: expected one CDS for each panel label")
    return targets, excluded


def _peptide_at_expected_position(target: dict, query: dict) -> bool:
    """Check the full peptide and its frame at the annotated contiguous locus."""
    if target["strand"] != query["strand"] or len(query["spans"]) != 1:
        return False
    parts = sorted(target["spans"])
    if any(left[1] != right[0] for left, right in zip(parts, parts[1:])):
        return False
    if target["strand"] == 1:
        nt_offset = parts[0][0] - query["spans"][0][0]
    else:
        nt_offset = query["spans"][0][1] - parts[-1][1]
    nt_offset += target["codon_offset"]
    aa_offset, remainder = divmod(nt_offset, 3)
    if nt_offset < 0 or remainder:
        return False
    expected = target["protein"]
    return query["protein"][aa_offset:aa_offset + len(expected)] == expected


def _measure(targets: list[dict], pool: list[dict], overlap: float) -> list[dict]:
    results = []
    for target in targets:
        locus = (target["strand"], target["spans"])
        first_locus = None
        first_peptide = None
        for rank, query in enumerate(pool, start=1):
            covering = _covered(locus, (query["strand"], query["spans"])) >= overlap
            if covering and first_locus is None:
                first_locus = rank
            if _peptide_at_expected_position(target, query) and first_peptide is None:
                first_peptide = rank
            if first_locus is not None and first_peptide is not None:
                break
        results.append({"first_locus_rank": first_locus, "first_peptide_rank": first_peptide})
    return results


def _run_record(record, targets: list[dict], finder, overlap: float, caps: list[int], ranking_followup=False) -> dict:
    dna = str(record.seq).upper()
    minimal = _candidates(dna, start_aware=False)
    start_ranked = _rank_by_start(minimal)
    start_trimmed = _candidates(dna, start_aware=True)
    called = _gene_caller_candidates(dna, finder)
    pools = {
        "minimal": minimal, "start_ranked": start_ranked,
        "start_trimmed": start_trimmed, "gene_caller": called,
    }
    if ranking_followup:
        from .ranking import Candidate, rank_indices

        for name, pool in list(pools.items()):
            geometry = [Candidate(unit["strand"], tuple(map(tuple, unit["spans"])), len(unit["protein"]))
                        for unit in pool]
            for policy in ("length_center_250", "defer_ceiling"):
                pools[f"{name}_{policy}"] = [pool[index] for index in rank_indices(geometry, policy)]
        pools["resolved_translation_segments"] = resolved_translation_segments(dna)
    return {
        method: {
            "pool_size": len(pool),
            "residues_selected": {cap: sum(len(row["protein"]) for row in pool[:cap]) for cap in caps},
            "targets": _measure(targets, pool, overlap),
        }
        for method, pool in pools.items()
    }


def _summary(rows: list[dict], caps: list[int]) -> dict:
    result = {}
    for method in rows[0]["methods"] if rows else ():
        columns = [row["methods"][method] for row in rows]
        targets = [entry for column in columns for entry in column["targets"]]
        result[method] = {
            "targets": len(targets),
            "pool_size": sum(column["pool_size"] for column in columns),
            "pool_locus_reach": sum(target["first_locus_rank"] is not None for target in targets),
            "pool_peptide_reach": sum(target["first_peptide_rank"] is not None for target in targets),
            "curve": {
                cap: {
                    "loci": sum(
                        target["first_locus_rank"] is not None and target["first_locus_rank"] <= cap
                        for target in targets
                    ),
                    "peptides": sum(
                        target["first_peptide_rank"] is not None and target["first_peptide_rank"] <= cap
                        for target in targets
                    ),
                    "residues": sum(column["residues_selected"][cap] for column in columns),
                }
                for cap in caps
            },
        }
    return result


def run(panel: dict, repo_paths: dict[str, Path], ranking_followup=False) -> dict:
    if version("pyrodigal") != PYRODIGAL_VERSION:
        raise ValueError(f"Install pyrodigal=={PYRODIGAL_VERSION} for this comparison")
    caps = panel["candidate_caps"]
    overlap = panel["locus_overlap"]
    finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    rows = []
    with TemporaryDirectory() as work:
        for source in panel["sources"]:
            for entry in source["records"]:
                raw = _source_record(repo_paths[source["id"]], source["revision"], entry["path"])
                content_hash = sha256(raw).hexdigest()
                if content_hash != entry["sha256"]:
                    raise ValueError(f"{entry['path']}: archive hash differs from the panel")
                path = Path(work) / Path(entry["path"]).name
                path.write_bytes(raw)
                record = SeqIO.read(path, "genbank")
                targets, excluded = _targets(record, path, set(entry["target_labels"]))
                rows.append({
                    "source": source["id"],
                    "path": entry["path"],
                    "sha256": content_hash,
                    "length_nt": len(record),
                    "n_ambiguous": sum(base not in "ACGT" for base in str(record.seq).upper()),
                    "valid_targets": [
                        {"label": target["label"], "length_aa": len(target["protein"]),
                         "location_parts": len(target["spans"]),
                         "fits_one_chunk": len(target["protein"]) <= MAX_RESIDUES}
                        for target in targets
                    ],
                    "excluded": excluded,
                    "methods": _run_record(record, targets, finder, overlap, caps, ranking_followup),
                })
    return {
        "schema": "archive_candidate_coverage_v1",
        "candidate_caps": caps,
        "locus_overlap": overlap,
        "versions": {name: version(name) for name in ("biopython", "pyrodigal")},
        "sources": {source["id"]: _summary(
            [row for row in rows if row["source"] == source["id"]], caps,
        ) for source in panel["sources"]},
        "records": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reclone", type=Path, required=True)
    parser.add_argument("--lampra", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranking-followup", action="store_true")
    args = parser.parse_args()
    panel = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    result = run(panel, {"open_yeast": args.reclone, "lampra": args.lampra}, args.ranking_followup)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
