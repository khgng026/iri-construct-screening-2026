"""Measure protein-to-DNA annotation for the public LAMPRA EGFP example."""
from __future__ import annotations

import argparse
from io import StringIO
import json
from pathlib import Path
import subprocess

from .archive_coverage import PROTOCOL, _source_record, _targets
from .translated_search import PROTEIN_REV, REFERENCES, _git_source, _sha, write_fasta


def models_from_gff(gff):
    models = {}
    for line in gff.splitlines():
        if not line or line.startswith("#"):
            continue
        _, _, kind, start, end, score, strand, phase, attrs = line.split("\t")
        fields = dict(field.split("=", 1) for field in attrs.split(";") if field)
        if kind == "mRNA":
            models[fields["ID"]] = {"strand": strand, "score": float(score), "exons": []}
        if kind == "CDS":
            models[fields["Parent"]]["exons"].append((int(start) - 1, int(end), int(phase)))
    return sorted(models.values(), key=lambda model: -model["score"])


def evaluate_model(record, model, target):
    strand = 1 if model["strand"] == "+" else -1
    spans = sorted(model["exons"], reverse=strand == -1)
    coding = [
        record.seq[left:right] if strand == 1 else record.seq[left:right].reverse_complement()
        for left, right, _ in spans
    ]
    dna = sum(coding, record.seq[:0])[spans[0][2]:] if spans else record.seq[:0]
    protein = str(dna.translate()).rstrip("*") if len(dna) % 3 == 0 else ""
    overlap = sum(
        max(0, min(right, stop) - max(left, start))
        for left, right, _ in spans
        for start, stop in target["spans"]
    ) if strand == target["strand"] else 0
    exact = sorted((left, right) for left, right, _ in spans) == sorted(target["spans"])
    return {
        "predicted_intervals": len(spans),
        "same_strand": strand == target["strand"],
        "coding_region_overlap": overlap / sum(right - left for left, right in target["spans"]),
        "exact_annotated_intervals": strand == target["strand"] and exact,
        "exact_protein": strand == target["strand"] and protein == target["protein"],
        "candidate_residues": len(protein),
    }


def run(args):
    from Bio import SeqIO

    panel = json.loads(PROTOCOL.read_bytes())
    source, = [entry for entry in panel["sources"] if entry["id"] == "lampra"]
    entry, = source["records"]
    raw = _source_record(args.lampra, source["revision"], entry["path"])
    if _sha(raw) != entry["sha256"]:
        raise ValueError("Archived plasmid map differs from the frozen source")
    args.work.mkdir(parents=True)
    record_file = args.work / "record.gbk"
    record_file.write_bytes(raw)
    record = SeqIO.read(record_file, "genbank")
    natural = _git_source(args.protein_repo, PROTEIN_REV, "data/alignments/extant_fluorescent_proteins.fasta")
    proteins = {item.description.split("|", 1)[0].split()[0]: str(item.seq)
                for item in SeqIO.parse(StringIO(natural.decode()), "fasta")}
    references = {name: proteins[name] for name in REFERENCES}
    dna, refs = args.work / "dna.fasta", args.work / "refs.fasta"
    write_fasta(dna, {"reporter_plasmid": str(record.seq).upper()})
    write_fasta(refs, references)
    completed = subprocess.run(
        [str(args.miniprot), "--gff", "-t", "1", str(dna), str(refs)],
        capture_output=True, text=True, check=True,
    )
    targets, excluded = _targets(record, record_file, set(entry["target_labels"]))
    if excluded or len(targets) != 1 or targets[0]["protein"] in references.values():
        raise ValueError("The documented reporter must be verified and absent from the references")
    target = targets[0]
    spans = sorted(target["spans"])
    models = [evaluate_model(record, model, target) for model in models_from_gff(completed.stdout)]
    return {
        "schema": "iri-compound-reporter-pilot-v1",
        "source_sha256": _sha(raw),
        "reference_sha256": _sha(refs.read_bytes()),
        "version": subprocess.check_output([str(args.miniprot), "--version"], text=True).strip(),
        "documented_intervals": len(spans),
        "intervening_bases": [start - end for (_, end), (start, _) in zip(spans, spans[1:])],
        "documented_residues": len(target["protein"]),
        "matches": len(models),
        "top_model": models[0] if models else None,
        "exact_protein_models": sum(model["exact_protein"] for model in models),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("lampra", "protein-repo", "miniprot", "work", "out"):
        parser.add_argument(f"--{option}", required=True, type=Path)
    args = parser.parse_args()
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(run(args), handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
