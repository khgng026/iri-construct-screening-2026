"""Prepare separated scorer inputs and evaluator files from reviewed archived GenBank records."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from experiments.paired_preprocessing.cohort_study import _aligner, similarity_tier
from experiments.paired_preprocessing.corpus_batch import protein_to_spans, transcript_positions
from experiments.paired_preprocessing.verify import _index
from experiments.preprocessing_exploration.translated_search import write_fasta

from .scoring import _references


PARTS = ("development", "holdout")


def _claim_split(claimed, key, partition, category):
    if claimed.setdefault(key, partition) != partition:
        raise ValueError(f"Development and holdout share {category}")


def _write_json(path, data):
    with path.open("x", encoding="utf-8") as output:
        json.dump(data, output, sort_keys=True, indent=2, allow_nan=False)
        output.write("\n")


def _corroborated_target(record, exclusion, rule):
    """Validate missing deposited translation against an independently archived exact peptide."""
    from Bio.SeqFeature import ExactPosition

    label = exclusion["label"]
    if exclusion["reason"] != "missing_translation" or not rule.get("corroborating_protein_sha256"):
        raise ValueError(f"{label}: frozen target failed annotation validation")
    feature, = [feature for feature in record.features
                if feature.type == "CDS" and label in feature.qualifiers.get("label", [])]
    location = feature.location
    exceptional = {"pseudo", "pseudogene", "exception", "transl_except", "ribosomal_slippage"}
    if (location is None or len(location.parts) != 1 or location.strand not in (-1, 1)
            or not isinstance(location.start, ExactPosition) or not isinstance(location.end, ExactPosition)
            or not 0 <= location.start < location.end <= len(record.seq)
            or set(feature.qualifiers) & exceptional):
        raise ValueError(f"{label}: expected a complete contiguous ordinary CDS")
    offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
    table = int(feature.qualifiers.get("transl_table", ["1"])[0])
    if offset not in (0, 1, 2):
        raise ValueError(f"{label}: unexpected codon offset")
    peptide = str(feature.extract(record.seq)[offset:].translate(table=table, cds=True))
    if sha256(peptide.encode()).hexdigest() != rule["corroborating_protein_sha256"]:
        raise ValueError(f"{label}: DNA translation differs from the corroborated peptide")
    return {"label": label, "strand": location.strand, "spans": [(int(location.start), int(location.end))],
            "protein": peptide, "codon_offset": offset}


def prepare(manifest_bytes, reference_bytes, sources, output):
    from Bio import SeqIO
    from experiments.preprocessing_exploration.archive_coverage import _targets

    if output.exists():
        raise ValueError("Choose a new output directory for the prepared archive")
    manifest = json.loads(manifest_bytes)
    reference_hash = sha256(reference_bytes).hexdigest()
    if (manifest["schema"] != "iri-whole-construct-archive-manifest-v1"
            or manifest["references_sha256"] != reference_hash):
        raise ValueError("The reference bytes differ from the archive manifest")
    references = json.loads(reference_bytes)
    ref_sequences, labels = _references(references)
    ref_proteins, possible_families = set(ref_sequences.values()), set(labels.values())
    aligner = _aligner()
    root = sources.resolve()
    claimed = {}
    sections = {part: {"records": {}, "truth": [], "metadata": [], "units": {}} for part in PARTS}
    entries = _index(manifest["records"], "id", "archive manifest")
    source_evidence = []
    with TemporaryDirectory() as tmp:
        copy = Path(tmp) / "source.gb"
        for name, entry in sorted(entries.items()):
            partition = entry["partition"]
            if partition not in sections or not entry["lineage_group"] or not entry["source_url"]:
                raise ValueError(f"{name}: partition, archive source and lineage evidence are required")
            for field in ("lineage_group", "protein_lineage", "pair_id"):
                if entry.get(field):
                    _claim_split(claimed, (field, entry[field]), partition, field)
            source = (root / entry["path"]).resolve()
            source.relative_to(root)
            contents = source.read_bytes()
            if sha256(contents).hexdigest() != entry["sha256"]:
                raise ValueError(f"{name}: archive bytes differ from the manifest")
            copy.write_bytes(contents)
            record = SeqIO.read(copy, "genbank")
            dna = str(record.seq).upper()
            _claim_split(claimed, ("DNA", sha256(dna.encode()).hexdigest()), partition, "identical DNA")
            target_rules = _index(entry["targets"], "label", f"{name} target labels")
            if not target_rules:
                raise ValueError(f"{name}: provide at least one independently supported target")
            targets, exclusions = _targets(record, copy, set(target_rules))
            targets.extend(_corroborated_target(record, row, target_rules[row["label"]]) for row in exclusions)
            data = sections[partition]
            key = f"{partition[0]}{len(data['records']) + 1:05d}"
            truth, units, meta = [], [], []
            for i, target in enumerate(targets, 1):
                family = target_rules[target["label"]]["family"]
                evidence = target_rules[target["label"]]["evidence"]
                if family not in possible_families or not evidence.strip():
                    raise ValueError(f"{name}: target lacks independent evidence or a scorable family")
                protein = target["protein"]
                protein_hash = sha256(protein.encode()).hexdigest()
                expected_hash = target_rules[target["label"]].get("corroborating_protein_sha256")
                if expected_hash is not None and protein_hash != expected_hash:
                    raise ValueError(f"{name}: the verified translation differs from the independent evidence")
                _claim_split(claimed, ("protein", protein_hash), partition, "an identical target protein")
                if protein in ref_proteins:
                    raise ValueError(f"{name}: an exact target protein occurs in the reference panel")
                cid, strand = f"{key}-{i}", target["strand"]
                positions = transcript_positions([{"start0": left, "end0": right, "strand": strand}
                                                 for left, right in target["spans"]], target["codon_offset"] + 1)
                spans = protein_to_spans(positions, 0, len(protein))
                if sum(end - start for start, end in spans) != 3 * len(protein):
                    raise ValueError(f"{name}: target residue coordinates disagree with the verified translation")
                truth.append({"id": cid, "family": family, "strand": strand, "spans": spans})
                units.append({"id": cid, "strand": strand, "spans": spans, "protein": protein})
                meta.append({"id": cid, "family": family, "protein_sha256": protein_hash, "evidence": evidence,
                             "similarity": similarity_tier(aligner, protein, family, references)})
            data["records"][key] = dna
            data["units"][key] = units
            data["truth"].append({"id": key, "length": len(dna), "components": truth,
                                  "annotation_complete": entry.get("annotation_complete", False) is True})
            data["metadata"].append({"id": key, "lineage_group": entry["lineage_group"], "components": meta,
                                     "protein_lineage": entry.get("protein_lineage"), "pair_id": entry.get("pair_id")})
            source_evidence.append({"id": key, "source_url": entry["source_url"], "source_sha256": entry["sha256"]})
    if any(not sections[part]["truth"] for part in PARTS):
        raise ValueError("The archive requires both development and holdout records")

    output.mkdir(parents=True)
    for partition, data in sections.items():
        method_dir, truth_dir = output / partition / "method_inputs", output / partition / "evaluator"
        for path in (method_dir, truth_dir):
            path.mkdir(parents=True)
        write_fasta(method_dir / "constructs.fasta", data["records"])
        _write_json(truth_dir / "answer-key.json", {"constructs": data["truth"]})
        _write_json(truth_dir / "construct-map.private.json", data["metadata"])
        _write_json(truth_dir / "annotation-units.json", data["units"])
    _write_json(output / "archive-audit.private.json", {
        "source_records": source_evidence, "references_sha256": reference_hash,
        "manifest_sha256": sha256(manifest_bytes).hexdigest(),
        "artifacts_sha256": {str(path.relative_to(output)): sha256(path.read_bytes()).hexdigest()
                             for path in sorted(output.rglob("*")) if path.is_file()},
        "documented_by_partition_and_tier": {
            part: Counter(unit["similarity"]["tier"] for row in data["metadata"] for unit in row["components"])
            for part, data in sections.items()
        },
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.manifest.read_bytes(), args.references.read_bytes(), args.sources, args.out)


if __name__ == "__main__":
    main()
