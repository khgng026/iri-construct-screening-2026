"""Audit independent reporter labels and candidate reach on published cloning vectors."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from importlib.metadata import version
from io import StringIO
import json
from pathlib import Path
import subprocess

from .evidence_audit import component_rows, read


PANEL = Path(__file__).with_name("open_reporters.freeze.json")
CAPS = (1, 2, 4)


def _git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args])


def _checkout_source(repo, source):
    if _git(repo, "rev-parse", "HEAD").decode().strip() != source["commit"]:
        raise ValueError(f"Expected the pinned commit for {source['repository']}")


def _sha(value):
    return hashlib.sha256(value).hexdigest()


def _protein_db(source_bytes):
    from Bio import SeqIO

    by_name = {}
    by_protein = defaultdict(list)
    for record in SeqIO.parse(StringIO(source_bytes.decode()), "fasta"):
        display_name, accession, _, _, species = record.description.split("|", 4)
        name = display_name.split()[0]
        protein = str(record.seq)
        if name in by_name:
            raise ValueError(f"Duplicate reference name: {name}")
        by_name[name] = protein
        by_protein[protein].append({"name": name, "accession": accession, "genus": species.split()[0]})
    return by_name, by_protein


def _reporter(source_bytes):
    from Bio import SeqIO
    from Bio.SeqFeature import ExactPosition

    record = SeqIO.read(StringIO(source_bytes.decode()), "genbank")
    if record.annotations.get("topology") != "circular":
        raise ValueError("Expected a circular whole-plasmid record")
    coding = [feature for feature in record.features if feature.type == "CDS"]
    markers = [feature for feature in coding if feature.qualifiers.get("label") == ["bla"]]
    coding = [feature for feature in coding if feature not in markers]
    if len(coding) != 1 or len(markers) != 1:
        raise ValueError("Expected one reporter and the documented backbone marker")
    feature = coding[0]
    if (
        len(feature.location.parts) != 1
        or feature.location.strand not in (-1, 1)
        or not isinstance(feature.location.start, ExactPosition)
        or not isinstance(feature.location.end, ExactPosition)
        or feature.qualifiers.get("codon_start", ["1"]) != ["1"]
    ):
        raise ValueError("Expected one complete stranded coding region")
    label = feature.qualifiers["label"][0]
    protein = str(feature.extract(record.seq).translate(table=1, cds=True))
    if "translation" in feature.qualifiers and feature.qualifiers["translation"] != [protein]:
        raise ValueError("The source translation differs from the reporter coding region")
    return record, {
        "name": label,
        "protein": protein,
        "strand": feature.location.strand,
        "spans": [(int(feature.location.start), int(feature.location.end))],
    }


def _reach(target, pool, overlap):
    from experiments.paired_preprocessing.verify import _covered

    locus = target["strand"], target["spans"]

    def covering(candidates):
        return [
            row for row in candidates if _covered(locus, (row["strand"], row["spans"])) >= overlap
        ]

    tiers = {}
    for cap in CAPS:
        reached = covering(pool[:cap])
        tiers[str(cap)] = {
            "reached": bool(reached),
            "exact": any(row["protein"] == target["protein"] for row in reached),
            "contains_documented_peptide": any(
                target["protein"] in row["protein"] for row in reached
            ),
            "inputs": min(cap, len(pool)),
        }
    return {"pool_reached": bool(covering(pool)), "caps": tiers}


def _counts(rows, cap):
    genera = [
        {match["genus"] for match in row["independent_matches"]}
        for row in rows
    ]
    return {
        "reporters": len(rows),
        "independent_exact_reference_matches": sum(bool(row["independent_matches"]) for row in rows),
        "genus_counts": dict(Counter(next(iter(group)) for group in genera if len(group) == 1)),
        "ambiguous_genus": sum(len(group) > 1 for group in genera),
        "cohort_v2_exact_matches": sum(row["cohort_v2_exact_match"] for row in rows),
        "primary_cap": {
            method: {
                field: sum(row["candidate_reach"][method]["caps"][str(cap)][field] for row in rows)
                for field in ("reached", "exact", "contains_documented_peptide", "inputs")
            }
            for method in ("minimal", "gene_caller")
        },
    }


def audit(construct_repo, protein_repo, project_root):
    from Bio.Data.CodonTable import TranslationError
    import pyrodigal

    from experiments.paired_preprocessing.distant_reporter_study import _identities
    from experiments.paired_preprocessing.reporter_pilot import (
        PYRODIGAL_VERSION, _candidates, _gene_caller_candidates,
    )

    panel = read(PANEL)
    source = panel["construct_source"]
    protein_source = panel["protein_source"]
    _checkout_source(construct_repo, source)
    _checkout_source(protein_repo, protein_source)
    if version("pyrodigal") != PYRODIGAL_VERSION:
        raise ValueError(f"Expected Pyrodigal {PYRODIGAL_VERSION}")
    finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    natural = _git(protein_repo, "show", f"{protein_source['commit']}:{protein_source['path']}")
    natural_by_name, natural_by_protein = _protein_db(natural)
    reference_sets = dict(panel["reference_sets"])
    reference_sets["all_groups"] = list(dict.fromkeys(
        name for names in reference_sets.values() for name in names
    ))
    cohort_hashes = {
        row["protein_sha256"]
        for split in ("development", "holdout")
        for row in component_rows(project_root, split)
    }

    rows, targets, rejected = {}, [], {}
    gfp_records = set(panel["gfp_like_record_numbers"])
    for number in sorted(gfp_records | set(panel["other_reporter_record_numbers"])):
        record_id = f"ODC_{number:04}"
        source_path = f"{source['directory']}/{record_id}.gb"
        contents = _git(construct_repo, "show", f"{source['commit']}:{source_path}")
        try:
            record, target = _reporter(contents)
        except (ValueError, TranslationError) as error:
            rejected[record_id] = {"source_sha256": _sha(contents), "reason_type": type(error).__name__}
            continue
        if record.id != record_id:
            raise ValueError("Source record ID differs from the frozen catalogue ID")
        group = "gfp_like" if number in gfp_records else "other_reporter"
        matches = natural_by_protein.get(target["protein"], []) if group == "gfp_like" else []
        dna = str(record.seq).upper()
        rows[record_id] = {
            "source_sha256": _sha(contents),
            "label": target["name"],
            "catalogue_group": group,
            "independent_matches": matches,
            "cohort_v2_exact_match": _sha(target["protein"].encode("ascii")) in cohort_hashes,
            "candidate_reach": {
                "minimal": _reach(target, _candidates(dna, start_aware=False), panel["overlap_at_least"]),
                "gene_caller": _reach(target, _gene_caller_candidates(dna, finder), panel["overlap_at_least"]),
            },
        }
        if matches:
            targets.append({"id": record_id, "family": "gfp_like", "protein": target["protein"]})

    cap = panel["primary_candidate_cap"]
    distances = {}
    for name, references in reference_sets.items():
        library = [{"accession": label, "family": "gfp_like", "sequence": natural_by_name[label]}
                   for label in references]
        measures = _identities(targets, library)
        for key, measure in measures.items():
            rows[key].setdefault("similarity", {})[name] = measure
        eligible = [
            rows[key] for key, measure in measures.items()
            if measure["max_same_family_global_identity"] < panel["identity_below"]
            and measure["max_covered_same_family_local_identity"] is None
        ]
        distances[name] = _counts(eligible, cap)

    return {
        "schema": panel["schema"],
        "panel_sha256": _sha(PANEL.read_bytes()),
        "natural_panel_sha256": _sha(natural),
        "records": rows,
        "rejected": rejected,
        "summary": {
            "all_translatable": _counts(list(rows.values()), cap),
            "independently_corroborated": _counts([rows[target["id"]] for target in targets], cap),
            "remote_by_reference_set": distances,
            "gfp_like_without_natural_exact_match": sum(
                row["catalogue_group"] == "gfp_like" and not row["independent_matches"] for row in rows.values()
            ),
            "environment": {name: version(name) for name in ("biopython", "pyrodigal")},
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construct-repo", type=Path, required=True)
    parser.add_argument("--protein-repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.construct_repo, args.protein_repo, Path("."))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
