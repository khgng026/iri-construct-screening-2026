"""Compare direct translated search and candidate preprocessing on archived reporters."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
import hashlib
from importlib.metadata import version
from io import StringIO
import json
from pathlib import Path
import subprocess
import time

from experiments.paired_preprocessing.reporter_pilot import _candidates, _gene_caller_candidates

CONSTRUCT_REV = "55d4b762194055c8027fd6c780068d8a78a03282"
PROTEIN_REV = "799bd76889918bf1551c60016fb23278024cac73"
REFERENCES = ("avGFP", "cgreGFP", "ppluGFP1", "ppluGFP2", "LanYFP", "LanFP1", "LanFP2")
CAPS = (1, 2, 4, 8)
FIELDS = ("qseqid", "sseqid", "qstart", "qend", "bitscore", "evalue")
E_VALUE = "1e-5"


@dataclass(frozen=True)
class Region:
    left: int
    right: int
    strand: int

    @property
    def length(self):
        return self.right - self.left

    @property
    def frame(self):
        return (self.left if self.strand == 1 else self.right) % 3

    def overlap(self, other):
        if self.strand != other.strand or self.frame != other.frame:
            return 0
        return max(0, min(self.right, other.right) - max(self.left, other.left))

    def protein_span(self, start, end):
        """Locate an inclusive amino-acid alignment on the construct."""
        if not 1 <= start <= end or 3 * end > self.length:
            raise ValueError("Protein alignment lies outside the candidate")
        if self.strand == 1:
            return Region(self.left + 3 * (start - 1), self.left + 3 * end, 1)
        return Region(self.right - 3 * end, self.right - 3 * (start - 1), -1)


@dataclass(frozen=True)
class Hit:
    region: Region
    score: float
    reference: str
    candidate_index: int | None


def _sha(contents):
    return hashlib.sha256(contents).hexdigest()


def _git_source(repo, revision, path):
    return subprocess.check_output(["git", "-C", str(repo), "show", f"{revision}:{path}"])


def _reporter(record):
    markers = [feature for feature in record.features if feature.type == "CDS"]
    reporters = [feature for feature in markers if feature.qualifiers.get("label") != ["bla"]]
    if len(markers) != 2 or len(reporters) != 1:
        raise ValueError("The archived vector should contain one reporter and its backbone marker")
    feature = reporters[0]
    if len(feature.location.parts) != 1 or feature.location.strand not in (-1, 1):
        raise ValueError("Expected one complete reporter coding region")
    protein = str(feature.extract(record.seq).translate(table=1, cds=True))
    left, right = int(feature.location.start), int(feature.location.end)
    # The stop codon is outside the peptide sequence used to evaluate hits.
    strand = feature.location.strand
    region = Region(left + (3 if strand == -1 else 0), right - (3 if strand == 1 else 0), strand)
    return region, protein


def prepare(args):
    from Bio import SeqIO

    audit_bytes = args.audit.read_bytes()
    audit = json.loads(audit_bytes)
    if audit["schema"] != "iri-q3-open-reporters-feasibility-v1":
        raise ValueError("Expected the frozen Open Reporters audit")
    source = _git_source(args.protein_repo, PROTEIN_REV, "data/alignments/extant_fluorescent_proteins.fasta")
    if _sha(source) != audit["natural_panel_sha256"]:
        raise ValueError("Natural reference archive differs from the audit")
    by_name = {
        row.description.split("|", 1)[0].split()[0]: str(row.seq)
        for row in SeqIO.parse(StringIO(source.decode()), "fasta")
    }
    references = {name: by_name[name] for name in REFERENCES}
    inputs, targets = {}, {}
    for construct, entry in audit["records"].items():
        if not entry["independent_matches"] and entry["catalogue_group"] != "other_reporter":
            continue
        path = f"Open Reporters Collection/Plasmids_Genbank/{construct}.gb"
        contents = _git_source(args.construct_repo, CONSTRUCT_REV, path)
        if _sha(contents) != entry["source_sha256"]:
            raise ValueError(f"{construct}: the archived plasmid differs from the audit")
        record = SeqIO.read(StringIO(contents.decode()), "genbank")
        if record.id != construct:
            raise ValueError(f"{construct}: source record ID changed")
        region, protein = _reporter(record)
        label = "corroborated_gfp_like" if entry["independent_matches"] else "other_catalogue_reporter"
        if entry["independent_matches"] and not any(
            by_name[match["name"]] == protein for match in entry["independent_matches"]
        ):
            raise ValueError(f"{construct}: independent protein match did not reproduce")
        if protein in references.values():
            raise ValueError(f"{construct}: exact reporter present in the search reference library")
        inputs[construct] = str(record.seq).upper()
        targets[construct] = (label, region)
    return inputs, targets, references, _sha(audit_bytes)


def write_fasta(path, records):
    with path.open("x", encoding="ascii") as file:
        for name, sequence in records.items():
            file.write(f">{name}\n{sequence}\n")


def _redundancy_filtered(hits):
    accepted = []
    for hit in sorted(hits, key=lambda item: (
        -item.score, item.region.left, item.region.right, item.region.strand, item.reference,
    )):
        if any(
            hit.region.overlap(old.region) >= 0.8 * min(hit.region.length, old.region.length)
            for old in accepted
        ):
            continue
        accepted.append(hit)
    return accepted


def parse_hits(text, candidate_index):
    by_construct = defaultdict(list)
    for line in text.splitlines():
        query, reference, start, end, score, _ = line.split("\t")
        start, end = int(start), int(end)
        if candidate_index is None:
            loc = Region(min(start, end) - 1, max(start, end), 1 if start <= end else -1)
            construct, index = query, None
        else:
            construct, index, region = candidate_index[query]
            loc = region.protein_span(start, end)
        if loc.length <= 0 or loc.length % 3:
            raise ValueError("Expected a non-empty translated hit spanning complete codons")
        by_construct[construct].append(Hit(loc, float(score), reference, index))
    return by_construct


def search(args, engine, mode, query, refs, database, candidate_index):
    if engine == "blast":
        command = [
            str(args.blast_bin / mode), "-query", str(query), "-subject", str(refs),
            "-evalue", E_VALUE, "-outfmt", "6 " + " ".join(FIELDS),
        ]
    else:
        command = [
            str(args.diamond), mode, "--query", str(query), "--db", str(database),
            "--threads", "1", "--evalue", E_VALUE, "--max-target-seqs", "0",
            "--max-hsps", "0", "--outfmt", "6", *FIELDS, "--quiet",
        ]
        if args.diamond_sensitivity == "very-sensitive":
            command.append("--very-sensitive")
    started = time.perf_counter()
    run = subprocess.run(command, capture_output=True, text=True, check=True)
    elapsed = time.perf_counter() - started
    return parse_hits(run.stdout, candidate_index), elapsed


def evaluate(targets, hits, candidates=None, cap=None):
    positives = {key for key, (label, _) in targets.items() if label == "corroborated_gfp_like"}
    negatives = set(targets) - positives
    selected, reportable = {}, {}
    for key in targets:
        selected[key] = candidates[key][:cap] if candidates is not None else ()
        usable = (
            hit for hit in hits.get(key, ())
            if hit.candidate_index is None or cap is None or hit.candidate_index < cap
        )
        reportable[key] = _redundancy_filtered(usable)

    def reaches(key, regions, fraction):
        target = targets[key][1]
        return any(region.overlap(target) >= fraction * target.length for region in regions)

    rows = {}
    for report_cap in (1, 2, None):
        active = {key: values[:report_cap] for key, values in reportable.items()}
        rows[str(report_cap or "all")] = {
            "positives_recovered_at_80pct": sum(
                reaches(key, (hit.region for hit in active[key]), 0.8) for key in positives
            ),
            "positives_recovered_at_50pct": sum(
                reaches(key, (hit.region for hit in active[key]), 0.5) for key in positives
            ),
            "catalogue_other_reporter_locus_hits_at_50pct": sum(
                reaches(key, (hit.region for hit in active[key]), 0.5) for key in negatives
            ),
            "reported_loci": sum(len(values) for values in active.values()),
        }

    output = {"positive_reporters": len(positives), "catalogue_other_reporters": len(negatives), "reports": rows}
    if candidates is not None:
        output.update({
            "candidates": sum(len(rows) for rows in selected.values()),
            "candidate_residues": sum(region.length // 3 for rows in selected.values() for region in rows),
            "positive_candidate_reach": sum(reaches(key, selected[key], 0.8) for key in positives),
        })
    return output


def make_candidates(inputs, method, ranking=None):
    if method == "gene_caller":
        import pyrodigal

        finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    indices, proteins, regions = {}, {}, {}
    for construct, dna in inputs.items():
        units = _gene_caller_candidates(dna, finder) if method == "gene_caller" else _candidates(dna, False)
        if ranking is not None:
            from .ranking import Candidate, rank_indices

            geometry = [Candidate(unit["strand"], tuple(map(tuple, unit["spans"])), len(unit["protein"]))
                        for unit in units]
            units = [units[index] for index in rank_indices(geometry, ranking)]
        regions[construct] = []
        for index, unit in enumerate(units):
            left, right = unit["spans"][0]
            location = Region(left, right, unit["strand"])
            name = f"{construct}_p{index + 1}"
            if location.length != 3 * len(unit["protein"]):
                raise ValueError("Candidate coordinates do not match the translation length")
            indices[name] = (construct, index, location)
            proteins[name] = unit["protein"]
            regions[construct].append(location)
    return indices, proteins, regions


def run(args):
    inputs, targets, refs, audit_hash = prepare(args)
    args.work.mkdir(parents=True)
    raw, ref_path, database = (args.work / name for name in ("dna.fasta", "refs.fasta", "refs"))
    write_fasta(raw, inputs)
    write_fasta(ref_path, refs)
    subprocess.run(
        [str(args.diamond), "makedb", "--in", str(ref_path), "--db", str(database), "--quiet"],
        check=True, capture_output=True,
    )
    output = {
        "schema": "iri-benign-translated-search-v1",
        "construct_revision": CONSTRUCT_REV,
        "protein_revision": PROTEIN_REV,
        "prior_audit_sha256": audit_hash,
        "dna_sha256": _sha(raw.read_bytes()),
        "reference_sha256": _sha(ref_path.read_bytes()),
        "reference_count": len(refs),
        "evalue_ceiling": float(E_VALUE),
        "input_bases": sum(len(sequence) for sequence in inputs.values()),
        "environment": {
            "blastx": subprocess.check_output([str(args.blast_bin / "blastx"), "-version"], text=True).splitlines()[0],
            "diamond": subprocess.check_output([str(args.diamond), "version"], text=True).strip(),
            "diamond_sensitivity": args.diamond_sensitivity,
            "biopython": version("biopython"),
            "pyrodigal": version("pyrodigal"),
        },
        "arms": {},
    }
    timing = {}
    for engine in ("blast", "diamond"):
        hits, seconds = search(args, engine, "blastx", raw, ref_path, database, None)
        output["arms"][f"{engine}_whole_dna"] = evaluate(targets, hits)
        timing[f"{engine}_whole_dna"] = {"search_seconds": seconds}
    arms = [(method, None) for method in ("minimal", "gene_caller")]
    if args.ranking_followup:
        arms += [(method, "length_center_250") for method in ("minimal", "gene_caller")]
    for method, ranking in arms:
        started = time.perf_counter()
        indices, proteins, regions = make_candidates(inputs, method, ranking)
        preprocess_seconds = time.perf_counter() - started
        key = f"{method}_{ranking}" if ranking else method
        for cap in (*CAPS, None):
            tier = str(cap or "all")
            selected = {key: seq for key, seq in proteins.items() if cap is None or indices[key][1] < cap}
            query_path = args.work / f"{key}-{tier}.fasta"
            write_fasta(query_path, selected)
            for engine in ("blast", "diamond"):
                hits, seconds = search(args, engine, "blastp", query_path, ref_path, database, indices)
                arm = f"{engine}_{key}_{tier}"
                output["arms"][arm] = evaluate(targets, hits, regions, cap)
                timing[arm] = {"prepare_candidates_seconds": preprocess_seconds, "search_seconds": seconds}
    return output, timing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("audit", "construct-repo", "protein-repo", "blast-bin", "diamond", "work", "out"):
        parser.add_argument(f"--{option}", required=True, type=Path)
    parser.add_argument("--diamond-sensitivity", choices=("default", "very-sensitive"), default="default")
    parser.add_argument("--ranking-followup", action="store_true")
    args = parser.parse_args()
    result, timing = run(args)
    (args.work / "timing.json").write_text(json.dumps(timing, indent=2, sort_keys=True) + "\n")
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
