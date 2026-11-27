"""Run uncapped protein-family characterization on archived benign research records."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from functools import lru_cache
from hashlib import sha256
from importlib.metadata import version
from io import StringIO
import json
from os import link
from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
from time import perf_counter

from experiments.paired_preprocessing.cohort_study import candidate_units, chunk_unit
from experiments.paired_preprocessing.reporter_pilot import (
    AMINO_ACIDS, MAX_RESIDUES, MIN_RESIDUES, MODEL_HASH, PYRODIGAL_VERSION,
)
from experiments.preprocessing_exploration.translated_search import write_fasta


PROTEIN_ARMS = ("minimal", "gene_caller", "annotation_reference")
BLAST_FIELDS = "qseqid sseqid qstart qend qframe bitscore evalue"


def _fasta(contents):
    from Bio import SeqIO

    records = {}
    for record in SeqIO.parse(StringIO(contents.decode("ascii")), "fasta"):
        if record.id in records:
            raise ValueError(f"Repeated DNA input identifier {record.id}")
        records[record.id] = str(record.seq).upper()
    return records


def _references(rows):
    by_id, families = {}, {}
    for row in rows:
        ref_id, protein, family = row["accession"], row["sequence"], row["family"]
        if (not ref_id or any(char.isspace() for char in ref_id) or ref_id in by_id
                or not family or not MIN_RESIDUES <= len(protein) <= MAX_RESIDUES
                or set(protein) - AMINO_ACIDS):
            raise ValueError(f"Invalid reference entry {ref_id!r}")
        by_id[ref_id] = protein
        families[ref_id] = family
    if not by_id:
        raise ValueError("The reference panel is empty")
    return by_id, families


def blastx_calls(records, references, binary, *, evalue=10.0):
    """Return every reported hit before family/locus consolidation and thresholding."""
    proteins, families = _references(references)
    calls = {name: [] for name in records}
    if not records:
        raise ValueError("The DNA input is empty")
    executable = str(binary.resolve())
    with TemporaryDirectory() as temp:
        query, subjects = Path(temp) / "constructs.fa", Path(temp) / "references.fa"
        write_fasta(query, records)
        write_fasta(subjects, proteins)
        command = [executable, "-query", str(query), "-subject", str(subjects),
                   "-evalue", str(evalue), "-max_target_seqs", str(max(5, len(proteins))),
                   "-outfmt", f"6 {BLAST_FIELDS}"]
        started = perf_counter()
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
        elapsed = perf_counter() - started
    for line in completed.stdout.splitlines():
        query_id, ref_id, start, end, frame, bitscore, expect = line.split("\t")
        start, end, frame = int(start), int(end), int(frame)
        left, right, strand = min(start, end) - 1, max(start, end), 1 if frame > 0 else -1
        if not (1 <= abs(frame) <= 3 and 0 <= left < right <= len(records[query_id])
                and (right - left) % 3 == 0 and (start <= end) == (strand == 1)):
            raise ValueError("BLASTx returned invalid translated coordinates")
        calls[query_id].append({"family": families[ref_id], "strand": strand,
                                "spans": [[left, right]], "score": float(bitscore),
                                "evalue": float(expect), "reference_id": ref_id})
    version_text = subprocess.run([executable, "-version"], check=True, capture_output=True, text=True)
    return calls, {"version": version_text.stdout.splitlines()[0], "wall_seconds": elapsed,
                   "constructs": len(records), "input_nt": sum(map(len, records.values())),
                   "reference_count": len(proteins), "reporting_evalue_ceiling": evalue}


@contextmanager
def esm2_family_scores(references):
    import torch
    from iri_constructs.metrics import cosine, mean_pool
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    model = LocalAdapter("esm2", prepare_device("cpu", 20261122), ESM2_REVISION,
                         MAX_RESIDUES, offline=True)
    try:
        if model.metadata["weight_sha256"] != MODEL_HASH:
            raise ValueError("The pinned ESM-2 weights differ")
        proteins, names = _references(references)
        vectors = [(names[ref_id], mean_pool(model.residues(protein))) for ref_id, protein in proteins.items()]

        @lru_cache(maxsize=16_384)
        def rank(protein):
            embedding = mean_pool(model.residues(protein))
            best = {}
            for family, vector in vectors:
                best[family] = max(best.get(family, -1.0), cosine(embedding, vector))
            return sorted(best.items(), key=lambda item: (-item[1], item[0]))

        yield rank, model.metadata
    finally:
        model.close()


def _minimal_parents(units):
    """Recover sibling IDs from adjacent model windows of a single six-frame translation."""
    def coding_start(unit):
        lo, hi = unit["spans"][0]
        return lo if unit["strand"] == 1 else hi

    by_frame = {}
    parents = {}
    for unit in sorted(units, key=lambda row: row["strand"] * coding_start(row)):
        strand = unit["strand"]
        start = coding_start(unit)
        key = strand, start % 3
        length = len(unit["protein"])
        end, identifier = by_frame.get(key, (None, None))
        parent = identifier if end == start else unit["id"]
        parents[unit["id"]] = parent
        if length == MAX_RESIDUES:
            by_frame[key] = (start + strand * 3 * length, parent)
        else:
            by_frame.pop(key, None)
    return [unit | {"parent_id": parents[unit["id"]]} for unit in units]


def _units_with_parents(method, original, finder):
    if method == "minimal":
        return _minimal_parents(candidate_units(method, original, finder))
    if method == "annotation_reference":
        return [chunk | {"parent_id": unit["id"]} for unit in original for chunk in chunk_unit(unit)]

    from experiments.preprocessing_exploration.validation import called_genes, grouped_windows

    groups = grouped_windows(called_genes(original, finder), balanced=False)
    chunks = [unit | {"parent_id": f"gene{i + 1}"} for i, group in enumerate(groups) for unit in group]
    return [unit | {"id": f"q{i + 1}"} for i, unit in enumerate(chunks)]


class Checkpoint:
    """Persist completed records and bind them to one immutable scoring run."""

    def __init__(self, path, identity):
        self.path = path
        path.mkdir(parents=True, exist_ok=True)
        manifest = path / "manifest.json"
        self.identity = {"schema": "iri-whole-construct-checkpoint-v1", **identity}
        if not manifest.exists():
            try:
                self._write(manifest, self.identity, exclusive=True)
            except FileExistsError:
                pass  # Another local runner bound the directory before us.
        if json.loads(manifest.read_bytes()) != self.identity:
            raise ValueError("Checkpoint belongs to a different scoring run")

    def _entry(self, record):
        return self.path / f"{sha256(record.encode()).hexdigest()}.json"

    @staticmethod
    def _write(path, value, *, exclusive=False):
        handle = NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False)
        temporary = Path(handle.name)
        try:
            with handle:
                json.dump(value, handle, sort_keys=True, allow_nan=False)
                handle.write("\n")
            if exclusive:
                link(temporary, path)
            else:
                temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)

    def get(self, record, metadata):
        entry = self._entry(record)
        if not entry.exists():
            return None
        cached = json.loads(entry.read_bytes())
        if cached["id"] != record or cached["model"] != metadata or cached["run"] != self.identity:
            raise ValueError("Checkpoint entry differs from the active scoring model or record")
        return cached["calls"], cached["work"]

    def save(self, record, metadata, calls, work):
        self._write(self._entry(record), {"id": record, "model": metadata, "run": self.identity,
                                         "calls": calls, "work": work})


def esm2_calls(method, records, references, *, scorer=esm2_family_scores, checkpoint=None):
    if method not in PROTEIN_ARMS or not records:
        raise ValueError("A protein scoring arm and non-empty inputs are required")
    finder = None
    if method == "gene_caller":
        import pyrodigal

        if version("pyrodigal") != PYRODIGAL_VERSION:
            raise ValueError(f"Install pyrodigal=={PYRODIGAL_VERSION}")
        finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    calls, accounting = {}, {}
    reused = 0
    started = perf_counter()
    with scorer(references) as (rank, metadata):
        for record, original in records.items():
            cached = checkpoint.get(record, metadata) if checkpoint is not None else None
            if cached is not None:
                calls[record], accounting[record] = cached
                reused += 1
                continue
            prep_start = perf_counter()
            units = [unit for unit in _units_with_parents(method, original, finder)
                     if MIN_RESIDUES <= len(unit["protein"]) <= MAX_RESIDUES
                     and not set(unit["protein"]) - AMINO_ACIDS]
            score_start = perf_counter()
            record_calls = []
            for unit in units:
                ranking = rank(unit["protein"])
                family, score = ranking[0]
                record_calls.append({"query_id": unit["id"], "family": family, "score": score,
                                     "parent_id": unit["parent_id"],
                                     "strand": unit["strand"], "spans": unit["spans"],
                                     "ranking": ranking})
            calls[record] = record_calls
            accounting[record] = {
                "inputs": len(units), "residues": sum(len(unit["protein"]) for unit in units),
                "preprocess_seconds": score_start - prep_start,
                "scoring_seconds": perf_counter() - score_start,
            }
            if checkpoint is not None:
                checkpoint.save(record, metadata, calls[record], accounting[record])
    return calls, {"wall_seconds": perf_counter() - started, "model": metadata,
                   "reused_record_checkpoints": reused, "per_record": accounting}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=("blastx", *PROTEIN_ARMS), required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--blastx", type=Path)
    parser.add_argument("--checkpoint-dir", type=Path, help="Save or resume per-record ESM-2 scores")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("--out must point to a new file")
    input_bytes, refs_bytes = args.inputs.read_bytes(), args.references.read_bytes()
    references = json.loads(refs_bytes)
    identity = {"arm": args.arm, "input_sha256": sha256(input_bytes).hexdigest(),
                "references_sha256": sha256(refs_bytes).hexdigest()}
    if args.arm == "blastx":
        if args.blastx is None:
            parser.error("--blastx is required for the BLASTx arm")
        if args.checkpoint_dir is not None:
            parser.error("--checkpoint-dir is for the ESM-2 arms")
        calls, work = blastx_calls(_fasta(input_bytes), references, args.blastx)
    else:
        inputs = json.loads(input_bytes) if args.arm == "annotation_reference" else _fasta(input_bytes)
        checkpoint = Checkpoint(args.checkpoint_dir, identity) if args.checkpoint_dir is not None else None
        calls, work = esm2_calls(args.arm, inputs, references, checkpoint=checkpoint)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump({"schema": "iri-whole-construct-calls-v1", **identity, "calls": calls, "work": work},
                  output, sort_keys=True, allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    main()
