"""Compatibility helpers from the upstream reviewed-pair fusion experiment.

The old input-validation and numerical comparison APIs are retained for existing
callers/tests. Scientific inference uses the shared reviewed-input CLI and its
T02 checksum gate; this module has no model loader, downloader or executable CLI.
"""
from __future__ import annotations
import csv
import hashlib
import statistics
from pathlib import Path
from .metrics import cosine as _cosine


def cosine(vector_a, vector_b):
    return _cosine([float(x) for x in vector_a], [float(x) for x in vector_b])


AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


REQUIRED = (
    "pair_id", "isolated_id", "fusion_id", "start_aa", "end_aa",
    "isolated_accession", "fusion_accession",
    "isolated_source_url", "fusion_source_url", "boundary_evidence",
    "review_status", "licensing_status", "safety_review",
)


OPTIONAL = ("reference_id", "decoy_id")


def parse_fasta(path: Path) -> dict[str, str]:
    """Read unaligned protein FASTA, with a unique first header token per record."""
    proteins: dict[str, str] = {}
    current: str | None = None
    parts: list[str] = []

    def save() -> None:
        if current is not None:
            sequence = "".join(parts).upper()
            if not sequence:
                raise ValueError(f"Empty FASTA entry: {current}")
            if any(aa not in AA_ALPHABET for aa in sequence):
                raise ValueError(f"Unsupported amino-acid symbol in entry: {current}")
            proteins[current] = sequence

    with path.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                current = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                if not current or current in proteins:
                    raise ValueError("Missing or repeated FASTA identifier")
                parts = []
            elif current is None:
                raise ValueError("FASTA must begin with a header")
            elif any(character.isspace() for character in line):
                raise ValueError(f"Whitespace within FASTA sequence: {current}")
            else:
                parts.append(line)
    save()
    if not proteins:
        raise ValueError("No proteins found in FASTA")
    return proteins


def segment_bounds(start_aa: int, end_aa: int, full_length: int) -> tuple[int, int]:
    """Convert 1-based inclusive amino-acid bounds to Python slices."""
    if not (1 <= start_aa <= end_aa <= full_length):
        raise ValueError("Invalid annotated fusion protein bounds")
    return start_aa - 1, end_aa


def validate_pair(row: dict[str, str], proteins: dict[str, str]) -> dict:
    missing = [k for k in REQUIRED if not (row.get(k) or "").strip()]
    if missing:
        raise ValueError(f"Missing manifest fields: {', '.join(missing)}")
    for gate in ("review_status", "licensing_status", "safety_review"):
        if row[gate].strip().lower() != "approved":
            raise ValueError(f"Pair {row['pair_id']} is not approved at gate {gate}")
    for name in ("isolated_source_url", "fusion_source_url"):
        if not row[name].startswith(("https://", "http://")):
            raise ValueError(f"Pair {row['pair_id']} needs a public source URL: {name}")

    isolated_id = row["isolated_id"].strip()
    fusion_id = row["fusion_id"].strip()
    if isolated_id not in proteins or fusion_id not in proteins:
        raise ValueError(f"Pair {row['pair_id']} has missing FASTA identifiers")
    if isolated_id == fusion_id:
        raise ValueError(f"Pair {row['pair_id']} must use distinct entries")
    isolated, fusion = proteins[isolated_id], proteins[fusion_id]
    try:
        start_aa, end_aa = int(row["start_aa"]), int(row["end_aa"])
    except ValueError as exc:
        raise ValueError(f"Pair {row['pair_id']} has non-integer boundaries") from exc
    start0, end0 = segment_bounds(start_aa, end_aa, len(fusion))
    if len(fusion) <= len(isolated):
        raise ValueError(f"Pair {row['pair_id']} must contain a longer fusion sequence")
    if fusion[start0:end0] != isolated:
        raise ValueError(
            f"Pair {row['pair_id']}: annotated target region differs from isolated "
            "protein; check boundaries or use a separately documented variant design"
        )
    reference = (row.get("reference_id") or "").strip()
    decoy = (row.get("decoy_id") or "").strip()
    if bool(reference) != bool(decoy):
        raise ValueError("Provide both reference_id and decoy_id, or neither")
    if reference:
        if reference not in proteins or decoy not in proteins:
            raise ValueError(f"Pair {row['pair_id']} has missing reference FASTA IDs")
        if len({isolated_id, reference, decoy}) != 3:
            raise ValueError("Isolated, reference and decoy IDs must differ")
        if proteins[reference] == isolated or proteins[decoy] == isolated:
            raise ValueError("Reference and decoy must not be exact isolated duplicates")
        if proteins[reference] == proteins[decoy]:
            raise ValueError("Reference and decoy sequences must differ")
    return {
        "pair_id": row["pair_id"].strip(),
        "isolated_id": isolated_id,
        "fusion_id": fusion_id,
        "start_aa": start_aa,
        "end_aa": end_aa,
        "start0": start0,
        "end0": end0,
        "reference_id": reference,
        "decoy_id": decoy,
        "isolated_length": len(isolated),
        "fusion_length": len(fusion),
    }


def load_pairs(path: Path, proteins: dict[str, str]) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not set(REQUIRED).issubset(reader.fieldnames):
            raise ValueError(f"Manifest requires columns: {', '.join(REQUIRED)}")
        pairs = [validate_pair(row, proteins) for row in reader]
    if not pairs:
        raise ValueError("Manifest has no approved pairs")
    names = [p["pair_id"] for p in pairs]
    if len(set(names)) != len(names):
        raise ValueError("Repeated pair_id in manifest")
    return pairs


def compare_pair(pair: dict, proteins: dict[str, str], encode) -> tuple[dict, dict]:
    isolate = encode(proteins[pair["isolated_id"]])
    fusion = encode(proteins[pair["fusion_id"]])
    isolated_mean = isolate.mean(axis=0)
    fusion_mean = fusion.mean(axis=0)
    region_mean = fusion[pair["start0"]:pair["end0"]].mean(axis=0)

    global_cos = cosine(isolated_mean, fusion_mean)
    region_cos = cosine(isolated_mean, region_mean)
    result = {
        "pair_id": pair["pair_id"],
        "isolated_id": pair["isolated_id"],
        "fusion_id": pair["fusion_id"],
        "isolated_length_aa": pair["isolated_length"],
        "fusion_length_aa": pair["fusion_length"],
        "region_start_aa": pair["start_aa"],
        "region_end_aa": pair["end_aa"],
        "cos_isolated_vs_whole_fusion": global_cos,
        "cos_isolated_vs_contextual_region": region_cos,
        "contextual_minus_whole": region_cos - global_cos,
        "margin_isolated": "",
        "margin_whole_fusion": "",
        "margin_contextual_region": "",
        "reference_preferred_isolated": "",
        "reference_preferred_whole_fusion": "",
        "reference_preferred_contextual_region": "",
    }
    vectors = {
        "isolated": isolated_mean,
        "whole_fusion": fusion_mean,
        "contextual_region": region_mean,
    }

    if pair["reference_id"]:
        reference = encode(proteins[pair["reference_id"]]).mean(axis=0)
        decoy = encode(proteins[pair["decoy_id"]]).mean(axis=0)
        for label, vec in (
            ("isolated", isolated_mean),
            ("whole_fusion", fusion_mean),
            ("contextual_region", region_mean),
        ):
            margin = cosine(vec, reference) - cosine(vec, decoy)
            result[f"margin_{label}"] = margin
            result[f"reference_preferred_{label}"] = int(margin > 0)
        vectors["reference"] = reference
        vectors["decoy"] = decoy

    return result, vectors


def average(records: list[dict], key: str):
    values = [float(record[key]) for record in records if record[key] != ""]
    return statistics.mean(values) if values else None


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            sha.update(chunk)
    return sha.hexdigest()
