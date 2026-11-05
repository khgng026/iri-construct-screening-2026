#!/usr/bin/env python3
"""T05: compare an isolated benign protein with its region inside an annotated fusion.

No sequences are downloaded or committed. Use only locally reviewed public records.
Coordinates are 1-based, inclusive amino-acid positions in the FUSION sequence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path

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


def cosine(vector_a, vector_b) -> float:
    """Cosine similarity; compatible with numpy vectors or Python lists."""
    a = [float(x) for x in vector_a]
    b = [float(x) for x in vector_b]
    if len(a) != len(b) or not a:
        raise ValueError("Incompatible embedding dimensions")
    ab = sum(x * y for x, y in zip(a, b))
    aa = math.sqrt(sum(x * x for x in a))
    bb = math.sqrt(sum(y * y for y in b))
    if aa == 0 or bb == 0:
        raise ValueError("Zero-norm embedding")
    return max(-1.0, min(1.0, ab / (aa * bb)))


def select_device(torch, requested: str):
    if requested != "auto":
        if requested == "cuda" and not torch.cuda.is_available():
            raise ValueError("CUDA requested but unavailable")
        if requested == "mps" and not torch.backends.mps.is_available():
            raise ValueError("MPS requested but unavailable")
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def create_encoder(model_name: str, requested_device: str):
    try:
        import torch
        import transformers
        from transformers import AutoTokenizer, EsmModel
    except ImportError as exc:
        raise SystemExit(
            "Install model dependencies: pip install -r "
            "experiments/t05_esm2_fusion/requirements.txt"
        ) from exc

    if not model_name.startswith("facebook/esm2_"):
        raise ValueError("Only Facebook ESM-2 model identifiers are supported")
    device = select_device(torch, requested_device)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name).to(device).eval()
    max_aa = min(1022, int(getattr(model.config, "max_position_embeddings", 1026)) - 2)
    cache: dict[str, object] = {}
    calls = 0

    def encode(seq: str):
        nonlocal calls
        if seq not in cache:
            if len(seq) > max_aa:
                raise ValueError(
                    f"One protein is {len(seq)} aa, exceeding the model's "
                    f"{max_aa}-residue limit; no truncation is performed"
                )
            tokens = tokenizer(seq, return_tensors="pt", add_special_tokens=True)
            if tokens["input_ids"].shape[1] != len(seq) + 2:
                raise ValueError("Unexpected tokenization length: residue alignment unsafe")
            tokens = {key: value.to(device) for key, value in tokens.items()}
            with torch.inference_mode():
                hidden = model(**tokens).last_hidden_state[0, 1:len(seq) + 1]
            cache[seq] = hidden.float().cpu().numpy().copy()
            calls += 1
        return cache[seq]

    def runtime_stats():
        return {
            "torch_version": torch.__version__,
            "transformers_version": transformers.__version__,
            "device": str(device),
            "max_residues": max_aa,
            "unique_model_inferences": calls,
            "cuda_peak_memory_bytes": (
                int(torch.cuda.max_memory_allocated(device))
                if device.type == "cuda" else None
            ),
        }

    return encode, runtime_stats


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True, help="Reviewed local protein FASTA")
    parser.add_argument("--pairs", type=Path, required=True, help="Reviewed local pair manifest CSV")
    parser.add_argument("--outdir", type=Path, default=Path("results/t05_esm2_fusion"))
    parser.add_argument("--model", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--device", choices=("auto", "cpu", "cuda", "mps"), default="auto")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    start = time.perf_counter()
    proteins = parse_fasta(args.fasta)
    pairs = load_pairs(args.pairs, proteins)
    print(f"Validated {len(pairs)} annotated benign pairs and {len(proteins)} FASTA records")
    if args.validate_only:
        return

    encode, runtime_stats = create_encoder(args.model, args.device)
    metrics: list[dict] = []
    arrays = {}
    index_rows = []
    for i, pair in enumerate(pairs):
        result, vectors = compare_pair(pair, proteins, encode)
        metrics.append(result)
        for label, vec in vectors.items():
            key = f"pair_{i:04d}_{label}"
            arrays[key] = vec
            index_rows.append({"pair_id": pair["pair_id"], "kind": label, "npz_key": key})
        print(f"Completed pair {i + 1}/{len(pairs)}: {pair['pair_id']}")

    import numpy as np

    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "pair_metrics.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    with (args.outdir / "embedding_index.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=("pair_id", "kind", "npz_key"))
        writer.writeheader()
        writer.writerows(index_rows)
    np.savez_compressed(args.outdir / "embeddings.npz", **arrays)
    metadata = {
        "experiment": "T05 annotated isolated protein versus fusion",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model": args.model,
        "pair_count": len(pairs),
        "pairs_with_reference_controls": sum(bool(p["reference_id"]) for p in pairs),
        "mean_cos_isolated_vs_whole_fusion": average(metrics, "cos_isolated_vs_whole_fusion"),
        "mean_cos_isolated_vs_contextual_region": average(metrics, "cos_isolated_vs_contextual_region"),
        "mean_contextual_minus_whole": average(metrics, "contextual_minus_whole"),
        "mean_margin_isolated": average(metrics, "margin_isolated"),
        "mean_margin_whole_fusion": average(metrics, "margin_whole_fusion"),
        "mean_margin_contextual_region": average(metrics, "margin_contextual_region"),
        "elapsed_seconds_including_model_load": round(time.perf_counter() - start, 3),
        "fasta_sha256": digest(args.fasta),
        "pairs_sha256": digest(args.pairs),
        "interpretation_limit": (
            "Cosine similarity is representational, not evidence of biological function. "
            "Reference margins compare two benign examples, not a validated function classifier. "
            "All fusion boundaries were taken from reviewed annotations."
        ),
        **runtime_stats(),
    }
    (args.outdir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote metrics, embeddings and metadata to {args.outdir}")


if __name__ == "__main__":
    main()
