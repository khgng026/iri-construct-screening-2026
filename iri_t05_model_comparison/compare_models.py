#!/usr/bin/env python3
"""IRI T05: benign construct annotation / protein embedding retrieval pilot.

This compares pretrained protein *representations*, not intrinsic function prediction.
No hazardous sequences, custom training, or generated biology are involved.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import sys
import time
import traceback
from collections import Counter, defaultdict
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

if __name__ == "__main__":
    # The historical functions below are retained for review. Executable use is
    # routed through the reviewed local-input workflow before optional imports.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from iri_constructs.experiment_cli import main as reviewed_main
    raise SystemExit(reviewed_main(default_experiment="comparison"))

import numpy as np

# Explicit benign, publicly documented GenBank vectors / reporters.
# The set is intentionally small: designed as a feasibility spike, not a validation set.
ACCESSIONS = [
    "U55762.1",  # pEGFP-N1
    "U55763.1",  # pEGFP-C1
    "U47295.2",  # pGL3-Basic
    "U47296",    # pGL3-Control (allow NCBI to resolve current version)
    "J01749.1",  # pBR322
    "L09137",    # pUC19
    "L08752.1",  # pUC18
    "U26267.1",  # pBluescript
    "U57608",    # pEGFP-N2
    "MK301207.1", # dsRed2 fluorescent reporter synthetic construct
    "FJ226077.1", # DsRed-Express2 fluorescent reporter synthetic construct
    "MW139745.1", # NanoLuc luciferase synthetic construct
    "U47122.2",   # pSP-luc+
    "AY738223.1", # pGL4.11 luc2P reporter vector
    "DQ188837.1", # pGL4.17 luc2/Neo reporter vector
    "DQ188844.1", # pGL4.80 Renilla luciferase reporter vector
    "MK484108.1", # pGL4Zc firefly luciferase reporter vector
]

# Annotation heuristics, not experimental evidence of biochemical function.
# Only known, benign reporter or vector-marker features are retained.
TRAIT_PATTERNS = [
    ("fluorescent_reporter", re.compile(r"\b(?:egfp|eyfp|ecfp|gfp|dsred(?:[a-z0-9-]*)?|mcherry|fluorescent protein|green fluorescent|red fluorescent|cyan fluorescent|yellow fluorescent)\b", re.I)),
    ("luminescent_reporter", re.compile(r"\b(?:luciferase|firefly luciferase|luc2p?|nluc|hrluc|rluc|luc\+?)\b", re.I)),
    ("beta_lactamase_marker", re.compile(r"\b(?:beta[- ]lactamase|bla|ampicillin resistance|ampicillin-resistance)\b", re.I)),
    ("aminoglycoside_marker", re.compile(r"\b(?:neomycin phosphotransferase|kanamycin resistance|kanamycin-resistance|neo\b|nptii\b|aminoglycoside phosphotransferase)\b", re.I)),
    ("gst_tag", re.compile(r"\b(?:glutathione[- ]s[- ]transferase|gst tag|gst-tag)\b", re.I)),
]

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYUOXB ZJ".replace(" ", ""))
MODEL_NAMES = ("esm2", "esmc", "esm3")
MODEL_IDS = {
    "esm2": "facebook/esm2_t30_150M_UR50D",
    "esmc": "esmc_300m (Biohub/ESMC via esm package)",
    "esm3": "esm3-sm-open-v1 (Biohub/esm3-sm-open-v1)",
}


def write_csv(path: Path, records: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(records[0]) if records else []
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def classify_annotation(feature) -> tuple[str, str]:
    annotations = feature.qualifiers
    product = " ".join(annotations.get("product", [])).strip()
    gene = " ".join(annotations.get("gene", [])).strip()
    label_text = f"{product} {gene}"
    hits = [trait for trait, pattern in TRAIT_PATTERNS if pattern.search(label_text)]
    if len(hits) != 1:
        return "", (product or gene)
    return hits[0], (product or gene)


def get_record(accession: str, data_dir: Path, email: str = ""):
    from Bio import SeqIO
    from io import StringIO
    path = data_dir / f"{accession}.gb"
    if not path.is_file():
        params = {
            "db": "nuccore", "id": accession, "rettype": "gbwithparts",
            "retmode": "text", "tool": "iri_t05_benign_comparison",
        }
        if email:
            params["email"] = email
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(params)
        err = None
        for attempt in range(3):
            try:
                req = Request(url, headers={"User-Agent": "iri-t05-model-comparison/0.1"})
                with urlopen(req, timeout=35) as reply:
                    data = reply.read().decode("utf-8")
                record = SeqIO.read(StringIO(data), "genbank")
                if len(record.seq) < 500 or accession.split(".")[0] not in record.id:
                    raise ValueError("unexpected GenBank record: " + record.id)
                path.write_text(data, encoding="utf-8")
                err = None
                break
            except (HTTPError, URLError, TimeoutError, ValueError, StopIteration) as e:
                err = e
                time.sleep(1.5 * (attempt + 1))
        if err is not None:
            raise RuntimeError(f"Could not download {accession}: {err}") from err
    return SeqIO.read(path, "genbank"), path


def prepare(data_dir: Path, email: str = "", allow_partial: bool = False) -> list[dict]:
    parts = []
    omitted = []
    status = []
    for accession in ACCESSIONS:
        try:
            record, source_path = get_record(accession, data_dir / "genbank", email)
        except Exception as exc:
            status.append({"accession": accession, "status": "download_failed", "detail": str(exc)})
            print(f"[WARN] Skipping {accession}: {exc}", file=sys.stderr)
            continue
        n_found = 0
        for n, f in enumerate(record.features):
            if f.type != "CDS":
                continue
            trait, product = classify_annotation(f)
            q = f.qualifiers
            sequence = "".join(q.get("translation", [])).strip().rstrip("*").upper()
            # No ORF inference or speculative translations: use annotated protein only.
            if not sequence or len(sequence) < 40 or len(sequence) > 900 or set(sequence) - AA_ALLOWED:
                omitted.append({"accession": record.id, "index": n, "product": product,
                                "why": "missing/short/long/unsupported annotated protein"})
                continue
            if not trait:
                omitted.append({"accession": record.id, "index": n, "product": product,
                                "why": "outside conservative benign trait vocabulary"})
                continue
            part_id = f"{record.id}_CDS{n}"
            parts.append({
                "part_id": part_id, "construct_id": record.id, "construct_name": record.description,
                "cds_feature_index": n, "annotated_product": product, "trait": trait,
                "protein": sequence, "aa_length": len(sequence),
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "source_url": "https://www.ncbi.nlm.nih.gov/nuccore/" + record.id,
                "source_genbank_path": str(source_path.relative_to(data_dir)),
                "label_source": "GenBank product/gene qualifier; heuristic mapping; not independently reviewed",
            })
            n_found += 1
        status.append({"accession": record.id, "status": "ok", "detail": f"{n_found} labeled CDS"})
    write_csv(data_dir / "parts.csv", parts, [
        "part_id", "construct_id", "construct_name", "cds_feature_index", "annotated_product",
        "trait", "protein", "aa_length", "sequence_sha256", "source_url",
        "source_genbank_path", "label_source",
    ])
    write_csv(data_dir / "omitted_features.csv", omitted, ["accession", "index", "product", "why"])
    write_csv(data_dir / "download_status.csv", status, ["accession", "status", "detail"])
    print(f"Prepared {len(parts)} CDS from {len(set(p['construct_id'] for p in parts))} constructs; "
          f"{len(set(p['trait'] for p in parts))} traits. Saved {data_dir / 'parts.csv'}")
    for trait, count in sorted(Counter(p["trait"] for p in parts).items()):
        print(f"  {trait}: {count}")
    if len(parts) < 2 and not allow_partial:
        raise RuntimeError("Insufficient data. Check NCBI connectivity / annotated CDS coverage. "
                           "You can add benign annotated GenBank files to the local cache or expand ACCESSIONS.")
    return parts


def load_model(name: str, device: str):
    import torch
    if name == "esm2":
        from transformers import AutoModel, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(MODEL_IDS[name])
        model = AutoModel.from_pretrained(MODEL_IDS[name]).to(device).eval()
        return model, tok
    if name == "esmc":
        from esm.models.esmc import ESMC
        model = ESMC.from_pretrained("esmc_300m").to(device).eval()
        return model, None
    if name == "esm3":
        from esm.models.esm3 import ESM3
        model = ESM3.from_pretrained("esm3-sm-open-v1").to(device).eval()
        return model, None
    raise ValueError(name)


def per_token_embedding(name: str, model, tokenizer, sequence: str, device: str):
    import torch
    with torch.inference_mode():
        if name == "esm2":
            data = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
            data = {k: v.to(device) for k, v in data.items()}
            states = model(**data).last_hidden_state[0]
            ids = data["input_ids"][0].tolist()
            keep = np.array([not v for v in tokenizer.get_special_tokens_mask(
                ids, already_has_special_tokens=True)], dtype=bool)
            keep &= data["attention_mask"][0].bool().cpu().numpy()
            res = states[torch.tensor(keep, device=states.device)]
        else:
            from esm.sdk.api import ESMProtein, LogitsConfig
            pt = model.encode(ESMProtein(sequence=sequence))
            output = model.logits(pt, LogitsConfig(sequence=True, return_embeddings=True))
            if output.embeddings is None:
                raise RuntimeError(f"{name} returned no embeddings; check esm package version")
            res = output.embeddings
            if res.ndim == 3:
                if res.shape[0] != 1:
                    raise ValueError("Expected single-sequence inference")
                res = res[0]
            if res.shape[0] == len(sequence) + 2:
                res = res[1:-1]  # BOS, EOS
            elif res.shape[0] != len(sequence):
                raise RuntimeError(f"Unexpected token count: {res.shape[0]} for {len(sequence)} residues")
        if res.shape[0] != len(sequence):
            raise RuntimeError(f"Embedding/token count mismatch: {res.shape[0]} vs {len(sequence)}")
        return res.float().mean(dim=0).cpu().numpy()


def normalize(v: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(v))
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError("Invalid embedding norm")
    return (v / norm).astype(np.float32)


def extract_embeddings(name: str, parts: list[dict], output: Path, device: str):
    import torch
    unique = {}
    for part in parts:
        unique.setdefault(part["sequence_sha256"], part["protein"])
    started = time.perf_counter()
    model, tokenizer = load_model(name, device)
    load_s = time.perf_counter() - started
    if device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats()
    embeddings = {}
    times = []
    started = time.perf_counter()
    for i, (key, sequence) in enumerate(unique.items(), 1):
        t0 = time.perf_counter()
        embedding = per_token_embedding(name, model, tokenizer, sequence, device)
        if device.startswith("cuda"):
            torch.cuda.synchronize()
        times.append(time.perf_counter() - t0)
        embeddings[key] = normalize(embedding)
        print(f"[{name}] {i}/{len(unique)} proteins; {times[-1]:.2f}s", flush=True)
    elapsed = time.perf_counter() - started
    gpu_gb = float(torch.cuda.max_memory_allocated() / 1e9) if device.startswith("cuda") else None
    # Store by digest, not by row order (and no guessed dimension across model families).
    np.savez_compressed(output / f"embeddings_{name}.npz", **embeddings)
    info = {
        "model": name, "checkpoint": MODEL_IDS[name], "status": "ok", "device": device,
        "model_load_seconds_excluding_earlier_downloads": round(load_s, 3),
        "inference_seconds_unique_proteins": round(elapsed, 3),
        "median_seconds_per_unique_protein": round(statistics.median(times), 3),
        "unique_proteins": len(unique), "peak_inference_gpu_allocated_gb": gpu_gb,
        "embedding_dimension": len(next(iter(embeddings.values()))),
        "torch_version": torch.__version__,
    }
    del model
    if device.startswith("cuda"):
        torch.cuda.empty_cache()
    return embeddings, info


def eligibility(parts: list[dict], avoid_duplicates: bool) -> dict[str, list[str]]:
    """All queries scored on same candidate references for all models."""
    refs = {}
    for part in parts:
        alternatives = [p for p in parts
                        if p["construct_id"] != part["construct_id"]
                        and (not avoid_duplicates or p["sequence_sha256"] != part["sequence_sha256"])]
        # The classification comparison requires at least 2 eligible ref traits, including the truth.
        if (len({p["trait"] for p in alternatives}) >= 2
                and any(p["trait"] == part["trait"] for p in alternatives)):
            refs[part["part_id"]] = [p["part_id"] for p in alternatives]
    return refs


def score_one_model(parts: list[dict], embeddings: dict[str, np.ndarray], mode: str) -> tuple[list[dict], dict]:
    pmap = {p["part_id"]: p for p in parts}
    elig = eligibility(parts, avoid_duplicates=(mode == "distinct_sequence"))
    predictions = []
    for qid, ref_ids in elig.items():
        q = pmap[qid]
        qvec = embeddings.get(q["sequence_sha256"])
        if qvec is None:
            continue
        ranked = []
        for rid in ref_ids:
            ref = pmap[rid]
            rvec = embeddings.get(ref["sequence_sha256"])
            if rvec is None:
                continue
            sim = float(np.dot(qvec, rvec))
            ranked.append((sim, ref))
        if len(ranked) < 2:
            continue
        ranked.sort(key=lambda x: (-x[0], x[1]["part_id"]))
        top_score, top = ranked[0]
        # This is a nearest-reference label, NOT the model's own function prediction.
        predictions.append({
            "mode": mode, "part_id": qid, "construct_id": q["construct_id"],
            "true_trait": q["trait"], "predicted_trait": top["trait"],
            "correct": int(q["trait"] == top["trait"]), "similarity": round(top_score, 6),
            "neighbor_part_id": top["part_id"], "neighbor_construct": top["construct_id"],
            "reference_n": len(ranked),
        })
    if not predictions:
        return [], {"mode": mode, "n": 0, "accuracy": None, "macro_f1": None,
                    "construct_recall": None, "construct_precision": None,
                    "construct_f1": None, "eligible_queries_in_manifest": len(elig)}
    true = [p["true_trait"] for p in predictions]
    pred = [p["predicted_trait"] for p in predictions]
    labels = sorted(set(true) | set(pred))
    f1_by_label = {}
    for lab in labels:
        tp = sum(a == lab and b == lab for a, b in zip(true, pred))
        fp = sum(a != lab and b == lab for a, b in zip(true, pred))
        fn = sum(a == lab and b != lab for a, b in zip(true, pred))
        f1_by_label[lab] = (2 * tp / (2 * tp + fp + fn)) if (2 * tp + fp + fn) else 0
    # Only use the set of *eligible* parts for each construct in retrieval evaluation.
    expected = defaultdict(set)
    inferred = defaultdict(set)
    for p in predictions:
        expected[p["construct_id"]].add(p["true_trait"])
        inferred[p["construct_id"]].add(p["predicted_trait"])
    c_recall, c_precision, c_f1 = [], [], []
    for cid, traits in expected.items():
        guess = inferred[cid]
        intersection = len(guess & traits)
        prec = intersection / len(guess) if guess else 0.0
        rec = intersection / len(traits) if traits else 0.0
        c_recall.append(rec)
        c_precision.append(prec)
        c_f1.append((2*prec*rec/(prec+rec)) if prec+rec else 0)
    metrics = {
        "mode": mode,
        "n": len(predictions),
        "n_constructs": len(expected),
        "n_true_traits": len(set(true)),
        "eligible_queries_in_manifest": len(elig),
        "accuracy": round(sum(p["correct"] for p in predictions) / len(predictions), 4),
        "macro_f1": round(float(statistics.mean(f1_by_label.values())), 4),
        "construct_recall": round(float(statistics.mean(c_recall)), 4),
        "construct_precision": round(float(statistics.mean(c_precision)), 4),
        "construct_f1": round(float(statistics.mean(c_f1)), 4),
    }
    return predictions, metrics


def whole_construct_baseline(parts: list[dict]) -> list[dict]:
    """A deterministic single-CDS baseline, independent of neural model predictions."""
    by_construct = defaultdict(list)
    for part in parts:
        by_construct[part["construct_id"]].append(part)
    rows = []
    for cid, group in by_construct.items():
        longest = max(group, key=lambda part: len(part["protein"]))
        true_traits = {p["trait"] for p in group}
        rows.append({
            "construct_id": cid, "n_labeled_cds": len(group), "n_true_traits": len(true_traits),
            "all_annotated_trait_set": ";".join(sorted(true_traits)),
            "longest_cds_trait": longest["trait"],
            "single_cds_trait_recall_upper_bound": round(1/len(true_traits), 4),
        })
    return rows


def recommend(stats: list[dict], parts: list[dict], info: list[dict], mode: str = "distinct_sequence") -> str:
    filtered = {s["model"]: s for s in stats if s["mode"] == mode and s["macro_f1"] is not None}
    good = {k: v for k, v in filtered.items()
            if v["n"] >= 12 and v["n_true_traits"] >= 3 and v["n_constructs"] >= 4}
    if len(good) < 2:
        return "INCONCLUSIVE: fewer than two models have >=12 comparable distinct-sequence test queries " \
               "spanning >=3 functional labels and >=4 constructs. Expand curated benign examples."
    # Require same query set for comparable model tests and choose cheapest sufficiently accurate.
    if len({(v["n"], v["n_constructs"], v["n_true_traits"]) for v in good.values()}) != 1:
        return "INCONCLUSIVE: model tests have different eligible sample counts. Do not compare their scores directly."
    best = max(v["macro_f1"] for v in good.values())
    top = [k for k, v in good.items() if v["macro_f1"] >= best - 0.05]
    # Compare median measured inference time over the same unique data; lower if near-best.
    timings = {p["model"]: p.get("median_seconds_per_unique_protein", float("inf"))
               for p in info if p.get("status") == "ok"}
    winner = min(top, key=lambda k: (timings.get(k, float("inf")),
                                     {"esm2": 150, "esmc": 300, "esm3": 1400}[k]))
    if len(top) >= 2:
        return f"PROVISIONAL: {winner} has the lowest measured median inference latency among models " \
               f"within 0.05 macro-F1 of the best observed result. Confirm with a larger independently " \
               f"curated benchmark before treating this as a scientific conclusion."
    return f"PROVISIONAL: {winner} has the highest distinct-sequence macro-F1 by >0.05; " \
           "confirm the difference on a larger, independently curated benchmark."


def run(data_dir: Path, output_dir: Path, model_list: list[str], device: str):
    import torch
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    parts = read_csv(data_dir / "parts.csv")
    if not parts:
        raise RuntimeError("No annotated CDS in input. Run --prepare first.")
    output_dir.mkdir(parents=True, exist_ok=True)
    info, scores, all_predictions = [], [], []
    for name in model_list:
        print(f"\n== {name} ({device}) ==", flush=True)
        try:
            embs, details = extract_embeddings(name, parts, output_dir, device)
            info.append(details)
            for mode in ("allow_reused_sequence", "distinct_sequence"):
                pred, metrics = score_one_model(parts, embs, mode)
                scores.append({"model": name, **metrics})
                all_predictions += [{"model": name, **p} for p in pred]
                print(f"  {mode}: N={metrics['n']}, macro-F1={metrics['macro_f1']}, "
                      f"construct recall={metrics['construct_recall']}", flush=True)
        except Exception as exc:
            print(f"[ERROR] {name}: {type(exc).__name__}: {exc}", file=sys.stderr)
            info.append({"model": name, "checkpoint": MODEL_IDS[name],
                         "status": "failed", "error": repr(exc),
                         "traceback": traceback.format_exc(limit=3)})
        finally:
            if device.startswith("cuda"):
                torch.cuda.empty_cache()
    write_csv(output_dir / "model_metrics.csv", scores, [
        "model", "mode", "n", "n_constructs", "n_true_traits", "eligible_queries_in_manifest",
        "accuracy", "macro_f1", "construct_recall", "construct_precision", "construct_f1",
    ])
    write_csv(output_dir / "predictions.csv", all_predictions, [
        "model", "mode", "part_id", "construct_id", "true_trait", "predicted_trait",
        "correct", "similarity", "neighbor_part_id", "neighbor_construct", "reference_n",
    ])
    write_csv(output_dir / "single_cds_baseline.csv", whole_construct_baseline(parts))
    (output_dir / "model_run_info.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    decision = recommend(scores, parts, info)
    report = [
        "# IRI T05: preliminary model comparison",
        "",
        f"**Decision:** {decision}",
        "",
        "## Evaluation protocol",
        "",
        "- Benign public plasmid records; translated proteins from annotated GenBank CDS only.",
        "- Labels mapped heuristically from GenBank gene/product text. Verify before publication.",
        "- Protein representations are residue-mean-pooled, L2-normalized, cosine compared.",
        "- Leave-one-CONSTRUCT-out nearest-reference classification (no same-construct references).",
        "- Primary score (`distinct_sequence`) excludes identical protein sequences in other constructs.",
        "- `allow_reused_sequence` is diagnostic; it can be severely inflated by shared vector genes.",
        "- Construct trait recall is evaluated over eligible labeled CDS only, not unannotated biology.",
        "- These are retrieval-based trait predictions, not direct biochemical function predictions.",
        "- No protein model was run on full plasmid DNA. The single-CDS baseline is non-neural.",
        "",
        "## Results",
        "",
        "| Model | Mode | Scored CDS | Macro-F1 | Construct recall |",
        "|---|---|---:|---:|---:|",
    ]
    for x in scores:
        report.append(f"| {x['model']} | {x['mode']} | {x['n']} | {x['macro_f1']} | {x['construct_recall']} |")
    report += ["", "## Caveats", "", "Small, non-random convenience panel; heuristic labels; "
               "closely related vectors; no wet-lab functionality; few or zero distinct homologs possible. "
               "An INCONCLUSIVE verdict is informative and should trigger sample expansion, not model selection.",
               "", "See `predictions.csv`, `model_metrics.csv`, `model_run_info.json`, "
               "`single_cds_baseline.csv` and `data/parts.csv` for auditing.", ""]
    (output_dir / "decision.md").write_text("\n".join(report), encoding="utf-8")
    print("\n" + decision)
    print("Results at:", output_dir.resolve())


def self_test():
    """Sanity test of labeling/evaluation only; does NOT test ESM inference."""
    from types import SimpleNamespace
    example = SimpleNamespace(qualifiers={"product": ["enhanced green fluorescent protein"], "gene": ["egfp"]})
    assert classify_annotation(example)[0] == "fluorescent_reporter"
    example2 = SimpleNamespace(qualifiers={"product": ["some unknown enzyme"]})
    assert classify_annotation(example2)[0] == ""
    def sample(name, construct, label, emb):
        key = hashlib.sha256(name.encode()).hexdigest()
        return {"part_id": name, "construct_id": construct, "trait": label,
                "sequence_sha256": key, "protein": "A" * 100}, key, normalize(np.asarray(emb, dtype=float))
    rows = [sample("a1", "A", "fluorescent_reporter", [1, 0, 0]),
            sample("b1", "B", "fluorescent_reporter", [1, 0.1, 0]),
            sample("c1", "C", "luminescent_reporter", [0, 1, 0]),
            sample("d1", "D", "luminescent_reporter", [0.1, 1, 0])]
    parts = [x[0] for x in rows]
    embs = {x[1]: x[2] for x in rows}
    predictions, result = score_one_model(parts, embs, "distinct_sequence")
    assert result["n"] == 4, result
    assert result["accuracy"] == 1 and result["macro_f1"] == 1.0, result
    assert len(predictions) == 4
    assert all(p["construct_id"] != p["neighbor_construct"] for p in predictions)
    # Exact cross-construct duplicates excluded from distinct-sequence mode.
    parts[1]["sequence_sha256"] = parts[0]["sequence_sha256"]
    assert len(eligibility(parts, avoid_duplicates=True)) <= 2
    print("SELF-TEST PASS: annotation mapping, nearest-neighbor scoring, leave-construct-out, duplicate exclusion")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true", help="Fetch annotated benign GenBank plasmids and create CSV")
    mode.add_argument("--run", action="store_true", help="Run selected ESM models and evaluate retrieval")
    mode.add_argument("--self-test", action="store_true", help="Offline scoring sanity test (no model weights)")
    parser.add_argument("--data", type=Path, default=Path("data"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--email", default="", help="Optional identifying email for NCBI e-utilities")
    parser.add_argument("--models", default="esm2,esmc,esm3", help="Comma-separated: esm2,esmc,esm3")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.prepare:
        prepare(args.data, args.email)
    elif args.run:
        names = [x.strip() for x in args.models.split(",") if x.strip()]
        if not names or any(x not in MODEL_NAMES for x in names) or len(names) != len(set(names)):
            parser.error("--models must list unique values among esm2,esmc,esm3")
        run(args.data, args.output, names, args.device)
