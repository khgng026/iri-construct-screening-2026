#!/usr/bin/env python3
"""Fetch public, versioned benign archive candidates. Nothing is automatically admitted.

An exact source is saved only under the ignored results/ tree. The accompanying
report distinguishes source-file checksums from source and family approval.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/manifests/open_archive_candidates_2026_11_06.json"

def retrieve(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "iri-benign-archive-review/0.1 (academic research; metadata only in CI)", "Accept": "text/plain"})
    with urlopen(request, timeout=35) as response:
        data = response.read(1_500_001)
    if len(data) > 1_500_000:
        raise ValueError("Response exceeds bounded size")
    if data.lstrip().lower().startswith((b"<!doctype", b"<html", b"<error")):
        raise ValueError("Provider returned an HTML/error response, not annotated sequence")
    return data

def acquire(accession: str):
    ncbi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(
        {"db": "nuccore", "id": accession, "rettype": "gbwithparts", "retmode": "text"})
    ena = f"https://www.ebi.ac.uk/ena/browser/api/embl/{accession}"
    failures = []
    for provider, url, fmt in [("NCBI GenBank", ncbi, "genbank"), ("ENA EMBL", ena, "embl")]:
        try:
            data = retrieve(url)
            text = data.decode("utf-8")
            if (fmt == "genbank" and not text.lstrip().startswith("LOCUS")) or (
                fmt == "embl" and not text.lstrip().startswith("ID")):
                raise ValueError("Not an expected flat-file record")
            record = SeqIO.read(StringIO(text), fmt)
            if record.id != accession:
                raise ValueError(f"Accession mismatch: {record.id!r} vs {accession!r}")
            if not (100 <= len(record.seq) <= 100_000):
                raise ValueError(f"Implausible length: {len(record.seq)}")
            if len(record.features) < 1:
                raise ValueError("Missing feature table")
            return data, record, provider, url, fmt, failures
        except (HTTPError, URLError, TimeoutError, UnicodeError, ValueError, OSError) as exc:
            failures.append({"provider": provider, "reason": str(exc)[:240]})
    return None, None, None, None, None, failures

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="results/open_archive_review_2026_11_06")
    ap.add_argument("--limit", type=int, default=0, help="For a bounded dry-run; 0 means all candidates")
    ap.add_argument("--pause", type=float, default=0.5, help="Seconds between distinct accession attempts")
    args = ap.parse_args()
    output = (ROOT / args.output).resolve()
    if not output.is_relative_to(ROOT / "results"):
        ap.error("Output must stay within the ignored repository results directory")
    output.mkdir(parents=True, exist_ok=True)
    originals = output / "raw_sources"
    originals.mkdir(exist_ok=True)
    inventory = json.loads(MANIFEST.read_text(encoding="utf-8"))
    candidates = inventory["candidates"][:args.limit or None]
    results = []
    for item in candidates:
        acc = item["accession_version"]
        if not re.fullmatch(r"[A-Z]{1,4}[0-9]{5,8}\.[0-9]+", acc):
            raise ValueError("Unexpected accession format")
        data, record, provider, url, fmt, failures = acquire(acc)
        if data is None:
            result = {"accession_version": acc, "status": "HOLD_access_failure", "attempts": failures}
        else:
            suffix = ".gb" if fmt == "genbank" else ".embl"
            raw_path = originals / (acc + suffix)
            raw_path.write_bytes(data)
            features = Counter(f.type for f in record.features)
            result = {
                "accession_version": acc, "status": "ACQUIRED_NOT_ADMITTED",
                "provider": provider, "retrieval_url": url, "format": fmt,
                "original_bytes": len(data), "original_sha256": sha256(data).hexdigest(),
                "sequence_length_nt": len(record.seq), "feature_count": len(record.features),
                "cds_count": features.get("CDS", 0), "feature_types": dict(sorted(features.items())),
                "title": record.description, "topology": record.annotations.get("topology"),
                "original_file_local_path": str(raw_path.relative_to(ROOT)),
                "license_and_exact_annotation_review": "PENDING",
                "homology_families": [], "reference_independence_groups": [],
                "model_scoring_approved": False, "attempts_before_success": failures
            }
        results.append(result)
        print(acc, result["status"], "CDS", result.get("cds_count", "unknown"), flush=True)
        time.sleep(max(0, args.pause))
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_level": "exact-archive-download-integrity_only; NOT identity-equivalence, annotation truth or scientific admission",
        "original_ten_approved_candidates_unchanged": True,
        "candidate_count": len(candidates),
        "acquired_count": sum(x["status"] == "ACQUIRED_NOT_ADMITTED" for x in results),
        "blocked_count": sum(x["status"] != "ACQUIRED_NOT_ADMITTED" for x in results),
        "records": results
    }
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Acquired {report['acquired_count']} / {len(candidates)} public candidate records; not admitted.")
    return 0 if report["acquired_count"] else 2

if __name__ == "__main__":
    raise SystemExit(main())
