"""Validate public, benign benchmark *metadata*. Does not inspect biological sequences."""
from __future__ import annotations
import argparse
import csv
from datetime import date
from pathlib import Path

REQUIRED_COLUMNS = ("record_id", "source", "source_url", "accession", "accession_version", "archetype", "sequence_length_nt", "feature_count", "annotation_status", "review_status", "licensing_status", "safety_review", "retrieved_on", "checksum_sha256", "notes")
ARCHETYPES = {"single_cds", "multi_cds", "tagged_fusion", "backbone_heavy"}
STATUS = {"pending", "approved", "rejected", "uncertain"}

def validate_manifest(path: str | Path) -> list[str]:
    """Return human-readable metadata errors; an empty list denotes valid structure."""
    errors: list[str] = []
    seen: set[str] = set()
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return ["Manifest is empty or missing a CSV header"]
        missing = set(REQUIRED_COLUMNS) - set(reader.fieldnames)
        if missing:
            return [f"Missing columns: {', '.join(sorted(missing))}"]
        for n, row in enumerate(reader, start=2):
            rid = (row["record_id"] or "").strip()
            if not rid:
                errors.append(f"Row {n}: record_id is blank")
            elif rid in seen:
                errors.append(f"Row {n}: repeated record_id {rid!r}")
            seen.add(rid)
            for field in ("source", "source_url", "accession", "accession_version"):
                if not (row[field] or "").strip():
                    errors.append(f"Row {n}: {field} is blank")
            if (row["archetype"] or "").strip() not in ARCHETYPES:
                errors.append(f"Row {n}: invalid archetype")
            for field in ("sequence_length_nt", "feature_count"):
                try:
                    number = int(row[field] or "")
                    if number < (1 if field == "sequence_length_nt" else 0):
                        raise ValueError
                except (ValueError, TypeError):
                    errors.append(f"Row {n}: invalid {field}")
            for field in ("annotation_status", "review_status", "licensing_status", "safety_review"):
                if (row[field] or "").strip() not in STATUS:
                    errors.append(f"Row {n}: invalid {field}")
            try:
                date.fromisoformat(row["retrieved_on"] or "")
            except (ValueError, TypeError):
                errors.append(f"Row {n}: retrieved_on must be YYYY-MM-DD")
            checksum = (row["checksum_sha256"] or "").strip()
            if checksum and (len(checksum) != 64 or any(x not in "0123456789abcdefABCDEF" for x in checksum)):
                errors.append(f"Row {n}: invalid checksum_sha256")
    return errors

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        errors = validate_manifest(args.manifest)
    except OSError as exc:
        parser.error(str(exc))
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("Manifest metadata validation passed (not biological validation).")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
