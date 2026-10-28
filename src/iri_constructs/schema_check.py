"""Structural checks for fictional JSON schema examples; not biological validation."""
from __future__ import annotations
import json
from pathlib import Path

ARCHETYPES = {"single_cds", "multi_cds", "tagged_fusion", "backbone_heavy"}
FEATURE_TYPES = {"CDS", "promoter", "terminator", "tag", "linker", "selection_marker", "backbone", "other", "unknown"}
ANNOTATION_STATUSES = {"known", "unknown", "ambiguous", "pending_review"}
INPUT_STATUSES = {"pending_model_selection", "ready", "unsupported", "needs_review"}
MOLECULE_TYPES = {"protein", "dna", "metadata_only"}

def _valid_span(loc: dict, max_length: int, *, allow_unknown: bool) -> bool:
    if not isinstance(loc, dict):
        return False
    start, end, strand = loc.get("start0"), loc.get("end0"), loc.get("strand")
    if start is None or end is None:
        return allow_unknown and start is None and end is None and strand is None
    return (isinstance(start, int) and not isinstance(start, bool)
            and isinstance(end, int) and not isinstance(end, bool)
            and 0 <= start < end <= max_length and strand in (1, -1, None))

def check_schema_example(path: Path) -> list[str]:
    """Return errors in a schematic example, if any; not full JSON Schema validation."""
    try:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return [f"Invalid example JSON: {exc}"]
    errs: list[str] = []
    if doc.get("schema_version") != "0.1":
        errs.append("Unrecognized schema version")
    parent = doc.get("construct", {})
    cid = parent.get("construct_id")
    length = parent.get("sequence_length_nt")
    if not isinstance(cid, str) or not cid:
        return ["Missing construct_id"]
    if not isinstance(length, int) or isinstance(length, bool) or length < 1:
        return ["Invalid sequence_length_nt"]
    if parent.get("primary_archetype") not in ARCHETYPES:
        errs.append("Invalid primary_archetype")
    if any(item not in ARCHETYPES for item in parent.get("secondary_archetypes", [])):
        errs.append("Invalid secondary_archetypes")
    if parent.get("annotation_status") not in ANNOTATION_STATUSES:
        errs.append("Invalid parent annotation_status")
    features = doc.get("features", [])
    units = doc.get("units", [])
    if not isinstance(features, list) or not isinstance(units, list):
        return errs + ["features and units must be lists"]
    fids: set[str] = set()
    for f in features:
        fid = f.get("feature_id")
        if not isinstance(fid, str) or not fid or fid in fids:
            errs.append("Duplicate/missing feature_id")
        else:
            fids.add(fid)
        if f.get("construct_id") != cid:
            errs.append("Feature parent ID mismatch")
        if f.get("feature_label") not in FEATURE_TYPES:
            errs.append("Invalid feature_label")
        if f.get("annotation_status") not in ANNOTATION_STATUSES:
            errs.append("Invalid feature annotation_status")
        if not _valid_span(f.get("location"), length, allow_unknown=True):
            errs.append("Invalid feature location")
        if not isinstance(f.get("annotation_source"), str) or not f["annotation_source"]:
            errs.append("Missing annotation_source")
    uids: set[str] = set()
    for u in units:
        uid = u.get("unit_id")
        if not isinstance(uid, str) or not uid or uid in uids:
            errs.append("Duplicate/missing unit_id")
        else:
            uids.add(uid)
        if u.get("construct_id") != cid:
            errs.append("Unit parent ID mismatch")
        if not isinstance(u.get("feature_ids"), list) or not u["feature_ids"]:
            errs.append("Unit feature_ids must be a nonempty list")
        elif any(fid not in fids for fid in u["feature_ids"]):
            errs.append("Unit references a nonexistent feature")
        if u.get("molecule_type") not in MOLECULE_TYPES:
            errs.append("Invalid molecule_type")
        if u.get("model_input_status") not in INPUT_STATUSES:
            errs.append("Invalid model_input_status")
        segs = u.get("source_segments", [])
        if not isinstance(segs, list) or any(not _valid_span(seg, length, allow_unknown=False) for seg in segs):
            errs.append("Invalid unit source_segments")
    return errs
