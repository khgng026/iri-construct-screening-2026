"""Structural checks for fictional JSON schema examples; not biological validation."""
from __future__ import annotations
import json
import re
from pathlib import Path, PurePosixPath

ARCHETYPES = {"single_cds", "multi_cds", "tagged_fusion", "backbone_heavy"}
FEATURE_TYPES = {"CDS", "promoter", "terminator", "tag", "linker", "selection_marker", "backbone", "other", "unknown"}
ANNOTATION_STATUSES = {"known", "unknown", "ambiguous", "pending_review"}
INPUT_STATUSES = {"pending_model_selection", "ready", "unsupported", "needs_review"}
MOLECULE_TYPES = {"protein", "dna", "metadata_only"}

def _valid_sequence_ref(ref: dict, molecule: str) -> bool:
    if not isinstance(ref, dict):
        return False
    path = PurePosixPath(ref.get("path", ""))
    return (not path.is_absolute() and ".." not in path.parts
            and path.parts[:2] in (("data", "raw"), ("data", "processed"))
            and ref.get("format") == "fasta"
            and isinstance(ref.get("record_id"), str) and bool(ref["record_id"])
            and bool(re.fullmatch(r"[0-9a-f]{64}", str(ref.get("sha256", ""))))
            and type(ref.get("length")) is int and ref["length"] > 0
            and ref.get("alphabet") == molecule)

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
    version = doc.get("schema_version")
    if version not in {"0.1", "0.2"}:
        errs.append("Unrecognized schema version")
    parent = doc.get("construct", {})
    cid = parent.get("construct_id")
    protein_parent = version == "0.2" and parent.get("molecule_type") == "protein"
    length = parent.get("sequence_length_aa" if protein_parent else "sequence_length_nt")
    coordinate_system = "protein_0_based_half_open" if protein_parent else "dna_0_based_half_open"
    if not isinstance(cid, str) or not cid:
        return ["Missing construct_id"]
    if not isinstance(length, int) or isinstance(length, bool) or length < 1:
        return ["Invalid parent sequence length"]
    if protein_parent:
        if parent.get("sequence_length_nt") is not None:
            errs.append("Protein-only parent must not invent nucleotide length")
        if parent.get("coordinate_system") != coordinate_system:
            errs.append("Invalid protein coordinate system")
        if not _valid_sequence_ref(parent.get("sequence_ref"), "protein"):
            errs.append("Invalid parent sequence_ref")
        elif parent["sequence_ref"]["length"] != length:
            errs.append("Parent reference length mismatch")
    allowed_archetypes = ARCHETYPES | ({"unknown"} if version == "0.2" else set())
    if parent.get("primary_archetype") not in allowed_archetypes:
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
        if protein_parent and (f.get("coordinate_system") != coordinate_system or f.get("location", {}).get("strand") is not None):
            errs.append("Protein feature coordinates must be explicit with null strand")
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
        if version == "0.2" and (u.get("model_input_status") == "ready" or u.get("sequence_ref") is not None):
            if not _valid_sequence_ref(u.get("sequence_ref"), u.get("molecule_type")):
                errs.append("Ready unit requires a valid sequence_ref")
        if protein_parent and u.get("coordinate_system") != coordinate_system:
            errs.append("Invalid unit coordinate system")
        segs = u.get("source_segments", [])
        if not isinstance(segs, list) or any(not _valid_span(seg, length, allow_unknown=False) for seg in segs):
            errs.append("Invalid unit source_segments")
    return errs
