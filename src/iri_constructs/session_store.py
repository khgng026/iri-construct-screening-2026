"""Durable local experiment units and append-only events; no biological logic."""
from __future__ import annotations

import datetime as dt
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import uuid


def utcnow():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_name(value):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value) or value in (".", ".."):
        raise ValueError("Unsafe artifact identifier")
    return value


def fsync_dir(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_bytes(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    try:
        with tmp.open("xb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
        fsync_dir(path.parent)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_json(path, value):
    atomic_bytes(path, canonical(value) + b"\n")


def append_event(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {"time_utc": utcnow(), **value}
    with path.open("ab") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(canonical(record) + b"\n")
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)
    return record


class UnitStore:
    """Immutable attempts: only a verified descriptor makes a unit complete.

    A crash before directory publication leaves a hidden staging directory. A
    corrupt completed attempt is retained for inspection, never overwritten.
    """

    def __init__(self, root, identity, ledger=None):
        self.root = Path(root)
        self.identity = identity
        self.identity_hash = digest(canonical(identity))
        self.root.mkdir(parents=True, exist_ok=True)
        self.ledger = Path(ledger) if ledger else self.root / "events.jsonl"
        manifest = self.root / "identity.json"
        if manifest.exists():
            if json.loads(manifest.read_text()) != identity:
                raise ValueError("Store identity differs: use a new job directory")
        else:
            atomic_json(manifest, identity)

    def event(self, **value):
        return append_event(self.ledger, {"store": str(self.root),
                                        "identity_hash": self.identity_hash, **value})

    def verify(self, path):
        path = Path(path)
        desc = json.loads((path / "complete.json").read_text())
        if desc["identity_hash"] != self.identity_hash:
            raise ValueError("Artifact identity mismatch")
        if not desc.get("files"):
            raise ValueError("Empty completion descriptor")
        for name, expected in desc["files"].items():
            safe_name(name)
            target = path / name
            if target.is_symlink() or target.stat().st_size != expected["bytes"] or file_hash(target) != expected["sha256"]:
                raise ValueError("Artifact bytes failed validation: " + name)
        return desc

    def find(self, key):
        base = self.root / "units" / safe_name(key)
        if not base.exists():
            return None
        for path in sorted(base.iterdir(), reverse=True):
            if not path.is_dir() or path.name.startswith("."):
                continue
            try:
                desc = self.verify(path)
                if desc["key"] != key:
                    raise ValueError("Artifact key mismatch")
                return path
            except (OSError, ValueError, KeyError, TypeError) as exc:
                self.event(event="corrupt_unit_detected", key=key, path=str(path), error=str(exc))
        return None

    def publish(self, key, record, arrays=None):
        import numpy as np
        key = safe_name(key)
        existing = self.find(key)
        if existing:
            raise FileExistsError("Verified unit is immutable: " + key)
        base = self.root / "units" / key
        base.mkdir(parents=True, exist_ok=True)
        attempt = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%f") + "-" + uuid.uuid4().hex[:8]
        stage = base / (".staging-" + attempt)
        stage.mkdir()
        atomic_json(stage / "record.json", record)
        if arrays is not None:
            for name, arr in arrays.items():
                safe_name(name)
                arr = np.asarray(arr)
                if arr.dtype.hasobject or not np.isfinite(arr).all():
                    raise ValueError("Only finite numeric arrays are supported")
            stream = io.BytesIO()
            np.savez_compressed(stream, **arrays)
            atomic_bytes(stage / "arrays.npz", stream.getvalue())
        files = {p.name: {"sha256": file_hash(p), "bytes": p.stat().st_size}
                 for p in sorted(stage.iterdir()) if p.is_file()}
        desc = {"version": 1, "key": key, "attempt": attempt,
                "identity_hash": self.identity_hash, "completed_utc": utcnow(), "files": files}
        atomic_json(stage / "complete.json", desc)
        self.verify(stage)
        final = base / attempt
        os.rename(stage, final)
        fsync_dir(base)
        self.event(event="unit_complete", key=key, path=str(final), descriptor_sha256=file_hash(final / "complete.json"))
        return final

    def read(self, key, arrays=False):
        path = self.find(key)
        if path is None:
            return None
        record = json.loads((path / "record.json").read_text())
        if not arrays:
            return record
        import numpy as np
        with np.load(path / "arrays.npz", allow_pickle=False) as z:
            data = {name: z[name].copy() for name in z.files}
        return record, data

    def complete(self, expected_keys, metadata=None):
        keys = list(expected_keys)
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate expected keys")
        found = {key: self.find(key) for key in keys}
        missing = [key for key, path in found.items() if path is None]
        if missing:
            raise ValueError(f"Incomplete outputs: {len(missing)} missing units")
        index = {key: {"path": str(path.relative_to(self.root)),
                       "descriptor_sha256": file_hash(path / "complete.json")}
                 for key, path in found.items()}
        result = {"identity_hash": self.identity_hash, "completed_utc": utcnow(),
                  "expected_count": len(keys), "units": index, "metadata": metadata or {}}
        # Completion itself is a fresh descriptor; prior descriptors remain intact.
        name = "completion-" + uuid.uuid4().hex + ".json"
        atomic_json(self.root / name, result)
        self.event(event="job_complete", descriptor=str(self.root / name), units=len(keys))
        return result


def snapshot_tree(source, destination):
    """Consistent snapshot of a quiescent job; caller must own the write slot."""
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns(".staging-*", "*.tmp-*"))
    index = {str(p.relative_to(destination)): {"sha256": file_hash(p), "bytes": p.stat().st_size}
             for p in sorted(destination.rglob("*")) if p.is_file()}
    atomic_json(destination / "snapshot_index.json", {"source": str(source), "files": index,
                "same_disk_limitation": "Protects against some accidental edits; not disk failure."})
    return index
