"""Bind an external scientific review to exact files; never issue that review."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .runtime import InputError, digest_file


def _hash(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode('utf-8')).hexdigest()


def _artifact(root: Path, path: str | Path) -> tuple[Path, str]:
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise InputError('Review artifacts require nonempty local file paths')
    root = Path(root).resolve()
    supplied = Path(path)
    target = (supplied if supplied.is_absolute() else root / supplied).resolve()
    if not target.is_relative_to(root) or target == root:
        raise InputError('Review artifacts must remain inside the repository')
    if not target.is_file():
        raise InputError('A required review artifact is absent; reacquire its exact reviewed bytes')
    return target, target.relative_to(root).as_posix()


def build_review_freeze(root: Path, artifacts: list[str | Path]) -> dict:
    """Hash every explicit artifact, including ignored source bytes when requested.

    The returned object is a prospective evidence binding, not source admission,
    an authenticity check on a review message, or experimental execution authority.
    """
    if not artifacts:
        raise InputError('Cannot freeze an empty review packet')
    rows = {}
    for path in artifacts:
        target, name = _artifact(root, path)
        if name in rows:
            raise InputError('Duplicate artifact identity in review packet')
        rows[name] = {'sha256': digest_file(target), 'bytes': target.stat().st_size}
    freeze = {
        'schema_version': 'iri-review-freeze-1',
        'scientific_clearance': 'NOT_GRANTED_BY_THIS_ARTIFACT',
        'artifacts': dict(sorted(rows.items())),
        'scope': 'Exact evidence/code/config/source bytes for separate ChatGPT delegated review. '
                 'Review authenticity and stage-specific scientific gates require the original review message.',
    }
    freeze['freeze_sha256'] = _hash(freeze)
    return freeze


def verify_review_freeze(root: Path, freeze: dict) -> dict:
    """Refuse changed or missing evidence without promoting a freeze to approval."""
    if not isinstance(freeze, dict) or freeze.get('schema_version') != 'iri-review-freeze-1':
        raise InputError('Expected an iri-review-freeze-1 evidence binding')
    if freeze.get('scientific_clearance') != 'NOT_GRANTED_BY_THIS_ARTIFACT':
        raise InputError('An evidence freeze cannot grant scientific clearance')
    value = {key: item for key, item in freeze.items() if key != 'freeze_sha256'}
    if freeze.get('freeze_sha256') != _hash(value):
        raise InputError('Review freeze content changed; refuse use of a stale review')
    artifacts = freeze.get('artifacts')
    if not isinstance(artifacts, dict) or not artifacts:
        raise InputError('Review freeze must bind nonempty file evidence')
    for name, expected in artifacts.items():
        target, canonical = _artifact(root, name)
        if canonical != name or not isinstance(expected, dict):
            raise InputError('Review artifact identity is not canonical')
        if (type(expected.get('bytes')) is not int or expected['bytes'] < 0
                or target.stat().st_size != expected['bytes']
                or digest_file(target) != expected.get('sha256')):
            raise InputError('Reviewed evidence bytes changed; obtain a review of the revised packet')
    return {'status': 'evidence_bytes_verified', 'freeze_sha256': freeze['freeze_sha256'],
            'artifact_count': len(artifacts), 'scientific_execution_authorized': False}
