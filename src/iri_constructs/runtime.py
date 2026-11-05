"""Sequence-free run records and resource measurements; optional imports are lazy."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class InputError(ValueError):
    """An actionable validation error whose message contains no sequence content."""


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    except (OSError, ValueError) as exc:
        raise InputError('Cannot read JSON configuration; check file existence and JSON syntax') from exc
    if not isinstance(value, dict):
        raise InputError('Configuration must be a JSON object')
    return value


def require_text(obj: dict, key: str) -> str:
    value = obj.get(key)
    if not isinstance(value, str) or not value.strip():
        raise InputError(f'Missing or invalid metadata field: {key}')
    return value


def local_output(path: Path) -> Path:
    """Require ignored repository output roots, independent of working directory."""
    path = Path(path).resolve()
    roots = [ROOT / 'results', ROOT / 'data' / 'processed']
    if not any(path.is_relative_to(root) and path != root for root in roots):
        raise InputError('Output must be a subdirectory of repository results/ or data/processed/')
    return path


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def run_metadata(kind: str, config: dict) -> dict:
    versions = {}
    for package in ('biopython', 'numpy', 'torch', 'transformers', 'esm', 'huggingface-hub', 'psutil'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    try:
        revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True,
                                           stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))
    except (OSError, subprocess.SubprocessError):
        revision, dirty = None, None
    return {
        'result_schema': 't05-results-1', 'experiment': kind,
        'started_at_utc': datetime.now(timezone.utc).isoformat(),
        'python': platform.python_version(), 'platform': platform.platform(),
        'cpu': platform.processor(), 'logical_cpu_count': os.cpu_count(),
        'packages': versions, 'git_revision': revision, 'working_tree_dirty': dirty,
        'source_code_sha256': {p.relative_to(ROOT).as_posix(): digest_file(p)
                               for folder in ('src', 'scripts', 'experiments')
                               for p in sorted((ROOT / folder).rglob('*.py'))},
        'installed_packages': {d.metadata['Name']: d.version for d in importlib.metadata.distributions()
                               if d.metadata.get('Name')},
        'configuration': config, 'status': 'started',
        'limitations': ['Research representations and annotation agreement only; no function assay or safety decision.',
                        'No training or fine-tuning. No cross-model embedding-space comparisons.',
                        'No approved-input or pretrained-model performance is implied by unit tests.'],
    }


def failure(exc: Exception) -> dict:
    # Dependency/network/runtime exceptions can contain sequence strings or credentials.
    # Persist only controlled validation text or an exception type, never arbitrary repr/tracebacks.
    return {'error_type': type(exc).__name__,
            'reason': str(exc) if isinstance(exc, InputError) else
            'Execution failed; check installed dependencies, model access/cache, device capacity and token support.'}


class Measure:
    """Wall time, sampled process RSS and CUDA allocator peaks (not total device RAM)."""
    def __init__(self, device: str = 'cpu'):
        self.device = device
        self.data = {}
        self._stop = threading.Event()

    def __enter__(self):
        self._started = time.perf_counter()
        self._peak = None
        self._thread = None
        try:
            import psutil
            process = psutil.Process()
            self._peak = process.memory_info().rss
            def sample():
                while not self._stop.wait(0.02):
                    try:
                        self._peak = max(self._peak, process.memory_info().rss)
                    except psutil.Error:
                        return
            self._process = process
            self._thread = threading.Thread(target=sample, daemon=True)
            self._thread.start()
        except ImportError:
            pass
        if self.device.startswith('cuda'):
            import torch
            torch.cuda.synchronize(self.device)
            torch.cuda.reset_peak_memory_stats(self.device)
        return self

    def __exit__(self, *_):
        if self.device.startswith('cuda'):
            import torch
            torch.cuda.synchronize(self.device)
            self.data['peak_cuda_allocated_bytes'] = torch.cuda.max_memory_allocated(self.device)
            self.data['peak_cuda_reserved_bytes'] = torch.cuda.max_memory_reserved(self.device)
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
            self._peak = max(self._peak, self._process.memory_info().rss)
        self.data.update(wall_seconds=time.perf_counter() - self._started,
                         sampled_peak_process_rss_bytes=self._peak,
                         rss_sampling_interval_seconds=0.02,
                         memory_scope='Absolute process RSS; CUDA allocator only. Sampling may miss brief peaks.')
