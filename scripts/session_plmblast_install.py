#!/usr/bin/env python3
"""Acquire hash-pinned official wheels into a session-only dependency overlay."""
import argparse
from email.parser import BytesParser
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import uuid
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import atomic_bytes, atomic_json, file_hash, append_event
from session_experiment import check_deadline

LOCK_SHA = '0c6d566c5de1c2f247f05ee5776e583bd9d3a0d0dc68eeb841b7e236b346350f'


def require(value, message):
    if not value:
        raise ValueError(message)


def main(args):
    session = args.session.resolve()
    check_deadline(session)
    lock = session / 'inputs/plmblast_overlay_requirements_v1.txt'
    require(file_hash(lock) == LOCK_SHA, 'Dependency lock differs from reviewed recipe')
    require(shutil.disk_usage(session).free > 26 * 2 ** 30, 'Dependency installation reserve unavailable')
    versions = {}
    for line in lock.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        match = re.fullmatch(r'([a-z]+)==([\d.]+) --hash=sha256:([0-9a-f]{64})', line)
        require(match is not None, 'Unexpected locked requirement syntax')
        name, version, sha = match.groups()
        require(name not in versions, 'Duplicate package lock')
        versions[name] = {'version': version, 'sha256': sha}
    require(set(versions) == {'numpy', 'pandas', 'scipy', 'numba', 'llvmlite', 'pytz', 'tzdata'}, 'Overlay package roster differs')
    wheelhouse = session / 'vendor/plmblast-wheelhouse-v1'
    overlay = session / 'vendor/plmblast-overlay-v1'
    destination = session / 'inputs/plmblast_overlay_installation_v1.json'
    if destination.exists():
        record = json.loads(destination.read_text())
        require(record['status'] == 'verified_installed' and record['lock_sha256'] == LOCK_SHA
                and Path(record['overlay']) == overlay, 'Saved installation identity differs')
        entries = list(overlay.rglob('*'))
        require(not any(path.is_symlink() for path in entries), 'Installed overlay contains a symbolic link')
        require({str(path.relative_to(overlay)) for path in entries if path.is_file()} == set(record['installed_files']),
                'Installed overlay file membership changed')
        for relative, expected in record['installed_files'].items():
            require(file_hash(overlay / relative) == expected, 'Installed overlay changed')
        print(json.dumps({'status': 'verified_existing', 'overlay': str(overlay)}), flush=True)
        return
    require(not overlay.exists(), 'Unverified destination exists; preserve it for review')
    attempt = session / 'vendor/plmblast-install-attempts' / uuid.uuid4().hex
    attempt.mkdir(parents=True, exist_ok=False)
    wheelhouse.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    staging = attempt / 'overlay'
    # pip's isolated mode still reads global/site configuration. Disable those
    # only for these owned subprocesses; no persistent user settings change.
    pip_environment = dict(os.environ, PIP_CONFIG_FILE=os.devnull)
    try:
        official = {}
        for name, expected in versions.items():
            check_deadline(session)
            url = f'https://pypi.org/pypi/{name}/{expected["version"]}/json'
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'IRI-local-research/1.0'}), timeout=30) as response:
                raw = response.read(4 * 2 ** 20 + 1)
            require(len(raw) <= 4 * 2 ** 20, 'Official package metadata exceeded the bounded payload')
            atomic_bytes(attempt / (name + '-pypi.json'), raw)
            payload = json.loads(raw)
            matches = [item for item in payload['urls'] if item['digests']['sha256'] == expected['sha256']]
            require(len(matches) == 1, 'Pinned wheel is not uniquely present in official PyPI metadata')
            item = matches[0]
            require(item['packagetype'] == 'bdist_wheel' and item['url'].startswith('https://files.pythonhosted.org/'),
                    'Pinned artifact is not an official binary wheel')
            official[name] = {'version': expected['version'], 'filename': item['filename'], 'sha256': expected['sha256'],
                              'bytes': item['size'], 'url': item['url'], 'pypi_url': url,
                              'pypi_response_sha256': file_hash(attempt / (name + '-pypi.json')),
                              'declared_license': payload['info'].get('license'),
                              'license_classifiers': [value for value in payload['info'].get('classifiers', []) if value.startswith('License ::')]}
        atomic_json(attempt / 'official_wheels.json', official)
        previous = {p.name: file_hash(p) for p in wheelhouse.glob('*.whl')}
        command = [sys.executable, '-m', 'pip', '--isolated', '--disable-pip-version-check', 'download',
                   '--index-url', 'https://pypi.org/simple', '--only-binary=:all:', '--no-deps', '--require-hashes',
                   '--retries', '2', '--timeout', '30', '--dest', str(wheelhouse), '-r', str(lock)]
        with (attempt / 'download.log').open('xb') as log:
            subprocess.run(command, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                           env=pip_environment, check=True, timeout=600)
        require({p.name for p in wheelhouse.glob('*.whl')} == {value['filename'] for value in official.values()},
                'Wheelhouse has unexpected or missing artifacts')
        for name, item in official.items():
            path = wheelhouse / item['filename']
            require(path.stat().st_size == item['bytes'] and file_hash(path) == item['sha256'], 'Wheel byte verification failed')
            with zipfile.ZipFile(path) as archive:
                names = [entry for entry in archive.namelist() if entry.endswith('.dist-info/METADATA')]
                require(len(names) == 1, 'Wheel package metadata is ambiguous')
                metadata = BytesParser().parsebytes(archive.read(names[0]))
                require(metadata['Name'].lower() == name and metadata['Version'] == item['version'], 'Wheel metadata identity differs')
            item['new_local_payload_bytes'] = 0 if previous.get(path.name) == item['sha256'] else item['bytes']
        atomic_json(attempt / 'verified_wheels.json', official)
        check_deadline(session)
        command = [sys.executable, '-m', 'pip', '--isolated', '--disable-pip-version-check', 'install',
                   '--no-index', '--only-binary=:all:', '--no-deps', '--no-compile', '--require-hashes',
                   '--find-links', str(wheelhouse), '--target', str(staging), '-r', str(lock)]
        with (attempt / 'install.log').open('xb') as log:
            subprocess.run(command, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                           env=pip_environment, check=True, timeout=180)
        installed = {}
        for path in staging.glob('*.dist-info/METADATA'):
            metadata = BytesParser().parsebytes(path.read_bytes())
            name = metadata['Name'].lower()
            require(name in versions and metadata['Version'] == versions[name]['version'], 'Installed distribution differs')
            installed[name] = metadata['Version']
        require(set(installed) == set(versions), 'Missing installed distribution')
        entries = sorted(staging.rglob('*'))
        require(not any(path.is_symlink() for path in entries), 'Staged overlay contains a symbolic link')
        hashes = {str(path.relative_to(staging)): file_hash(path) for path in entries if path.is_file()}
        require(hashes, 'Empty overlay')
        staging.rename(overlay)
        record = {'status': 'verified_installed', 'overlay': str(overlay), 'lock_sha256': LOCK_SHA,
                  'script_sha256': file_hash(__file__), 'attempt': str(attempt), 'wheels': official,
                  'installed_versions': installed, 'installed_files': hashes,
                  'elapsed_seconds': time.monotonic() - started,
                  'scope': 'Session-only overlay; existing model environments unchanged. Import/ABI/JIT compatibility remains a separate runtime gate.'}
        atomic_json(destination, record)
        append_event(session / 'job_ledger.jsonl', {'event': 'isolated_dependency_installation', 'state': 'complete',
                     'record': str(destination), 'sha256': file_hash(destination), 'elapsed_seconds': record['elapsed_seconds']})
        print(json.dumps({'status': 'verified_installed', 'overlay': str(overlay), 'packages': installed,
                          'new_local_payload_bytes': sum(value['new_local_payload_bytes'] for value in official.values())}), flush=True)
    except BaseException as error:
        atomic_json(attempt / 'failure.json', {'status': 'failed', 'exception': f'{type(error).__name__}: {error}',
                    'elapsed_seconds': time.monotonic() - started, 'partial_artifacts_preserved': True})
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    main(parser.parse_args())
