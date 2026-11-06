#!/usr/bin/env python3
"""Portable entry point for core input verification and fresh CPU experiments."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
PACKET = ROOT / 'config/reproduction/packet'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def install_bundle(path):
    expected = json.loads((ROOT / 'config/reproduction/bundle.json').read_text())
    if sha(path.read_bytes()) != expected['sha256']:
        raise ValueError('Input archive checksum mismatch; no files installed')
    manifest = json.loads((ROOT / 'config/reproduction/data_bundle_manifest.json').read_text())
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)) or set(names) != set(manifest) | {'DATA_NOTICE.txt'}:
            raise ValueError('Unexpected archive membership')
        payloads = {}
        for name, digest in manifest.items():
            target = (ROOT / name).resolve()
            if not target.is_relative_to(ROOT.resolve() / 'data'):
                raise ValueError('Unsafe archive path')
            data = archive.read(name)
            if sha(data) != digest:
                raise ValueError('Input member checksum mismatch: ' + name)
            if target.exists() and target.read_bytes() != data:
                raise ValueError('Refusing to replace different existing input: ' + name)
            payloads[target] = data
        # Validate every member and existing destination before any writes.
        for target, data in payloads.items():
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
    return verify_inputs()


def verify_inputs():
    from iri_constructs.session_benchmark import load_suite
    suite = load_suite(ROOT, PACKET)
    print(json.dumps(suite['verification'], indent=2))
    return suite


def initialize(session, hours):
    if not 0 < hours <= 12:
        raise ValueError('Session duration must be between zero and 12 hours')
    verify_inputs()
    session.mkdir(parents=True, exist_ok=False)
    shutil.copytree(PACKET, session / 'inputs')
    end = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=hours)
    (session / 'launch.json').write_text(json.dumps({
        'hard_stop_utc': end.isoformat(), 'resource_guard': str(ROOT / 'scripts/resource_guard.py'),
        'note': 'Fresh reproduction session; historical launch/approvals are not reused.'}, indent=2))
    print('Created', session)


def model_config(name):
    directory = ROOT / 'config/reproduction/models'
    path = directory / (name + '.json')
    if path.parent != directory or not path.is_file():
        raise ValueError('Choose a configuration filename stem from config/reproduction/models')
    return path


def download(name):
    from huggingface_hub import hf_hub_download
    config = json.loads(model_config(name).read_text())
    if 'model_id' not in config:
        raise ValueError('Alignment does not require a checkpoint')
    files = ['model.safetensors', 'config.json']
    if config['family'] == 'esm2':
        files += ['vocab.txt', 'tokenizer_config.json', 'special_tokens_map.json']
    for filename in files:
        path = Path(hf_hub_download(config['model_id'], filename, revision=config['revision'], token=False))
        expected = config.get('expected_file_sha256', {}).get(filename)
        if expected and sha(path.read_bytes()) != expected:
            raise ValueError('Checkpoint checksum mismatch')
    print('Pinned checkpoint files available:', name)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest='command', required=True)
    install = commands.add_parser('install-inputs'); install.add_argument('archive', type=Path)
    commands.add_parser('verify-inputs')
    init = commands.add_parser('init'); init.add_argument('session', type=Path); init.add_argument('--hours', type=float, default=2)
    acquire = commands.add_parser('download-model'); acquire.add_argument('model')
    run = commands.add_parser('run'); run.add_argument('session', type=Path); run.add_argument('--model', default='alignment')
    run.add_argument('--limit', type=int, default=4); run.add_argument('--timeout', type=int, default=1800)
    args = p.parse_args()
    if args.command == 'install-inputs': install_bundle(args.archive)
    elif args.command == 'verify-inputs': verify_inputs()
    elif args.command == 'init': initialize(args.session.resolve(), args.hours)
    elif args.command == 'download-model': download(args.model)
    else:
        if args.limit < 0 or not 0 < args.timeout <= 43200:
            p.error('limit must be nonnegative (0 means full suite); timeout must be 1..43200 seconds')
        config = model_config(args.model)
        cmd = [sys.executable, str(ROOT / 'scripts/session_experiment.py'), '--session', str(args.session.resolve()),
               '--config', str(config), '--name', args.model, '--mode', 'alignment' if args.model == 'alignment' else 'model',
               '--limit', str(args.limit)]
        subprocess.run([sys.executable, str(ROOT / 'scripts/session_guard_entry.py'), '--timeout', str(args.timeout), '--', *cmd], check=True, cwd=ROOT)


if __name__ == '__main__':
    main()
