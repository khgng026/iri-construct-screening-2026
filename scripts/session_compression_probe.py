#!/usr/bin/env python3
"""Opt-in macOS transparent-compression probe using ONLY newly created toy JSON.

No scientific input path is accepted. Run only in the lead's granted lightweight
slot. This is a feasibility probe, not a production artifact compressor.
"""
from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile

UF_COMPRESSED = 0x20  # Local SDK sys/stat.h.
F_FULLFSYNC = 51      # Local SDK sys/fcntl.h; documented by fsync(2).
PEAK_SCRATCH_BYTES = 16 * 1024 * 1024
STORAGE_XATTRS = {'com.apple.decmpfs', 'com.apple.ResourceFork'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def sync_file(path):
    with path.open('rb') as stream:
        os.fsync(stream.fileno())
        fcntl.fcntl(stream.fileno(), F_FULLFSYNC)


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def describe(path):
    info = path.lstat()
    require(stat.S_ISREG(info.st_mode) and not path.is_symlink(), 'Probe requires a regular owned file')
    names = subprocess.run(['/usr/bin/xattr', str(path)], text=True, capture_output=True,
                           timeout=5, check=True).stdout.splitlines()
    attributes = {}
    for name in sorted(set(names)-STORAGE_XATTRS):
        encoded = subprocess.run(['/usr/bin/xattr', '-px', name, str(path)], text=True,
                                 capture_output=True, timeout=5, check=True).stdout
        attributes[name] = hashlib.sha256(bytes.fromhex(encoded)).hexdigest()
    return {'path': str(path), 'sha256': sha(path), 'logical_bytes': info.st_size,
            'allocated_bytes': info.st_blocks * 512, 'inode': info.st_ino, 'device': info.st_dev,
            'mode': stat.S_IMODE(info.st_mode), 'uid': info.st_uid, 'gid': info.st_gid,
            'mtime_ns': info.st_mtime_ns, 'flags': info.st_flags,
            'compressed_flag': bool(info.st_flags & UF_COMPRESSED),
            'nonstorage_xattrs': attributes}


def verify_copy(before, after):
    for field in ('sha256', 'logical_bytes', 'mode', 'uid', 'gid', 'mtime_ns', 'nonstorage_xattrs'):
        require(before[field] == after[field], 'Toy copy differs in ' + field)
    require(before['flags'] & ~UF_COMPRESSED == after['flags'] & ~UF_COMPRESSED,
            'Toy copy changed unrelated file flags')


def probe(args):
    require(sys.platform == 'darwin', 'This probe requires native macOS')
    require(args.reserve_bytes >= 0, 'Disk reserve must be nonnegative and inherited from the existing guard')
    parent = args.scratch_parent.resolve(strict=True)
    require(parent.is_dir(), 'Scratch parent must be an existing owned directory')
    require(shutil.disk_usage(parent).free >= args.reserve_bytes + PEAK_SCRATCH_BYTES,
            'Existing disk reserve plus bounded scratch headroom is unavailable')
    root = Path(tempfile.mkdtemp(prefix='iri-toy-compression-', dir=parent))
    ledger = root / 'ledger.jsonl'

    def event(name, **fields):
        value = {'time_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'event': name, **fields}
        with ledger.open('ab') as stream:
            stream.write(json.dumps(value, sort_keys=True, allow_nan=False).encode() + b'\n')
            stream.flush(); os.fsync(stream.fileno()); fcntl.fcntl(stream.fileno(), F_FULLFSYNC)
        sync_directory(root)

    original, compressed, rollback, decompressed = [root/name for name in (
        'record.json', 'compressed-copy.json', 'rollback-original.json', 'decompressed-copy.json')]
    version = subprocess.run(['/usr/bin/sw_vers'], text=True, capture_output=True, timeout=5, check=True)
    tool_evidence = {'ditto_path': '/usr/bin/ditto', 'ditto_sha256': sha(Path('/usr/bin/ditto')),
                     'manual_path': '/usr/share/man/man1/ditto.1',
                     'manual_sha256': sha(Path('/usr/share/man/man1/ditto.1')),
                     'macos_version': version.stdout.strip(), 'kernel_release': os.uname().release,
                     'python_version': sys.version, 'python_executable': sys.executable}
    # Fictional repetitive metadata, bounded below 4 MiB; no source/result data.
    payload = json.dumps({'fictional_metadata_only': True,
        'rows': [{'id': i, 'category': 'toy-metadata', 'description': 'repeated benign fictional metadata'}
                 for i in range(16384)]}, sort_keys=True, separators=(',', ':')).encode()+b'\n'
    require(len(payload) < 4*1024*1024, 'Toy file exceeds prespecified byte cap')
    with original.open('xb') as stream:
        stream.write(payload); stream.flush(); os.fsync(stream.fileno())
    os.chmod(original, 0o640)
    subprocess.run(['/usr/bin/xattr', '-w', 'com.openai.iri.compression_probe',
                    'fictional-toy-only', str(original)], timeout=5, check=True)
    os.utime(original, ns=(1700000000000000000, 1700000000000000000))
    sync_file(original); sync_directory(root)
    before = describe(original)
    require(not before['compressed_flag'], 'New toy input is unexpectedly already compressed')
    event('toy_created', before=before, reserve_bytes=args.reserve_bytes, scratch_limit_bytes=PEAK_SCRATCH_BYTES,
          tool_evidence=tool_evidence)
    env = {key: value for key, value in os.environ.items() if not key.startswith('DITTO')}

    def copy(options, source, target):
        require(not target.exists(), 'Probe copy target already exists')
        require(shutil.disk_usage(root).free >= args.reserve_bytes + 2*len(payload), 'Disk reserve/headroom changed')
        command = ['/usr/bin/ditto', '--noclone', '--rsrc', '--extattr', '--acl', *options,
                   str(source), str(target)]
        event('copy_intent', command=command)
        result = subprocess.run(command, env=env, text=True, capture_output=True, timeout=30)
        event('copy_exit', returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
        require(result.returncode == 0, 'Native ditto copy failed; inspect toy ledger')
        sync_file(target); sync_directory(root)

    try:
        copy(['--hfsCompression'], original, compressed)
        packed = describe(compressed); verify_copy(before, packed)
        require(describe(original)['sha256'] == before['sha256'], 'Toy source changed while copying')
        event('compressed_copy_verified', after=packed)
        if not packed['compressed_flag'] or packed['allocated_bytes'] >= before['allocated_bytes']:
            outcome = {'status': 'no_effective_transparent_compression', 'before': before, 'compressed': packed,
                       'replacement_performed': False, 'root': str(root), 'tool_evidence': tool_evidence}
        else:
            # Test atomic replacement only on this owned toy file. The hard link
            # retains the original inode until rollback is independently checked.
            event('replace_intent', original=str(original), temporary=str(compressed), rollback=str(rollback),
                  expected_sha256=before['sha256'], expected_bytes=before['logical_bytes'])
            os.link(original, rollback); sync_directory(root)
            os.replace(compressed, original); sync_directory(root)
            installed = describe(original); verify_copy(before, installed)
            command = [sys.executable, '-c',
                'import hashlib,json,pathlib,sys; p=pathlib.Path(sys.argv[1]); b=p.read_bytes(); '
                'assert len(b)==int(sys.argv[2]); assert hashlib.sha256(b).hexdigest()==sys.argv[3]; '
                'v=json.loads(b); assert v["fictional_metadata_only"] is True; assert len(v["rows"])==16384; '
                'print("fresh_process_logical_bytes_and_json_verified")', str(original), str(len(payload)), before['sha256']]
            result = subprocess.run(command, text=True, capture_output=True, timeout=15)
            require(result.returncode == 0, 'Fresh-process verification failed')
            event('replace_verified', after=installed, independent_process_output=result.stdout.strip())
            copy(['--nohfsCompression', '--nopreserveHFSCompression'], original, decompressed)
            unpacked = describe(decompressed); verify_copy(before, unpacked)
            require(not unpacked['compressed_flag'], 'Documented decompression copy left compression enabled')
            event('decompressed_copy_verified', after=unpacked)
            # Roll back the toy pathname to its original inode, not a new source.
            event('rollback_intent', original=str(original), rollback=str(rollback))
            os.replace(rollback, original); sync_directory(root)
            restored = describe(original); verify_copy(before, restored)
            require(restored['inode'] == before['inode'] and not restored['compressed_flag'],
                    'Toy rollback did not restore original inode/storage state')
            event('rollback_verified', after=restored)
            outcome = {'status': 'toy_round_trip_passed', 'before': before, 'compressed': packed,
                       'decompressed': unpacked, 'restored': restored, 'replacement_performed': True,
                       'scientific_artifacts_touched': False, 'root': str(root),
                       'tool_evidence': tool_evidence,
                       'allocated_byte_reduction_on_toy': before['allocated_bytes']-packed['allocated_bytes'],
                       'limitations': ['Toy support only; actual record effectiveness is untested.',
                           'APFS snapshots/shared extents can prevent equal immediate free-space reclamation.',
                           'No scientific artifact or completion descriptor was modified.']}
        (root/'result.json').write_text(json.dumps(outcome, indent=2, sort_keys=True)+'\n')
        sync_file(root/'result.json'); sync_directory(root)
        event('probe_complete', status=outcome['status'], result_sha256=sha(root/'result.json'))
        print(json.dumps({'status': outcome['status'], 'root': str(root), 'result_sha256': sha(root/'result.json')}), flush=True)
    except BaseException as exc:
        # Preserve all owned toy evidence. If the original backup exists, restore
        # it at the original toy pathname before recording the failure.
        if rollback.exists():
            os.replace(rollback, original); sync_directory(root)
        event('probe_failed', error=f'{type(exc).__name__}: {exc}', original=describe(original))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch-parent', type=Path, required=True)
    parser.add_argument('--reserve-bytes', type=int, required=True,
                        help='Unchanged minimum free bytes from the existing resource guard')
    probe(parser.parse_args())


if __name__ == '__main__':
    main()
