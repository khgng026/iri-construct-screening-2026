#!/usr/bin/env python3
"""Journalled native compression of an explicitly frozen completed-store allowlist.

Stages: freeze, apply --scope trial, verify, finalize; then separately authorized
apply --scope full, verify, finalize. Scientific logical bytes never change.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import datetime as dt
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import stat
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from session_independent_audit import verify_unit

INITIAL_STORES = ('expanded1000-8m-rank-v1-f61a400723f8',
                  'expanded1000-150m-rank-v1-81ed0141a739',
                  'expanded1000-650m-rank-v1-ea5b8e31cca5')
RESERVE = 25*2**30
SCRATCH = 16*2**20
MEMORY = 512*2**20
MAX_FILE = 64*2**20
UF_COMPRESSED = 0x20
STORAGE_XATTRS = {'com.apple.decmpfs', 'com.apple.ResourceFork'}
STABLE_METADATA = ('mode', 'uid', 'gid', 'mtime_ns', 'xattrs', 'acl')


def require(value, message):
    if not value: raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def read(path): return json.loads(Path(path).read_text())
def digest(value): return hashlib.sha256(value).hexdigest()
def now(): return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def sync_file(path):
    with Path(path).open('rb') as stream:
        os.fsync(stream.fileno())
        if sys.platform == 'darwin': fcntl.fcntl(stream.fileno(), 51)


def sync_dir(path):
    fd = os.open(path, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


def atomic(path, value):
    path = Path(path); temporary = path.with_name('.'+path.name+'.'+uuid.uuid4().hex)
    with temporary.open('xb') as stream:
        stream.write(canonical(value)+b'\n'); stream.flush(); os.fsync(stream.fileno())
        if sys.platform == 'darwin': fcntl.fcntl(stream.fileno(), 51)
    os.replace(temporary, path); sync_dir(path.parent)


def immutable(path, value):
    require(not Path(path).exists(), 'Refuse to overwrite immutable operational evidence')
    atomic(path, value)


def event(state, name, **fields):
    value = {'time_utc': now(), 'event': name, **fields}
    with (state/'ledger.jsonl').open('ab') as stream:
        stream.write(canonical(value)+b'\n'); stream.flush(); os.fsync(stream.fileno())
        if sys.platform == 'darwin': fcntl.fcntl(stream.fileno(), 51)
    sync_dir(state)


def no_symlinks(path, base):
    path, base = Path(path).absolute(), Path(base).absolute()
    require(path.is_relative_to(base), 'Path escapes approved hierarchy')
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        require(not current.is_symlink(), 'Symlink component is not eligible')
    return path


def command_output(command):
    return subprocess.run(command, text=True, capture_output=True, timeout=10, check=True).stdout


def metadata(path):
    info = Path(path).lstat()
    require(stat.S_ISREG(info.st_mode) and not Path(path).is_symlink(), 'Only regular files qualify')
    names = command_output(['/usr/bin/xattr', str(path)]).splitlines()
    xattrs = {}
    for name in sorted(set(names)-STORAGE_XATTRS):
        data = bytes.fromhex(command_output(['/usr/bin/xattr', '-px', name, str(path)]))
        require(len(data) <= 2**20, 'Unexpectedly large non-storage extended attribute')
        xattrs[name] = digest(data)
    lines = command_output(['/bin/ls', '-lde', str(path)]).splitlines()
    acl = [line.strip() for line in lines[1:] if line.lstrip()[:1].isdigit()]
    return {'inode': info.st_ino, 'device': info.st_dev, 'nlink': info.st_nlink,
        'mode': stat.S_IMODE(info.st_mode), 'uid': info.st_uid, 'gid': info.st_gid,
        'size': info.st_size, 'allocated_bytes': info.st_blocks*512, 'flags': info.st_flags,
        'mtime_ns': info.st_mtime_ns, 'atime_ns': info.st_atime_ns, 'ctime_ns': info.st_ctime_ns,
        'birthtime': info.st_birthtime, 'xattrs': xattrs, 'acl': acl,
        'storage_xattrs': sorted(set(names)&STORAGE_XATTRS)}


def stable(before, after):
    require(all(before[k] == after[k] for k in STABLE_METADATA), 'Logical file metadata changed')
    require(before['flags'] & ~UF_COMPRESSED == after['flags'] & ~UF_COMPRESSED,
            'Unrelated file flags changed')


def exact(path, entry):
    require(Path(path).stat().st_size == entry['bytes'] and sha(path) == entry['sha256'],
            'Logical file size/hash changed')


def safe_record(path, entry, *, allow_link=False):
    value = metadata(path)
    require(value['uid'] == os.getuid() and value['flags'] in (0, UF_COMPRESSED), 'Owner/flags not eligible')
    require(allow_link or value['nlink'] == 1, 'Unexpected hard links; refuse replacement')
    require(0 < value['size'] <= MAX_FILE, 'File is outside bounded compression size')
    if not value['flags'] & UF_COMPRESSED:
        require(not value['storage_xattrs'], 'Uncompressed source has unexpected storage/resource-fork attributes')
    stable(entry['metadata'], value); exact(path, entry)
    return value


class Guard:
    """Import the prepared checks unchanged; add stronger local finite limits."""
    def __init__(self, session, state, seconds, expected=None):
        require(0 < seconds <= 1800, 'Overall timeout must be finite and at most 1800 seconds')
        self.started = time.monotonic(); self.seconds = seconds; self.state = state
        path = Path(read(session/'launch.json')['resource_guard']).resolve()
        if expected: require(str(path) == expected['path'] and sha(path) == expected['sha256'], 'Prepared guard identity changed')
        self.identity = {'path': str(path), 'sha256': sha(path)}
        spec = importlib.util.spec_from_file_location('unchanged_prepared_guard', path)
        self.module = importlib.util.module_from_spec(spec); spec.loader.exec_module(self.module)
        self.baseline = self.module.snapshot(); self.last = 0
        self.check(force=True)

    def check(self, *, force=False, headroom=SCRATCH, child=None):
        require(time.monotonic()-self.started < self.seconds, 'Overall compression/verification time limit reached')
        own_peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        child_rss = 0
        if child is not None and child.poll() is None:
            result = subprocess.run(['/bin/ps', '-o', 'rss=', '-p', str(child.pid)], text=True,
                                    capture_output=True, timeout=5)
            require(result.returncode == 0 or child.poll() is not None, 'Could not inspect live ditto memory')
            child_rss = int(result.stdout.strip() or 0)*1024
        require(own_peak+child_rss < MEMORY, '512 MiB aggregate worker/ditto memory limit reached')
        if not force and time.monotonic()-self.last < 2: return
        sample = self.module.snapshot()
        reasons = self.module.reasons(sample, self.baseline)
        free = shutil.disk_usage(self.state).free
        if free < RESERVE+headroom: reasons.append('25 GiB reserve plus scratch headroom unavailable')
        event(self.state, 'resource_check', sample=sample, blocked_by=reasons,
              worker_peak_bytes=own_peak, child_rss_bytes=child_rss, available_bytes=free, required_headroom=headroom)
        self.last = time.monotonic()
        require(not reasons, 'Resource guard blocked: '+', '.join(reasons))


@contextmanager
def state_lock(state):
    # A transaction-journal lock only; never acquire the global compute.lock.
    with (state/'journal.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def quiet_jobs(plan):
    for job in plan['jobs']:
        result = subprocess.run(['/usr/sbin/lsof', '-nP', '-t', '+D', job['path']],
                                text=True, capture_output=True, timeout=20)
        require(result.returncode == 1 and not result.stdout.strip() and not result.stderr.strip(),
                'Store has open files, or quiescence check failed: '+job['name'])


def expected_units(identity):
    require(identity['mode'] == 'expanded_model_retrieval_v1', 'Only completed expanded ranking stores qualify')
    ids = identity['query_ids']
    require(ids and len(ids) == len(set(ids)), 'Invalid completed query coverage')
    return {'query-'+qid for qid in ids} | {'summary'}


def freeze(args):
    session = args.session.resolve(strict=True); state = args.state.absolute()
    require(state.is_relative_to(session/'audit') and not state.exists(), 'Use a new owned audit state directory')
    no_symlinks(state, session); state.mkdir(); (state/'transactions').mkdir(); sync_dir(state.parent)
    with state_lock(state):
        guard = Guard(session, state, args.timeout)
        approvals = read(args.approved_stores) if args.approved_stores else None
        names = [x['job'] for x in approvals['stores']] if approvals else list(INITIAL_STORES)
        require(1 <= len(names) <= 5 and len(names) == len(set(names)), 'Allowlist must contain 1–5 explicit unique completed stores')
        if approvals:
            require(approvals['authorization_note'] and all(Path(n).name == n and n not in ('.','..') for n in names),
                    'Future explicit allowlist requires an authorization note and safe names')
        jobs, lists = [], []
        for index, name in enumerate(names):
            path = no_symlinks(session/'jobs'/name, session/'jobs')
            identity = read(path/'identity.json'); ih = digest(canonical(identity))
            completions = sorted(path.glob('completion-*.json'))
            require(completions, 'Allowlisted store is not completed')
            selected = path/approvals['stores'][index]['completion_file'] if approvals else max(completions, key=lambda p:(read(p)['completed_utc'], p.name))
            require(selected in completions, 'Approved completion is absent')
            completion = read(selected); expected = expected_units(identity)
            require(completion['identity_hash'] == ih and completion['expected_count'] == len(expected)
                    and set(completion['units']) == expected, 'Completion unit coverage differs')
            if approvals:
                approved = approvals['stores'][index]
                require(sha(path/'identity.json') == approved['identity_file_sha256']
                        and sha(selected) == approved['completion_sha256'], 'Future allowlist pins do not match')
            job = {'name': name, 'path': str(path), 'identity_hash': ih,
                'identity_file_sha256': sha(path/'identity.json'), 'completion_file': selected.name,
                'completion_sha256': sha(selected), 'completion_files': {p.name:sha(p) for p in completions},
                'units': completion['units']}
            jobs.append(job); entries = []
            for key, reference in sorted(completion['units'].items()):
                guard.check()
                unit = verify_unit(path, key, reference, ih)
                target = no_symlinks(unit/'record.json', path)
                descriptor = read(unit/'complete.json'); value = descriptor['files']['record.json']
                entry = {'id': digest(str(target.relative_to(session)).encode())[:24], 'job_index': index,
                    'key': key, 'path': str(target), 'bytes': value['bytes'], 'sha256': value['sha256'],
                    'unit_descriptor_sha256': reference['descriptor_sha256'], 'metadata': metadata(target)}
                safe_record(target, entry); entries.append(entry)
            lists.append(entries)
        files = [values[i] for i in range(max(map(len, lists))) for values in lists if i < len(values)]
        plan = {'version': 'native-storage-compression-v1', 'created_utc': now(), 'session': str(session),
            'state': str(state), 'jobs': jobs, 'files': files, 'guard': guard.identity,
            'script_sha256': sha(Path(__file__)), 'independent_verifier_sha256': sha(ROOT/'scripts/session_independent_audit.py'),
            'ditto_sha256': sha(Path('/usr/bin/ditto')), 'macos': command_output(['/usr/bin/sw_vers']).strip(),
            'allowlist_source': {'path': str(args.approved_stores.resolve()), 'sha256': sha(args.approved_stores)} if approvals else 'fixed_original_three',
            'reserve_bytes': RESERVE, 'scratch_bytes': SCRATCH, 'memory_limit_bytes': MEMORY,
            'max_file_bytes': MAX_FILE, 'trial_file_ids': [x['id'] for x in files[:10]],
            'quiescence_authority': 'Explicit lead allowlist; runtime lsof and unchanged completion checks also required.'}
        quiet_jobs(plan); immutable(state/'plan.json', plan)
        event(state, 'plan_frozen', plan_sha256=sha(state/'plan.json'), files=len(files),
              logical_bytes=sum(x['bytes'] for x in files), allocated_bytes=sum(x['metadata']['allocated_bytes'] for x in files))
        print(json.dumps({'status':'frozen_read_only','state':str(state),'plan_sha256':sha(state/'plan.json'),
                          'files':len(files),'logical_bytes':sum(x['bytes'] for x in files)}), flush=True)


def load_plan(state, expected_sha256):
    state = state.resolve(strict=True); plan = read(state/'plan.json')
    require(expected_sha256 and sha(state/'plan.json') == expected_sha256,
            'Root-reviewed plan SHA256 is required and must match')
    require(plan['version'] == 'native-storage-compression-v1' and Path(plan['state']) == state, 'Wrong operational state')
    require(plan['script_sha256'] == sha(Path(__file__)) and plan['ditto_sha256'] == sha(Path('/usr/bin/ditto')),
            'Compression implementation/tool changed; do not reuse frozen plan')
    require(plan['independent_verifier_sha256'] == sha(ROOT/'scripts/session_independent_audit.py'), 'Independent verifier changed')
    if plan['allowlist_source'] == 'fixed_original_three':
        require([j['name'] for j in plan['jobs']] == list(INITIAL_STORES), 'Frozen initial allowlist expanded')
    else:
        approval = plan['allowlist_source']
        require(sha(Path(approval['path'])) == approval['sha256'], 'Explicit future allowlist approval changed')
        require([x['job'] for x in read(approval['path'])['stores']] == [j['name'] for j in plan['jobs']],
                'Explicit allowlist membership differs')
    require(len({e['id'] for e in plan['files']}) == len(plan['files']), 'Duplicate frozen file IDs')
    session=Path(plan['session']); seen=set()
    for index,job in enumerate(plan['jobs']):
        require(Path(job['path']) == session/'jobs'/job['name'], 'Job root is outside the explicit allowlist')
    for entry in plan['files']:
        job=plan['jobs'][entry['job_index']]; reference=job['units'][entry['key']]
        expected=Path(job['path'])/reference['path']/'record.json'
        require(Path(entry['path']) == expected and entry['id']==digest(str(expected.relative_to(session)).encode())[:24],
                'Only descriptor-bound record.json paths may be compressed')
        descriptor=read(expected.parent/'complete.json')
        require(descriptor['files']['record.json']=={'bytes':entry['bytes'],'sha256':entry['sha256']}
                and sha(expected.parent/'complete.json')==entry['unit_descriptor_sha256']==reference['descriptor_sha256'],
                'Frozen file is not bound to its completed descriptor')
        seen.add((entry['job_index'],entry['key']))
    require(seen=={(i,key) for i,j in enumerate(plan['jobs']) for key in j['units']}, 'Frozen plan omitted or added units')
    require(plan['trial_file_ids']==[e['id'] for e in plan['files'][:10]], 'First10 trial selection changed')
    return plan


def validate_chain(plan, guard, *, all_bytes=False):
    count = total = 0
    for job in plan['jobs']:
        root = no_symlinks(Path(job['path']), Path(plan['session'])/'jobs')
        require(sha(root/'identity.json') == job['identity_file_sha256'], 'Frozen identity bytes changed')
        require({p.name:sha(p) for p in root.glob('completion-*.json')} == job['completion_files'], 'Completion set/bytes changed')
        completion = read(root/job['completion_file'])
        require(sha(root/job['completion_file']) == job['completion_sha256'] and completion['units'] == job['units'], 'Pinned completion changed')
        require(expected_units(read(root/'identity.json')) == set(job['units']), 'Expected coverage changed')
        for key, reference in job['units'].items():
            guard.check()
            unit = no_symlinks(root/reference['path'], root)
            require(sha(unit/'complete.json') == reference['descriptor_sha256'], 'Descriptor bytes changed')
            if all_bytes:
                verified = verify_unit(root, key, reference, job['identity_hash'])
                descriptor = read(verified/'complete.json')
                total += sum(v['bytes'] for v in descriptor['files'].values())
            count += 1
    return {'verified_units':count, 'verified_logical_bytes':total}


def transactions(state, plan):
    entries = {x['id']:x for x in plan['files']}; result = []
    for directory in sorted((state/'transactions').iterdir()):
        require(directory.is_dir() and not directory.is_symlink(), 'Unexpected transaction path')
        if not (directory/'intent.json').exists():
            contents=list(directory.iterdir())
            require(any(directory.name.startswith(file_id+'-') for file_id in entries)
                    and all(p.is_file() and not p.is_symlink() and p.name.startswith('.intent.json.') for p in contents),
                    'Unpublished transaction contains unexpected data; preserve for review')
            # A crash while publishing the first intent cannot have reached the
            # copier/link/rename. Retain partial intent bytes without adopting them.
            continue
        intent = read(directory/'intent.json'); entry = entries.get(intent['file_id'])
        require(entry is not None and intent['plan_sha256'] == sha(state/'plan.json')
                and intent['source'] == entry['path'], 'Transaction is outside frozen plan')
        require(intent['backup'] == str(directory/'original.link')
                and intent['temporary'] == str(directory/'compressed.tmp'), 'Transaction-owned paths changed')
        result.append((directory, intent, entry))
    return result


def terminal(directory):
    names = [n for n in ('verified.json','rolled_back.json','skipped.json') if (directory/n).exists()]
    require(len(names) <= 1, 'Conflicting transaction outcomes')
    return names[0] if names else None


def state_digest(state, plan):
    rows = []
    for directory, intent, entry in transactions(state, plan):
        ending = terminal(directory)
        require(ending is not None, 'Unresolved intent; recover before independent verification')
        rows.append({'transaction':directory.name, 'intent':sha(directory/'intent.json'), 'outcome':ending,
                     'outcome_sha256':sha(directory/ending)})
    return digest(canonical(rows))


def rollback(state, directory, intent, entry, *, reason):
    source, backup = Path(intent['source']), Path(intent['backup'])
    if backup.exists():
        require(backup.lstat().st_ino == intent['before']['inode'], 'Rollback link no longer refers to original inode')
        exact(backup, entry); stable(intent['before'], metadata(backup))
        immutable(directory/'rollback_intent.json', {'time_utc':now(),'reason':reason}) if not (directory/'rollback_intent.json').exists() else None
        # Retain the installed copy as owned failure evidence before restoring.
        if source.exists() and source.lstat().st_ino != backup.lstat().st_ino and not (directory/'replaced-copy.link').exists():
            event(state,'rollback_copy_link_intent',source=str(source),destination=str(directory/'replaced-copy.link'))
            os.link(source, directory/'replaced-copy.link'); sync_dir(directory)
        if source.exists() and source.lstat().st_ino == backup.lstat().st_ino:
            # rename(backup,source) would be a no-op for two links to one inode.
            event(state,'rollback_redundant_backup_unlink_intent',path=str(backup))
            backup.unlink(); sync_dir(directory)
        else:
            event(state,'rollback_restore_intent',source=str(backup),destination=str(source))
            os.replace(backup, source); sync_dir(source.parent); sync_dir(directory)
    else:
        require(source.exists() and source.lstat().st_ino == intent['before']['inode'],
                'Missing rollback link and source inode differs; preserve evidence for manual review')
    exact(source, entry); stable(intent['before'], metadata(source))
    immutable(directory/'rolled_back.json', {'time_utc':now(),'reason':reason,'source':metadata(source)})
    event(state, 'transaction_rolled_back', transaction=directory.name, file_id=entry['id'], reason=reason)


def recover_pending(state, plan):
    count = 0
    for directory, intent, entry in transactions(state, plan):
        if terminal(directory) is None:
            rollback(state, directory, intent, entry, reason='recover durable incomplete intent'); count += 1
    return count


def native_copy(source, target, guard, seconds):
    require(0 < seconds <= 60, 'Per-file timeout must be finite and at most 60 seconds')
    env = {k:v for k,v in os.environ.items() if not k.startswith('DITTO')}
    output = target.parent/'ditto.log'
    with output.open('xb') as stream:
        child = subprocess.Popen(['/usr/bin/ditto','--noclone','--rsrc','--extattr','--acl','--hfsCompression',
                                  str(source),str(target)], stdout=stream, stderr=subprocess.STDOUT,
                                  env=env)
        started = time.monotonic()
        try:
            while child.poll() is None:
                guard.check(child=child)
                require(time.monotonic()-started < seconds, 'Per-file ditto timeout')
                time.sleep(.1)
            require(child.returncode == 0, 'Native ditto failed; retained log contains details')
        finally:
            if child.poll() is None:
                child.terminate()
                try: child.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    child.kill(); child.wait(timeout=3)
            stream.flush(); os.fsync(stream.fileno())
    sync_file(target); sync_dir(target.parent)


def compress_one(state, plan, entry, guard, *, per_file_seconds=30, copier=native_copy, fault=lambda stage: None):
    source = no_symlinks(Path(entry['path']), Path(plan['jobs'][entry['job_index']]['path']))
    guard.check(force=True, headroom=SCRATCH+2*entry['bytes'])
    before = safe_record(source, entry)
    if before['flags'] & UF_COMPRESSED:
        event(state,'already_compressed_exact_skip',file_id=entry['id'],current=before); return 'already_compressed'
    directory = state/'transactions'/(entry['id']+'-'+uuid.uuid4().hex[:12]); directory.mkdir(); sync_dir(directory.parent)
    require(directory.stat().st_dev == before['device'], 'Staging and source must share one filesystem')
    temporary, backup = directory/'compressed.tmp', directory/'original.link'
    intent = {'time_utc':now(),'file_id':entry['id'],'plan_sha256':sha(state/'plan.json'),'source':str(source),
              'temporary':str(temporary),'backup':str(backup),'before':before,
              'expected_sha256':entry['sha256'],'expected_bytes':entry['bytes'],'writer_pid':os.getpid()}
    immutable(directory/'intent.json',intent); event(state,'transaction_intent',transaction=directory.name,file_id=entry['id'])
    try:
        fault('intent'); copier(source,temporary,guard,per_file_seconds)
        exact(temporary,entry); packed = metadata(temporary); stable(before,packed)
        immutable(directory/'copy_verified.json', {'time_utc':now(),'temporary':packed,'sha256':sha(temporary)})
        fault('copy_verified')
        if not packed['flags'] & UF_COMPRESSED or packed['allocated_bytes'] >= before['allocated_bytes']:
            immutable(directory/'skipped.json',{'time_utc':now(),'reason':'no effective storage reduction','temporary':packed})
            event(state,'transaction_skipped',transaction=directory.name); return 'no_savings'
        guard.check(force=True,headroom=SCRATCH)
        current = safe_record(source,entry)
        require(current['inode'] == before['inode'] and current['ctime_ns'] == before['ctime_ns'], 'Source changed during staging')
        immutable(directory/'backup_link_intent.json',{'time_utc':now(),'source_inode':before['inode'],'backup':str(backup)})
        os.link(source,backup); sync_dir(directory); sync_file(backup)
        require(source.lstat().st_ino == backup.lstat().st_ino == before['inode'], 'Source inode changed before replacement')
        exact(backup,entry); stable(before,metadata(backup))
        immutable(directory/'replace_intent.json',{'time_utc':now(),'source_inode':before['inode'],'temporary_inode':packed['inode']})
        fault('backup_retained'); os.replace(temporary,source); sync_dir(source.parent); sync_dir(directory)
        fault('replaced'); exact(source,entry); installed = metadata(source); stable(before,installed)
        job=plan['jobs'][entry['job_index']]
        verify_unit(Path(job['path']),entry['key'],job['units'][entry['key']],job['identity_hash'])
        fault('unit_verified')
        immutable(directory/'verified.json',{'time_utc':now(),'installed':installed,'logical_sha256':sha(source),
                  'allocated_byte_reduction':before['allocated_bytes']-installed['allocated_bytes']})
        event(state,'transaction_verified_backups_retained',transaction=directory.name,file_id=entry['id'])
        return 'compressed'
    except BaseException as exc:
        # Restore the currently incomplete transaction even when a resource stop
        # caused the exception. Previous committed backups remain untouched.
        if terminal(directory) is None:
            rollback(state,directory,intent,entry,reason=f'{type(exc).__name__}: {exc}')
        raise


def apply(args, state, plan, guard):
    require(args.quiescent, 'Explicit quiescent-store assertion is required at apply time')
    quiet_jobs(plan); recover_pending(state,plan); validate_chain(plan,guard,all_bytes=True)
    if args.scope == 'full':
        marker=state/'trial_finalized.json'
        require(marker.exists(), 'Full scope requires separately verified/finalized first10 trial')
        trial=read(marker)
        require(trial['plan_sha256']==sha(state/'plan.json') and trial['trial_file_ids']==plan['trial_file_ids'],
                'Finalized trial belongs to another plan or file selection')
    selected = plan['files'] if args.scope == 'full' else [e for e in plan['files'] if e['id'] in plan['trial_file_ids']]
    run_id=uuid.uuid4().hex; event(state,'apply_started',scope=args.scope,run_id=run_id,selected_files=len(selected))
    counts={}; last_quiet=time.monotonic()
    for entry in selected:
        if time.monotonic()-last_quiet >= 30: quiet_jobs(plan); last_quiet=time.monotonic()
        result=compress_one(state,plan,entry,guard,per_file_seconds=args.per_file_seconds)
        counts[result]=counts.get(result,0)+1
    summary={'time_utc':now(),'scope':args.scope,'run_id':run_id,'counts':counts,'selected_file_ids':[e['id'] for e in selected],
             'state_digest':state_digest(state,plan),'available_bytes':shutil.disk_usage(state).free,
             'backups_retained':True,'awaiting_independent_whole_store_verification':True}
    immutable(state/('apply-'+run_id+'.json'),summary)
    event(state,'apply_completed_backups_retained',**summary)
    print(json.dumps({'status':'applied_backups_retained','scope':args.scope,'counts':counts,'state':str(state)}),flush=True)


def verify(args,state,plan,guard):
    quiet_jobs(plan)
    require(all(intent.get('writer_pid') != os.getpid() for directory,intent,entry in transactions(state,plan)),
            'Whole-store verification must run in a different process from every mutation')
    token=state_digest(state,plan); evidence=validate_chain(plan,guard,all_bytes=True)
    require(token == state_digest(state,plan), 'Transactions changed during independent verification')
    for directory,intent,entry in transactions(state,plan):
        if (directory/'original.link').exists():
            exact(directory/'original.link',entry)
            require((directory/'original.link').stat().st_ino == intent['before']['inode'], 'Backup identity changed')
    report={'status':'passed','time_utc':now(),'pid':os.getpid(),'plan_sha256':sha(state/'plan.json'),
            'state_digest':token,'verification_implementation':'session_independent_audit.verify_unit, separate command/process',
            'verification_implementation_sha256':sha(ROOT/'scripts/session_independent_audit.py'),
            'jobs':[{'name':j['name'],'identity_file_sha256':j['identity_file_sha256'],'completion_sha256':j['completion_sha256']} for j in plan['jobs']],
            **evidence,'available_bytes':shutil.disk_usage(state).free}
    path=state/('verification-'+uuid.uuid4().hex+'.json'); immutable(path,report)
    event(state,'whole_store_verification_passed',report=str(path),report_sha256=sha(path))
    print(json.dumps({'status':'verified','report':str(path),'sha256':sha(path)}),flush=True)


def finalize(args,state,plan,guard):
    require(args.verification is not None, 'Finalize requires a pinned independent verification report')
    report_path=no_symlinks(args.verification.resolve(strict=True),state); report=read(report_path)
    require(args.verification_sha256 and sha(report_path)==args.verification_sha256,
            'Root-reviewed verification SHA256 is required and must match')
    require(report['status']=='passed' and report['pid']!=os.getpid() and report['plan_sha256']==sha(state/'plan.json')
            and report['state_digest']==state_digest(state,plan)
            and report['verification_implementation_sha256']==plan['independent_verifier_sha256'],
            'Verification is not a fresh-process complete pass for this exact state')
    expected_jobs=[{'name':j['name'],'identity_file_sha256':j['identity_file_sha256'],
                    'completion_sha256':j['completion_sha256']} for j in plan['jobs']]
    expected_bytes=sum(value['bytes'] for job in plan['jobs'] for reference in job['units'].values()
        for value in read(Path(job['path'])/reference['path']/'complete.json')['files'].values())
    require(report['jobs']==expected_jobs and report['verified_units']==sum(len(j['units']) for j in plan['jobs'])
            and report['verified_logical_bytes']==expected_bytes
            and all(intent.get('writer_pid')!=report['pid'] for directory,intent,entry in transactions(state,plan)),
            'Verification report omits store/files or was written by a mutating process')
    quiet_jobs(plan); validate_chain(plan,guard)
    removed=0
    for directory,intent,entry in transactions(state,plan):
        guard.check(); source=Path(entry['path']); exact(source,entry); stable(entry['metadata'],metadata(source))
        for name in ('original.link','compressed.tmp','replaced-copy.link'):
            path=directory/name
            if not path.exists(): continue
            require(not path.is_symlink() and stat.S_ISREG(path.lstat().st_mode),'Unexpected owned cleanup path')
            # Never discard a mismatching copy retained as failure evidence.
            if path.stat().st_size != entry['bytes'] or sha(path)!=entry['sha256']:
                event(state,'nonredundant_failure_evidence_retained',path=str(path)); continue
            observed=metadata(path)
            try: stable(entry['metadata'],observed)
            except ValueError:
                event(state,'metadata_failure_evidence_retained',path=str(path),observed=observed)
                continue
            if name=='original.link': require(path.stat().st_ino==intent['before']['inode'],'Original backup inode differs')
            event(state,'owned_redundant_link_cleanup_intent',path=str(path),verification_sha256=sha(report_path))
            path.unlink(); sync_dir(directory); removed+=1
            event(state,'owned_redundant_link_cleanup_complete',path=str(path))
    applications=[read(p) for p in state.glob('apply-*.json')]
    trial=[r for r in applications if r['scope']=='trial' and set(r['selected_file_ids'])==set(plan['trial_file_ids'])]
    if trial and not (state/'trial_finalized.json').exists():
        immutable(state/'trial_finalized.json',{'time_utc':now(),'verification_sha256':sha(report_path),
                  'trial_file_ids':plan['trial_file_ids'],'plan_sha256':sha(state/'plan.json')})
    evidence={'time_utc':now(),'verification_sha256':sha(report_path),'removed_owned_redundant_links':removed,
              'available_bytes':shutil.disk_usage(state).free,'logical_result_files_deleted':0}
    immutable(state/('finalization-'+uuid.uuid4().hex+'.json'),evidence)
    event(state,'finalization_complete',**evidence); print(json.dumps({'status':'finalized',**evidence}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['freeze','apply','recover','verify','finalize'])
    parser.add_argument('--session',type=Path)
    parser.add_argument('--state',type=Path,required=True)
    parser.add_argument('--approved-stores',type=Path,help='Optional new explicit pinned allowlist; never expands an existing plan')
    parser.add_argument('--scope',choices=['trial','full'],default='trial')
    parser.add_argument('--quiescent',action='store_true')
    parser.add_argument('--verification',type=Path)
    parser.add_argument('--plan-sha256',help='Exact frozen plan hash reviewed by root; required after freeze')
    parser.add_argument('--verification-sha256',help='Exact independent verification report hash; required for finalize')
    parser.add_argument('--timeout',type=float,default=1800)
    parser.add_argument('--per-file-seconds',type=float,default=30)
    args=parser.parse_args()
    require(sys.platform=='darwin','Native macOS is required')
    signal.signal(signal.SIGTERM, lambda signum, frame: (_ for _ in ()).throw(SystemExit(128+signum)))
    # This macOS/Python rejects attempts to reduce RLIMIT_DATA. Admission and
    # every active ditto poll enforce the aggregate 512 MiB observed-RSS gate;
    # hashes stream 1 MiB blocks and eligible files are bounded at 64 MiB.
    if args.action=='freeze':
        require(args.session is not None,'Freeze requires the session path'); freeze(args); return
    state=args.state.resolve(strict=True)
    with state_lock(state):
        plan=load_plan(state,args.plan_sha256); guard=Guard(Path(plan['session']),state,args.timeout,plan['guard'])
        if args.action=='apply': apply(args,state,plan,guard)
        elif args.action=='recover':
            quiet_jobs(plan); print(json.dumps({'recovered':recover_pending(state,plan)}))
        elif args.action=='verify': verify(args,state,plan,guard)
        else: finalize(args,state,plan,guard)


if __name__=='__main__': main()
