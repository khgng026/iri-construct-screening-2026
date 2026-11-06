"""Bounded job execution and explicit optional-failure isolation."""
from __future__ import annotations
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import traceback
import uuid

from .session_store import append_event, atomic_json, utcnow


class BlockedJob(Exception):
    pass


class OwnedCleanupFailure(BaseException):
    """Stop the queue when owned-process cleanup cannot be established."""


def cleanup_owned_group(pgid, *, grace_seconds=5.0):
    """Terminate only a group recorded when this supervisor created a session.

    The group can outlive its leader. Checking only Popen.poll()/wait() misses a
    descendant that ignores SIGTERM after the leader has already exited.
    """
    if type(pgid) is not int or pgid <= 1 or pgid == os.getpgrp():
        raise OwnedCleanupFailure('Refusing invalid or supervisor process group')
    result={'pgid':pgid,'term_sent':False,'kill_sent':False,'live_pids_after':[]}
    def exists():
        try: os.killpg(pgid,0);return True
        except ProcessLookupError:return False
    if not exists():
        result['quiescent']=True;return result
    try: os.killpg(pgid,signal.SIGTERM);result['term_sent']=True
    except ProcessLookupError:
        result['quiescent']=True;return result
    limit=time.monotonic()+grace_seconds
    while exists() and time.monotonic()<limit:time.sleep(.1)
    if exists():
        try:os.killpg(pgid,signal.SIGKILL);result['kill_sent']=True
        except ProcessLookupError:pass
    limit=time.monotonic()+2
    while exists() and time.monotonic()<limit:time.sleep(.1)
    if exists():
        # An already dead orphan can briefly remain a zombie. Inspect only group
        # membership/status; do not persist unrelated process commands or data.
        snapshot=subprocess.check_output(['/bin/ps','-axo','pid=,pgid=,stat='],text=True,timeout=5)
        for line in snapshot.splitlines():
            fields=line.split()
            if len(fields)==3 and int(fields[1])==pgid and fields[2][0] not in 'ZX':
                result['live_pids_after'].append(int(fields[0]))
    result['quiescent']=not result['live_pids_after']
    if not result['quiescent']:
        raise OwnedCleanupFailure('Owned group survived bounded cleanup: '+str(result))
    return result


def _finish_guard(process, ownership_path, token):
    """Allow wrapper cleanup, then check its recorded child groups independently."""
    forced=False
    if process.poll() is None:
        process.terminate()
        try:process.wait(timeout=20)
        except subprocess.TimeoutExpired:forced=True
    outer=cleanup_owned_group(process.pid,grace_seconds=0)
    if process.poll() is None:process.wait(timeout=5)
    recorded=[]
    if ownership_path.exists():
        ownership=json.loads(ownership_path.read_text())
        if ownership.get('token')!=token or ownership.get('wrapper_pid')!=process.pid:
            raise OwnedCleanupFailure('Guard ownership record identity mismatch')
        for entry in ownership['child_groups']:
            if entry.get('pid')!=entry.get('pgid') or entry.get('start_new_session') is not True:
                raise OwnedCleanupFailure('Invalid recorded child-group ownership')
            recorded.append(cleanup_owned_group(entry['pgid'],grace_seconds=2))
    else:
        raise OwnedCleanupFailure('Guard ownership record absent; cannot establish descendant cleanup')
    return {'guard_forced':forced,'guard_group':outer,'recorded_child_groups':recorded}


def _finish_awake(process):
    if process.poll() is None:
        process.terminate()
        try:process.wait(timeout=2)
        except subprocess.TimeoutExpired:pass
    result=cleanup_owned_group(process.pid,grace_seconds=0)
    if process.poll() is None:process.wait(timeout=5)
    return result


def attempt_job(job_id, action, ledger, *, optional=True):
    attempt = uuid.uuid4().hex
    start = time.monotonic()
    append_event(ledger, {"event": "job_attempt", "job_id": job_id, "attempt": attempt, "state": "running"})
    try:
        outputs = action()
        result = {"state": "complete", "outputs": outputs}
    except BlockedJob as exc:
        result = {"state": "blocked", "exception": str(exc)}
    except Exception as exc:
        result = {"state": "failed", "exception": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}
    except BaseException as exc:
        append_event(ledger, {"event": "job_attempt_interrupted", "job_id": job_id,
                     "attempt": attempt, "state": "failed", "exception": f"{type(exc).__name__}: {exc}",
                     "elapsed_seconds": time.monotonic()-start})
        raise
    result.update(job_id=job_id, attempt=attempt, elapsed_seconds=time.monotonic() - start)
    append_event(ledger, {"event": "job_attempt_finished", **result})
    if not optional and result["state"] != "complete":
        raise RuntimeError(result)
    return result


def remaining_session_seconds(deadline_utc, anchor_path=None):
    remaining = (dt.datetime.fromisoformat(deadline_utc) - dt.datetime.now(dt.timezone.utc)).total_seconds()
    if anchor_path is not None and Path(anchor_path).exists():
        anchor=json.loads(Path(anchor_path).read_text())
        if anchor['hard_stop_utc'] != deadline_utc:
            raise BlockedJob('Session deadline and immutable monotonic anchor disagree')
        boot=hashlib.sha256(subprocess.check_output(['/usr/sbin/sysctl','-n','kern.bootsessionuuid'],timeout=10)).hexdigest()
        if boot != anchor['boot_session_sha256']:
            raise BlockedJob('Machine reboot invalidated session monotonic anchor; preserve work for a fresh readiness handoff')
        remaining=min(remaining,anchor['monotonic_hard_stop']-time.monotonic())
    return remaining


def run_guarded(command, *, interpreter, guard, timeout, deadline_utc, log_path):
    """The prepared guard retains power/memory/disk ownership and thresholds."""
    path = Path(log_path)
    remaining = remaining_session_seconds(deadline_utc,path.parent.parent/'inputs/deadline_anchor.json')
    if remaining <= 0:
        raise BlockedJob("Session hard deadline reached")
    allowed = min(float(timeout), remaining, 43200)
    args = [str(interpreter), str(guard), "--timeout", str(allowed), "--", *map(str, command)]
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path=path.with_name(path.stem+'.attempt-'+uuid.uuid4().hex+path.suffix)
    ownership_path=path.with_name(path.name+'.ownership.json')
    cleanup_path=path.with_name(path.name+'.cleanup.json')
    token=uuid.uuid4().hex
    environment=dict(os.environ,IRI_GUARD_OWNERSHIP_PATH=str(ownership_path.resolve()),IRI_GUARD_OWNERSHIP_TOKEN=token)
    with path.open("xb") as log:
        process = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True,env=environment)
        # This assertion is finite and ends if the actual guarded runner exits.
        awake = None
        cleanup={}
        try:
            awake = subprocess.Popen(["/usr/bin/caffeinate", "-i", "-s", "-t", str(max(1, int(allowed))), "-w", str(process.pid)],start_new_session=True)
            code = process.wait(timeout=allowed + 30)
        except subprocess.TimeoutExpired:
            raise RuntimeError('Guard exceeded its finite supervisory wait')
        finally:
            try:
                cleanup['guard']=_finish_guard(process,ownership_path,token)
            except BaseException as exc:
                cleanup['guard_error']=f'{type(exc).__name__}: {exc}'
                raise OwnedCleanupFailure(cleanup['guard_error']) from exc
            finally:
                # Guard cleanup failures must not skip power-assertion cleanup.
                try:
                    if awake is not None:cleanup['caffeinate']=_finish_awake(awake)
                except BaseException as exc:
                    cleanup['caffeinate_error']=f'{type(exc).__name__}: {exc}'
                    raise OwnedCleanupFailure(cleanup['caffeinate_error']) from exc
                finally:atomic_json(cleanup_path,cleanup)
    if code != 0:
        raise RuntimeError(f"Guarded job exit {code}; inspect {path}")
    return {"returncode": code, "guard_log": str(path), "timeout_seconds": allowed,
            "ownership_record":str(ownership_path),'cleanup_record':str(cleanup_path)}
