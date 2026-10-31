"""Bounded job execution and explicit optional-failure isolation."""
from __future__ import annotations
import datetime as dt
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
    result.update(job_id=job_id, attempt=attempt, elapsed_seconds=time.monotonic() - start)
    append_event(ledger, {"event": "job_attempt_finished", **result})
    if not optional and result["state"] != "complete":
        raise RuntimeError(result)
    return result


def run_guarded(command, *, interpreter, guard, timeout, deadline_utc, log_path):
    """The prepared guard retains power/memory/disk ownership and thresholds."""
    remaining = (dt.datetime.fromisoformat(deadline_utc) - dt.datetime.now(dt.timezone.utc)).total_seconds()
    if remaining <= 0:
        raise BlockedJob("Session hard deadline reached")
    allowed = min(float(timeout), remaining, 43200)
    args = [str(interpreter), str(guard), "--timeout", str(allowed), "--", *map(str, command)]
    path = Path(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as log:
        process = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT)
        # This assertion is finite and ends if the actual guarded runner exits.
        awake = None
        try:
            awake = subprocess.Popen(["/usr/bin/caffeinate", "-i", "-s", "-t", str(max(1, int(allowed))), "-w", str(process.pid)])
            code = process.wait(timeout=allowed + 30)
        except subprocess.TimeoutExpired:
            # The guard gets its own cleanup opportunity; never target other jobs.
            process.terminate()
            code = process.wait(timeout=20)
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=20)
            if awake is not None and awake.poll() is None:
                awake.terminate()
                awake.wait(timeout=5)
    if code != 0:
        raise RuntimeError(f"Guarded job exit {code}; inspect {path}")
    return {"returncode": code, "guard_log": str(path), "timeout_seconds": allowed}
