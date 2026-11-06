#!/usr/bin/env python3
"""Finite actual-dispatcher termination/restart trial; execute only in a lead slot.

Two toy cases check cooperative and SIGTERM-resistant descendants. Every output,
queue and numeric checkpoint belongs to a fresh ignored directory. Existing guard
evidence and the prepared guard source are hashed before and after, never edited.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
PREPARED_GUARD = ROOT / "scripts/resource_guard.py"


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    temporary = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    with temporary.open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def ps_snapshot():
    text = subprocess.check_output(
        ["/bin/ps", "-ww", "-axo", "pid=,ppid=,pgid=,stat=,lstart=,command="], text=True, timeout=5)
    rows = {}
    for line in text.splitlines():
        fields = line.strip().split(None, 9)
        if len(fields) != 10:
            continue
        pid, ppid, pgid = map(int, fields[:3])
        rows[pid] = {"pid": pid, "ppid": ppid, "pgid": pgid, "state": fields[3],
                     "started": " ".join(fields[4:9]), "command": fields[9]}
    return rows


def owned_snapshot(dispatcher_pid, rows):
    owned = {dispatcher_pid} if dispatcher_pid in rows else set()
    while True:
        added = {pid for pid, row in rows.items() if row["ppid"] in owned} - owned
        if not added:
            break
        owned.update(added)
    return {pid: rows[pid] for pid in owned}


def same_process(left, right):
    # argv is display evidence, not a birth identity: macOS briefly reports
    # '(caffeinate)' while its state can still be S during exit. Preserve the
    # original argv and match the recorded PID/group/start-time identity.
    return right is not None and all(left[k] == right[k] for k in ("pid", "pgid", "started"))


def live_owned(tracked, rows=None):
    rows = ps_snapshot() if rows is None else rows
    return [pid for pid, row in tracked.items() if same_process(row, rows.get(pid))
            and rows[pid]["state"][0] not in "ZX"]


def collect(dispatcher, tracked):
    rows = ps_snapshot()
    for pid, row in owned_snapshot(dispatcher.pid, rows).items():
        # macOS can replace argv with a defunct-process label after exit. It is
        # no longer a live identity to refresh or an emergency-cleanup target.
        if row["state"][0] in "ZX":
            continue
        if pid in tracked:
            require(same_process(tracked[pid], row), "Live tracked owned PID changed identity: " +
                    json.dumps({"before": tracked[pid], "after": row}, sort_keys=True))
        else:
            tracked[pid] = row
    return rows


def wait_exit(process, tracked, timeout):
    end = time.monotonic() + timeout
    while process.poll() is None and time.monotonic() < end:
        collect(process, tracked)
        time.sleep(.1)
    require(process.poll() is not None, "Owned dispatcher exceeded finite exit wait")
    return process.returncode


def wait_quiescent(tracked, timeout):
    end = time.monotonic() + timeout
    while True:
        live = live_owned(tracked)
        if not live or time.monotonic() >= end:
            return live
        time.sleep(.1)


def emergency_cleanup(tracked):
    """Failure cleanup only, using freshly matched recorded identities, never names."""
    actions = []
    for sig in (signal.SIGTERM, signal.SIGKILL):
        rows = ps_snapshot()
        for pid in live_owned(tracked, rows):
            require(pid != os.getpid(), "Refusing to signal trial controller")
            try:
                os.kill(pid, sig)
                actions.append({"pid": pid, "signal": int(sig)})
            except ProcessLookupError:
                pass
        if sig == signal.SIGTERM:
            wait_quiescent(tracked, 2)
    return {"actions": actions, "live_after": wait_quiescent(tracked, 3)}


def caffeinate_rows(tracked):
    result = []
    for row in tracked.values():
        words = row["command"].split()
        if words and Path(words[0]).name == "caffeinate":
            require("-t" in words and "-w" in words, "Caffeinate missing finite duration/PID scope")
            seconds = int(words[words.index("-t") + 1])
            watched = int(words[words.index("-w") + 1])
            require(0 < seconds <= 90 and watched in tracked, "Caffeinate scope/duration differs from owned case")
            result.append({**row, "finite_seconds": seconds, "watched_pid": watched})
    return result


def toy(args):
    root, token = args.case_root.resolve(), args.token
    require(len(token) == 32 and all(c in "0123456789abcdef" for c in token), "Invalid owned toy token")
    require(read(root / "toy-identity.json")["token"] == token, "Toy identity mismatch")
    if args.toy_after:
        write_json(root / "after-marker.json", {"token": token, "completed": True})
        return 0
    if args.toy_grandchild:
        if args.resistant:
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
        write_json(root / "grandchild-ready.json", {"pid": os.getpid(), "pgid": os.getpgrp(),
                   "ppid": os.getppid(), "token": token, "resistant": args.resistant})
        end = time.monotonic() + 60
        while time.monotonic() < end:
            time.sleep(.1)
        return 0
    checkpoint = root / "numeric-prefix.json"
    if not checkpoint.exists():
        write_json(checkpoint, {"token": token, "completed_numeric_prefix": [1, 4, 9]})
    require(read(checkpoint) == {"token": token, "completed_numeric_prefix": [1, 4, 9]}, "Numeric checkpoint changed")
    if (root / "resume-requested.json").exists():
        require(read(root / "resume-requested.json")["token"] == token, "Resume identity mismatch")
        write_json(root / "resumed.json", {"token": token, "prefix_sha256": sha(checkpoint), "completed": True})
        return 0
    command = [sys.executable, str(SCRIPT), "--toy-grandchild", "--case-root", str(root), "--token", token]
    if args.resistant:
        command.append("--resistant")
    child = subprocess.Popen(command)  # inherit the actual guard-created worker group
    end = time.monotonic() + 10
    while not (root / "grandchild-ready.json").exists():
        require(child.poll() is None and time.monotonic() < end, "Toy descendant failed readiness")
        time.sleep(.05)
    write_json(root / "worker-ready.json", {"pid": os.getpid(), "pgid": os.getpgrp(),
               "guard_pid": os.getppid(), "grandchild_pid": child.pid, "token": token})
    # Default SIGTERM intentionally exits the leader promptly. The resistant
    # descendant checks whether the supervisors clean its recorded group too.
    try:
        child.wait(timeout=60)
    except subprocess.TimeoutExpired:
        child.kill()
        child.wait(timeout=3)
    return 0


def execute_case(directory, name, parent_deadline, design):
    directory.mkdir(exist_ok=False)
    for child in ("inputs", "logs", "audit", "jobs"):
        (directory / child).mkdir()
    token = uuid.uuid4().hex
    resistant = name == "sigterm_resistant_descendant"
    write_json(directory / "toy-identity.json", {"token": token, "case": name})
    deadline = min(parent_deadline, dt.datetime.now(dt.timezone.utc) + dt.timedelta(seconds=90))
    remaining = (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds()
    require(remaining >= 75, "Insufficient session time for a finite dispatcher case")
    write_json(directory / "launch.json", {"hard_stop_utc": deadline.isoformat()})
    boot = hashlib.sha256(subprocess.check_output(["/usr/sbin/sysctl", "-n", "kern.bootsessionuuid"], timeout=5)).hexdigest()
    write_json(directory / "inputs/deadline_anchor.json", {"hard_stop_utc": deadline.isoformat(),
               "boot_session_sha256": boot, "monotonic_hard_stop": time.monotonic() + remaining})
    first_command = [sys.executable, str(SCRIPT), "--toy-worker", "--case-root", str(directory), "--token", token]
    if resistant:
        first_command.append("--resistant")
    queue = {"jobs": [
        {"id": "toy-block", "command": first_command, "timeout": 60, "required_file": "resumed.json"},
        {"id": "toy-after", "command": [sys.executable, str(SCRIPT), "--toy-after", "--case-root", str(directory), "--token", token],
         "timeout": 60, "requires": ["toy-block"], "required_file": "after-marker.json"}]}
    queue_path = directory / "inputs/toy-queue.json"
    write_json(queue_path, queue)
    command = [sys.executable, str(ROOT / "scripts/session_dispatch.py"), "--session", str(directory), "--queue", str(queue_path)]
    record = {"case": name, "status": "running", "directory": str(directory), "command": command,
              "tracked": {}, "controller_cleanup_used": False}
    tracked, processes = {}, []
    try:
        with (directory / "dispatcher-first.log").open("xb") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
            processes.append(process)
            ready_limit = time.monotonic() + design["readiness_timeout_seconds"]
            while True:
                collect(process, tracked)
                require(process.poll() is None, "Dispatcher exited before owned-worker readiness; inspect guard refusal/logs")
                if (directory / "worker-ready.json").exists() and len(caffeinate_rows(tracked)) >= 2:
                    break
                require(time.monotonic() < ready_limit, "Dispatcher/worker/caffeinate readiness timed out")
                time.sleep(.1)
            ready = read(directory / "worker-ready.json")
            grandchild = read(directory / "grandchild-ready.json")
            require(ready["token"] == grandchild["token"] == token and ready["pid"] == ready["pgid"] == grandchild["pgid"], "Toy group ownership mismatch")
            require(ready["grandchild_pid"] == grandchild["pid"] and all(pid in tracked for pid in
                    (process.pid, ready["pid"], ready["guard_pid"], grandchild["pid"])), "Ready process missing from actual dispatcher descendants")
            record["ready"] = ready
            record["finite_caffeinate_before_signal"] = caffeinate_rows(tracked)
            prefix_hash = sha(directory / "numeric-prefix.json")
            record["numeric_prefix_sha256"] = prefix_hash
            write_json(directory / "owned-before-signal.json", list(tracked.values()))
            signalled = time.monotonic()
            process.send_signal(signal.SIGTERM)
            record["dispatcher_returncode"] = wait_exit(process, tracked, design["dispatcher_exit_timeout_seconds"])
            record["termination_seconds"] = time.monotonic() - signalled
            record["live_owned_after_dispatcher_exit"] = wait_quiescent(tracked, design["post_exit_quiescence_timeout_seconds"])
            require(record["dispatcher_returncode"] == 143, "Dispatcher did not preserve SIGTERM exit143")
            require(not record["live_owned_after_dispatcher_exit"], "Owned descendant or caffeinate survived dispatcher cleanup")
        ledger_path = directory / "job_ledger.jsonl"
        before = ledger_path.read_bytes()
        events = [json.loads(line) for line in before.splitlines()]
        attempts = [e for e in events if e["event"] == "job_attempt"]
        interrupted = [e for e in events if e["event"] == "job_attempt_interrupted"]
        require(len(attempts) == len(interrupted) == 1 and attempts[0]["attempt"] == interrupted[0]["attempt"]
                and interrupted[0]["state"] == "failed", "Ledger lacks exactly one matching interrupted attempt")
        require(not any(e["event"] == "job_attempt_finished" for e in events), "Interrupted action recorded a finished result")
        require(not any(e.get("job_id") == "toy-after" for e in events) and not (directory / "after-marker.json").exists(),
                "Following queue job started after termination")
        record["interruption_ledger_passed"] = True
        record["ledger_before_resume_sha256"] = hashlib.sha256(before).hexdigest()
        write_json(directory / "ledger-before-resume.json", events)
        write_json(directory / "resume-requested.json", {"token": token})
        with (directory / "dispatcher-resume.log").open("xb") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
            processes.append(process)
            record["resume_returncode"] = wait_exit(process, tracked, design["resume_timeout_seconds"])
            require(record["resume_returncode"] == 0, "Actual dispatcher restart failed")
            record["live_owned_after_resume"] = wait_quiescent(tracked, design["post_exit_quiescence_timeout_seconds"])
            require(not record["live_owned_after_resume"], "Owned process survived completed restart")
        after = ledger_path.read_bytes()
        require(after.startswith(before), "Restart rewrote interruption ledger history")
        events = [json.loads(line) for line in after.splitlines()]
        finished = [e for e in events if e["event"] == "job_attempt_finished"]
        require([(e["job_id"], e["state"]) for e in finished] == [("toy-block", "complete"), ("toy-after", "complete")],
                "Restart lacks exactly the expected successful queue attempts")
        require(sha(directory / "numeric-prefix.json") == prefix_hash and read(directory / "resumed.json")["prefix_sha256"] == prefix_hash,
                "Restart changed the durable numeric prefix")
        require(read(directory / "after-marker.json")["token"] == token, "Following toy marker identity mismatch")
        record["finite_caffeinate_all_observed"] = caffeinate_rows(tracked)
        record["resume_preserved_prefix_and_ledger"] = True
        record["status"] = "passed"
    except BaseException as error:
        record["status"] = "interrupted" if isinstance(error, (SystemExit, KeyboardInterrupt)) else "failed"
        record["error"] = f"{type(error).__name__}: {error}"
        if isinstance(error, (SystemExit, KeyboardInterrupt)):
            raise
    finally:
        for process in processes:
            if process.poll() is None:
                collect(process, tracked)
        record["tracked"] = list(tracked.values())
        if live_owned(tracked):
            record["controller_cleanup_used"] = True
            record["controller_cleanup"] = emergency_cleanup(tracked)
            if record["status"] == "passed":
                record["status"] = "failed"
        for process in processes:
            if process.poll() is None:
                process.wait(timeout=5)
        record["evidence_files"] = {str(path.relative_to(directory)): sha(path) for path in directory.rglob("*") if path.is_file()
                                    and path.name != "case-report.json" and ".tmp-" not in path.name}
        write_json(directory / "case-report.json", record)
    return record


def run(args):
    session = args.session.resolve()
    session.relative_to((ROOT / "results").resolve())
    design_path = session / "inputs/dispatcher_termination_trial_design_v1.json"
    design = read(design_path)
    require(design["version"] == "actual-dispatcher-termination-v1" and design["status"] == "frozen_before_execution", "Unfrozen trial design")
    require(sha(PREPARED_GUARD) == design["prepared_guard_sha256"], "Prepared guard changed")
    output = session / "audit" / ("independent_dispatch_termination_" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:10])
    output.mkdir(exist_ok=False)
    before = {relative: sha(session / relative) for relative in design["preserve_prior_evidence"]}
    report = {"version": design["version"], "status": "running", "output": str(output),
              "design_sha256": sha(design_path), "prior_evidence_sha256_before": before,
              "source_sha256": {str(path): sha(path) for path in (SCRIPT, PREPARED_GUARD,
                  ROOT / "scripts/session_dispatch.py", ROOT / "scripts/session_guard_entry.py", ROOT / "src/iri_constructs/session_jobs.py")},
              "cases": []}
    write_json(output / "frozen-design.json", design)
    write_json(output / "report.json", report)
    def interrupted(signum, frame):
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    signal.signal(signal.SIGALRM, interrupted)
    signal.setitimer(signal.ITIMER_REAL, design["controller_total_cap_seconds"])
    exit_code = 1
    try:
        deadline = dt.datetime.fromisoformat(read(session / "launch.json")["hard_stop_utc"])
        for name in design["cases"]:
            report["cases"].append(execute_case(output / name, name, deadline, design))
            write_json(output / "report.json", report)
            print(json.dumps({"case": name, "status": report["cases"][-1]["status"]}), flush=True)
        report["status"] = "passed" if all(case["status"] == "passed" for case in report["cases"]) else "failed"
        exit_code = 0 if report["status"] == "passed" else 1
    except BaseException as error:
        report["status"] = "interrupted" if isinstance(error, (SystemExit, KeyboardInterrupt)) else "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        exit_code = int(error.code) if isinstance(error, SystemExit) and isinstance(error.code, int) else 1
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        report["prior_evidence_sha256_after"] = {relative: sha(session / relative) for relative in before}
        report["prior_evidence_preserved"] = report["prior_evidence_sha256_after"] == before
        report["prepared_guard_unchanged"] = sha(PREPARED_GUARD) == design["prepared_guard_sha256"]
        if not report["prior_evidence_preserved"] or not report["prepared_guard_unchanged"]:
            report["status"] = "failed"
            exit_code = 1
        write_json(output / "report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output / "report.json")}))
    return exit_code


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", type=Path)
    parser.add_argument("--case-root", type=Path)
    parser.add_argument("--token")
    parser.add_argument("--resistant", action="store_true")
    roles = parser.add_mutually_exclusive_group()
    roles.add_argument("--toy-worker", action="store_true")
    roles.add_argument("--toy-grandchild", action="store_true")
    roles.add_argument("--toy-after", action="store_true")
    args = parser.parse_args()
    if args.toy_worker or args.toy_grandchild or args.toy_after:
        require(args.case_root is not None and args.token is not None, "Toy role requires owned root/token")
        sys.exit(toy(args))
    require(args.session is not None, "Trial requires --session")
    sys.exit(run(args))
