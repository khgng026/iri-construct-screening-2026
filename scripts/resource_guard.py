"""Read-only preflight, or monitor ONLY a command explicitly supplied by its caller.
No experiment, downloads, power assertions, or system setting changes by default.
"""
import argparse
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent

def capture(*args):
    return subprocess.check_output(args, text=True, timeout=10).strip()

def snapshot():
    if sys.platform != "darwin":
        import psutil
        memory = psutil.virtual_memory()
        battery = psutil.sensors_battery()
        return dict(time=datetime.datetime.now().astimezone().isoformat(),
                    pressure=1 if memory.available >= max(2 * 2**30, memory.total * .15) else 2,
                    swap_mib=psutil.swap_memory().used / 2**20,
                    disk_free_gib=shutil.disk_usage(ROOT).free / 2**30,
                    ac_power=battery is None or battery.power_plugged,
                    power="psutil battery status; no battery treated as stationary host")
    swap = capture('/usr/sbin/sysctl', '-n', 'vm.swapusage')
    power = capture('/usr/bin/pmset', '-g', 'batt')
    return dict(time=datetime.datetime.now().astimezone().isoformat(),
                pressure=int(capture('/usr/sbin/sysctl', '-n', 'kern.memorystatus_vm_pressure_level')),
                swap_mib=float(re.search(r'used =\s*([\d.]+)M', swap)[1]),
                disk_free_gib=shutil.disk_usage(ROOT).free / 2**30,
                ac_power="'AC Power'" in power, power=power)

def reasons(s, baseline=None, require_ac=True):
    r = []
    if s['pressure'] != 1:
        r.append('memory pressure is not normal')
    if s['disk_free_gib'] < 20:
        r.append('less than 20 GiB disk free')
    if require_ac and not s['ac_power']:
        r.append('AC power disconnected')
    if baseline and s['swap_mib'] - baseline['swap_mib'] >= 512:
        r.append('swap grew by at least 512 MiB')
    return r

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--timeout', type=float, default=300, help='seconds; finite, max 12h')
    p.add_argument('--interval', type=float, default=2)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    if not (0 < a.timeout <= 43200 and 0.5 <= a.interval <= 30):
        p.error('timeout must be 0..43200s and interval 0.5..30s')
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    s = snapshot()
    print(json.dumps({'snapshot': s, 'blocked_by': reasons(s)}), flush=True)
    if not command:
        return 2 if reasons(s) else 0
    # Advisory lock coordinates users of THIS wrapper, not unrelated processes.
    with (ROOT / 'compute.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('Another guarded job holds the compute lock.'); return 2
        if reasons(s):
            return 2
        run = ROOT / 'logs' / datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')
        run.mkdir(parents=True)
        (run / 'command.json').write_text(json.dumps(command))
        with (run / 'output.log').open('w') as output, (run / 'resources.jsonl').open('w') as log:
            child = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
            start = time.monotonic()
            stop = []
            try:
                while child.poll() is None:
                    current = snapshot()
                    log.write(json.dumps(current) + '\n'); log.flush()
                    stop = reasons(current, s)
                    if time.monotonic() - start >= a.timeout:
                        stop.append('wall-time limit')
                    if stop:
                        break
                    time.sleep(a.interval)
            finally:
                # Never signal any process not in the child group created here.
                if child.poll() is None:
                    os.killpg(child.pid, signal.SIGTERM)
                    try:
                        child.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL); child.wait()
            result = dict(exit_code=child.returncode, stopped_for=stop, log_dir=str(run))
            (run / 'status.json').write_text(json.dumps(result, indent=2))
            print(json.dumps(result))
            return 2 if stop else child.returncode

if __name__ == '__main__':
    sys.exit(main())
