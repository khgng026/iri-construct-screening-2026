"""Run one cleared local stage with a wall/RSS ceiling and source-free costs."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import psutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or args.output.exists():
        parser.error('Supply a command and a new resource report path')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    sampled_peak, stop_reason = 0, None
    with args.output.with_suffix('.log').open('w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        watched = psutil.Process(process.pid)
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                process.returncode = os.waitstatus_to_exitcode(status)
                break
            try:
                rss = watched.memory_info().rss
                sampled_peak = max(sampled_peak, rss)
            except psutil.NoSuchProcess:
                rss = 0
            if stop_reason is None and (rss > 6 * 1024**3 or time.monotonic() - started > 1800):
                stop_reason = 'resource_ceiling_exceeded'
                os.killpg(process.pid, signal.SIGKILL)
            time.sleep(0.02)
    peak = usage.ru_maxrss * 1024  # Linux wait4 reports KiB.
    if peak > 6 * 1024**3:
        stop_reason = 'resource_ceiling_exceeded'
    report = {'command': command, 'exit_code': process.returncode, 'stop_reason': stop_reason,
              'wall_seconds': time.monotonic() - started,
              'process_cpu_seconds': usage.ru_utime + usage.ru_stime,
              'os_peak_process_rss_bytes': peak, 'sampled_peak_process_rss_bytes': sampled_peak,
              'limits': {'wall_seconds': 1800, 'process_rss_bytes': 6 * 1024**3},
              'rss_policy': '20ms sampling with termination; OS high-water mark checked after exit',
              'gpu_seconds': None, 'monetary_cost': None}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: report[k] for k in ['exit_code', 'stop_reason', 'wall_seconds', 'os_peak_process_rss_bytes']}))
    return process.returncode if process.returncode else (1 if stop_reason else 0)


if __name__ == '__main__':
    raise SystemExit(main())
