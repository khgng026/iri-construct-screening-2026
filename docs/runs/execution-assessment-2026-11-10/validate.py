"""Reproduce software/environment checks only; this never grants run clearance."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.runtime import local_output, run_metadata, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = local_output(args.output)
    output.mkdir(parents=True, exist_ok=False)
    report = run_metadata('execution_request_technical_validation', {
        'scientific_experiments_requested': False,
        'biological_inference_requested': False,
        'downloads': False,
    })
    report['validation_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    commands = [
        ('smoke', ['scripts/smoke_test.py']),
        ('unit-tests', ['-m', 'unittest', 'discover', '-s', 'tests', '-v']),
        ('historical-release', ['scripts/verify_pilot_release.py']),
        ('dependency-consistency', ['-m', 'pip', 'check']),
        ('environment', ['scripts/t05_environment.py', '--load-esm2',
                         '--output', str(output / 'environment')]),
    ]
    report['checks'] = []
    for name, arguments in commands:
        argv = [sys.executable, *arguments]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.perf_counter()
        with (output / f'{name}.txt').open('w', encoding='utf-8') as log:
            result = subprocess.run(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        report['checks'].append({
            'name': name, 'argv': argv, 'exit_code': result.returncode,
            'wall_seconds': time.perf_counter() - started,
            'child_user_cpu_seconds': after.ru_utime - before.ru_utime,
            'child_system_cpu_seconds': after.ru_stime - before.ru_stime,
            'log': f'{name}.txt',
        })
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    report['validation_cost_scope'] = 'Software checks only; not E2/E3/E4 scientific model costs.'
    report['child_peak_rss_bytes'] = peak * 1024 if platform.system() == 'Linux' else None
    report['child_peak_rss_scope'] = 'OS high-water mark across child processes; not additive or a per-model measurement.'
    report['gpu_execution_seconds'] = None
    report['monetary_cost'] = None
    report['inference_executed_on_biological_inputs'] = False
    report['new_scientific_scores_produced'] = False
    report['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
    report['status'] = ('software_checks_passed' if all(c['exit_code'] == 0 for c in report['checks'])
                        else 'software_check_failed')
    report['limitations'].append('Environment preflight may exit zero with unavailable model dependencies; inspect environment/report.json.')
    write_json(output / 'validation.json', report)
    paths = sorted(p for p in output.rglob('*') if p.is_file())
    write_json(output / 'checksums.json', {
        'sha256': {p.relative_to(output).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in paths}})
    print(report['status'], output)
    return 0 if report['status'] == 'software_checks_passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
