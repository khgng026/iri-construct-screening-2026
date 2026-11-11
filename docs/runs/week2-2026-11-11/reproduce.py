"""Reproduce technical handoff checks; never run scored biological experiments."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.runtime import local_output, read_json, run_metadata, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--curated-inputs', type=Path)
    parser.add_argument('--supplemental-inputs', type=Path)
    args = parser.parse_args()
    output = local_output(args.output)
    output.mkdir(parents=True, exist_ok=False)
    checks = [('smoke', ['scripts/smoke_test.py']),
              ('tests', ['-m', 'unittest', 'discover', '-s', 'tests', '-v']),
              ('historical-release', ['scripts/verify_pilot_release.py']),
              ('pip-check', ['-m', 'pip', 'check']),
              ('demo', ['scripts/validate_week2.py', '--inputs', 'data/manifests/week2/demo.preparatory.json',
                        '--output', str(output / 'demo')])]
    paired = []
    for kind, manifest in [('curated', args.curated_inputs), ('supplemental', args.supplemental_inputs)]:
        if manifest:
            for number in (1, 2):
                name = f'{kind}-{number}'
                checks.append((name, ['scripts/validate_week2.py', '--inputs', str(manifest.resolve()),
                                      '--output', str(output / name)]))
            paired.append(kind)
    report = run_metadata('week2_technical_handoff_validation', {
        'scientific_execution_requested': False,
        'curated_inputs': str(args.curated_inputs) if args.curated_inputs else None,
        'supplemental_inputs': str(args.supplemental_inputs) if args.supplemental_inputs else None})
    report.update(evidence_kind='technical_validation_only', biological_inference_executed=False,
                  scientific_scores_produced=False, checks=[])
    report['validation_driver_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for name, arguments in checks:
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.perf_counter()
        with (output / (name + '.txt')).open('w', encoding='utf-8') as log:
            result = subprocess.run([sys.executable, *arguments], cwd=ROOT,
                                    stdout=log, stderr=subprocess.STDOUT)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        report['checks'].append({'name': name, 'arguments': arguments, 'exit_code': result.returncode,
            'wall_seconds': time.perf_counter() - started,
            'child_cpu_seconds': after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime})
    report['determinism_checks'] = {}
    for kind in paired:
        paths = [output / f'{kind}-{n}' / 'report.json' for n in (1, 2)]
        values = [read_json(p) if p.exists() else {} for p in paths]
        report['determinism_checks'][kind] = bool(values[0] and values[0] == values[1]
            and values[0].get('blocked_records') == 0
            and values[0].get('scientific_scores_produced') is False)
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    report.update(child_peak_rss_bytes=peak * 1024 if platform.system() == 'Linux' else None,
                  cost_scope='Technical subprocess checks; OS child high-water RSS, not additive/per-model RSS. '
                             'Excludes acquisition, prior model validation, inspection and tool overhead.',
                  gpu_seconds=None, monetary_cost=None)
    passed = (all(c['exit_code'] == 0 for c in report['checks'])
              and all(report['determinism_checks'].values()))
    report['status'] = 'technical_checks_passed' if passed else 'technical_checks_failed'
    write_json(output / 'verification.json', report)
    print(report['status'], output)
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
