"""Release gate: build the wheel, install it in a fresh environment outside the
checkout and exercise the installed CLI and API there.

    python scripts/verify_wheel_install.py [--keep DIRECTORY]

Needs ordinary package-index access to resolve the build backend and Biopython.
Preprocessing itself is run with sockets disabled by the copied test-suite.
Nothing from the checkout is on the import path of the fresh environment: the
tests and fictional fixtures are copied to a temporary directory first.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import venv
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def run(command, **kwargs):
    done = subprocess.run([str(c) for c in command], capture_output=True, text=True, **kwargs)
    if done.returncode != 0:
        sys.stderr.write(done.stdout + done.stderr)
        raise SystemExit(f'FAILED: {" ".join(str(c) for c in command)}')
    return done


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--keep', type=Path, help='keep the working directory here instead of deleting it')
    args = parser.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix='iri-wheel-check-'))
    if ROOT in work.resolve().parents:
        raise SystemExit('Temporary directory resolves inside the checkout; set TMPDIR elsewhere')
    evidence = {'python': sys.version.split()[0], 'workdir_outside_checkout': True}
    try:
        run([sys.executable, '-m', 'pip', 'wheel', '--no-deps', '--wheel-dir', work / 'dist', ROOT])
        wheel, = (work / 'dist').glob('iri_constructs-*.whl')
        names = zipfile.ZipFile(wheel).namelist()
        stray = [n for n in names if not n.startswith(('iri_constructs/', 'iri_constructs-'))]
        if stray:
            raise SystemExit(f'Wheel contains files outside the package: {stray[:5]}')
        evidence.update(wheel=wheel.name, wheel_files=len(names))
        venv.EnvBuilder(with_pip=True, clear=True).create(work / 'env')
        python = work / 'env' / ('Scripts' if os.name == 'nt' else 'bin') / 'python'
        command = python.parent / 'iri-constructs'
        run([python, '-m', 'pip', 'install', '--quiet', wheel])
        installed = json.loads(run([python, '-m', 'pip', 'list', '--format', 'json']).stdout)
        packages = {p['name'].lower() for p in installed}
        if packages & {'torch', 'transformers', 'esm', 'huggingface-hub'}:
            raise SystemExit('Model packages were installed by the core wheel')
        evidence['installed_packages'] = sorted(f"{p['name']}=={p['version']}" for p in installed if p['name'].lower() != 'pip')
        site = work / 'site'
        (site / 'tests').mkdir(parents=True)
        shutil.copy(ROOT / 'tests' / 'test_preprocess.py', site / 'tests')
        shutil.copytree(ROOT / 'examples' / 'synthetic', site / 'examples' / 'synthetic')
        clean = {k: v for k, v in os.environ.items() if k not in ('PYTHONPATH', 'PYTHONHOME')}
        clean['IRI_TEST_INSTALLED'] = '1'
        where = run([python, '-c', 'import iri_constructs, sys; print(iri_constructs.__file__); print(iri_constructs.__version__)'],
                    cwd=site, env=clean).stdout.split()
        if str(ROOT) in where[0]:
            raise SystemExit('Fresh environment imported the checkout, not the installed wheel')
        evidence.update(imported_from_site_packages=True, version=where[1])
        evidence['cli_version'] = run([command, '--version'], cwd=site, env=clean).stdout.strip()
        run([command, 'preprocess', '--help'], cwd=site, env=clean)
        single = subprocess.run([str(command), 'preprocess', 'examples/synthetic/single_cds.gb', '--out', 'out-single'],
                                cwd=site, env=clean, capture_output=True, text=True)
        batch = subprocess.run([str(command), 'preprocess', 'examples/synthetic', '--out', 'out-batch', '--max-model-residues', '1022'],
                               cwd=site, env=clean, capture_output=True, text=True)
        if (single.returncode, batch.returncode) != (0, 2):
            sys.stderr.write(single.stderr + batch.stderr)
            raise SystemExit(f'Unexpected CLI exit codes {(single.returncode, batch.returncode)}; expected (0, 2)')
        report = json.loads((site / 'out-batch' / 'report.json').read_text())
        evidence['cli_single_exit'], evidence['cli_batch_exit'] = single.returncode, batch.returncode
        evidence['batch_summary'] = {k: report['summary'][k] for k in ('records_ok', 'records_failed', 'proteins_exported', 'cds_withheld')}
        tests = subprocess.run([str(python), '-m', 'unittest', 'discover', '-s', 'tests', '-v'], cwd=site, env=clean,
                               capture_output=True, text=True)
        tail = tests.stderr.strip().splitlines()[-3:]
        evidence['installed_test_result'] = ' '.join(tail)
        if tests.returncode != 0:
            sys.stderr.write(tests.stderr)
            raise SystemExit('Installed-package tests failed')
        print(json.dumps(evidence, indent=2, sort_keys=True))
        print('Wheel install verification passed.')
    finally:
        if args.keep:
            shutil.copytree(work, args.keep, dirs_exist_ok=True, ignore=shutil.ignore_patterns('env'))
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
