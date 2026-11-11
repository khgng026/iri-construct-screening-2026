"""Validate pinned local model interfaces on a fictional alphabet fixture only.

This is software validation, not biological inference or family retrieval. No
source sequence files are accepted by this command, and it cannot issue data or
scientific execution clearance.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.metrics import mean_pool
from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device
from iri_constructs.runtime import InputError, Measure, failure, local_output, run_metadata, write_json

# One occurrence of each canonical residue in alphabet order; no biological
# record, annotation, family label or target function is associated with it.
FIXTURE = 'ACDEFGHIKLMNPQRSTVWY'


def representation_contract(rows, length):
    """Check finite per-residue rows and verify the existing mean-pooling rule."""
    if (not isinstance(rows, list) or len(rows) != length or not rows[0]
            or any(not isinstance(row, list) or len(row) != len(rows[0]) for row in rows)
            or not all(math.isfinite(x) for row in rows for x in row)):
        raise InputError('Invalid fixture residue representation')
    pooled = mean_pool(rows)
    expected = [math.fsum(row[i] for row in rows) / length for i in range(len(rows[0]))]
    if pooled != expected:
        raise InputError('Residue mean pooling disagrees with the declared contract')
    # Only fingerprints are persisted. These checks are not retrieval metrics.
    fingerprint = hashlib.sha256(json.dumps(rows, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    return {'residue_count': len(rows), 'embedding_width': len(rows[0]),
            'finite': True, 'mean_pool_verified': True,
            'representation_sha256': fingerprint}


def check_model(name, revision, offline, threads, seed):
    adapter = None
    load = None
    load_cpu_start = None
    result = {'name': name, 'requested_revision': revision, 'status': 'started',
              'biological_inference_executed': False, 'fixture_forward_calls': 0}
    try:
        import torch
        torch.set_num_threads(threads)
        device = prepare_device('cpu', seed)
        cpu_start = load_cpu_start = time.process_time()
        with Measure() as load:
            adapter = LocalAdapter(name, device, revision, offline=offline)
        result['load'] = dict(load.data, process_cpu_seconds=time.process_time() - cpu_start)
        result['model'] = adapter.metadata
        dtypes = sorted({str(p.dtype) for p in adapter.model.parameters() if p.is_floating_point()})
        result['floating_parameter_dtypes'] = dtypes
        if dtypes != ['torch.float32']:
            raise InputError('Model parameters do not satisfy the float32 policy')
        cpu_start = time.process_time()
        with Measure() as first:
            states = adapter.residues(FIXTURE)
        result['fixture_forward_calls'] += 1
        result['first_fixture_forward'] = dict(first.data, process_cpu_seconds=time.process_time() - cpu_start)
        result['contract'] = representation_contract(states, len(FIXTURE))
        cpu_start = time.process_time()
        with Measure() as repeat:
            repeated = adapter.residues(FIXTURE)
        result['fixture_forward_calls'] += 1
        result['repeat_fixture_forward'] = dict(repeat.data, process_cpu_seconds=time.process_time() - cpu_start)
        result['repeat_contract'] = representation_contract(repeated, len(FIXTURE))
        if states != repeated:
            raise InputError('Repeated deterministic fixture forward changed representations')
        result['deterministic_repeat'] = True
        # All models reject ambiguous/cap-exceeding inputs before a forward.
        rejection_checks = {}
        for label, value in [('empty', ''), ('ambiguous', 'X'), ('over_cap', 'A' * 1023)]:
            try:
                adapter.residues(value)
            except InputError:
                rejection_checks[label] = True
            else:
                raise InputError('Adapter accepted an unsupported software fixture')
        result['unsupported_input_rejection'] = rejection_checks
        result['status'] = 'software_fixture_validated'
    except Exception as exc:
        result.update(status='blocked', **failure(exc))
        if load is not None and 'load' not in result:
            result['load'] = dict(load.data, process_cpu_seconds=time.process_time() - load_cpu_start)
    finally:
        if adapter is not None:
            adapter.close()
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--models', default='esm2')
    parser.add_argument('--esm2-revision', default=ESM2_REVISION)
    parser.add_argument('--esmc-revision', default='main')
    parser.add_argument('--allow-download', action='store_true')
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args(argv)
    names = args.models.split(',')
    if len(set(names)) != len(names) or not names or any(n not in ('esm2', 'esmc') for n in names):
        parser.error('--models must be distinct members of esm2,esmc')
    if not 1 <= args.threads <= 3:
        parser.error('--threads must be in 1..3 for this bounded CPU validation')
    try:
        output = local_output(args.output)
        output.mkdir(parents=True, exist_ok=False)
    except (InputError, FileExistsError) as exc:
        parser.error(str(exc))
    report = run_metadata('model_interface_software_validation', {
        'models': names, 'downloads_allowed': args.allow_download,
        'device': 'cpu', 'threads': args.threads, 'seed': 20261105,
        'fixture_kind': 'canonical alphabet in order, no biological source',
        'fixture_length': len(FIXTURE), 'fixture_sha256': hashlib.sha256(FIXTURE.encode()).hexdigest(),
    })
    report.update(evidence_kind='software_fixture_only', biological_inference_executed=False,
                  scientific_scores_produced=False, models=[])
    write_json(output / 'report.json', report)
    for name in names:
        report['active_model'] = name
        write_json(output / 'report.json', report)
        revision = args.esm2_revision if name == 'esm2' else args.esmc_revision
        result = check_model(name, revision, not args.allow_download, args.threads, 20261105)
        report['models'].append(result)
        write_json(output / 'report.json', report)
    report.pop('active_model', None)
    report['status'] = ('software_models_validated' if all(r['status'] == 'software_fixture_validated'
                                                        for r in report['models']) else 'partial_or_blocked')
    report['limitations'].extend([
        'No biological source sequences, labels, retrieval comparisons or E1-E4 outcomes were used.',
        'Load timing includes snapshot/cache resolution, hash verification and CPU encoder loading.',
        'Second fixture pass is a deterministic repeat, not an independent biological observation.',
        'CPU time is process time; RSS is sampled absolute process memory and may miss brief peaks.',
        'Successful software validation does not confer dataset admission or scientific clearance.',
    ])
    write_json(output / 'report.json', report)
    distributions = sorted(f'{d.metadata["Name"]}=={d.version}' for d in importlib.metadata.distributions()
                           if d.metadata.get('Name'))
    (output / 'environment-lock.txt').write_text('\n'.join(distributions) + '\n')
    print(report['status'], output / 'report.json')
    return 0 if report['status'] == 'software_models_validated' else 1


if __name__ == '__main__':
    raise SystemExit(main())
