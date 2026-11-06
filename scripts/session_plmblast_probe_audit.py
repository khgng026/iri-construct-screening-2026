#!/usr/bin/env python3
"""Independent read-only audit of the frozen pLM-BLAST resource probe.

Standard library only. No runner/upstream/numerical/model imports, no alignment,
no threshold modification and no admission writing. Retrieval scores are never
used in the forecast or resource decision.
"""
import argparse
import json
import math
from pathlib import Path
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import atomic_json, file_hash, canonical, digest
from session_context_profile_audit import completed_job

DESIGN_SHA = '98528bfad9d28a6b2c9de7a1c8d50c050cae5c9759bbb3dbd104e12fcfcd83b1'


def require(condition, message):
    if not condition: raise ValueError(message)


def read(path): return json.loads(Path(path).read_text())


def near(actual, expected, label):
    require(math.isfinite(actual) and math.isfinite(expected) and abs(actual-expected) <= 1e-9 + abs(expected)*1e-12,
            label + ' arithmetic mismatch')


def unit(job, key): return read(job['units'][key] / 'record.json')


def audit(args):
    started = time.monotonic(); session = args.session.resolve(); job_path = args.job.resolve()
    design_path = session / 'inputs/plmblast_design_v1.json'
    require(file_hash(design_path) == DESIGN_SHA, 'Original frozen v1 design changed')
    design = read(design_path); queries = design['probe_query_ids']; references = design['probe_reference_accessions']
    pairs = [(q, a) for q in queries for a in references]
    def pair_key(q, a): return 'pair-' + q + '--' + a
    expected = {pair_key(q, a) for q, a in pairs} | {'timing-' + pair_key(q, a) for q, a in pairs} | {'evaluation-' + q for q in queries} | {'startup', 'summary'}
    probe = completed_job(job_path, expected); identity = probe['identity']
    require(identity['mode'] == 'plmblast_alignment_v1' and identity['panel'] == 'probe12' and identity['design_sha256'] == DESIGN_SHA and
            identity['query_ids'] == queries and identity['reference_accessions'] == references and set(identity['unit_keys']) == expected,
            'Probe panel/design/coverage differs')
    for path, sha in identity['code_hashes'].items():
        require(file_hash(ROOT / path) == sha, 'Recorded runner dependency changed: ' + path)
    source_manifest = session / 'audit/plmblast_source_acquisition.json'
    require(file_hash(source_manifest) == design['source_manifest_sha256'] == identity['source']['manifest_sha256'], 'Pinned upstream acquisition differs')
    source = read(source_manifest)
    require(source['revision'] == design['source_revision'] == identity['source']['revision'] and source['licence'] == 'MIT', 'Upstream revision/license differs')
    vendor = session / 'vendor' / ('pLM-BLAST-' + source['revision'])
    require(identity['source']['files'] == {r['path']: r['sha256'] for r in source['files']}, 'Upstream source roster differs')
    for record in source['files']:
        path = (vendor / record['path']).resolve()
        require(path.is_relative_to(vendor) and path.stat().st_size == record['bytes'] and file_hash(path) == record['sha256'], 'Upstream source bytes differ')
    challenge_path = session / 'inputs/challenge_manifest.json'
    require(file_hash(challenge_path) == design['challenge_sha256'], 'Frozen challenge changed')
    challenge = read(challenge_path); by_query = {q['query_id']: q for q in challenge['queries']}
    base_path = ROOT / 'data/manifests/t05_expanded_catalogue.json'
    require(file_hash(base_path) == challenge['source_files_sha256']['data/manifests/t05_expanded_catalogue.json'], 'Original reference metadata changed')
    base = read(base_path); lengths = {r['accession']: r['sequence_length_aa'] for r in base['records']}
    require(len(lengths) == 105 and [by_query[q]['length_aa'] for q in queries] == design['probe_query_lengths'] and
            [lengths[a] for a in references] == design['probe_reference_lengths'], 'Probe canonical lengths differ')
    encoding_path = Path(identity['encoding_source']['path']); encoding_identity = read(encoding_path / 'identity.json')
    encoding_expected = {'reference-' + a for a in references} | {'query-' + q for q in queries} | {'summary'}
    encoding = completed_job(encoding_path, encoding_expected, identity['encoding_source']['completion_sha256'])
    require(encoding['identity_hash'] == identity['encoding_source']['identity_sha256'] and encoding_identity['panel'] == 'probe12' and
            encoding_identity['method'] == identity['method'] and encoding_identity['design_sha256'] == DESIGN_SHA and
            encoding_identity['code_hashes'] == identity['code_hashes'], 'Encoding source binding differs')
    model_path = Path(encoding_identity['model_probe']['path']); model_identity = read(model_path / 'identity.json')
    model_probe = completed_job(model_path, set(queries), encoding_identity['model_probe']['completion_sha256'])
    require(model_probe['identity_hash'] == encoding_identity['model_probe']['identity_sha256'] and model_identity['mode'] == 'admission_probe_v1' and
            model_identity['query_ids'] == queries and model_identity['repeats'] == 3 and model_identity['model'] == identity['method'], 'Original model probe binding differs')
    for name in ('model_id', 'revision', 'layer', 'num_layers', 'dtype', 'device'):
        require(identity['method'][name] == design['model'][name], 'Frozen model representation differs')
    require(identity['method']['checkpoint_files']['model.safetensors']['sha256'] == design['model']['checkpoint_sha256'], 'Frozen model checkpoint differs')
    runtime = identity['runtime']; overlay = Path(runtime['overlay']).resolve()
    require(overlay == session / 'vendor/plmblast-overlay-v1' and runtime['torch_threads'] == 4 and runtime['numba_threads'] == 1,
            'Runtime overlay/thread policy differs')
    installation_path = session / 'inputs/plmblast_overlay_installation_v1.json'; installation = read(installation_path)
    require(installation['status'] == 'verified_installed' and Path(installation['overlay']) == overlay, 'Overlay installation record differs')
    entries = list(overlay.rglob('*')); require(not any(p.is_symlink() for p in entries), 'Overlay contains symlink')
    actual_files = {str(p.relative_to(overlay)) for p in entries if p.is_file()}
    recorded_files = set(installation['installed_files'])
    require(not recorded_files - actual_files, 'Recorded overlay file is missing')
    extra_files = []
    for relative in sorted(actual_files - recorded_files):
        path = overlay / relative
        require(path.suffix == '.pyc' and path.parent.name == '__pycache__', 'Unexpected non-cache overlay file')
        source_path = path.parent.parent / (path.name.split('.cpython-')[0] + '.py')
        require(str(source_path.relative_to(overlay)) in recorded_files, 'Extra cache lacks an installed source')
        extra_files.append({'path': relative, 'sha256': file_hash(path), 'bytes': path.stat().st_size,
                            'mtime_ns': path.stat().st_mtime_ns, 'installed_source': str(source_path.relative_to(overlay))})
    for relative, sha in installation['installed_files'].items():
        require(file_hash(overlay / relative) == sha, 'Installed overlay bytes changed')
    for name, record in runtime['packages'].items():
        require(file_hash(record['module_path']) == record['module_sha256'] and
                file_hash(record['distribution_metadata_path']) == record['distribution_metadata_sha256'], 'Recorded runtime module metadata changed')
        if name in design['alignment_environment']:
            require(record['version'] == design['alignment_environment'][name] and Path(record['module_path']).is_relative_to(overlay), 'Pinned runtime package differs')
    binary_files = {str(p.relative_to(overlay)) for p in entries if p.is_file() and p.suffix in ('.so', '.dylib')}
    require(binary_files == set(runtime['overlay_binary_sha256']), 'Runtime binary coverage differs')
    require(all(installation['installed_files'][p] == sha for p, sha in runtime['overlay_binary_sha256'].items()), 'Runtime binary hashes differ from verified install')
    summary = unit(probe, 'summary'); starts = []
    for record in summary['startup_evidence']:
        ready_path = Path(record['ready_path']); ready = record['ready']
        require(file_hash(ready_path) == record['ready_sha256'] and read(ready_path) == ready and ready['runtime'] == runtime, 'Startup runtime evidence changed')
        require(ready['fixtures']['status'] == 'passed' and all(ready['fixtures'][k] is True for k in
            ('identity_pipeline_equivalence', 'native_no_hit', 'cutoff_before_rounding', 'rounding_tie', 'coordinate_validation')),
            'Numeric startup fixture did not pass')
        require(ready['fixtures']['orientation_cosine_max_abs'] <= 2e-6 and
                ready['peak_process_rss_mib'] <= design['limits']['maximum_worker_peak_mib'], 'Startup numeric/resource gate failed')
        starts.append(ready)
    require(starts, 'Missing startup observations')
    observations = []
    for q, a in pairs:
        key = pair_key(q, a); row = unit(probe, key); timing = unit(probe, 'timing-' + key)
        require(row['query_id'] == q and row['accession'] == a and row['status'] in ('eligible_alignment', 'no_eligible_alignment') and
                timing['status'] == 'passed' and timing['attempted'] is True and timing['seconds'] > 0 and math.isfinite(timing['seconds']),
                'Pair/source/resource outcome failed')
        require(row['query_sequence_sha256'] == encoding_identity['sequence_hashes']['query-' + q] and
                row['reference_sequence_sha256'] == encoding_identity['sequence_hashes']['reference-' + a], 'Pair input hashes differ from encoded states')
        near(summary['pair_timings'][key], timing['seconds'], 'Saved pair timing')
        n, m = by_query[q]['length_aa'], lengths[a]; work = n*m*(n+m)
        require(row['peak_process_rss_mib'] <= design['limits']['maximum_worker_peak_mib'], 'Observed pair memory gate failed')
        observations.append({'pair_key': key, 'query_length': n, 'reference_length': m, 'work': work,
            'wall_seconds': timing['seconds'], 'native_seconds': row['native_seconds'],
            'rate_seconds_per_work': timing['seconds']/work, 'peak_process_rss_mib': row['peak_process_rss_mib']})
    require(len(observations) == summary['expected_pair_count'] == 12 and summary['failed_pair_count'] == 0 and
            summary['failed_or_unadmitted_pair_count'] == 0 and summary['resource_gate_failures'] == [] and summary['summaries'] is None,
            'Probe coverage/failure/performance scope differs')
    maximum_rate = max(observations, key=lambda r: r['rate_seconds_per_work'])
    repeat_work = sum(by_query[q]['length_aa']*m*(by_query[q]['length_aa']+m) for q in design['repeat_query_ids'] for m in lengths.values())
    startup = max(r['startup_seconds'] for r in starts)
    factor = design['forecast']['safety_factor']; require(factor == 1.5, 'Frozen safety factor changed')
    alignment_seconds = factor * (startup + maximum_rate['rate_seconds_per_work']*repeat_work)
    model_rows = []
    for q in queries:
        row = unit(model_probe, q); times = row['times_seconds']
        require(row['query_id'] == q and row['length'] == by_query[q]['length_aa'] and len(times) == 3 and
                all(math.isfinite(t) and t > 0 for t in times), 'Model-probe timing coverage differs')
        model_rows.append({'query_id': q, 'length': row['length'], 'times_seconds': times,
                           'max_seconds_per_length_squared': max(times)/row['length']**2,
                           'peak_process_rss_mib': row['peak_process_rss_mib']})
    encoding_rate = max(model_rows, key=lambda r: r['max_seconds_per_length_squared'])
    encoding_work = sum(m*m for m in lengths.values()) + sum(by_query[q]['length_aa']**2 for q in design['repeat_query_ids'])
    encoding_seconds = factor * encoding_rate['max_seconds_per_length_squared'] * encoding_work
    saved = summary['forecast']
    near(saved['alignment_seconds'], alignment_seconds, 'Alignment forecast')
    near(saved['encoding_seconds'], encoding_seconds, 'Encoding forecast')
    near(saved['maximum_observed_pair_seconds_per_work'], maximum_rate['rate_seconds_per_work'], 'Maximum pair rate')
    require(saved['repeat_work'] == repeat_work and saved['safety_factor'] == factor and saved['is_guarantee'] is False, 'Forecast work/interpretation differs')
    within = alignment_seconds <= design['limits']['repeat_alignment_seconds'] and encoding_seconds <= design['limits']['repeat_encoding_seconds']
    require(saved['within_frozen_budget'] is within and not within and summary['resource_recommendation'] == 'not_admitted',
            'Original negative v1 resource decision differs')
    return {'status': 'passed_resource_arithmetic_with_overlay_membership_finding' if extra_files else 'passed_original_negative_gate_confirmed',
        'auditor_sha256': file_hash(__file__), 'design_sha256': DESIGN_SHA,
        'probe_job': str(job_path), 'probe_identity_sha256': probe['identity_hash'], 'probe_completion_sha256': probe['completion_sha256'],
        'encoding_identity_sha256': encoding['identity_hash'], 'encoding_completion_sha256': encoding['completion_sha256'],
        'model_probe_identity_sha256': model_probe['identity_hash'], 'model_probe_completion_sha256': model_probe['completion_sha256'],
        'installation_record_sha256': file_hash(installation_path), 'verified_installed_file_count': len(installation['installed_files']),
        'overlay_exact_membership_matches_installation': not extra_files, 'unrecorded_overlay_cache_files': extra_files,
        'overlay_finding': 'Recorded installed bytes all match, but extra bytecode caches violate exact membership; cache/source semantic equality not checked.' if extra_files else None,
        'verified_runtime_binary_count': len(binary_files), 'pair_observations': observations, 'model_probe_observations': model_rows,
        'largest_normalized_pair_rate_observation': maximum_rate, 'encoding_rate_observation': encoding_rate,
        'startup_seconds': startup, 'startup_peak_rss_mib': max(r['peak_process_rss_mib'] for r in starts),
        'observed_pair_peak_rss_mib': max(r['peak_process_rss_mib'] for r in observations),
        'observed_model_probe_peak_rss_mib': max(r['peak_process_rss_mib'] for r in model_rows),
        'pair_wall_seconds_sum': math.fsum(r['wall_seconds'] for r in observations),
        'pair_wall_seconds_range': [min(r['wall_seconds'] for r in observations), max(r['wall_seconds'] for r in observations)],
        'repeat_work': repeat_work, 'repeat_pair_count': len(design['repeat_query_ids'])*len(lengths), 'encoding_work': encoding_work,
        'recomputed_alignment_seconds': alignment_seconds, 'recomputed_encoding_seconds': encoding_seconds,
        'combined_forecast_seconds': alignment_seconds+encoding_seconds,
        'original_alignment_cap_seconds': design['limits']['repeat_alignment_seconds'], 'original_encoding_cap_seconds': design['limits']['repeat_encoding_seconds'],
        'original_v1_admission': 'not_admitted', 'admission_written': False, 'original_design_or_runner_modified': False,
        'retrieval_scores_used_for_resource_decision': False, 'fresh_alignment_or_model_inference': False,
        'scope': 'Resource/provenance/forecast audit only; observations do not prove full-panel runtime or memory upper bounds.',
        'elapsed_seconds': time.monotonic()-started}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--session', type=Path, required=True); p.add_argument('--job', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); require(not args.output.exists(), 'Preserve the existing immutable audit output')
    def expired(signum, frame): raise TimeoutError('Resource probe audit exceeded60seconds')
    signal.signal(signal.SIGALRM, expired); signal.alarm(60)
    try:
        report = audit(args)
    except BaseException as exc:
        atomic_json(args.output, {'status': 'failed', 'auditor_sha256': file_hash(__file__), 'exception_type': type(exc).__name__, 'exception': str(exc)})
        raise
    finally:
        signal.alarm(0)
    atomic_json(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('pair_observations','model_probe_observations','unrecorded_overlay_cache_files')}), flush=True)
