#!/usr/bin/env python3
"""Independent saved-array audit of profile bins and matched contrasts.

Does not import production profile helpers, model adapters or inference code.
Discarded joined residue tensors prevent a fresh scalar cosine/L2 audit; that
limit is explicitly retained in the report, even when all available checks pass.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import atomic_json, canonical, digest, file_hash

PROBABILITIES = (0., .05, .25, .5, .75, .95, 1.)
BIN_LABELS = {-1: 'clean_no_artificial_join', 0: '0-7', 1: '8-31', 2: '32-63', 3: '64-127', 4: '128+'}


def require(condition, message):
    if not condition: raise ValueError(message)


def read(path): return json.loads(Path(path).read_text())


def check_close(actual, expected, label, absolute=1e-10, relative=1e-12):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and set(actual) == set(expected), label + ': dictionary keys differ')
        for key, value in expected.items(): check_close(actual[key], value, label + '/' + key, absolute, relative)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), label + ': list length differs')
        for i, value in enumerate(expected): check_close(actual[i], value, label + '/' + str(i), absolute, relative)
    elif isinstance(expected, float):
        require(isinstance(actual, (int, float)) and math.isfinite(actual) and math.isfinite(expected) and
                abs(actual - expected) <= absolute + relative * abs(expected), label + ': numeric value differs')
    else:
        require(type(actual) is type(expected) and actual == expected, label + ': exact value differs')


def direct_bins(length, left, right):
    require(type(length) is int and length > 0 and type(left) is bool and type(right) is bool, 'Invalid bin geometry')
    result = defaultdict(list)
    for offset in range(length):
        if not left and not right:
            key = -1
        else:
            distances = ([offset] if left else []) + ([length - offset - 1] if right else [])
            distance = min(distances)
            key = sum(distance >= boundary for boundary in (8, 32, 64, 128))
        result[key].append(offset)
    return dict(sorted(result.items()))


def independent_mean(matrix, indices=None):
    import numpy as np
    require(matrix.ndim == 2 and matrix.shape[0] and matrix.shape[1] and np.isfinite(matrix).all(), 'Invalid mean source')
    indices = list(range(len(matrix))) if indices is None else indices
    require(indices, 'Empty mean support')
    return np.asarray([math.fsum(float(matrix[i, j]) for i in indices) / len(indices)
                       for j in range(matrix.shape[1])], dtype=np.float64)


def weighted_mean(vectors, weights):
    import numpy as np
    require(len(vectors) == len(weights) and len(vectors) > 0 and all(w > 0 for w in weights), 'Invalid reconstruction weights')
    require(len({len(v) for v in vectors}) == 1, 'Reconstruction width differs')
    denominator = sum(weights)
    return np.asarray([math.fsum(w * float(v[j]) for w, v in zip(weights, vectors)) / denominator
                       for j in range(len(vectors[0]))], dtype=np.float64)


def max_abs(left, right):
    import numpy as np
    require(left.shape == right.shape and left.size and np.isfinite(left).all() and np.isfinite(right).all(), 'Invalid comparison vectors')
    return max(abs(float(a) - float(b)) for a, b in zip(left.flat, right.flat))


def direct_distribution(values):
    values = sorted(float(v) for v in values)
    require(values and all(math.isfinite(v) for v in values), 'Invalid scalar distribution')
    quantiles = []
    for p in PROBABILITIES:
        position = (len(values) - 1) * p
        lower, upper = math.floor(position), math.ceil(position)
        fraction = position - lower
        quantiles.append(values[lower] + fraction * (values[upper] - values[lower]))
    total = math.fsum(values)
    return {'count': len(values), 'sum': total, 'mean': total / len(values),
            'quantile_probabilities': list(PROBABILITIES), 'quantiles': quantiles}


def independent_summaries(records):
    groups = defaultdict(list); paired = defaultdict(list)
    for record in records:
        query = record['query']; components = defaultdict(dict)
        for row in record['bins']:
            for suite in ('all', query['suite']):
                groups[(record['scope'], suite, row['bin'])].append((query, row))
            require(row['bin_index'] not in components[row['accession']], 'Repeated target bin')
            components[row['accession']][row['bin_index']] = row
        if record['scope'] == 'combined':
            for target, bins in components.items():
                if 0 in bins and 4 in bins:
                    paired[(query['group_id'], target)].append((query['query_id'],
                        bins[4]['cosine']['mean'] - bins[0]['cosine']['mean'],
                        bins[4]['relative_l2']['mean'] - bins[0]['relative_l2']['mean']))
    output = []
    for (scope, suite, bin_name), observations in sorted(groups.items()):
        count = sum(row['residue_count'] for _, row in observations)
        targets = sorted({row['accession'] for _, row in observations})
        summary = {'scope': scope, 'suite': suite, 'bin': bin_name, 'residue_observations': count,
            'component_bin_observations': len(observations), 'queries': len({q['query_id'] for q, _ in observations}),
            'groups': len({q['group_id'] for q, _ in observations}), 'targets': targets, 'target_count': len(targets),
            'group_target_pairs': len({(q['group_id'], r['accession']) for q, r in observations})}
        for metric in ('cosine', 'relative_l2'):
            target_means = defaultdict(list); group_target_means = defaultdict(list)
            for query, row in observations:
                target_means[row['accession']].append(row[metric]['mean'])
                group_target_means[(query['group_id'], row['accession'])].append(row[metric]['mean'])
            def macro(grouped):
                return math.fsum(math.fsum(v) / len(v) for v in grouped.values()) / len(grouped)
            summary[metric] = {'residue_weighted_mean': math.fsum(row[metric]['sum'] for _, row in observations) / count,
                'component_bin_mean_distribution': direct_distribution(row[metric]['mean'] for _, row in observations),
                'target_macro_mean': macro(target_means), 'group_target_macro_mean': macro(group_target_means)}
        output.append(summary)
    matched = []
    for (group, target), observations in sorted(paired.items()):
        ids = [v[0] for v in observations]
        require(len(ids) == len(set(ids)), 'Repeated matched query')
        matched.append({'group_id': group, 'accession': target, 'matched_query_count': len(ids), 'query_ids': sorted(ids),
                        'cosine_far_minus_near': math.fsum(v[1] for v in observations) / len(ids),
                        'relative_l2_far_minus_near': math.fsum(v[2] for v in observations) / len(ids)})
    return output, matched


def completed_job(path, expected, pinned_completion=None):
    path = Path(path).resolve(); identity = read(path / 'identity.json'); identity_hash = digest(canonical(identity))
    candidates = list(path.glob('completion-*.json'))
    if pinned_completion:
        candidates = [p for p in candidates if file_hash(p) == pinned_completion]
        require(len(candidates) == 1, 'Pinned completion is missing/ambiguous')
    else:
        require(candidates, 'No completion descriptor')
        candidates = [max(candidates, key=lambda p: (read(p)['completed_utc'], p.name))]
    completion_path = candidates[0]; completion = read(completion_path)
    require(completion['identity_hash'] == identity_hash and completion['expected_count'] == len(expected) and
            set(completion['units']) == set(expected), 'Completion identity/coverage differs')
    units = {}; unit_hashes = {}
    for key, entry in completion['units'].items():
        unit = (path / entry['path']).resolve()
        require(unit.is_relative_to(path), 'Unit path escaped source job')
        require(file_hash(unit / 'complete.json') == entry['descriptor_sha256'], 'Unit descriptor checksum differs')
        descriptor = read(unit / 'complete.json')
        require(descriptor['key'] == key and descriptor['identity_hash'] == identity_hash and 'record.json' in descriptor['files'],
                'Unit identity differs')
        for name, details in descriptor['files'].items():
            require(Path(name).name == name and name not in ('.', '..'), 'Unsafe unit filename')
            artifact = unit / name
            require(not artifact.is_symlink() and artifact.stat().st_size == details['bytes'] and file_hash(artifact) == details['sha256'],
                    'Unit bytes/checksum differ')
        units[key] = unit; unit_hashes[key] = entry['descriptor_sha256']
    return {'path': path, 'identity': identity, 'identity_hash': identity_hash, 'completion_sha256': file_hash(completion_path),
            'units': units, 'unit_descriptor_sha256': unit_hashes}


def arrays_for(job, key):
    import numpy as np
    unit = job['units'][key]
    require('arrays.npz' in read(unit / 'complete.json')['files'], 'Arrays are not checksum-bound')
    with np.load(unit / 'arrays.npz', allow_pickle=False) as data:
        return {name: data[name].copy() for name in data.files}


def precision(record, length, width):
    require(record['dtype'] == 'torch.float32' and record['device_type'] == 'cpu' and record['device'] == 'cpu' and
            record['shape'] == [1, length + 2, width] and record['requires_grad'] is False and record['model_training'] is False and
            record['autocast_enabled'] is False and record['adapter_return_dtype'] == 'float32', 'Saved raw precision evidence differs')


def audit(args):
    import numpy as np
    from iri_constructs.session_benchmark import load_suite
    started = time.monotonic(); job_path = args.job.resolve(); session = args.session.resolve()
    deadline = started + args.maximum_seconds
    def bounded(): require(time.monotonic() < deadline, 'Independent saved-array audit deadline reached')
    identity = read(job_path / 'identity.json')
    require(identity['mode'] == 'context_profile_v1' and identity['subset'] == 'full', 'Audit expects the completed full profile')
    design_path = session / 'inputs/context_profile_independent_audit_design_v1.json'
    audit_design = read(design_path)
    require(audit_design['version'] == 'context-profile-independent-audit-v1' and
            audit_design['vector_absolute_tolerance'] == 1e-6, 'Independent audit design differs')
    source_design = session / 'inputs/context_profile_design_v1.json'
    require(file_hash(source_design) == identity['design_sha256'] and read(source_design) == identity['design'], 'Profile design bytes differ')
    require(identity['design']['bins_distance_to_nearest_artificial_join_residues'] == [[0, 8], [8, 32], [32, 64], [64, 128], [128, None]] and
            identity['quantiles'] == list(PROBABILITIES) and identity['tolerance'] == 1e-6, 'Frozen bins/quantiles/tolerance changed')
    suite = load_suite(ROOT, session / 'inputs')
    lookup = {q['query_id']: q for q in suite['queries']}; ids = identity['query_ids']; targets = identity['targets']
    require(len(ids) == len(set(ids)) == 1146 and set(ids) == set(lookup) and
            targets == suite['manifest']['targets'] and len(targets) == 16, 'Full query/target coverage differs')
    require(identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'] and
            identity['sequence_sha256'] == {a: digest(suite['sequences'][a].encode()) for a in targets}, 'Canonical profile sources differ')
    expected = {'query-' + q for q in ids} | {'standalone-' + a for a in targets} | {'summary'}
    profile = completed_job(job_path, expected)
    parent_path = Path(identity['parent_job']).resolve()
    require(parent_path.is_relative_to(session / 'jobs') and job_path.is_relative_to(session / 'jobs'), 'Source outside owned session jobs')
    parent_identity = read(parent_path / 'identity.json')
    parent_ids = parent_identity['query_ids']
    parent_expected = {'query-' + q for q in parent_ids} | {'catalogue-' + a for a in suite['sequences']} | {'summary'}
    parent = completed_job(parent_path, parent_expected, identity['parent_completion_sha256'])
    require(parent['identity_hash'] == identity['parent_identity_hash'] and parent_identity['mode'] == 'model' and
            len(parent_ids) == len(set(parent_ids)) == 1146 and set(parent_ids) == set(ids) and
            parent_identity['method'] == identity['method'] and parent_identity['config'] == identity['config'], 'Profile/parent representation or query identity differs')
    require(parent_identity['challenge_sha256'] == identity['challenge_sha256'] and
            parent_identity['catalogue_accessions'] == sorted(suite['sequences']) and
            parent_identity['catalogue_sequence_sha256'] == {a: digest(s.encode()) for a, s in suite['sequences'].items()}, 'Parent catalogue source differs')
    errors = defaultdict(float); standalone = {}; standalone_means = {}; width = None
    def vector_check(name, left, right):
        error = max_abs(left, right)
        require(error <= 1e-6, name + ' exceeds frozen absolute tolerance: ' + str(error))
        errors[name] = max(errors[name], error)
        return error
    for accession in targets:
        bounded(); key = 'standalone-' + accession; record = read(profile['units'][key] / 'record.json'); arrays = arrays_for(profile, key)
        sequence = suite['sequences'][accession]
        require(record['accession'] == accession and record['sequence_sha256'] == digest(sequence.encode()) and record['length'] == len(sequence), 'Standalone source identity differs')
        require(set(arrays) == {'residues', 'mean'}, 'Standalone array keys differ')
        states, mean = arrays['residues'], arrays['mean']
        require(states.ndim == 2 and states.shape[0] == len(sequence) and states.dtype == np.float32 and
                mean.dtype == np.float32 and mean.shape == (states.shape[1],) and np.isfinite(states).all() and np.isfinite(mean).all(), 'Invalid standalone arrays')
        width = states.shape[1] if width is None else width
        require(states.shape[1] == width and (np.linalg.norm(states.astype(np.float64), axis=1) > 0).all(), 'Standalone width/norm mismatch')
        reconstructed = independent_mean(states)
        vector_check('standalone_mean', reconstructed, mean)
        core_record = read(parent['units']['catalogue-' + accession] / 'record.json')
        require(core_record['sequence_sha256'] == record['sequence_sha256'] and core_record['length'] == len(sequence), 'Parent standalone metadata differs')
        core = arrays_for(parent, 'catalogue-' + accession)
        error = vector_check('standalone_parent_mean', reconstructed, core['mean'])
        stored_error = max_abs(reconstructed.astype(np.float32), core['mean'])
        check_close(record['parent_mean_max_abs'], stored_error, 'Stored standalone check', 1e-12, 0.)
        precision(record['precision'], len(sequence), width)
        standalone[accession] = states
        standalone_means[accession] = reconstructed
    isolated_cache = {}; scalar_records = []; total_bins = total_components = total_residues = 0
    for counter, qid in enumerate(ids):
        bounded(); query = lookup[qid]; key = 'query-' + qid
        record = read(profile['units'][key] / 'record.json'); arrays = arrays_for(profile, key)
        public = {k: v for k, v in query.items() if k != 'sequence'}
        require(record['query'] == public and record['scope'] == ('clean' if query['component_count'] == 1 else 'combined') and
                record['tolerance'] == 1e-6, 'Query identity/scope/tolerance differs')
        core_record = read(parent['units'][key] / 'record.json'); core = arrays_for(parent, key)
        require(core_record['query'] == public, 'Parent query identity differs')
        rows = record['bins']; require(rows and [r['row_index'] for r in rows] == list(range(len(rows))), 'Bin row indices differ')
        require(set(arrays) == {'contextual_bin_means', 'isolated_bin_means'}, 'Profile array keys differ')
        require(all(v.shape == (len(rows), width) and v.dtype == np.float32 and np.isfinite(v).all() for v in arrays.values()), 'Invalid bin vectors')
        precision(record['precision'], query['length_aa'], width)
        expected_order = []; reconstructed_components = []
        for index, span in enumerate(query['spans']):
            accession = span['accession']; length = span['end0'] - span['start0']; left = index > 0; right = index + 1 < query['component_count']
            bins = direct_bins(length, left, right); selected = [r for r in rows if r['component_index'] == index]
            require([r['bin_index'] for r in selected] == list(bins), 'Missing/extra distance bins')
            for row in selected:
                bin_index = row['bin_index']; support = bins[bin_index]; expected_order.append((index, bin_index))
                require(row['accession'] == accession and row['component_length'] == length and row['residue_count'] == len(support) and
                        row['left_artificial_join'] is left and row['right_artificial_join'] is right and row['bin'] == BIN_LABELS[bin_index],
                        'Bin membership/count/label differs')
                cache_key = (accession, left, right, bin_index)
                if cache_key not in isolated_cache:
                    isolated_cache[cache_key] = independent_mean(standalone[accession], support)
                vector_check('isolated_bin_mean', isolated_cache[cache_key], arrays['isolated_bin_means'][row['row_index']])
                for metric in ('cosine', 'relative_l2'):
                    stats = row[metric]
                    require(set(stats) == {'count', 'sum', 'mean', 'quantile_probabilities', 'quantiles'} and stats['count'] == len(support) and
                            stats['quantile_probabilities'] == list(PROBABILITIES) and len(stats['quantiles']) == len(PROBABILITIES), 'Scalar support/quantile schema differs')
                    values = [stats['sum'], stats['mean'], *stats['quantiles']]
                    require(all(math.isfinite(v) for v in values) and stats['quantiles'] == sorted(stats['quantiles']), 'Invalid scalar values/quantiles')
                    check_close(stats['mean'], stats['sum'] / len(support), 'Scalar sum/count', 1e-12, 0.)
                    lo, hi = stats['quantiles'][0], stats['quantiles'][-1]
                    require(lo - 1e-12 <= stats['mean'] <= hi + 1e-12 and
                            (lo >= -1 and hi <= 1 if metric == 'cosine' else lo >= 0), 'Scalar statistic range differs')
            require(sum(r['residue_count'] for r in selected) == length, 'Bins do not partition component')
            reconstructed = {}
            for kind in ('contextual', 'isolated'):
                vectors = [arrays[kind + '_bin_means'][r['row_index']] for r in selected]
                reconstructed[kind] = weighted_mean(vectors, [r['residue_count'] for r in selected])
                error = vector_check(kind + '_component_parent', reconstructed[kind], core[kind + '_' + accession])
                check_close(record['checks'][accession][kind + '_bin_parent_max_abs'], error, 'Stored component reconstruction', 1e-12, 0.)
                for name in (kind + '_parent_max_abs', kind + '_bin_fresh_mean_max_abs'):
                    require(math.isfinite(record['checks'][accession][name]) and 0 <= record['checks'][accession][name] <= 1e-6,
                            'Invalid saved nonreconstructible fresh-check evidence')
            vector_check('isolated_component_residue_mean', reconstructed['isolated'], standalone_means[accession])
            reconstructed_components.append((length, reconstructed['contextual']))
        require([(r['component_index'], r['bin_index']) for r in rows] == expected_order, 'Unexpected bin row order')
        whole = weighted_mean([v for _, v in reconstructed_components], [n for n, _ in reconstructed_components])
        error = vector_check('whole_parent', whole, core['whole'])
        check_close(record['checks']['binned_whole_parent_max_abs'], error, 'Stored whole reconstruction', 1e-12, 0.)
        if query['component_count'] == 1:
            require(math.isfinite(record['checks']['clean_residue_repeat_max_abs']) and
                    0 <= record['checks']['clean_residue_repeat_max_abs'] <= 1e-6, 'Invalid clean fresh-residue repeat evidence')
            vector_check('clean_bin_mean_repeat', arrays['contextual_bin_means'], arrays['isolated_bin_means'])
        total_bins += len(rows); total_components += query['component_count']; total_residues += query['length_aa']
        scalar_records.append({'query': public, 'scope': record['scope'], 'bins': rows})
        if (counter + 1) % 200 == 0: print(json.dumps({'audited_queries': counter + 1, 'expected_queries': len(ids)}), flush=True)
    bounded(); expected_bins, expected_matched = independent_summaries(scalar_records)
    summary = read(profile['units']['summary'] / 'record.json')
    require(summary['query_count'] == 1146 and summary['standalone_target_count'] == 16, 'Summary population differs')
    check_close(summary['summary']['by_scope_suite_bin'], expected_bins, 'Aggregate bin summaries')
    check_close(summary['summary']['matched_far128plus_vs_near0to7_by_group_target'], expected_matched, 'Matched near/far summaries')
    require(summary['summary']['quantile_scope'] == 'Per-unit quantiles are over residues; aggregate quantiles are over component-bin means, never pooled residue quantiles.',
            'Quantile scope claim differs')
    require(total_components == 3536, 'Component observation total differs')
    return {'status': 'passed', 'auditor_sha256': file_hash(__file__), 'audit_design_sha256': file_hash(design_path),
        'job': str(job_path), 'profile_identity_sha256': profile['identity_hash'], 'profile_completion_sha256': profile['completion_sha256'],
        'parent_job': str(parent_path), 'parent_identity_sha256': parent['identity_hash'], 'parent_completion_sha256': parent['completion_sha256'],
        'unit_descriptor_sha256': {'profile': profile['unit_descriptor_sha256'], 'parent': parent['unit_descriptor_sha256']},
        'query_count': len(ids), 'standalone_target_count': len(targets), 'component_count': total_components,
        'component_bin_count': total_bins, 'residue_observation_count': total_residues,
        'aggregate_bin_strata': len(expected_bins), 'matched_group_target_count': len(expected_matched),
        'matched_query_component_observations': sum(v['matched_query_count'] for v in expected_matched),
        'vector_dimension': width, 'maximum_absolute_errors': dict(errors),
        'all_saved_bin_geometry_reconstructions_and_aggregate_contrasts_checked': True,
        'joined_per_residue_cosine_l2_recomputed': False,
        'limitation': audit_design['unrecoverable_from_saved_outputs'],
        'production_profile_helpers_imported': False, 'fresh_model_inference': False,
        'elapsed_seconds': time.monotonic() - started}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--session', type=Path, required=True); p.add_argument('--job', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--maximum-seconds', type=int, default=180)
    args = p.parse_args()
    require(1 <= args.maximum_seconds <= 600 and not args.output.exists(), 'Invalid deadline or existing immutable audit output')
    def expired(signum, frame): raise TimeoutError('Saved-array audit exceeded fixed wall deadline')
    signal.signal(signal.SIGALRM, expired); signal.alarm(args.maximum_seconds)
    try:
        report = audit(args)
    except BaseException as exc:
        report = {'status': 'failed', 'auditor_sha256': file_hash(__file__), 'job': str(args.job.resolve()),
                  'exception_type': type(exc).__name__, 'exception': str(exc), 'fresh_model_inference': False}
        atomic_json(args.output, report)
        raise
    finally:
        signal.alarm(0)
    atomic_json(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k != 'unit_descriptor_sha256'}), flush=True)
