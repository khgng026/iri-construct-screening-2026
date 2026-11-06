#!/usr/bin/env python3
"""Optional guarded three-query fresh replay of context-profile scalar statistics.

Static authoring does not grant execution. Explicit queue admission and the
prepared guard's recorded process-group ownership are mandatory. No production
profile arithmetic is imported. Only the already approved model adapter is reused.
"""
from __future__ import annotations
import argparse
from collections import Counter
import datetime as dt
import importlib.metadata
import math
import os
from pathlib import Path
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from session_context_profile_audit import (read, require, completed_job, arrays_for, direct_bins,
    independent_mean, weighted_mean, direct_distribution, max_abs, check_close, precision)

GUARD_SHA = 'a5435e14ce302c8d7a78844baed1e015ce526931965a33f960552716d439673d'


def selected_queries(queries):
    return [min((q for q in queries if q['component_count'] == count), key=lambda q: q['query_id']) for count in (2, 3, 4)]


def independent_residue_changes(contextual, isolated):
    """Coordinatewise compensated scalar sums, independent of einsum/norm."""
    import numpy as np
    require(contextual.ndim == 2 and contextual.shape == isolated.shape and np.isfinite(contextual).all() and
            np.isfinite(isolated).all() and min(contextual.shape) > 0, 'Invalid paired residue states')
    cosine, relative = [], []
    for left, right in zip(contextual, isolated):
        a, b = [float(v) for v in left], [float(v) for v in right]
        a2 = math.fsum(v * v for v in a); b2 = math.fsum(v * v for v in b)
        require(a2 > 0 and b2 > 0, 'Undefined zero-norm residue')
        value = math.fsum(x * y for x, y in zip(a, b)) / (math.sqrt(a2) * math.sqrt(b2))
        require(abs(value) <= 1 + 1e-12, 'Cosine arithmetic outside roundoff range')
        cosine.append(min(1., max(-1., value)))
        relative.append(math.sqrt(math.fsum((x-y) ** 2 for x, y in zip(a, b))) / math.sqrt(b2))
    return cosine, relative


def guard_admission():
    path = os.environ.get('IRI_GUARD_OWNERSHIP_PATH'); token = os.environ.get('IRI_GUARD_OWNERSHIP_TOKEN')
    require(path and token, 'Fresh replay requires explicit lead execution under the prepared guard wrapper')
    deadline = time.monotonic() + 5
    while True:
        try:
            record = read(path)
        except FileNotFoundError:
            record = None
        if record is not None:
            require(record['token'] == token and record['guard_source_sha256'] == GUARD_SHA, 'Prepared guard identity differs')
            if any(entry['pgid'] == os.getpgrp() and entry['pid'] == entry['pgid'] and entry['start_new_session']
                   for entry in record['child_groups']):
                break
        # The child can start before its supervisor finishes the atomic
        # ownership publication. Await that record briefly before loading.
        require(time.monotonic() < deadline, 'Current process group is not recorded by the unchanged prepared guard')
        time.sleep(.02)
    return {'ownership_path': path, 'guard_source_sha256': GUARD_SHA, 'owned_pgid': os.getpgrp()}


class RawForwardCheck:
    """Independent observation of actual FP32 CPU parameters and raw outputs."""
    def __init__(self, encoder):
        self.encoder = encoder; self.observations = []
        torch = encoder.torch; counts = Counter()
        require(encoder.config['family'] == 'esm2' and encoder.config['device'] == 'cpu' and
                encoder.config['dtype'] == 'float32' and encoder.layer == encoder.num_layers, 'Replay is restricted to final-layer FP32 CPU ESM2')
        for name, tensor in encoder.model.named_parameters():
            if tensor.is_floating_point():
                require(tensor.dtype == torch.float32 and tensor.device.type == 'cpu', 'Actual model parameter precision differs')
                counts[str(tensor.dtype) + '@' + str(tensor.device)] += tensor.numel()
        require(counts, 'No floating model parameters')
        self.parameter_counts = dict(counts)
        def observe(module, inputs, output):
            states = output.last_hidden_state
            require(states.dtype == torch.float32 and states.device.type == 'cpu' and
                    not states.requires_grad and not module.training, 'Raw forward dtype/device/training differs')
            self.observations.append({'dtype': str(states.dtype), 'device': str(states.device), 'device_type': states.device.type,
                                     'shape': list(states.shape), 'requires_grad': states.requires_grad, 'model_training': module.training,
                                     'autocast_enabled': torch.is_autocast_enabled('cpu'), 'adapter_return_dtype': 'float32'})
        self.handle = encoder.model.register_forward_hook(observe)

    def encode(self, sequence):
        import numpy as np
        self.observations.clear()
        require(not self.encoder.torch.is_autocast_enabled('cpu'), 'Autocast enabled during replay')
        states = self.encoder.encode(sequence)
        require(len(self.observations) == 1 and states.dtype == np.float32 and states.ndim == 2 and
                states.shape[0] == len(sequence) and np.isfinite(states).all(), 'Invalid fresh adapter output')
        precision(self.observations[0], len(sequence), states.shape[1])
        return states, self.observations[0]

    def close(self): self.handle.remove()


def compare_replay(query, states, standalone, saved, saved_arrays, core_arrays):
    import numpy as np
    errors = {}; rows = []; contextual_means = []; isolated_means = []; components = []
    def compare(name, actual, expected):
        error = max_abs(actual, expected)
        require(error <= 1e-6, name + ' differs beyond frozen tolerance')
        errors[name] = max(errors.get(name, 0.), error)
        return error
    compare('whole_vs_parent', independent_mean(states), core_arrays['whole'])
    require(saved['query'] == {k: v for k, v in query.items() if k != 'sequence'}, 'Replay query metadata differs')
    for index, span in enumerate(query['spans']):
        a = span['accession']; contextual = states[span['start0']:span['end0']]; isolated = standalone[a]
        require(contextual.shape == isolated.shape, 'Component residue shape differs')
        cosine, relative = independent_residue_changes(contextual, isolated)
        bins = direct_bins(len(contextual), index > 0, index + 1 < query['component_count'])
        indices = []
        for bin_index, support in bins.items():
            row_index = len(rows); indices.append(row_index)
            expected = saved['bins'][row_index]
            require(expected['row_index'] == row_index and expected['accession'] == a and expected['component_index'] == index and
                    expected['bin_index'] == bin_index and expected['residue_count'] == len(support), 'Replay bin membership differs')
            context_mean = independent_mean(contextual, support).astype(np.float32)
            isolated_mean = independent_mean(isolated, support).astype(np.float32)
            compare('contextual_bin_mean', context_mean, saved_arrays['contextual_bin_means'][row_index])
            compare('isolated_bin_mean', isolated_mean, saved_arrays['isolated_bin_means'][row_index])
            row = {'row_index': row_index, 'accession': a, 'component_index': index, 'bin_index': bin_index,
                   'residue_count': len(support), 'cosine': direct_distribution(cosine[i] for i in support),
                   'relative_l2': direct_distribution(relative[i] for i in support)}
            for metric in ('cosine', 'relative_l2'):
                check_close(expected[metric], row[metric], 'Fresh residue ' + metric)
                errors[metric + '_mean_max_abs'] = max(errors.get(metric + '_mean_max_abs', 0.), abs(expected[metric]['mean'] - row[metric]['mean']))
                errors[metric + '_quantile_max_abs'] = max(errors.get(metric + '_quantile_max_abs', 0.),
                    max(abs(x-y) for x, y in zip(expected[metric]['quantiles'], row[metric]['quantiles'])))
            rows.append(row); contextual_means.append(context_mean); isolated_means.append(isolated_mean)
        length = len(contextual)
        contextual_reconstructed = weighted_mean([contextual_means[i] for i in indices], [rows[i]['residue_count'] for i in indices])
        isolated_reconstructed = weighted_mean([isolated_means[i] for i in indices], [rows[i]['residue_count'] for i in indices])
        compare('contextual_component_vs_parent', independent_mean(contextual), core_arrays['contextual_' + a])
        compare('isolated_component_vs_parent', independent_mean(isolated), core_arrays['isolated_' + a])
        compare('contextual_component_bin_reconstruction', contextual_reconstructed, core_arrays['contextual_' + a])
        compare('isolated_component_bin_reconstruction', isolated_reconstructed, core_arrays['isolated_' + a])
        components.append((length, contextual_reconstructed))
    require(len(rows) == len(saved['bins']), 'Replay bin coverage differs')
    compare('whole_bin_reconstruction', weighted_mean([v for _, v in components], [n for n, _ in components]), core_arrays['whole'])
    return {'status': 'passed', 'query_id': query['query_id'], 'query_sequence_sha256': digest(query['sequence'].encode()),
            'component_count': query['component_count'], 'length': len(states), 'bin_count': len(rows), 'bins': rows,
            'maximum_absolute_errors': errors}, {'contextual_bin_means': np.stack(contextual_means), 'isolated_bin_means': np.stack(isolated_means)}


def run(args):
    require(args.execute_fresh_replay, 'Fresh inference was not explicitly selected')
    guard = guard_admission(); started = time.monotonic()
    require(file_hash(__file__) == args.expected_code_sha256, 'Replay code differs from explicit invocation freeze')
    session = args.session.resolve(); profile_path = args.profile_job.resolve()
    design_path = session / 'inputs/context_profile_replay_audit_design_v1.json'
    require(file_hash(design_path) == args.expected_design_sha256, 'Replay design differs from explicit invocation freeze')
    design = read(design_path)
    require(design['version'] == 'context-profile-fresh-replay-independent-v1' and design['status'] == 'frozen_before_any_fresh_replay' and
            design['vector_absolute_tolerance'] == 1e-6 and design['maximum_seconds'] == 600, 'Replay design changed')
    require(file_hash(ROOT / 'scripts/session_context_profile_audit.py') == design['stored_auditor_sha256'], 'Independent arithmetic helpers changed')
    saved_audit = read(args.saved_audit)
    require(saved_audit['status'] == 'passed' and saved_audit['auditor_sha256'] == design['stored_auditor_sha256'] and
            saved_audit['job'] == str(profile_path) and saved_audit['query_count'] == 1146, 'Completed all-query stored audit required')
    from iri_constructs.session_benchmark import load_suite
    suite = load_suite(ROOT, session / 'inputs'); chosen = selected_queries(suite['queries'])
    require(suite['verification']['challenge_manifest_sha256'] == design['challenge_sha256'] and
            [q['query_id'] for q in chosen] == design['query_ids'] and [q['length_aa'] for q in chosen] == design['query_lengths'] and
            [q['component_count'] for q in chosen] == design['component_counts'], 'Frozen lexicographic replay selection differs')
    profile_identity = read(profile_path / 'identity.json')
    require(profile_path.is_relative_to(session / 'jobs') and profile_identity['mode'] == 'context_profile_v1' and
            profile_identity['subset'] == 'full', 'Wrong complete profile source')
    profile_expected = {'query-' + q for q in profile_identity['query_ids']} | {'standalone-' + a for a in profile_identity['targets']} | {'summary'}
    profile = completed_job(profile_path, profile_expected, saved_audit['profile_completion_sha256'])
    require(profile['identity_hash'] == saved_audit['profile_identity_sha256'], 'Saved profile audit identity differs')
    parent_path = Path(profile_identity['parent_job']).resolve(); parent_identity = read(parent_path / 'identity.json')
    parent_expected = {'query-' + q for q in parent_identity['query_ids']} | {'catalogue-' + a for a in suite['sequences']} | {'summary'}
    parent = completed_job(parent_path, parent_expected, profile_identity['parent_completion_sha256'])
    require(parent['identity_hash'] == profile_identity['parent_identity_hash'] == saved_audit['parent_identity_sha256'] and
            parent_identity['method'] == profile_identity['method'] and parent_identity['config'] == profile_identity['config'], 'Parent profile/model binding differs')
    require(parent_identity['config']['model_id'] in design['parent_models_allowed'], 'Unplanned replay representation')
    for name, version in parent_identity['environment'].items():
        require(importlib.metadata.version(name) == version, 'Parent runtime package differs: ' + name)
    for name in ('src/iri_constructs/session_models.py', 'src/iri_constructs/session_benchmark.py', 'src/iri_constructs/session_store.py'):
        require(parent_identity['code_hashes'][name] == file_hash(ROOT / name), 'Parent representation dependency changed')
    code = {name: file_hash(ROOT / name) for name in ('scripts/session_context_profile_replay_audit.py',
            'scripts/session_context_profile_audit.py', 'src/iri_constructs/session_models.py',
            'src/iri_constructs/session_benchmark.py', 'src/iri_constructs/session_store.py')}
    identity = {'mode': 'independent_context_profile_fresh_replay_v1', 'design_sha256': file_hash(design_path),
        'profile_identity_sha256': profile['identity_hash'], 'profile_completion_sha256': profile['completion_sha256'],
        'parent_identity_sha256': parent['identity_hash'], 'parent_completion_sha256': parent['completion_sha256'],
        'saved_audit_sha256': file_hash(args.saved_audit), 'method': profile_identity['method'], 'config': profile_identity['config'],
        'query_ids': design['query_ids'], 'query_sequence_sha256': {q['query_id']: digest(q['sequence'].encode()) for q in chosen},
        'code_hashes': code, 'information_scope': 'Three fixed fresh joined-query replays; independent residue arithmetic; saved approved standalone matrices'}
    store = UnitStore(session / 'jobs' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    deadline = min(dt.datetime.fromisoformat(read(session / 'launch.json')['hard_stop_utc']).timestamp(), time.time() + 600)
    def bounded(): require(time.time() < deadline, 'Fresh replay session/budget deadline reached')
    from iri_constructs.session_models import load_encoder
    encoder = observer = None; expected = []; records = []
    try:
        bounded(); encoder = load_encoder(parent_identity['config'])
        require({k: v for k, v in encoder.metadata.items() if k != 'load_seconds'} == profile_identity['method'], 'Loaded exact model identity differs')
        observer = RawForwardCheck(encoder)
        standalone = {}
        for accession in sorted({a for q in chosen for a in q['components']}):
            source = read(profile['units']['standalone-' + accession] / 'record.json')
            require(source['sequence_sha256'] == digest(suite['sequences'][accession].encode()) and source['length'] == len(suite['sequences'][accession]),
                    'Standalone replay source differs')
            standalone[accession] = arrays_for(profile, 'standalone-' + accession)['residues']
        store.event(event='job_started', state='running', guard=guard, query_count=3)
        for query in chosen:
            bounded(); key = 'query-' + query['query_id']; expected.append(key)
            old = store.read(key, arrays=True)
            if old is not None:
                record, replay_arrays = old
                require(record['status'] == 'passed' and record['query_sequence_sha256'] == identity['query_sequence_sha256'][query['query_id']],
                        'Saved replay unit identity/status differs')
                saved = read(profile['units'][key] / 'record.json'); saved_arrays = arrays_for(profile, key)
                require(record['bin_count'] == len(saved['bins']) and len(record['bins']) == len(saved['bins']), 'Saved replay bin coverage differs')
                for actual, previous in zip(record['bins'], saved['bins']):
                    require(all(actual[name] == previous[name] for name in ('row_index', 'accession', 'component_index', 'bin_index', 'residue_count')),
                            'Saved replay geometry differs')
                    for metric in ('cosine', 'relative_l2'): check_close(actual[metric], previous[metric], 'Saved replay scalar')
                for kind in ('contextual_bin_means', 'isolated_bin_means'):
                    require(max_abs(replay_arrays[kind], saved_arrays[kind]) <= 1e-6, 'Saved replay means differ')
            else:
                tick = time.monotonic(); states, raw = observer.encode(query['sequence'])
                record, replay_arrays = compare_replay(query, states, standalone, read(profile['units'][key] / 'record.json'),
                                                      arrays_for(profile, key), arrays_for(parent, key))
                del states
                record.update(precision=raw, parameter_dtype_counts=observer.parameter_counts, seconds=time.monotonic() - tick,
                              source_profile_unit_descriptor_sha256=file_hash(profile['units'][key] / 'complete.json'))
                store.publish(key, record, replay_arrays)
            records.append(record)
        summary = {'status': 'passed', 'query_count': 3, 'query_ids': design['query_ids'],
            'total_residue_observations_replayed': sum(r['length'] for r in records),
            'total_component_bins_replayed': sum(r['bin_count'] for r in records),
            'maximum_absolute_errors': {name: max(r['maximum_absolute_errors'].get(name, 0.) for r in records)
                                        for name in sorted(set().union(*(set(r['maximum_absolute_errors']) for r in records)))},
            'scope': design['interpretation']}
        old = store.read('summary')
        if old is None: store.publish('summary', summary)
        else: require(old == summary, 'Recomputed replay summary differs')
        store.complete(expected + ['summary'], {'elapsed_seconds': time.monotonic() - started, 'guard': guard})
        import json
        print(json.dumps({'status': 'passed', 'job': str(store.root), **summary}), flush=True)
    except BaseException as exc:
        store.event(event='job_failed', state='failed', exception_type=type(exc).__name__, exception=str(exc))
        raise
    finally:
        if observer is not None: observer.close()
        if encoder is not None: encoder.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute-fresh-replay', action='store_true', required=True)
    p.add_argument('--session', type=Path, required=True); p.add_argument('--profile-job', type=Path, required=True)
    p.add_argument('--saved-audit', type=Path, required=True); p.add_argument('--name', required=True)
    p.add_argument('--expected-code-sha256', required=True); p.add_argument('--expected-design-sha256', required=True)
    def interrupted(signum, frame): raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted); signal.signal(signal.SIGINT, interrupted)
    run(p.parse_args())
