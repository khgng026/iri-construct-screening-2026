#!/usr/bin/env python3
"""Guarded official-API check of a frozen intermediate-layer probe's residue states.

The public model hidden_states tuple is checked on the shortest frozen clean
query only. Other admission probes cover three longer inputs and repeats. No
whole-layer tensors are saved, and no model generation or training is used.
"""
import argparse
import importlib.metadata
import json
import math
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from iri_constructs.session_benchmark import load_suite
from session_experiment import check_deadline, model_identity
from session_context_profile_audit import completed_job, arrays_for, read, require
from session_context_profile_replay_audit import guard_admission


def validate_probe_record(record, query, config):
    require(record['query_id'] == query['query_id'] and record['length'] == query['length_aa'], 'Probe unit query or length differs')
    times = record['times_seconds']
    require(len(times) == 3 and all(isinstance(t, (int, float)) and math.isfinite(t) and t > 0 for t in times),
            'Probe needs three finite positive timing observations')
    outputs = record['actual_forward_outputs']
    require(len(outputs) == 3 and all(v == {'dtype': 'torch.float32', 'device': 'cpu'} for v in outputs),
            'Probe needs three actual FP32 CPU raw-output observations')
    difference = record['within_config_max_abs']
    require(isinstance(difference, (int, float)) and math.isfinite(difference) and 0 <= difference <= 1e-6 and
            record['actual_parameter_dtypes'] == ['torch.float32'], 'Source intermediate probe failed numerical admission')
    meta = record['encode']
    require(meta['layer'] == config['layer'] and meta['device'] == 'cpu' and meta['dtype'] == 'float32' and
            meta['truncated'] is False and meta['model_training'] is False and meta['inference_mode'] is True and
            meta['residue_count'] == query['length_aa'] and meta['token_count'] == query['length_aa'] + 2,
            'Source intermediate probe mapping/layer/precision metadata differs')


def run(args):
    session = args.session.resolve()
    check_deadline(session)
    require(file_hash(__file__) == args.expected_code_sha256, 'Layer audit implementation changed after admission')
    guard = guard_admission()
    suite = load_suite(ROOT, session / 'inputs')
    design_path = session / 'inputs/intermediate_layers_design_v1.json'
    require(file_hash(design_path) == args.expected_design_sha256, 'Frozen midpoint design changed')
    design = read(design_path)
    require(design['version'] == 'iri-midpoint-layer-sensitivity-v1' and
            design['status'] == 'frozen_before_intermediate_layer_scoring', 'Unexpected layer design')
    probe_path = args.probe_job.resolve()
    require(probe_path.is_relative_to(session / 'jobs'), 'Probe is outside this session job directory')
    probe_identity = read(probe_path / 'identity.json')
    require(probe_path.name.endswith('-' + digest(canonical(probe_identity))[:12]) and
            probe_identity['probe_code_sha256'] == file_hash(ROOT / 'scripts/session_model_probe.py'),
            'Probe directory identity or frozen admission implementation differs')
    require(probe_identity['mode'] == 'admission_probe_v1' and probe_identity['repeats'] == 3 and
            probe_identity['whole_residue_inputs_only'] is True, 'Not a complete standard admission probe')
    config = probe_identity['config']
    require(config['family'] == 'esm2' and config['device'] == 'cpu' and config['dtype'] == 'float32',
            'Only frozen ESM2 CPU FP32 midpoint configurations are admitted')
    matching = [row for row in design['configs'] if row['model'] == config['model_id']]
    require(len(matching) == 1, 'Model is not a unique frozen midpoint candidate')
    entry = matching[0]
    require(read(ROOT / entry['config']) == config and file_hash(ROOT / entry['config']) == entry['config_sha256'],
            'Probe configuration differs from the prospective design')
    require(config['layer'] == entry['layer'] == probe_identity['model']['num_layers'] // 2,
            'Probe layer is not the fixed architecture midpoint')
    clean = sorted([q for q in suite['queries'] if q['suite'] == 'clean'], key=lambda q: (q['length_aa'], q['query_id']))
    longest = max(suite['queries'], key=lambda q: (q['length_aa'], q['query_id']))
    panel = [clean[0], clean[len(clean) // 2], clean[-1], longest]
    require(probe_identity['query_ids'] == [q['query_id'] for q in panel] and
            probe_identity['query_hashes'] == {q['query_id']: q['joined_sequence_sha256'] for q in panel},
            'Probe examples differ from the original admission panel')
    source = completed_job(probe_path, set(probe_identity['query_ids']))
    for query in panel:
        record = read(source['units'][query['query_id']] / 'record.json')
        validate_probe_record(record, query, config)
    query = clean[0]
    saved = arrays_for(source, query['query_id'])['residues']
    capturing_path = Path(importlib.metadata.distribution('transformers').locate_file('transformers/utils/output_capturing.py')).resolve()
    capturing_sha256 = file_hash(capturing_path)
    require(capturing_sha256 == 'a188aeb10cca3a3238d7df1d4e3c08b44dac283ba5b8b7a7cafd2d2bbc532e7c',
            'Official hidden-state tuple capture implementation differs from static review')
    identity = {'mode': 'official_intermediate_layer_probe_audit_v1', 'design_sha256': file_hash(design_path),
                'source_probe': str(probe_path), 'source_identity_hash': source['identity_hash'],
                'source_completion_sha256': source['completion_sha256'], 'query_id': query['query_id'],
                'query_sequence_sha256': query['joined_sequence_sha256'], 'config': config,
                'comparison': 'Official EsmModel(..., output_hidden_states=True).hidden_states[fixed_midpoint][0, 1:-1]',
                'absolute_tolerance': 1e-6,
                'official_tuple_capture': {'path': str(capturing_path), 'sha256': capturing_sha256},
                'code': {name: file_hash(ROOT / name) for name in ('scripts/session_layer_probe_audit.py',
                    'scripts/session_context_profile_audit.py', 'scripts/session_context_profile_replay_audit.py',
                    'scripts/session_model_probe.py', 'scripts/session_experiment.py', 'src/iri_constructs/session_benchmark.py',
                    'src/iri_constructs/session_models.py', 'src/iri_constructs/session_store.py')}}
    store = UnitStore(session / 'jobs' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    old = store.read('audit')
    if old is not None:
        require(old['status'] == 'passed', 'Saved audit did not pass')
        store.complete(['audit'], {'resumed_verified': True, 'guard': guard})
        print(json.dumps({'status': 'passed', 'job': str(store.root), 'resumed': True}), flush=True)
        return
    from iri_constructs.session_models import load_encoder, validate_token_mapping
    import numpy as np
    started = time.monotonic()
    encoder = load_encoder(config)
    try:
        require(model_identity(encoder) == probe_identity['model'], 'Official verification loaded a different model/environment')
        torch = encoder.torch
        require(not torch.is_autocast_enabled('cpu') and not encoder.model.training, 'Autocast or training is enabled')
        require(all(p.dtype == torch.float32 and p.device.type == 'cpu' for p in encoder.model.parameters() if p.is_floating_point()),
                'Actual model parameter precision differs')
        tokens = encoder.tokenizer(query['sequence'], return_tensors='pt', truncation=False, padding=False, return_special_tokens_mask=True)
        mask = validate_token_mapping(query['sequence'], tokens['input_ids'][0].numpy(), tokens['attention_mask'][0].numpy(),
                                      tokens.pop('special_tokens_mask')[0].numpy(), encoder.tokenizer)
        inputs = {k: v for k, v in tokens.items() if k in ('input_ids', 'attention_mask')}
        with torch.inference_mode():
            output = encoder.model(**inputs, output_hidden_states=True, output_attentions=False, return_dict=True)
        require(len(output.hidden_states) == encoder.num_layers + 1, 'Official hidden-state tuple depth changed')
        selected = output.hidden_states[config['layer']]
        require(selected.dtype == torch.float32 and selected.device.type == 'cpu' and not selected.requires_grad and
                tuple(selected.shape[:2]) == (1, query['length_aa'] + 2), 'Official selected hidden state precision or shape differs')
        actual = selected[0, torch.as_tensor(mask)].numpy()
        require(actual.shape == saved.shape and saved.dtype == np.float32 and np.isfinite(actual).all() and np.isfinite(saved).all(),
                'Official selected residues do not match saved probe shape or finiteness')
        error = float(np.max(np.abs(actual.astype(np.float64) - saved.astype(np.float64))))
        require(error <= 1e-6, 'Captured midpoint states differ from the official API tuple')
        record = {'status': 'passed', 'query_id': query['query_id'], 'length': query['length_aa'], 'layer': config['layer'],
                  'official_tuple_length': len(output.hidden_states), 'residue_shape': list(actual.shape),
                  'max_absolute_error': error, 'tolerance': 1e-6, 'actual_dtype': str(selected.dtype), 'actual_device': str(selected.device),
                  'elapsed_seconds': time.monotonic() - started,
                  'peak_process_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20 if sys.platform == 'darwin' else 1024),
                  'scope': 'One fixed short-query intermediate-state interface audit; standard source probe separately checks all four lengths and repeatability. No exhaustive independent layer verification.'}
        store.publish('audit', record)
        del output, selected, actual
        store.complete(['audit'], {'guard': guard})
        print(json.dumps({'status': 'passed', 'job': str(store.root), **record}), flush=True)
    finally:
        encoder.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--probe-job', type=Path, required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('--expected-code-sha256', required=True)
    parser.add_argument('--expected-design-sha256', required=True)
    run(parser.parse_args())
