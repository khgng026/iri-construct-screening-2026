#!/usr/bin/env python3
"""Reproducible local T05 experiment. Run --help; no sequence contents in logs."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
os.environ.setdefault('HF_HUB_DISABLE_TELEMETRY', '1')
from iri_constructs.t05_data import build_query, load_catalogue, local_data_path, sha256
from iri_constructs.schema_check import check_schema_example

MANIFEST = ROOT / 'data/manifests/t05_benign_catalogue.json'
MODEL = ROOT / 'config/t05_model.json'
QUERIES = ROOT / 'config/t05_queries.json'
RESULTS = ROOT / 'results/t05'


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def save_sequence(identifier, sequence):
    relative = f'data/processed/t05/queries/{identifier}.fasta'
    path = local_data_path(ROOT, relative)
    data = f'>{identifier}\n{sequence}\n'.encode('ascii')
    if path.exists() and path.read_bytes() != data:
        raise ValueError('Existing artificial query differs; refusing overwrite')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return {'path': relative, 'format': 'fasta', 'record_id': identifier,
            'sha256': sha256(sequence.encode('ascii')), 'length': len(sequence), 'alphabet': 'protein'}


def normalized_query(spec, sequences):
    query, oracle, (start, end) = build_query(spec['components'], spec['target'], sequences)
    qid = spec['id']
    coordinate = 'protein_0_based_half_open'
    whole_ref = save_sequence(qid + '-whole', query)
    oracle_ref = save_sequence(qid + '-oracle', oracle)
    features, spans = [], []
    offset = 0
    for i, accession in enumerate(spec['components']):
        loc = {'start0': offset, 'end0': offset + len(sequences[accession]), 'strand': None}
        features.append({'feature_id': f'{qid}-f{i+1}', 'construct_id': qid,
                         'feature_label': 'other', 'display_name': accession,
                         'location': loc, 'coordinate_system': coordinate,
                         'annotation_source': f't05-fixed-v1 known construction; canonical UniProt {accession}',
                         'source_record_id': accession, 'annotation_status': 'known'})
        spans.append(loc)
        offset = loc['end0']
    target_index = spec['components'].index(spec['target'])
    units = []
    for arm, ref, fids, segments in [
        ('whole', whole_ref, [f['feature_id'] for f in features], spans),
        ('oracle', oracle_ref, [features[target_index]['feature_id']], [spans[target_index]])]:
        units.append({'unit_id': f'{qid}-{arm}', 'construct_id': qid, 'feature_ids': fids,
                      'source_segments': segments, 'molecule_type': 'protein',
                      'coordinate_system': coordinate, 'model_input_status': 'ready',
                      'sequence_ref': ref, 'representation_note':
                      'Direct digital concatenation' if arm == 'whole' else 'Known target interval, sliced then re-encoded; oracle control'})
    document = {'schema_version': '0.2', 'construct': {
        'construct_id': qid, 'source': 'artificial_construction' if spec['artificial'] else 'clean_reference_control',
        'source_url': None, 'accession': None, 'accession_version': 't05-fixed-v1',
        'molecule_type': 'protein', 'sequence_length_nt': None, 'sequence_length_aa': len(query),
        'coordinate_system': coordinate, 'sequence_ref': whole_ref, 'topology': 'linear',
        'primary_archetype': 'unknown', 'secondary_archetypes': [], 'annotation_status': 'known',
        'artificial': spec['artificial'], 'license_note': 'UniProt CC BY 4.0; source manifest supplies attribution. Local sequences only.',
        'review_note': 'Protein-only digital test; no observed DNA construct, CDS coordinates or fusion function inferred.'},
        'features': features, 'units': units}
    path = ROOT / f'data/processed/t05/normalized/{qid}.json'
    write_json(path, document)
    if check_schema_example(path):
        raise ValueError('Normalized T05 schema validation failed')
    return query, oracle, (start, end), document


def prepare(download):
    manifest = read_json(MANIFEST)
    sequences = load_catalogue(ROOT, manifest, download=download)
    for spec in read_json(QUERIES)['queries']:
        normalized_query(spec, sequences)
    if download:
        from huggingface_hub import snapshot_download
        config = read_json(MODEL)
        started = perf_counter()
        snapshot_download(config['model_id'], revision=config['revision'], token=False,
                          allow_patterns=['model.safetensors', 'config.json', 'tokenizer_config.json',
                                          'special_tokens_map.json', 'vocab.txt', 'README.md'])
        print(f'Model preparation: {perf_counter()-started:.3f} seconds (may use existing cache)')
    print(f'Prepared {len(sequences)} reviewed references and four fixed queries; no sequence contents logged.')


def run(run_id):
    import numpy as np
    import torch
    from huggingface_hub import hf_hub_download
    from iri_constructs.t05_embeddings import Encoder, rank_cosine

    if not re.fullmatch(r'[A-Za-z0-9_-]+', run_id):
        raise ValueError('Use a simple alphanumeric run ID')
    out = RESULTS / run_id
    out.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    result = {'run_id': run_id, 'started_at_utc': datetime.now(timezone.utc).isoformat(),
              'failures': [], 'queries': [], 'catalogue': [], 'network_policy': 'local-only; no downloads'}
    stage = 'input_validation'
    try:
        manifest = read_json(MANIFEST)
        config = read_json(MODEL)
        sequences = load_catalogue(ROOT, manifest)
        specs = read_json(QUERIES)['queries']
        result.update({'catalogue_size': len(sequences), 'excluded_records': manifest['excluded_records'],
                       'model': config, 'input_hashes': {str(p.relative_to(ROOT)): sha256(p.read_bytes())
                       for p in (MANIFEST, MODEL, QUERIES, ROOT / 'requirements-t05.lock.txt')},
                       'code_hashes': {str(p.relative_to(ROOT)): sha256(p.read_bytes()) for p in
                       [Path(__file__), *sorted((ROOT / 'src/iri_constructs').glob('*.py'))]},
                       'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                       'git_worktree_dirty': bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()),
                       'environment': {'python': platform.python_version(), 'os': platform.platform(),
                       'machine': platform.machine(),
                       'cpu': subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip() if sys.platform == 'darwin' else platform.processor(),
                       'memory_gb': int(subprocess.check_output(['sysctl','-n','hw.memsize'],text=True)) / 2**30 if sys.platform == 'darwin' else None,
                       'packages': dict(sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions())),
                       'mps_available': torch.backends.mps.is_available(), 'mps_used': False}})
        stage = 'model_load'
        encoder = Encoder(config, offline=True)
        result['model_load_seconds'] = encoder.load_seconds
        result['model_details'] = {'layers': encoder.model.config.num_hidden_layers,
            'hidden_size': encoder.model.config.hidden_size,
            'loaded_parameter_count': sum(p.numel() for p in encoder.model.parameters()),
            'model_config_max_position_embeddings': encoder.model.config.max_position_embeddings,
            'position_embedding_type': encoder.model.config.position_embedding_type,
            'tokenizer_model_max_length': encoder.tokenizer.model_max_length,
            'operational_residue_limit': config['max_residues'],
            'dtype': str(next(encoder.model.parameters()).dtype),
            'load_note': 'Masked-language head intentionally unused; no model training or pooler.'}
        result['model_file_sha256'] = {}
        for name in ['model.safetensors','config.json','tokenizer_config.json','special_tokens_map.json','vocab.txt']:
            cached = hf_hub_download(config['model_id'], name, revision=config['revision'], local_files_only=True)
            result['model_file_sha256'][name] = sha256(Path(cached).read_bytes())
        stage = 'catalogue_embedding'
        vectors, accessions = [], sorted(sequences)
        arrays = {}
        for accession in accessions:
            pooled, residues, info = encoder.encode([sequences[accession]])
            vectors.append(pooled[0])
            arrays[f'catalogue_{accession}'] = pooled[0].numpy()
            result['catalogue'].append({'accession': accession, **info})
        references = torch.stack(vectors)
        stage = 'padding_integration_check'
        # A real short protein paired with a longer one must retain its embedding.
        batch, _, batch_info = encoder.encode([sequences['P42212'], sequences['P00560']])
        clean_index = accessions.index('P42212')
        max_diff = (batch[0]-references[clean_index]).abs().max().item()
        padding_passed = torch.allclose(batch[0], references[clean_index], atol=1e-5, rtol=1e-5)
        result['padding_check'] = {'passed': padding_passed, 'max_abs_difference': max_diff,
                                   'atol': 1e-5, 'rtol': 1e-5, **batch_info}
        if not padding_passed:
            raise ValueError('Single-versus-padded embedding check failed')
        for spec in specs:
            stage = f'query_{spec["id"]}'
            construction_start = perf_counter()
            query, oracle, (start, end), document = normalized_query(spec, sequences)
            construction_seconds = perf_counter() - construction_start
            item = {**spec, 'whole_length_aa': len(query), 'target_interval_aa': [start, end],
                    'target_fraction': len(oracle)/len(query), 'construction_seconds': construction_seconds,
                    'arms': {}}
            for arm, sequence in [('whole',query),('oracle',oracle)]:
                pooled, residues, info = encoder.encode([sequence])
                rank_start = perf_counter()
                ranking = rank_cosine(pooled[0], references, accessions)
                rank_seconds = perf_counter() - rank_start
                target = next(row for row in ranking if row['accession']==spec['target'])
                item['arms'][arm] = {**info, 'sequence_length_aa': len(sequence),
                    'sequence_sha256': sha256(sequence.encode('ascii')), 'target_rank': target['rank'],
                    'target_score': target['cosine_similarity'], 'catalogue_size': len(accessions),
                    'rank_seconds': rank_seconds, 'ranking': ranking}
                arrays[f'{spec["id"]}_{arm}_pooled'] = pooled[0].numpy()
                arrays[f'{spec["id"]}_{arm}_residues'] = residues[0].numpy()
            if spec['id'] == 'clean_gfp' and (item['arms']['whole']['target_rank'] != 1 or
                                               abs(item['arms']['whole']['target_score']-1) > 1e-5):
                raise ValueError('Clean self-retrieval sanity check failed')
            result['queries'].append(item)
            print(f'{spec["id"]}: whole rank {item["arms"]["whole"]["target_rank"]}/{len(accessions)}, '
                  f'oracle rank {item["arms"]["oracle"]["target_rank"]}/{len(accessions)}')
        np.savez_compressed(out / 'embeddings.npz', **arrays)
        with (out / 'rankings.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['query_id','arm','expected_target','catalogue_size',
                                                       'rank','accession','cosine_similarity','exact_score_tie_count'])
            writer.writeheader()
            for query in result['queries']:
                for arm, info in query['arms'].items():
                    for row in info['ranking']:
                        writer.writerow({'query_id':query['id'],'arm':arm,'expected_target':query['target'],
                                         'catalogue_size':len(accessions), **row})
        result['status'] = 'passed'
    except Exception as exc:
        result['status'] = 'failed'
        # Do not echo exceptions from third-party parsers: they may contain raw data.
        result['failures'].append({'stage':stage,'error_type':type(exc).__name__})
        raise
    finally:
        result['total_seconds'] = perf_counter() - started
        write_json(out / 'results.json', result)
    print(f'Run {run_id}: passed in {result["total_seconds"]:.3f} seconds')


def compare(first_id, second_id):
    import numpy as np
    for identifier in [first_id, second_id]:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', identifier):
            raise ValueError('Invalid run ID')
    first_dir, second_dir = RESULTS / first_id, RESULTS / second_id
    first, second = read_json(first_dir/'results.json'), read_json(second_dir/'results.json')
    for key in ['input_hashes','code_hashes','model','model_file_sha256','catalogue_size','environment']:
        if first[key] != second[key]:
            raise ValueError('Runs have mismatched provenance/configuration')
    if first['status'] != 'passed' or second['status'] != 'passed':
        raise ValueError('Cannot compare failed runs')
    if [q['id'] for q in first['queries']] != [q['id'] for q in second['queries']] or not first['queries']:
        raise ValueError('Empty or mismatched query sets')
    score_diff, count = 0.0, 0
    for a,b in zip(first['queries'],second['queries'],strict=True):
        for arm in ['whole','oracle']:
            x,y = a['arms'][arm],b['arms'][arm]
            if x['sequence_sha256'] != y['sequence_sha256']:
                raise ValueError('Query checksums differ')
            if [r['accession'] for r in x['ranking']] != [r['accession'] for r in y['ranking']]:
                raise ValueError('Ranking differs')
            for rx,ry in zip(x['ranking'],y['ranking'],strict=True):
                score_diff = max(score_diff,abs(rx['cosine_similarity']-ry['cosine_similarity']))
                count += 1
    a,b = np.load(first_dir/'embeddings.npz'),np.load(second_dir/'embeddings.npz')
    if not a.files or set(a.files) != set(b.files):
        raise ValueError('Empty or mismatched embedding arrays')
    vector_diff = 0.0
    for name in a.files:
        if a[name].shape != b[name].shape:
            raise ValueError('Embedding shapes differ')
        vector_diff = max(vector_diff,float(np.max(np.abs(a[name]-b[name]))))
        if not np.allclose(a[name], b[name], atol=1e-6, rtol=0):
            raise ValueError('Embedding reproducibility check failed')
    if score_diff > 1e-6:
        raise ValueError('Score reproducibility check failed')
    comparison = {'passed':True,'runs':[first_id,second_id],'score_count':count,
                  'embedding_array_count':len(a.files),'atol':1e-6,'rtol':0,
                  'max_score_abs_difference':score_diff,'max_embedding_abs_difference':vector_diff,
                  'all_rankings_identical':True,'all_provenance_hashes_identical':True}
    write_json(RESULTS / f'comparison-{first_id}-{second_id}.json', comparison)
    print(json.dumps(comparison,indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    p = sub.add_parser('prepare',help='Verify frozen data, construct local inputs; optionally download missing files')
    p.add_argument('--download',action='store_true')
    p = sub.add_parser('run',help='Run CPU inference entirely offline; refuses to overwrite results')
    p.add_argument('--run-id',required=True)
    p = sub.add_parser('compare',help='Compare full outputs of two completed runs')
    p.add_argument('first'); p.add_argument('second')
    args = parser.parse_args()
    try:
        if args.command == 'prepare': prepare(args.download)
        elif args.command == 'run': run(args.run_id)
        else: compare(args.first,args.second)
    except Exception as exc:
        print(f'T05 failed ({type(exc).__name__}); check recorded failure stage and input files.',file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
