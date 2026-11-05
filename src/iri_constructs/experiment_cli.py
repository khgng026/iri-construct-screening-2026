"""Consistent T05 CLI and sequence-free JSON reports."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import logging
import statistics
from pathlib import Path

from .bench_inputs import load_proteins
from .config import configure_logging
from .fusion import compare_case, validate_case
from .metrics import mean_pool, retrieve, reviewed_labels
from .model_adapters import ESM2_REVISION, MODEL_REPOS, LocalAdapter, prepare_device
from .runtime import InputError, Measure, ROOT, failure, local_output, run_metadata, write_json


def execute(args, report: dict, adapter_factory=LocalAdapter) -> int:
    config, proteins, provenance, exclusions = load_proteins(args.inputs)
    report.update(provenance=provenance, inputs=[p.provenance() for p in proteins.values()],
                  input_exclusions=exclusions)
    if not proteins:
        raise InputError('No validated annotated protein units are available')
    labels, label_exclusions = reviewed_labels(config, proteins)
    report['label_exclusions'] = label_exclusions
    if args.experiment == 'comparison' and not labels:
        raise InputError('Comparison requires independently reviewed known family annotations')
    cases = config.get('fusion_cases', [])
    if args.experiment == 'fusion':
        if not isinstance(cases, list) or not cases:
            raise InputError('Fusion experiment requires a nonempty fusion_cases list')
        for case in cases:
            validate_case(case, proteins)
        ids = [case['case_id'] for case in cases]
        if len(ids) != len(set(ids)):
            raise InputError('Duplicate fusion case_id')
    if args.validate_only:
        report.update(status='validated_only', inference_executed=False)
        return 0
    device = prepare_device(args.device, args.seed)
    report['resolved_device'] = device
    if device.startswith('cuda'):
        import torch
        report['cuda_device_name'] = torch.cuda.get_device_name(device)
        report['cuda_total_memory_bytes'] = torch.cuda.get_device_properties(device).total_memory
    model_results, all_vectors = {}, {}
    for name in args.models:
        adapter = None
        revision = getattr(args, name + '_revision')
        load_measure = Measure(device)
        entry = {'status': 'started', 'inference_executed': False,
                 'requested_model': {'name': name, 'repository': MODEL_REPOS[name], 'revision': revision}}
        model_results[name] = entry
        logging.getLogger('iri_constructs').info('event=model_start model=%s', name)
        try:
            with load_measure:
                adapter = adapter_factory(name, device, revision, args.max_residues,
                                          offline=not args.allow_download)
            entry.update(model=adapter.metadata, loading_resources=load_measure.data)
            if args.experiment == 'fusion':
                rows = []
                for case in cases:
                    try:
                        rows.append(compare_case(case, proteins, adapter))
                    except Exception as exc:
                        rows.append({'case_id': case['case_id'], 'status': 'failed', **failure(exc)})
                entry.update(cases=rows, inference_executed=True,
                             status='ok' if all(r['status'] == 'ok' for r in rows) else 'partial',
                             coverage={'requested_cases': len(rows),
                                       'complete_cases': sum(r['status'] == 'ok' for r in rows)})
            else:
                vectors, details = {}, []
                selected = proteins if args.experiment == 'probe' else {k: proteins[k] for k in labels}
                cache = {}
                for iid, protein in selected.items():
                    item = {'input_id': iid}
                    with Measure(device) as measured:
                        try:
                            key = protein.sequence_sha256
                            cached = key in cache
                            if not cached:
                                cache[key] = mean_pool(adapter.residues(protein.sequence))
                            vectors[iid] = cache[key]
                            item.update(status='ok', embedding_dimension=len(vectors[iid]),
                                        cached_exact_sequence=cached)
                        except Exception as exc:
                            item.update(status='failed', **failure(exc))
                    item['resources'] = measured.data
                    details.append(item)
                entry.update(inputs=details, inference_executed=True,
                             status='ok' if len(vectors) == len(selected) else 'partial',
                             coverage={'requested_units': len(selected), 'embedded_units': len(vectors),
                                       'fraction': len(vectors)/len(selected) if selected else 0.0})
                unique_times = [i['resources']['wall_seconds'] for i in details
                                if i['status'] == 'ok' and not i['cached_exact_sequence']]
                entry['median_seconds_per_unique_protein'] = statistics.median(unique_times) if unique_times else None
                entry['timing_policy'] = 'Sequential first-pass inference, pooling and host transfer; no warm-up. Exact cache hits excluded from median.'
                if args.experiment == 'comparison':
                    entry['retrieval_full_panel'] = retrieve(vectors, labels, proteins)
                    archetypes = {p['record_id']: p['archetype'] for p in provenance['sources']}
                    entry['retrieval_by_archetype'] = {
                        archetype: retrieve(vectors, labels, proteins,
                            [iid for iid in labels if archetypes[proteins[iid].construct_id] == archetype])
                        for archetype in sorted(set(archetypes.values()))}
                    all_vectors[name] = vectors
        except Exception as exc:
            entry.update(status='blocked', **failure(exc))
        finally:
            entry['loading_resources'] = load_measure.data
            if adapter is not None:
                adapter.close()
        logging.getLogger('iri_constructs').info('event=model_end model=%s status=%s', name, entry['status'])
    report['models'] = model_results
    if args.experiment == 'comparison':
        common = set(labels)
        for name in args.models:
            common &= all_vectors.get(name, {}).keys()
        report['common_panel_ids'] = sorted(common)
        report['common_panel_retrieval'] = {
            name: retrieve(all_vectors.get(name, {}), labels, proteins, common, common) for name in args.models}
        report['limitations'].extend([
            'Exact duplicates, same construct and reviewer-linked independence groups are excluded as references.',
            'Residual homology, pretraining overlap and small convenience panels remain confounders.',
            'No automatic model selection; compare full-panel coverage and the exact common panel together.'])
    statuses = [entry['status'] for entry in model_results.values()]
    report['status'] = 'ok' if all(s == 'ok' for s in statuses) else 'blocked' if all(s == 'blocked' for s in statuses) else 'partial'
    report['inference_executed'] = any(e['inference_executed'] for e in model_results.values())
    return 0 if report['status'] == 'ok' else 1


def main(default_experiment: str | None = None, argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', choices=('fusion', 'comparison', 'probe'), default=default_experiment or 'comparison')
    parser.add_argument('--inputs', type=Path, required=True, help='Reviewed local t05-inputs-1 JSON; paths relative to this file')
    parser.add_argument('--output', type=Path, help='New subdirectory under repository results/ or data/processed/')
    parser.add_argument('--models', default=None, help='Comma-separated esm2,esmc,esm3; fusion/probe require esm2')
    parser.add_argument('--device', choices=('auto', 'cpu', 'cuda'), default='cpu')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--max-residues', type=int, default=1022)
    parser.add_argument('--allow-download', action='store_true', help='Allow pretrained weights downloads; biological inputs are always local')
    parser.add_argument('--validate-only', action='store_true', help='Validate reviewed inputs without loading model weights')
    parser.add_argument('--esm2-revision', default=ESM2_REVISION)
    parser.add_argument('--esmc-revision', default='main', help='Use resolved SHA from a prior report to freeze a rerun')
    parser.add_argument('--esm3-revision', default='main', help='Use resolved SHA from a prior report to freeze a rerun')
    args = parser.parse_args(argv)
    args.models = (args.models or ('esm2,esmc,esm3' if args.experiment == 'comparison' else 'esm2')).split(',')
    if (not args.models or len(args.models) != len(set(args.models)) or
            any(m not in ('esm2', 'esmc', 'esm3') for m in args.models)):
        parser.error('--models must contain unique names esm2,esmc,esm3')
    if args.experiment != 'comparison' and args.models != ['esm2']:
        parser.error('fusion and probe experiments require --models esm2')
    if not 1 <= args.max_residues <= 1022:
        parser.error('--max-residues must be in 1..1022')
    if not 0 <= args.seed < 2**32:
        parser.error('--seed must be in 0..2**32-1')
    configure_logging('INFO')
    try:
        output = local_output(args.output or ROOT / 'results' / ('t05_' + args.experiment) /
                              datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        if (output / 'report.json').exists():
            raise InputError('Output report already exists; choose a new output directory')
    except InputError as exc:
        parser.error(str(exc))
    settings = {k: v for k,v in vars(args).items() if k not in ('inputs', 'output')}
    report = run_metadata(args.experiment, settings)
    with Measure() as measured:
        try:
            code = execute(args, report)
        except Exception as exc:
            report.update(status='blocked', inference_executed=False, **failure(exc))
            code = 1
    report['resources_total'] = measured.data
    report['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    write_json(output / 'report.json', report)
    logging.getLogger('iri_constructs').info('event=experiment_end experiment=%s status=%s', args.experiment, report['status'])
    print('Report:', output / 'report.json')
    return code
