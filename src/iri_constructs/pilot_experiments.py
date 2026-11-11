"""A priori, reviewed-input pilot execution built on the existing T05 interfaces.

No record admission, source downloading, invented biological inputs or training.
Scientific results require the same reviewed inputs and frozen reference panel.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
from pathlib import Path

from .bench_inputs import ProteinInput, load_proteins
from .config import configure_logging
from .fusion import compare_case, validate_case
from .metrics import cosine, mean_pool, reviewed_labels
from .model_adapters import ESM2_REVISION, LocalAdapter, prepare_device
from .pilot_scoring import (freeze_panel, load_scoring_config, paired_cluster_bootstrap,
                            query_independence_clusters, rank_scores, score_rankings)
from .runtime import (InputError, Measure, digest_file, failure, local_output,
                      read_json, run_metadata, write_json)


class _CachedAdapter:
    """Reuse exact forward passes, recording cache use separately from biological units."""
    def __init__(self, adapter):
        self.adapter, self.device = adapter, adapter.device
        self.cache, self.events = {}, []

    def residues(self, sequence):
        digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        cached = digest in self.cache
        event = {'sequence_sha256': digest, 'length_aa': len(sequence), 'cache_hit': cached,
                 'status': 'attempted'}
        self.events.append(event)
        try:
            if not cached:
                self.cache[digest] = self.adapter.residues(sequence)
            event['status'] = 'ok'
        except Exception:
            event['status'] = 'failed'
            raise
        return self.cache[digest]


def _vectors(adapter, proteins, ids):
    pooled, cache, records = {}, {}, []
    for iid in sorted(ids):
        protein = proteins[iid]
        item = {'input_id': iid, 'length_aa': len(protein.sequence)}
        with Measure(adapter.device) as measured:
            try:
                cached = protein.sequence_sha256 in cache
                if not cached:
                    cache[protein.sequence_sha256] = mean_pool(adapter.residues(protein.sequence))
                pooled[iid] = cache[protein.sequence_sha256]
                item.update(status='ok', cache_hit=cached, embedding_dimension=len(pooled[iid]))
            except Exception as exc:
                item.update(status='failed', **failure(exc))
        item['resources'] = measured.data
        records.append(item)
    return pooled, records


def _cosine_scores(queries, references, panel):
    scores = {}
    for qid in panel['query_ids']:
        allowed = panel['eligible_reference_ids'][qid]
        if qid in panel['query_exclusions'] or qid not in queries or any(r not in references for r in allowed):
            continue
        try:
            scores[qid] = {rid: cosine(queries[qid], references[rid]) for rid in allowed}
        except InputError:
            # No partial reference ranking after invalid representations.
            continue
    return scores


def _score(scores, panel, labels):
    result = score_rankings(rank_scores(scores, panel), panel, labels)
    result['status'] = ('ok' if result['scored_queries'] == result['eligible_queries'] and result['scored_queries']
                        else 'partial' if result['scored_queries'] else 'blocked_no_scorable_queries')
    return result


def _fusion_panel(inputs, proteins, labels, query_ids, reference_ids, archetypes):
    cases = inputs.get('fusion_cases')
    if not isinstance(cases, list) or not cases:
        raise InputError('Experiment B requires reviewed documented fusion_cases')
    seen, used, lineage = set(), set(), []
    selected, groups, components, digests = [], set(), [], set()
    for case in cases:
        fusion, target, _, _ = validate_case(case, proteins)
        cid = case['case_id']
        if cid in seen:
            raise InputError('Duplicate fusion case_id')
        seen.add(cid)
        if case['target_id'] not in labels:
            raise InputError('Every fusion target requires a known independently reviewed family annotation')
        if query_ids is not None and case['target_id'] not in query_ids:
            continue
        qid = 'case/' + cid
        if qid in proteins:
            raise InputError('Fusion case query ID collides with a source input ID')
        tokens = {'parent:' + target.construct_id, 'parent:' + fusion.construct_id,
                  'group:' + labels[case['target_id']]['independence_group']}
        groups.add(labels[case['target_id']]['independence_group'])
        if case['fusion_id'] in labels:
            tokens.add('group:' + labels[case['fusion_id']]['independence_group'])
            groups.add(labels[case['fusion_id']]['independence_group'])
        for component in list(components):
            if component & tokens:
                tokens |= component
                components.remove(component)
        components.append(tokens)
        selected.append((qid, case, fusion, target, tokens))
        used.update((target.construct_id, fusion.construct_id))
        digests.update((target.sequence_sha256, fusion.sequence_sha256))
    virtual_queries, selected_cases = [], []
    for qid, case, fusion, target, tokens in selected:
        component = next(c for c in components if c & tokens)
        cluster = 'fusion-cluster-' + hashlib.sha256('|'.join(sorted(component)).encode()).hexdigest()[:20]
        proteins[qid] = ProteinInput(qid, fusion.construct_id, target.sequence, fusion.source_sha256, fusion.feature_id)
        labels[qid] = {**labels[case['target_id']], 'input_id': qid, 'independence_group': cluster,
                       'archetype': archetypes.get(fusion.construct_id, 'unknown')}
        virtual_queries.append(qid)
        selected_cases.append((qid, case))
        lineage.append({'query_id': qid, 'case_id': case['case_id'], 'target_id': target.input_id,
                        'fusion_id': fusion.input_id, 'target_parent': target.construct_id,
                        'fusion_parent': fusion.construct_id, 'independence_cluster': cluster})
    refs, excluded = [], []
    for rid in reference_ids:
        reasons = []
        if proteins[rid].construct_id in used:
            reasons.append('Same target/fusion source parent')
        if proteins[rid].sequence_sha256 in digests:
            reasons.append('Exact target or complete-fusion sequence copy')
        if labels[rid]['independence_group'] in groups:
            reasons.append('Reviewer-linked target/fusion independence group')
        if reasons:
            excluded.append({'input_id': rid, 'reasons': reasons})
        else:
            refs.append(rid)
    return selected_cases, virtual_queries, refs, lineage, excluded


def _compare_methods(methods, panel, config):
    comparisons = []
    clusters = query_independence_clusters(panel)
    for left, right in itertools.combinations(methods, 2):
        a = {r['query_id']: r['documented_family_recall_at_k'] for r in methods[left].get('rows', []) if not r.get('budget_saturated')}
        b = {r['query_id']: r['documented_family_recall_at_k'] for r in methods[right].get('rows', []) if not r.get('budget_saturated')}
        comparison = paired_cluster_bootstrap(a, b, clusters, seed=config['seed'],
            replicates=config['bootstrap_replicates'], min_clusters=config['min_independent_clusters'])
        comparisons.append({'left': left, 'right': right,
                            'excluded_budget_saturated_query_ids': panel.get('budget_saturated_query_ids', []),
                            **comparison})
    return comparisons


def run_pilot_experiment(protocol_path, inputs_path, output_path, *, device='cpu',
                         allow_download=False, adapter_factory=LocalAdapter):
    """Save a design freeze before admission, and a panel freeze before producing scores."""
    output = local_output(Path(output_path))
    if output.exists() and any(output.iterdir()):
        raise InputError('Experiment output must be a new empty ignored directory')
    output.mkdir(parents=True, exist_ok=True)
    report = run_metadata('reviewed_construct_pilot', {'device': device, 'allow_download': allow_download})
    report.update(status='blocked', inference_executed=False, scientific_scores_produced=False,
                  requested_experiment=None, methods={})
    with Measure() as total:
        try:
            config = load_scoring_config(protocol_path)
            report.update(requested_experiment=config['experiment'], protocol=config,
                          protocol_file_sha256=digest_file(Path(protocol_path)))
            # This is written even if no input has been admitted. No results were observed first.
            write_json(output / 'protocol-freeze.json', {'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
                'protocol': config, 'protocol_file_sha256': report['protocol_file_sha256'],
                'stage': 'before source admission and all scientific scoring'})
            if config['experiment'] == 'A':
                from .pilot import prepare_pilot
                audit = prepare_pilot(Path(inputs_path), output / 'annotation')
                report.update(annotation_preservation=audit, status=audit['status'],
                              scientific_scores_produced=bool(audit.get('scored_records', 0)))
                if audit['status'] == 'ok' and not audit.get('scored_records', 0):
                    report.update(status='blocked_no_observable_annotations',
                                  reason='Parsing succeeded but no source annotation evidence was available to score')
            else:
                raw_config = read_json(Path(inputs_path))
                if not raw_config.get('records'):
                    report.update(status='blocked_no_reviewed_inputs',
                                  reason='No admitted local records; candidate metadata and review queues do not authorize scoring')
                else:
                    _execute(config, inputs_path, output, report, device, allow_download, adapter_factory)
        except Exception as exc:
            report.update(status='blocked', **failure(exc))
    report['resources_total'] = total.data
    report['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    write_json(output / 'report.json', report)
    return report


def _execute(config, inputs_path, output, report, device, allow_download, adapter_factory):
    inputs, proteins, provenance, exclusions = load_proteins(Path(inputs_path))
    evidence_kind = inputs.get('evidence_kind', 'unclassified')
    if evidence_kind not in ('admitted_research', 'fictional_software_fixture', 'unclassified'):
        raise InputError('Unknown evidence_kind; distinguish admitted research from software fixtures')
    report.update(provenance=provenance, input_exclusions=exclusions,
                  inputs=[p.provenance() for p in proteins.values()], evidence_kind=evidence_kind)
    labels, label_exclusions = reviewed_labels(inputs, proteins)
    report['label_exclusions'] = label_exclusions
    archetypes = {source['record_id']: source['archetype'] for source in provenance['sources']}
    labels = {iid: {**value, 'archetype': archetypes.get(proteins[iid].construct_id, 'unknown')}
              for iid, value in labels.items()}
    queries = inputs.get('query_ids')
    refs = inputs.get('reference_ids', sorted(labels))
    if queries is not None and (not isinstance(queries, list) or len(queries) != len(set(queries))
                                 or any(q not in labels for q in queries)):
        raise InputError('query_ids must be a list of unique reviewed family input IDs')
    if not isinstance(refs, list) or any(r not in labels for r in refs):
        raise InputError('reference_ids must select reviewed known family annotations')
    cases = []
    if config['experiment'] == 'B':
        cases, queries, refs, lineage, removed = _fusion_panel(inputs, proteins, labels, queries, refs, archetypes)
        report['fusion_lineage'] = lineage
        report['fusion_reference_exclusions'] = removed
    panel = freeze_panel(config, proteins, labels, queries, refs)
    report['panel'] = panel
    write_json(output / 'panel-freeze.json', panel)
    if panel['status'] != 'frozen':
        report.update(status='blocked_insufficient_panel', reason='No query has an independent positive and competing reference family')
        return
    methods = report['methods']
    models = [name for name in config['models'] if name != 'alignment']
    if 'alignment' in config['models']:
        from .pilot_alignment import alignment_scores
        try:
            alignment = alignment_scores(proteins, panel, max_residues=config['max_residues'])
            for kind, key in [('alignment_raw', 'raw_scores'), ('alignment_reference_normalized', 'normalized_scores')]:
                methods[kind] = {**_score(alignment[key], panel, labels),
                    'resources': alignment['resources'], 'alignment_configuration': alignment['configuration'],
                    'pairwise_score_calls': alignment['pairwise_score_calls'],
                    'reference_self_score_calls': alignment['reference_self_score_calls'],
                    'failures': alignment['failures'], 'shared_execution': 'one raw-alignment pass; normalization reuses scores'}
        except Exception as exc:
            for kind in ('alignment_raw', 'alignment_reference_normalized'):
                methods[kind] = {'status': 'blocked', **failure(exc)}
    for name in models:
        adapter = None
        cached = None
        model = {'status': 'blocked', 'inference_executed': False}
        report.setdefault('models', {})[name] = model
        measured = Measure()
        try:
            resolved_device = prepare_device(device, config['seed'])
            measured = Measure(resolved_device)
            revision = config.get('model_revisions', {}).get(name, ESM2_REVISION if name == 'esm2' else 'main')
            with measured:
                adapter = adapter_factory(name, resolved_device, revision, config['max_residues'], offline=not allow_download)
            model.update(metadata=adapter.metadata, loading_resources=measured.data)
            # Pin the resolved weight revision before passing any sequence to the model.
            write_json(output / ('model-freeze-' + name + '.json'), adapter.metadata)
            cached = _CachedAdapter(adapter)
            if config['experiment'] == 'B':
                condition_vectors = {condition: {} for condition in ('isolated_target', 'complete_fusion', 'contextual_target')}
                case_results = []
                for qid, case in cases:
                    pooled = {}
                    result = compare_case(case, proteins, cached, pooled_sink=pooled)
                    case_results.append({'query_id': qid, **result})
                    for condition, vector in pooled.items():
                        condition_vectors[condition][qid] = vector
                references, details = _vectors(cached, proteins, set(refs))
                for condition, query_vectors in condition_vectors.items():
                    methods[name + '/' + condition] = _score(_cosine_scores(query_vectors, references, panel), panel, labels)
                model.update(cases=case_results, reference_inputs=details)
            else:
                wanted = set(panel['query_ids']) | set(panel['reference_ids'])
                vectors, details = _vectors(cached, proteins, wanted)
                methods[name] = _score(_cosine_scores(vectors, vectors, panel), panel, labels)
                model.update(inputs=details)
            model.update(status='attempted')
        except Exception as exc:
            model.update(status='blocked', loading_resources=measured.data, **failure(exc))
            names = [name + '/' + c for c in ('isolated_target', 'complete_fusion', 'contextual_target')] if config['experiment'] == 'B' else [name]
            for key in names:
                methods.setdefault(key, {'status': 'blocked', **failure(exc)})
        finally:
            if cached is not None:
                successful = sum(not e['cache_hit'] and e['status'] == 'ok' for e in cached.events)
                model.update(forward_requests=cached.events, unique_forward_calls=successful,
                             inference_attempted=bool(cached.events), inference_executed=bool(successful))
                report['inference_executed'] |= bool(successful)
                cached.cache.clear()
            if adapter is not None:
                adapter.close()
    report['paired_comparisons'] = _compare_methods(methods, panel, config)
    common = set(panel['query_ids'])
    for method in methods.values():
        common &= {row['query_id'] for row in method.get('rows', [])}
    report['common_scored_query_ids'] = sorted(common)
    report['scores_produced'] = any(method.get('scored_queries', 0) for method in methods.values())
    report['scientific_scores_produced'] = report['scores_produced'] and evidence_kind == 'admitted_research'
    report['status'] = ('ok' if methods and all(method['status'] == 'ok' for method in methods.values())
                        else 'partial' if report['scores_produced'] else 'blocked_no_scientific_scores')
    eligible = {q for q in panel['query_ids'] if q not in panel['query_exclusions']}
    if (report['scores_produced'] and eligible and
            set(panel.get('budget_saturated_query_ids', [])) == eligible):
        report['status'] = 'descriptive_only_budget_saturated_panel'
    report['limitations'].extend([
        'Known fusion boundaries are annotation-guided controls, never automatic window predictions.',
        'No raw DNA, artificial concatenations, exact-self retrieval, training or biological function assay.',
        'Reference budget is fixed for every condition; model calls/lengths and reusable caches are reported separately.',
        'Paired intervals use matched independence groups; absent pairs and fewer than eight groups limit interpretation.'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--device', choices=('cpu', 'cuda', 'auto'), default='cpu')
    parser.add_argument('--allow-download', action='store_true', help='Allow model weights only; no biological record download')
    args = parser.parse_args(argv)
    log = configure_logging('INFO')
    try:
        report = run_pilot_experiment(args.protocol, args.inputs, args.output, device=args.device,
                                      allow_download=args.allow_download)
    except InputError as exc:
        parser.error(str(exc))
    log.info('event=pilot_experiment_end experiment=%s status=%s', report['requested_experiment'], report['status'])
    print('Experiment report:', args.output / 'report.json')
    return 0 if report['status'] == 'ok' else 1
