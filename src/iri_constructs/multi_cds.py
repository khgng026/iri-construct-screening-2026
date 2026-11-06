"""Prospective T07 v1.1 E2: documented CDS units, one parent retrieval budget.

The low-level freeze/scoring functions accept already admitted metadata. The
executable path enforces the existing T02 source gates before reading sequences.
No DNA model inputs, concatenated proteins, or inferred translations are made.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

from .annotations import normalize_record
from .bench_inputs import ProteinInput, read_source, reviewed_sources, validate_sequence
from .metrics import cosine, mean_pool, reviewed_labels
from .model_adapters import ESM2_REVISION, LocalAdapter, prepare_device
from .pilot_scoring import paired_cluster_bootstrap
from .runtime import (InputError, Measure, digest_file, failure, local_output,
                      read_json, require_text, run_metadata, write_json)


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode('utf-8')).hexdigest()


def _unique_ids(ids, name):
    if (not isinstance(ids, list) or any(not isinstance(i, str) or not i for i in ids)
            or len(ids) != len(set(ids))):
        raise InputError(name + ' must be a list of unique nonempty metadata IDs')
    return ids


def load_multi_cds_inputs(inputs_path, *, max_residues=1022):
    """Reuse T05 admission and T03 normalization, retaining withheld CDS denominators."""
    inputs, sources, provenance = reviewed_sources(Path(inputs_path))
    if any(source.review['annotation_status'] != 'approved' for source in sources):
        raise InputError('E2 requires approved T02 annotation review before any parsing')
    dna_parents = [s.review['record_id'] for s in sources if s.molecule == 'dna']
    selected = _unique_ids(inputs.get('parent_ids', dna_parents), 'parent_ids')
    if len(dna_parents) != len(set(dna_parents)) or any(p not in dna_parents for p in selected):
        raise InputError('E2 parent roles must identify each reviewed DNA parent exactly once')
    proteins, parents, all_unit_ids = {}, {}, set()
    for source in sources:
        record = read_source(source)
        checksum = source.review['checksum_sha256'].lower()
        if source.molecule == 'protein':
            proteins[source.input_id] = ProteinInput(source.input_id, source.review['record_id'],
                                                     str(record.seq).upper(), checksum)
            continue
        document, audit, extracted = normalize_record(record, source.review, checksum)
        features = {f['feature_id']: f for f in document['features']}
        units = [u for u in document['units'] if features[u['feature_ids'][0]]['feature_label'] == 'CDS']
        unit_ids = [source.input_id + '/' + u['unit_id'] for u in units]
        all_unit_ids.update(unit_ids)
        for unit in units:
            iid = source.input_id + '/' + unit['unit_id']
            if unit['unit_id'] in extracted:
                proteins[iid] = ProteinInput(iid, source.review['record_id'], extracted[unit['unit_id']],
                                            checksum, unit['feature_ids'][0])
        if source.review['record_id'] not in selected:
            continue
        digests, withheld = set(), []
        for unit in units:
            iid = source.input_id + '/' + unit['unit_id']
            original = audit['features'][unit['feature_ids'][0]]
            for translation in original['qualifiers'].get('translation', []):
                existing = ''.join(translation.split()).upper()
                if existing.endswith('*'):
                    existing = existing[:-1]
                if existing:
                    digests.add(hashlib.sha256(existing.encode('utf-8')).hexdigest())
            if unit['unit_id'] not in extracted:
                withheld.append({'unit_id': iid, 'reason': unit['representation_note']})
        parents[source.review['record_id']] = {
            'unit_ids': sorted(unit_ids), 'source_cds_count': len(units),
            'archetype': source.review['archetype'], 'source_sha256': checksum,
            'all_component_sequence_sha256': sorted(digests), 'unit_exclusions': withheld,
            'documented_family_by_unit': {},
        }
    raw_labels = inputs.get('families', [])
    if not isinstance(raw_labels, list) or any(not isinstance(e, dict) for e in raw_labels):
        raise InputError('families must contain separately reviewed annotation objects')
    seen, known = set(), {}
    for entry in raw_labels:
        iid = require_text(entry, 'input_id')
        if iid in seen:
            raise InputError('Duplicate E2 unit family annotation')
        seen.add(iid)
        if entry.get('annotation_status') == 'known' and entry.get('review_status') == 'approved':
            for field in ('family', 'evidence', 'independence_group'):
                require_text(entry, field)
            if iid not in proteins and iid not in all_unit_ids:
                raise InputError('A known E2 family annotation references no documented source unit')
            known[iid] = entry
    eligible_labels, label_exclusions = reviewed_labels(
        {'families': [entry for entry in raw_labels if entry['input_id'] in proteins]}, proteins)
    for parent in parents.values():
        parent['documented_family_by_unit'] = {iid: known[iid]['family'] for iid in parent['unit_ids'] if iid in known}
        parent['reviewed_groups'] = sorted({known[iid]['independence_group'] for iid in parent['unit_ids'] if iid in known})
        for iid in parent['unit_ids']:
            if iid in proteins and iid not in eligible_labels:
                parent['unit_exclusions'].append({'unit_id': iid, 'reason': 'Family annotation is missing, unknown or unapproved'})
            elif iid in proteins and len(proteins[iid].sequence) > max_residues:
                parent['unit_exclusions'].append({'unit_id': iid, 'reason': 'Annotated protein exceeds the frozen common residue cap'})
    references = _unique_ids(inputs.get('reference_ids'), 'reference_ids')
    if any(r not in eligible_labels for r in references):
        raise InputError('E2 references require supported proteins and known approved family annotations')
    provenance['label_exclusions'] = label_exclusions
    return inputs, proteins, eligible_labels, parents, provenance


def freeze_multi_cds_panel(config, proteins, labels, parents, reference_ids):
    """Freeze parent units, gold families and all-component leakage exclusions.

    Parent dictionaries include every documented source CDS in ``unit_ids`` and
    separately reviewed ``documented_family_by_unit``; absent annotations stay
    unknown. No method result changes eligibility or family denominators.
    """
    from .extension_config import validate_extension_config
    validate_extension_config(config,expected='E2')
    if config.get('extension', config.get('experiment')) != 'E2' or config.get('protocol_version') != '1.1':
        raise InputError('E2 requires a separately versioned T07 v1.1 protocol')
    if config.get('top_k') != 5 or config.get('retrieval_budget') != 5:
        raise InputError('E2 requires exactly five total parent hits, never five per CDS')
    cap = config.get('max_residues')
    if type(cap) is not int or not 1 <= cap <= 1022:
        raise InputError('E2 max_residues must be in 1..1022')
    refs = sorted(_unique_ids(reference_ids, 'reference_ids'))
    admitted_labels, excluded = reviewed_labels({'families': list(labels.values())}, proteins)
    if (excluded or set(admitted_labels) != set(labels) or
            any(iid != entry.get('input_id') for iid,entry in labels.items())):
        raise InputError('E2 panel contains missing, unknown or unapproved unit labels')
    if any(r not in labels for r in refs):
        raise InputError('E2 reference IDs lack reviewed family labels')
    snapshots, fingerprints = {}, {}
    for pid, parent in sorted(parents.items()):
        all_ids = sorted(_unique_ids(parent.get('unit_ids'), 'parent unit_ids'))
        count = parent.get('source_cds_count')
        if type(count) is not int or count != len(all_ids):
            raise InputError('Every documented CDS must remain in the parent source-count denominator')
        family_by_unit = parent.get('documented_family_by_unit', {})
        if (not isinstance(family_by_unit, dict) or set(family_by_unit) - set(all_ids)
                or any(not isinstance(f, str) or not f for f in family_by_unit.values())):
            raise InputError('Documented parent family expectations must map original CDS IDs to known families')
        eligible, missing = [], list(parent.get('unit_exclusions', []))
        digests = set(parent.get('all_component_sequence_sha256', []))
        groups = set(parent.get('reviewed_groups', []))
        for uid in all_ids:
            if uid in proteins:
                if proteins[uid].construct_id != pid:
                    raise InputError('A CDS input belongs to a different parent construct')
                digests.add(proteins[uid].sequence_sha256)
            if uid in labels:
                groups.add(labels[uid]['independence_group'])
                if family_by_unit.get(uid) != labels[uid]['family']:
                    raise InputError('Frozen documented-family expectation differs from reviewed source label')
            if uid not in proteins or uid not in labels:
                continue
            try:
                validate_sequence(proteins[uid].sequence, 'protein')
                if len(proteins[uid].sequence) > cap:
                    raise InputError('CDS exceeds common experimental cap')
                eligible.append(uid)
            except InputError as exc:
                missing.append({'unit_id': uid, 'reason': str(exc)})
        baseline = min(eligible, key=lambda u: (-len(proteins[u].sequence), u)) if eligible else None
        allowed, removals = [], []
        for rid in refs:
            reasons = []
            if proteins[rid].construct_id == pid:
                reasons.append('same_parent')
            if proteins[rid].sequence_sha256 in digests:
                reasons.append('exact_copy_of_parent_component')
            if labels[rid]['independence_group'] in groups:
                reasons.append('reviewer_linked_parent_component_group')
            if reasons:
                removals.append({'reference_id': rid, 'reasons': reasons})
            else:
                allowed.append(rid)
        all_families = sorted(set(family_by_unit.values()))
        observed = {labels[r]['family'] for r in allowed}
        unit_families = {labels[u]['family'] for u in eligible}
        supported = sorted(set(all_families) & observed & unit_families)
        unknown = sorted(set(all_ids) - set(family_by_unit))
        full_reason = ('no_eligible_documented_cds' if not eligible else
                       'unknown_documented_family_annotations' if unknown else
                       'no_independent_reference_for_every_documented_family' if not set(all_families) <= observed else
                       'no_competing_reference_family' if not observed - set(all_families) else None)
        matched_reason = ('no_independently_supported_eligible_families' if not supported else
                          'no_competing_reference_family' if not observed - set(supported) else None)
        snapshots[pid] = {
            'unit_ids': all_ids, 'eligible_unit_ids': sorted(eligible), 'baseline_unit_id': baseline,
            'source_cds_count': count, 'eligible_cds_count': len(eligible),
            'parsed_cds_count': sum(u in proteins for u in all_ids),
            'parsed_cds_coverage': sum(u in proteins for u in all_ids)/count if count else None,
            'model_label_eligible_cds_coverage': len(eligible)/count if count else None,
            'archetype': parent.get('archetype', 'unknown'),
            'source_sha256': parent.get('source_sha256'),
            'documented_family_by_unit': dict(sorted(family_by_unit.items())),
            'full_expected_families': all_families, 'matched_supported_families': supported,
            'unsupported_documented_families': sorted(set(all_families)-set(supported)),
            'unknown_family_unit_ids': unknown, 'unit_exclusions': missing,
            'reviewed_groups': sorted(groups), 'all_component_sequence_sha256': sorted(digests),
            'eligible_reference_ids': allowed, 'reference_exclusions': removals,
            'full_abstention_reason': full_reason, 'matched_abstention_reason': matched_reason,
            'budget_saturated': len(allowed) <= 5,
        }
        for uid in eligible:
            fingerprints[uid] = {**proteins[uid].provenance(), 'family': labels[uid]['family'],
                                 'independence_group': labels[uid]['independence_group']}
    for rid in refs:
        fingerprints[rid] = {**proteins[rid].provenance(), 'family': labels[rid]['family'],
                             'independence_group': labels[rid]['independence_group']}
    # Multiple parents sharing a documented component/group are correlated. Bind
    # conservative connected clusters before any method produces scores.
    components = []
    tokens_by_parent = {}
    for pid,parent in snapshots.items():
        tokens = {'parent:'+pid} | {'group:'+g for g in parent['reviewed_groups']} | {
            'component:'+d for d in parent['all_component_sequence_sha256']}
        tokens_by_parent[pid] = tokens.copy()
        for component in list(components):
            if tokens & component:
                tokens |= component
                components.remove(component)
        components.append(tokens)
    for pid,parent in snapshots.items():
        component = next(c for c in components if c & tokens_by_parent[pid])
        parent['independence_cluster'] = 'parent-cluster-'+_hash(sorted(component))[:20]
    panel = {'schema_version': 'multi-cds-panel-1', 'protocol_sha256': _hash(config),
             'protocol_version': '1.1', 'top_k': 5, 'retrieval_budget': 5,
             'parent_ids': sorted(snapshots), 'reference_ids': refs, 'parents': snapshots,
             'input_fingerprints': fingerprints,
             'status': 'frozen' if any(p['matched_abstention_reason'] is None for p in snapshots.values())
                       else 'blocked_insufficient_parent_panel',
             'merge_policy': 'maximum per reference across all eligible CDS scores, then one stable top-five parent ranking'}
    panel['panel_sha256'] = _hash(panel)
    return panel


def _validate_panel(panel, labels):
    if panel.get('schema_version') != 'multi-cds-panel-1' or panel.get('panel_sha256') != _hash(
            {k:v for k,v in panel.items() if k != 'panel_sha256'}):
        raise InputError('E2 frozen parent panel changed after admission')
    for iid, fingerprint in panel['input_fingerprints'].items():
        entry = labels.get(iid, {})
        if (entry.get('family') != fingerprint['family'] or
                entry.get('independence_group') != fingerprint['independence_group'] or
                entry.get('annotation_status') != 'known' or entry.get('review_status') != 'approved'):
            raise InputError('E2 labels changed after parent panel freeze')


def _ranking(scores, needed_units, allowed):
    if any(uid not in scores for uid in needed_units):
        return None, 'missing_unit_score_panel'
    merged = {}
    for uid in needed_units:
        row = scores[uid]
        if not isinstance(row, dict):
            raise InputError('Each CDS score row must be keyed by frozen reference ID')
        if set(row)-set(allowed):
            raise InputError('E2 score row contains excluded or out-of-panel reference leakage')
        if set(row) != set(allowed):
            return None, 'incomplete_reference_scores'
        if any(isinstance(s, bool) or not isinstance(s, (int,float)) or not math.isfinite(s) for s in row.values()):
            return None, 'nonfinite_reference_scores'
        for rid, score in row.items():
            merged[rid] = max(merged.get(rid, -math.inf), score)
    return sorted(merged, key=lambda r: (-merged[r], r))[:5], None


def _family_row(pid, parent, hits, expected, labels):
    recovered = set(expected) & {labels[r]['family'] for r in hits}
    return {'parent_id': pid, 'archetype': parent['archetype'], 'reference_ids': hits,
            'expected_families': expected, 'recovered_families': sorted(recovered),
            'distinct_expected_families': len(expected), 'distinct_recovered_families': len(recovered),
            'documented_family_recall_at_k': len(recovered)/len(expected),
            'all_documented_families_recovered': int(set(expected) <= recovered),
            'returned_hits': len(hits), 'budget_saturated': parent['budget_saturated']}


def _aggregate(rows, requested):
    return {'requested_parents': requested, 'scored_parents': len(rows),
            'coverage': len(rows)/requested if requested else 0.0,
            'documented_family_recall_at_k': sum(r['documented_family_recall_at_k'] for r in rows)/len(rows) if rows else None,
            'complete_family_recovery_rate': sum(r['all_documented_families_recovered'] for r in rows)/len(rows) if rows else None,
            'rows': rows}


def score_multi_cds(unit_reference_scores, panel, labels):
    """Compare one-CDS versus all-CDS representations on identical frozen references."""
    _validate_panel(panel, labels)
    query_units = {u for p in panel['parents'].values() for u in p['eligible_unit_ids']}
    if not isinstance(unit_reference_scores, dict) or set(unit_reference_scores)-query_units:
        raise InputError('E2 scores contain units outside the frozen parent panel')
    outputs = {condition: {'full': [], 'matched': [], 'abstentions': []}
               for condition in ('longest_cds', 'all_cds')}
    controls, control_abstentions = [], []
    for pid, parent in panel['parents'].items():
        allowed = parent['eligible_reference_ids']
        for condition in outputs:
            units = [parent['baseline_unit_id']] if condition == 'longest_cds' and parent['baseline_unit_id'] else parent['eligible_unit_ids']
            hits, failure_reason = _ranking(unit_reference_scores, units, allowed) if units else (None, 'no_eligible_documented_cds')
            for subset, families in [('full', parent['full_expected_families']), ('matched', parent['matched_supported_families'])]:
                reason = parent[subset+'_abstention_reason'] or failure_reason
                if reason:
                    outputs[condition]['abstentions'].append({'parent_id': pid, 'subset': subset, 'reason': reason})
                else:
                    outputs[condition][subset].append(_family_row(pid, parent, hits, families, labels))
        # Same CDS scores, query IDs and reference IDs: this control isolates no
        # per-unit accuracy effect; parent gains must not be relabeled unit gains.
        baseline = parent['baseline_unit_id']
        if baseline is not None:
            hits, reason = _ranking(unit_reference_scores, [baseline], allowed)
            observed = {labels[r]['family'] for r in allowed}
            if labels[baseline]['family'] not in observed or not observed-{labels[baseline]['family']}:
                reason = 'matched_unit_lacks_independent_positive_or_competing_family'
            if reason is None:
                controls.append({'parent_id': pid, 'unit_id': baseline, 'reference_ids': hits,
                    'expected_family': labels[baseline]['family'],
                    'family_retrieved': labels[baseline]['family'] in {labels[r]['family'] for r in hits},
                    'comparison': 'identical unit scores in longest and all-CDS conditions; not additional independent evidence'})
            else:
                control_abstentions.append({'parent_id':pid,'unit_id':baseline,'reason':reason})
    for condition, item in outputs.items():
        item['full'] = _aggregate(item['full'], len(panel['parent_ids']))
        item['matched'] = _aggregate(item['matched'], len(panel['parent_ids']))
        item['status'] = 'scored' if item['matched']['scored_parents'] else 'blocked_no_scorable_parents'
    matched = sorted(set(r['parent_id'] for r in outputs['longest_cds']['matched']['rows']) &
                     set(r['parent_id'] for r in outputs['all_cds']['matched']['rows']))
    return {'schema_version': 'multi-cds-results-1', 'panel_sha256': panel['panel_sha256'],
            'conditions': outputs, 'common_scored_parent_ids': matched, 'matched_per_unit_control': controls,
            'matched_per_unit_abstentions':control_abstentions,
            'claim_scope': 'Distinct documented-family recovery per parent; more represented CDS units is not higher per-unit accuracy'}


def compare_multi_cds_conditions(result, panel, config):
    """Paired parent uncertainty; saturated catalogues cannot establish a gain."""
    comparisons = {}
    clusters = {pid:p['independence_cluster'] for pid,p in panel['parents'].items()}
    saturated = {pid for pid,p in panel['parents'].items() if p['budget_saturated']}
    for subset in ('full', 'matched'):
        endpoints = {}
        for metric in ('documented_family_recall_at_k', 'all_documented_families_recovered'):
            values = [{r['parent_id']:r[metric] for r in result['conditions'][condition][subset]['rows']
                       if r['parent_id'] not in saturated} for condition in ('longest_cds', 'all_cds')]
            endpoints[metric] = paired_cluster_bootstrap(*values, clusters, seed=config['seed'],
                replicates=config['bootstrap_replicates'], min_clusters=config['min_independent_clusters'])
        comparisons[subset] = endpoints
    return {'direction':'all-CDS minus longest-CDS', 'endpoints':comparisons,
            'excluded_budget_saturated_parent_ids':sorted(saturated),
            'limitation':'Equal-weight reviewed parent clusters; family coverage gains do not establish per-unit accuracy gains'}


def produce_multi_cds_scores(method, proteins, panel, config, *, device='cpu',
                             allow_download=False, adapter_factory=LocalAdapter,
                             freeze_model=None):
    """Produce one score per eligible CDS/reference pair; never change the panel.

    ``freeze_model(metadata)`` is invoked after loading but before a forward pass.
    The producer is injectable for numerical tests and E3/E4 interoperability.
    """
    for iid, frozen in panel['input_fingerprints'].items():
        if iid not in proteins or proteins[iid].provenance() != {k:frozen[k] for k in proteins[iid].provenance()}:
            raise InputError('E2 protein/source fingerprint changed after panel freeze')
    baseline = sorted({p['baseline_unit_id'] for p in panel['parents'].values() if p['baseline_unit_id']})
    units = sorted({u for p in panel['parents'].values() for u in p['eligible_unit_ids']})
    references = sorted({r for p in panel['parents'].values() for r in p['eligible_reference_ids']})
    allowed = {u:p['eligible_reference_ids'] for p in panel['parents'].values() for u in p['eligible_unit_ids']}
    if method == 'alignment':
        from .pilot_alignment import alignment_scores
        result = alignment_scores(proteins, {'query_ids':units,'eligible_reference_ids':allowed},
                                  max_residues=config['max_residues'])
        return {'unit_scores':result['raw_scores'],'unit_failures':result['failures'],
                'method_metadata':result['configuration'],'stage_resources':{'alignment':result['resources']},
                'costs':{'cold_load_seconds':None,'inference_pooling_seconds':result['resources']['wall_seconds'],
                    'total_wall_seconds':result['resources']['wall_seconds'],'model_calls':None,
                    'alignment_calls':result['pairwise_score_calls']+result['reference_self_score_calls'],
                    'unique_encoded_sequences':None,'processed_residues':None,'cpu_seconds':None,'gpu_seconds':None,
                    'peak_process_rss_bytes':result['resources']['sampled_peak_process_rss_bytes'],
                    'peak_cuda_allocated_bytes':None,'cache_hits':None},
                'inference_executed':False,'physical_execution_note':'Raw BLOSUM62 local alignment; no protein encoder calls'}
    adapter, vectors, cache, events, failures, stages = None, {}, {}, [], [], {}
    resolved, gpu_name, gpu_total_memory = None, None, None
    load = Measure()
    metadata = {'model':method,'requested_revision':config.get('model_revisions',{}).get(method)}
    with Measure() as total:
        try:
            resolved = prepare_device(device, config['seed'])
            if resolved.startswith('cuda'):
                import torch
                gpu_name = torch.cuda.get_device_name(resolved)
                gpu_total_memory = torch.cuda.get_device_properties(resolved).total_memory
            load = Measure(resolved)
            with load:
                revision = config.get('model_revisions',{}).get(method, ESM2_REVISION if method=='esm2' else 'main')
                adapter = adapter_factory(method,resolved,revision,config['max_residues'],offline=not allow_download)
            metadata = adapter.metadata
            if freeze_model is not None:
                freeze_model(metadata)
            stage_ids = [('references',references),('baseline_cds',baseline),
                         ('additional_cds',sorted(set(units)-set(baseline)))]
            for stage, ids in stage_ids:
                with Measure(resolved) as measured:
                    for iid in ids:
                        protein = proteins[iid]
                        event = {'input_id':iid,'stage':stage,'length_aa':len(protein.sequence),
                                 'cache_hit':protein.sequence_sha256 in cache,'status':'attempted'}
                        events.append(event)
                        try:
                            if protein.sequence_sha256 not in cache:
                                cache[protein.sequence_sha256] = mean_pool(adapter.residues(protein.sequence))
                            vectors[iid] = cache[protein.sequence_sha256]
                            event['status'] = 'ok'
                        except Exception as exc:
                            event.update(status='failed',**failure(exc))
                            failures.append({'input_id':iid,**failure(exc)})
                stages[stage] = measured.data
        except Exception as exc:
            failures.append({'stage':'model_loading_or_setup',**failure(exc)})
        finally:
            if adapter is not None:
                adapter.close()
    scores = {}
    for uid in units:
        refs = allowed[uid]
        if uid not in vectors or any(r not in vectors for r in refs):
            continue
        try:
            scores[uid] = {r:cosine(vectors[uid],vectors[r]) for r in refs}
        except InputError as exc:
            failures.append({'input_id':uid,**failure(exc)})
    successful = [e for e in events if e['status']=='ok' and not e['cache_hit']]
    costs = {'cold_load_seconds':load.data.get('wall_seconds'),
             'inference_pooling_seconds':sum(s['wall_seconds'] for s in stages.values()) if stages else None,
             'total_wall_seconds':total.data['wall_seconds'],'model_calls':len(successful),
             'attempted_model_calls':sum(not e['cache_hit'] for e in events),
             'alignment_calls':None,'unique_encoded_sequences':len(successful),
             'processed_residues':sum(e['length_aa'] for e in successful),
             'cpu_seconds':None,'gpu_seconds':None,
             'peak_process_rss_bytes':total.data.get('sampled_peak_process_rss_bytes'),
             'peak_cuda_allocated_bytes':max((s.get('peak_cuda_allocated_bytes',0) for s in stages.values()),default=0)
                                       if resolved is not None and resolved.startswith('cuda') else None,
             'cache_hits':sum(e['status']=='ok' and e['cache_hit'] for e in events)}
    return {'unit_scores':scores,'unit_failures':failures,'method_metadata':metadata,
            'costs':costs,'stage_resources':{'cold_load':load.data,**stages},'events':events,
            'hardware':{'device':resolved,'gpu':gpu_name,'gpu_total_memory_bytes':gpu_total_memory},
            'inference_executed':bool(successful),
            'physical_execution_note':'References once, baseline units, then extra CDSs; exact-sequence cache shared, no separate-condition reruns'}


def _cost_records(result, execution, panel, config, inputs, provenance, environment, method, run_id):
    records = []
    kinds = inputs.get('evidence_kind','unclassified')
    if kinds not in ('admitted_research','fictional_software_fixture','unclassified'):
        raise InputError('Unknown E2 evidence_kind declaration')
    runtime_policy = {'cache':'exact-sequence-within-run','precision':'float32',
                      'execution':'reference then baseline then additional CDS sequential shared pass'}
    for condition,item in result['conditions'].items():
        summary = item['matched']
        condition_units = {u for p in panel['parents'].values() for u in
            ([p['baseline_unit_id']] if condition=='longest_cds' and p['baseline_unit_id'] else p['eligible_unit_ids'])}
        scored_ids = {row['parent_id'] for row in summary['rows']}
        clusters = {panel['parents'][p]['independence_cluster'] for p in scored_ids}
        # Shared cache work is measured once. Do not copy the combined measured
        # cost onto both conditions as if two separate model runs occurred.
        costs = dict(execution['costs']) if condition=='all_cds' else {key:None for key in execution['costs']}
        if condition=='longest_cds' and execution['stage_resources'].get('baseline_cds'):
            stages = execution['stage_resources']
            costs['cold_load_seconds'] = execution['costs'].get('cold_load_seconds')
            costs['inference_pooling_seconds'] = sum(stages[s]['wall_seconds'] for s in ('references','baseline_cds'))
            events = [e for e in execution.get('events',[]) if e['stage'] in ('references','baseline_cds')]
            successes = [e for e in events if e['status']=='ok' and not e['cache_hit']]
            costs.update(model_calls=len(successes),unique_encoded_sequences=len(successes),
                         processed_residues=sum(e['length_aa'] for e in successes),
                         cache_hits=sum(e['cache_hit'] and e['status']=='ok' for e in events))
        record = {'schema_version':'cost-quality-record-1','run_id':run_id,'method':method+'/'+condition,
            'evidence_kind':kinds,'status':item['status'],'actual_measurements':bool(execution.get('stage_resources')),
            'source_admission':{'all_sources_approved':True,'review_manifest_sha256':provenance['review_manifest_sha256']},
            'panel_sha256':panel['panel_sha256'],'protocol_sha256':panel['protocol_sha256'],
            'input_config_sha256':provenance['input_config_sha256'],'budget':5,
            'budget_scope':'five total references per parent','comparison_task':'E2_parent_family_recovery',
            'reference_tier':'same-parent-all-component-exact-copy-reviewed-group-excluded',
            'runtime_policy_sha256':_hash(runtime_policy),'cache_policy':runtime_policy['cache'],'precision':'float32',
            'cache_state':'fresh within-run cache, one shared reference/baseline/additional execution',
            'shared_execution_id':run_id+'/'+method,
            'scored_observation_ids':sorted(scored_ids),
            'method_configuration_sha256':_hash({'method':method,'metadata':execution.get('method_metadata',{})}),
            'model_identity':execution.get('method_metadata',{}),
            'hardware':{**{k:environment.get(k) for k in ('platform','cpu','logical_cpu_count')},
                        **execution.get('hardware',{}),
                        'device':execution.get('hardware',{}).get('device',execution.get('method_metadata',{}).get('device','unknown'))},
            'software':{k:environment.get(k) for k in ('python','packages')},
            'metrics':{'coverage':summary['coverage'],
                'documented_family_recall_at_k':summary['documented_family_recall_at_k'],
                'complete_family_recovery_rate':summary['complete_family_recovery_rate'],
                'independent_groups':len(clusters),'scored_observations':summary['scored_parents'],
                'requested_observations':summary['requested_parents'],
                'budget_saturated_observations':sum(panel['parents'][p]['budget_saturated'] for p in scored_ids)},
            'costs':costs,'logical_query_units':len(condition_units),
            'cost_attribution':'Shared measured execution; longest condition has measured reference+baseline stage subset, not a cold rerun'}
        records.append(record)
    return records


def run_multi_cds(inputs_path, protocol_path, output_path, *, device='cpu', allow_download=False,
                  adapter_factory=LocalAdapter, score_producer=produce_multi_cds_scores):
    """Execute E2 on a reviewed, pre-frozen parent panel; empty data stays blocked."""
    from .extension_config import load_extension_config
    output = local_output(Path(output_path))
    if output.exists() and any(output.iterdir()):
        raise InputError('E2 output directory must be new or empty')
    output.mkdir(parents=True,exist_ok=True)
    report = run_metadata('E2_multi_cds',{'device':device,'allow_download':allow_download})
    report.update(status='blocked',inference_executed=False,scientific_scores_produced=False,methods={})
    cost_records = []
    with Measure() as measured:
        try:
            config = load_extension_config(protocol_path,expected='E2')
            report['protocol'] = config
            report['protocol_sha256'] = _hash(config)
            write_json(output/'protocol-freeze.json',{'protocol':config,'protocol_sha256':_hash(config),
                'protocol_file_sha256':digest_file(Path(protocol_path)),'stage':'before source admission and any scoring'})
            raw = read_json(Path(inputs_path))
            if not raw.get('records'):
                report.update(status='blocked_no_reviewed_inputs',reason='No exact-record-admitted source inputs; shortlist approval does not admit raw data')
            else:
                inputs,proteins,labels,parents,provenance = load_multi_cds_inputs(inputs_path,max_residues=config['max_residues'])
                panel = freeze_multi_cds_panel(config,proteins,labels,parents,inputs['reference_ids'])
                report.update(provenance=provenance,panel=panel)
                write_json(output/'panel-freeze.json',panel)
                if panel['status'] != 'frozen':
                    report['status'] = 'blocked_insufficient_parent_panel'
                else:
                    for method in config['models']:
                        try:
                            execution = score_producer(method,proteins,panel,config,device=device,
                                allow_download=allow_download,adapter_factory=adapter_factory,
                                freeze_model=lambda metadata,m=method:write_json(output/('model-freeze-'+m+'.json'),metadata))
                            result = score_multi_cds(execution['unit_scores'],panel,labels)
                            result['paired_condition_comparison'] = compare_multi_cds_conditions(result,panel,config)
                            expected = sum(p['matched_abstention_reason'] is None for p in panel['parents'].values())
                            counts = [c['matched']['scored_parents'] for c in result['conditions'].values()]
                            method_status = 'ok' if counts and all(c==expected and c>0 for c in counts) else 'partial' if any(counts) else 'blocked_no_scorable_parents'
                            report['methods'][method] = {'status':method_status,'results':result,
                                **{k:v for k,v in execution.items() if k != 'unit_scores'}}
                            report['inference_executed'] |= execution.get('inference_executed',False)
                            cost_records.extend(_cost_records(result,execution,panel,config,inputs,provenance,report,method,output.name))
                        except Exception as exc:
                            report['methods'][method] = {'status':'blocked',**failure(exc)}
                    counts = [condition.get('matched',{}).get('scored_parents',0)
                              for m in report['methods'].values()
                              for condition in m.get('results',{}).get('conditions',{}).values()]
                    report['scores_produced'] = any(counts)
                    report['scientific_scores_produced'] = any(counts) and inputs.get('evidence_kind')=='admitted_research'
                    report['status'] = ('ok' if report['methods'] and all(m['status']=='ok' for m in report['methods'].values())
                                        else 'partial' if any(counts) else 'blocked_no_scores')
                    if any(counts) and all(p['budget_saturated'] for p in panel['parents'].values() if p['matched_abstention_reason'] is None):
                        report['status'] = 'descriptive_only_budget_saturated_panel'
        except Exception as exc:
            report.update(status='blocked',**failure(exc))
    report['cost_quality_records'] = cost_records
    report['resources_total'] = measured.data
    report['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    report['limitations'].extend(['One-CDS reduction deliberately loses represented-family coverage; parent gains are not per-unit accuracy gains.',
        'Matched supported families and full documented expectations have separate frozen denominators; unknown families remain unscored.',
        'All-CDS inference reuses baseline/reference cache; only actual observed stage costs are attributed, with null unmeasured costs.'])
    write_json(output/'cost-quality-records.json',{'schema_version':'cost-quality-records-1','records':cost_records})
    write_json(output/'report.json',report)
    return report
