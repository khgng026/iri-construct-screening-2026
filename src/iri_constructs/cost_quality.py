"""Measured quality/compute records with explicit comparability and missingness.

This module aggregates recorded observations only. It never runs a model,
approves a source, fills unmeasured costs with zero or selects a winning model.
Legacy records lacking explicit evidence/admission context remain descriptive.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

from .runtime import InputError, local_output, require_text, run_metadata, write_json

METRICS = ('coverage', 'documented_family_recall_at_k', 'complete_family_recovery_rate',
           'independent_groups', 'scored_observations', 'requested_observations',
           'false_family_hit_fraction', 'budget_saturated_observations')
COUNTS = {'independent_groups', 'scored_observations', 'requested_observations', 'budget_saturated_observations',
          'model_calls', 'alignment_calls', 'unique_encoded_sequences', 'processed_residues',
          'peak_process_rss_bytes', 'peak_cuda_allocated_bytes', 'cache_hits'}
COSTS = ('cold_load_seconds', 'inference_pooling_seconds', 'total_wall_seconds',
         'model_calls', 'alignment_calls', 'unique_encoded_sequences', 'processed_residues',
         'cpu_seconds', 'gpu_seconds', 'peak_process_rss_bytes', 'peak_cuda_allocated_bytes', 'cache_hits')
IDENTITY = ('panel_sha256', 'protocol_sha256', 'input_config_sha256', 'budget', 'budget_scope',
            'comparison_task', 'reference_tier', 'runtime_policy_sha256', 'cache_policy', 'cache_state',
            'precision', 'hardware', 'software')
HASH_FIELDS = {'panel_sha256', 'protocol_sha256', 'input_config_sha256', 'runtime_policy_sha256',
               'method_configuration_sha256'}
EVIDENCE = {'admitted_research', 'fictional_software_fixture', 'unclassified'}


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode('utf-8')).hexdigest()


def _number(obj, key):
    value = obj.get(key)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int,float)) or not math.isfinite(value) or value < 0:
        raise InputError('Recorded metrics/costs must be finite nonnegative numbers or explicit null: ' + key)
    if key in COUNTS and type(value) is not int:
        raise InputError('Recorded observation/work/memory counts must be integers: ' + key)
    if key in METRICS and key not in COUNTS and value > 1:
        raise InputError('Recorded quality fractions must be in 0..1: ' + key)
    return value


def _valid_hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def records_from_report(value):
    """Read explicit normalized records, or conservatively expose legacy A-D rows.

    The legacy adapter does not turn unclassified historical records into admitted
    measurements. Missing cache/model identity and resource fields remain null.
    """
    if not isinstance(value, dict):
        raise InputError('Measurement input must be a JSON object')
    if value.get('schema_version') == 'cost-quality-record-1':
        return [value]
    if 'records' in value:
        if not isinstance(value['records'], list):
            raise InputError('Measurement records must be an explicit list')
        return value['records']
    methods = value.get('methods')
    if not isinstance(methods, dict):
        return []
    result = []
    for name, method in methods.items():
        if not isinstance(method, dict):
            raise InputError('Legacy method results must be objects')
        model = value.get('models', {}).get(name.split('/')[0], {})
        resources = method.get('resources', {})
        panel = value.get('panel', {})
        config = value.get('protocol', value.get('configuration', {}))
        provenance = value.get('provenance', {})
        source_rows = provenance.get('sources', [])
        approved = bool(source_rows) and all(all(s.get(k) == 'approved' for k in
            ('review_status','licensing_status','safety_review','annotation_status')) for s in source_rows)
        rows = method.get('rows', [])
        overall = method.get('overall', {})
        result.append({'schema_version': 'cost-quality-record-1',
            'run_id': value.get('run_id', value.get('started_at_utc', 'legacy-unspecified')),
            'method': name, 'status': method.get('status', value.get('status', 'unknown')),
            'evidence_kind': value.get('evidence_kind', 'unclassified'),
            'actual_measurements': value.get('actual_measurements') is True,
            'source_admission': {'all_sources_approved': approved,
                'review_manifest_sha256': provenance.get('review_manifest_sha256')},
            'panel_sha256': method.get('panel_sha256', panel.get('panel_sha256')),
            'protocol_sha256': panel.get('protocol_sha256'), 'input_config_sha256': provenance.get('input_config_sha256'),
            'budget': panel.get('retrieval_budget', config.get('retrieval_budget')),
            'budget_scope': config.get('budget_scope'), 'comparison_task': value.get('experiment'),
            'reference_tier': panel.get('reference_tier', 'original_t07'),
            'runtime_policy_sha256': value.get('runtime_policy_sha256'),
            'cache_policy': value.get('cache_policy'), 'cache_state': value.get('cache_state'),
            'precision': model.get('metadata', {}).get('precision_policy'),
            'hardware': {k:value[k] for k in ('platform','cpu','logical_cpu_count','resolved_device') if k in value},
            'software': {'python':value.get('python'), 'packages':value.get('packages')},
            'metrics': {'coverage': method.get('coverage'),
                'documented_family_recall_at_k': overall.get('documented_family_recall_at_k'),
                'complete_family_recovery_rate': overall.get('all_documented_families_recovered_rate'),
                'independent_groups': len({r['independence_group'] for r in rows if 'independence_group' in r}) if rows else None,
                'scored_observations': method.get('scored_queries'), 'requested_observations': method.get('requested_queries')},
            'scored_observation_ids': [r['query_id'] for r in rows if 'query_id' in r],
            'costs': {'cold_load_seconds':model.get('loading_resources', {}).get('wall_seconds'),
                'total_wall_seconds':resources.get('wall_seconds'),
                'model_calls': model.get('unique_forward_calls'),
                'alignment_calls': method.get('pairwise_score_calls'),
                'peak_process_rss_bytes':resources.get('sampled_peak_process_rss_bytes'),
                'peak_cuda_allocated_bytes':resources.get('peak_cuda_allocated_bytes')},
            'legacy_adapter_note':'Unrecorded context/costs remain null; no source or measurement approval is inferred'})
    return result


def _normalize(record, min_independent_groups):
    if not isinstance(record, dict) or record.get('schema_version') != 'cost-quality-record-1':
        raise InputError('Expected explicit cost-quality-record-1 per-method measurement record')
    row = {'run_id':require_text(record,'run_id'), 'method':require_text(record,'method'),
           'status':record.get('status', 'unknown'), 'evidence_kind':record.get('evidence_kind','unclassified'),
           'actual_measurements': record.get('actual_measurements') is True}
    if row['evidence_kind'] not in EVIDENCE:
        raise InputError('Unknown evidence_kind; distinguish research, fixtures and unclassified records')
    metric_values, costs = record.get('metrics', {}), record.get('costs', {})
    if not isinstance(metric_values, dict) or not isinstance(costs, dict):
        raise InputError('Recorded metrics and costs must be objects')
    row['metrics'] = {key:_number(metric_values,key) for key in METRICS}
    row['costs'] = {key:_number(costs,key) for key in COSTS}
    for key in IDENTITY:
        row[key] = record.get(key)
    for key in ('method_configuration_sha256','model_identity','shared_execution_id'):
        row[key] = record.get(key)
    row['source_admission'] = record.get('source_admission', {})
    if not isinstance(row['source_admission'], dict):
        raise InputError('Source admission evidence must be an object')
    metrics = row['metrics']
    scored, requested = metrics['scored_observations'], metrics['requested_observations']
    if scored is not None and requested is not None and scored > requested:
        raise InputError('Scored observations exceed the recorded requested denominator')
    if scored == 0 and any(metrics[k] is not None for k in
        ('documented_family_recall_at_k','complete_family_recovery_rate')):
        raise InputError('No-observation agreement metrics must remain null')
    ids = record.get('scored_observation_ids')
    if ids is not None:
        if (not isinstance(ids,list) or any(not isinstance(i,str) or not i for i in ids)
                or len(ids) != len(set(ids)) or scored is not None and len(ids) != scored):
            raise InputError('Scored observation IDs must be unique and match the recorded count')
        support = sorted(ids)
    elif metrics['coverage'] == 1 and scored is not None and scored == requested and scored > 0:
        support = 'all_frozen_panel_queries'
    else:
        support = None
    row['scored_observation_ids'] = sorted(ids) if ids is not None else None
    row['support_fingerprint'] = _hash(support) if support is not None else None
    reasons = []
    if row['evidence_kind'] != 'admitted_research':
        reasons.append('Not an explicitly admitted research measurement')
    if not row['actual_measurements']:
        reasons.append('Actual method execution/measurement is not attested')
    if row['status'] not in ('ok','completed','scored','partial','descriptive_only_budget_saturated_panel'):
        reasons.append('Method record is not a completed or measured partial execution')
    admission = row['source_admission']
    if admission.get('all_sources_approved') is not True or not _valid_hash(admission.get('review_manifest_sha256')):
        reasons.append('Exact source review/checksum admission evidence is missing')
    for key in HASH_FIELDS:
        if not _valid_hash(row.get(key)):
            reasons.append('Missing immutable identity: ' + key)
    for key in IDENTITY:
        value = row[key]
        if key in HASH_FIELDS:
            continue
        if key == 'budget':
            if type(value) is not int or value < 1:
                reasons.append('Missing/invalid total return budget')
        elif key in ('hardware','software'):
            if not isinstance(value,dict) or not value or not any(v is not None for v in value.values()):
                reasons.append('Missing measured context: ' + key)
        elif not isinstance(value,str) or not value.strip():
            reasons.append('Missing measured context: ' + key)
    hardware = row['hardware'] if isinstance(row['hardware'],dict) else {}
    if not hardware.get('device') or not hardware.get('platform') or not hardware.get('cpu'):
        reasons.append('Device/platform/CPU identity is incomplete')
    if str(hardware.get('device','')).startswith('cuda') and not hardware.get('gpu'):
        reasons.append('GPU identity is unavailable')
    software = row['software'] if isinstance(row['software'],dict) else {}
    if not software.get('python'):
        reasons.append('Python/runtime version is unavailable')
    if support is None:
        reasons.append('Exact common scored observation identities are unavailable')
    if not scored or metrics['documented_family_recall_at_k'] is None or metrics['coverage'] is None:
        reasons.append('No measured quality endpoint with eligible denominators')
    if metrics['independent_groups'] is None or metrics['independent_groups'] < min_independent_groups:
        reasons.append('Insufficient independently reviewed groups for quality comparison')
    if metrics['budget_saturated_observations'] and scored and metrics['budget_saturated_observations'] >= scored:
        reasons.append('Returning the entire eligible catalogue is descriptive, not a method-quality gain')
    row['comparison_exclusions'] = reasons
    row['direct_comparison_eligible'] = not reasons
    row['record_sha256'] = _hash(record)
    if not reasons:
        context = {key:row[key] for key in IDENTITY}
        context.update(support_fingerprint=row['support_fingerprint'],
                       review_manifest_sha256=admission['review_manifest_sha256'])
        row['comparison_group'] = _hash(context)
    else:
        row['comparison_group'] = None
    # Throughput is reported only with actual positive measured time and numerator.
    wall = costs.get('inference_pooling_seconds')
    row['scored_observations_per_inference_second'] = scored/wall if row['actual_measurements'] and scored and isinstance(wall,(int,float)) and wall > 0 else None
    residues = costs.get('processed_residues')
    row['processed_residues_per_inference_second'] = residues/wall if row['actual_measurements'] and residues is not None and isinstance(wall,(int,float)) and wall > 0 else None
    return row


def summarize_cost_quality(records, *, min_independent_groups=8):
    """Retain all actual/null fields; reject incompatible rows as direct comparisons."""
    if not isinstance(records,list) or type(min_independent_groups) is not int or min_independent_groups < 2:
        raise InputError('Cost/quality aggregation requires a record list and explicit cluster threshold')
    rows, seen = [], set()
    for record in records:
        row = _normalize(record,min_independent_groups)
        key = (row['run_id'],row['method'],row['reference_tier'])
        if key in seen:
            raise InputError('Duplicate/conflicting measurement identity; repeated files cannot inflate the summary')
        seen.add(key)
        rows.append(row)
    groups = {}
    for index,row in enumerate(rows):
        if row['comparison_group'] is not None:
            groups.setdefault(row['comparison_group'],[]).append(index)
    # A shared display label must not silently pool two different checkpoints or
    # method configurations alongside another comparator. Exclude ambiguous rows
    # before forming any pair/group/plot, while preserving their descriptive data.
    for indices in groups.values():
        variants = {}
        for index in indices:
            row = rows[index]
            variants.setdefault(row['method'],set()).add(row['method_configuration_sha256'])
        ambiguous = {name for name,values in variants.items() if len(values)>1}
        for index in indices:
            if rows[index]['method'] in ambiguous:
                rows[index]['comparison_exclusions'].append('Ambiguous same-method label pools different underlying model/configuration variants')
                rows[index]['direct_comparison_eligible'] = False
                rows[index]['comparison_group'] = None
    groups = {}
    for index,row in enumerate(rows):
        if row['comparison_group'] is not None:
            groups.setdefault(row['comparison_group'],[]).append(index)
    pair_checks = []
    for a,b in itertools.combinations(range(len(rows)),2):
        left,right = rows[a],rows[b]
        changed = [key for key in IDENTITY if left[key] != right[key]]
        if left['support_fingerprint'] != right['support_fingerprint']:
            changed.append('scored_observation_ids')
        if left['source_admission'].get('review_manifest_sha256') != right['source_admission'].get('review_manifest_sha256'):
            changed.append('review_manifest_sha256')
        if left['method'] == right['method'] and left['method_configuration_sha256'] != right['method_configuration_sha256']:
            changed.append('same method label uses a different underlying method/model configuration')
        accepted = left['direct_comparison_eligible'] and right['direct_comparison_eligible'] and not changed
        pair_checks.append({'left_row':a,'right_row':b,
            'status':'matched_measured_context' if accepted else 'rejected_direct_comparison',
            'mismatched_fields':changed,
            'incomplete_or_unadmitted_rows':[i for i in (a,b) if not rows[i]['direct_comparison_eligible']]})
    group_records = []
    for group_id,indices in sorted(groups.items()):
        accepted_pairs=[p for p in pair_checks if p['left_row'] in indices and p['right_row'] in indices
                        and p['status']=='matched_measured_context']
        methods=sorted({rows[i]['method'] for i in indices})
        group_records.append({'comparison_group':group_id,'row_indices':indices,
            'status':'matched_methods_available' if accepted_pairs and len(methods)>1 else 'single_or_incompatible_method_only',
            'methods':methods,'matched_pair_count':len(accepted_pairs)})
    actual = [r for r in rows if r['evidence_kind']=='admitted_research' and r['actual_measurements']
              and r['source_admission'].get('all_sources_approved') is True
              and _valid_hash(r['source_admission'].get('review_manifest_sha256'))]
    status = ('blocked_no_measurements' if not rows else 'descriptive_only_no_admitted_research'
              if not actual else 'ok' if any(p['status']=='matched_measured_context' for p in pair_checks)
              else 'descriptive_only_no_matched_comparison')
    return {'schema_version':'cost-quality-summary-1','status':status,'rows':rows,
        'input_records':len(rows),'actual_admitted_research_records':len(actual),
        'direct_comparison_groups':group_records,'pair_comparability':pair_checks,
        'min_independent_groups':min_independent_groups,'automatic_winner':None,
        'shared_execution_note':'Rows reusing physical inference/alignment across tiers or conditions are never summed; inspect shared_execution_id.',
        'limitations':['Tables describe recorded observations only; missing metrics/costs remain null.',
            'Matched hardware/software/protocol/cache/panel context does not establish causal performance differences.',
            'No model winner or cost benefit is inferred automatically.']}


def plot_cost_quality(summary, output_path, *, comparison_group=None, _allow_fictional_test_fixture=False):
    """Optional standalone matplotlib plot; production plots need admitted measurements."""
    output_path=Path(output_path)
    local_output(output_path.parent)
    if output_path.suffix.lower() not in ('.svg','.png'):
        raise InputError('Cost/quality plot must be a standalone SVG or PNG')
    rows=summary.get('rows',[])
    if _allow_fictional_test_fixture:
        selected=[r for r in rows if r['evidence_kind']=='fictional_software_fixture']
        title='FICTIONAL SOFTWARE FIXTURE - NOT RESEARCH'
    else:
        groups=[g for g in summary.get('direct_comparison_groups',[]) if g['status']=='matched_methods_available']
        if comparison_group is not None:
            groups=[g for g in groups if g['comparison_group']==comparison_group]
        if not groups:
            return {'status':'blocked_no_comparable_admitted_measurements','path':None}
        if len(groups)!=1:
            return {'status':'blocked_select_one_comparison_group','path':None,
                    'comparison_groups':[g['comparison_group'] for g in groups]}
        selected=[rows[i] for i in groups[0]['row_indices'] if rows[i]['direct_comparison_eligible']]
        title='Exploratory measured quality / compute on one frozen matched panel'
    points=[r for r in selected if r['costs'].get('inference_pooling_seconds') is not None
            and r['costs']['inference_pooling_seconds']>0 and r['metrics'].get('documented_family_recall_at_k') is not None]
    if len(points)<2:
        return {'status':'blocked_missing_measured_cost_quality_points','path':None}
    try:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
    except ImportError:
        return {'status':'unavailable_plot_dependency','dependency':'matplotlib','path':None}
    figure,axis=plt.subplots(figsize=(8,5),layout='constrained')
    try:
        for row in points:
            x=row['costs']['inference_pooling_seconds']
            y=row['metrics']['documented_family_recall_at_k']
            axis.scatter([x],[y])
            axis.annotate(row['method'],(x,y),xytext=(5,5),textcoords='offset points',fontsize=8)
        axis.set(xlabel='Measured inference / pooling wall seconds',
                 ylabel='Documented-family recall at fixed total K',title=title,ylim=(-.02,1.05))
        axis.grid(alpha=.25)
        output_path.parent.mkdir(parents=True,exist_ok=True)
        figure.savefig(output_path)
    finally:
        plt.close(figure)
    return {'status':'written_fictional_software_fixture' if _allow_fictional_test_fixture else 'written_admitted_measurements',
            'path':output_path.name,'point_count':len(points),'automatic_winner':None}


def write_cost_quality_report(records, output_path, *, protocol=None):
    """Write E4 JSON/table, and a plot only when admitted matching observations exist."""
    output=local_output(Path(output_path))
    if output.exists() and any(output.iterdir()):
        raise InputError('Cost/quality output must be a new empty ignored directory')
    output.mkdir(parents=True,exist_ok=True)
    if protocol is not None:
        from .extension_config import validate_extension_config
        validate_extension_config(protocol,expected='E4')
        write_json(output/'protocol-freeze.json',{'protocol':protocol,'protocol_sha256':_hash(protocol),
            'stage':'before aggregating supplied recorded measurements'})
    minimum=(protocol or {}).get('min_independent_clusters',8)
    summary=summarize_cost_quality(records,min_independent_groups=minimum)
    summary['plot']=plot_cost_quality(summary,output/'cost-quality.svg')
    write_json(output/'cost-quality.json',summary)
    columns=['run_id','method','reference_tier','evidence_kind','status','actual_measurements','direct_comparison_eligible',
             *METRICS,*COSTS,'scored_observations_per_inference_second','processed_residues_per_inference_second']
    with (output/'cost-quality.csv').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns)
        writer.writeheader()
        for row in summary['rows']:
            flat={key:row.get(key) for key in columns}
            flat.update(row['metrics']);flat.update(row['costs'])
            writer.writerow({key:('NULL' if flat.get(key) is None else flat[key]) for key in columns})
    report=run_metadata('T07_v1.1_E4',protocol or {})
    report.update(status=summary['status'],inference_executed=False,scientific_scores_produced=False,
                  actual_admitted_research_records=summary['actual_admitted_research_records'],
                  cost_quality=summary,plot=summary['plot'])
    write_json(output/'report.json',report)
    return report
