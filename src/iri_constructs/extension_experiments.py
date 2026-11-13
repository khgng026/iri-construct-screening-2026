"""T07 v1.1 shared execution and tier integration; no source admission decisions."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

from .bench_inputs import load_proteins, reviewed_sources
from .config import configure_logging
from .extension_config import load_extension_config
from .metrics import cosine, reviewed_labels
from .model_adapters import ESM2_REVISION, LocalAdapter, prepare_device
from .pilot_experiments import _CachedAdapter, _compare_methods, _vectors
from .runtime import InputError, Measure, digest_file, failure, local_output, read_json, run_metadata, write_json


def _hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',', ':'),allow_nan=False).encode()).hexdigest()


def _policy(inputs, key, inputs_path):
    value = inputs.get(key)
    if value is None:
        return None
    if not isinstance(value,str) or not value.strip():
        raise InputError(key + ' must reference a local reviewed policy JSON file')
    return read_json((Path(inputs_path).resolve().parent/value).resolve())


def _tier_scores(scores, panel):
    return {q:{r:value for r,value in row.items() if r in panel['eligible_reference_ids'][q]}
            for q,row in scores.items() if q in panel['query_ids'] and q not in panel['query_exclusions']}


def _vector_scores(vectors, panel):
    scores = {}
    for q in panel['query_ids']:
        if q not in vectors:
            continue
        row = {}
        for r in panel['eligible_reference_ids'][q]:
            if r in vectors:
                try:
                    row[r] = cosine(vectors[q],vectors[r])
                except InputError:
                    pass
        scores[q] = row
    return scores


def _cost_record(report, method, tier_name, panel, result, costs, config):
    sources = report['provenance']['sources']
    all_approved = all(all(s[k]=='approved' for k in ('review_status','licensing_status','safety_review','annotation_status')) for s in sources)
    hardware = {'platform':report['platform'],'cpu':report['cpu'],'logical_cpu_count':report['logical_cpu_count'],
                'device':'cpu' if method.startswith('alignment') else report['models'][method]['resolved_device']}
    if str(hardware['device']).startswith('cuda'):
        hardware['gpu']=report['models'][method].get('gpu')
    runtime_policy = {'scope':'Complete frozen source catalogue once; tier scores reuse outputs, costs are shared and not additive across tiers',
                      'load':'separate cold model initialization; no warm-up', 'cache':'exact-sequence-within-run',
                      'input_cap':config['max_residues'],'precision':'float32 embeddings; native alignment scoring',
                      'resource_scope':'process RSS and process CPU; CUDA allocator peaks when measured'}
    rows = result.get('rows',[])
    identity = ({'implementation':'Bio.Align.PairwiseAligner','method':method,
                 'configuration':report.get('alignment',{}).get('configuration')}
                if method.startswith('alignment') else report['models'][method].get('metadata'))
    return {'schema_version':'cost-quality-record-1','run_id':report['run_id'],'method':method,
        'evidence_kind':report['evidence_kind'],'status':'completed' if rows else 'blocked',
        'actual_measurements':bool(costs.get('model_calls') or costs.get('alignment_calls')),
        'source_admission':{'all_sources_approved':bool(sources) and all_approved,
                            'review_manifest_sha256':report['provenance']['review_manifest_sha256']},
        'panel_sha256':panel['panel_sha256'],'protocol_sha256':panel['protocol_sha256'],
        'input_config_sha256':report['provenance']['input_config_sha256'],
        'budget':5,'budget_scope':'five total hits per matched protein query','comparison_task':'E3 independent-family retrieval',
        'reference_tier':tier_name,'runtime_policy_sha256':_hash(runtime_policy),
        'runtime_policy':runtime_policy,'shared_execution_id':report['run_id']+('/alignment' if method.startswith('alignment') else '/'+method),
        'cache_policy':'exact-sequence-within-run',
        'cache_state':'new-model-instance; empty-sequence-cache; first-pass-no-warmup; cached-weights unless download explicitly allowed',
        'weight_cache_observation':'offline cached files' if not report.get('downloads_allowed') else 'downloads allowed; filesystem/OS cache not experimentally controlled',
        'precision':runtime_policy['precision'],
        'model_identity':identity,'method_configuration_sha256':_hash(identity),
        'scored_observation_ids':sorted(r['query_id'] for r in rows),
        'hardware':hardware,'software':{'python':report['python'],'packages':report['packages']},
        'method_thread_settings':{k:v for k,v in report['models'].get(method,{}).get('metadata',{}).items() if 'threads' in k},
        'metrics':{'coverage':result.get('coverage'),'documented_family_recall_at_k':result.get('overall',{}).get('documented_family_recall_at_k'),
                   'complete_family_recovery_rate':result.get('overall',{}).get('all_documented_families_recovered_rate'),
                   'independent_groups':len({r['independence_cluster'] for r in rows}),
                   'scored_observations':len(rows),'requested_observations':len(panel['query_ids']),
                   'budget_saturated_observations':sum(r.get('budget_saturated',False) for r in rows)},
        'costs':costs}


def run_reference_generalization(inputs_path, protocol_path, output_path, *, device='cpu',
                                 allow_download=False, adapter_factory=LocalAdapter):
    """Freeze every tier before producing method scores; retain all original T07 guards."""
    from .reference_tiers import freeze_reference_tiers, score_reference_tier
    config = load_extension_config(protocol_path,expected='E3')
    output = local_output(Path(output_path))
    if output.exists() and any(output.iterdir()):
        raise InputError('Extension output must be a new empty ignored directory')
    output.mkdir(parents=True,exist_ok=True)
    report = run_metadata('T07_v1.1_E3',config)
    report.update(run_id=output.name,status='blocked',inference_executed=False,
                  scientific_scores_produced=False,evidence_kind='unclassified',tiers={},models={},downloads_allowed=allow_download)
    records = []
    with Measure() as total:
        try:
            write_json(output/'protocol-freeze.json',{'protocol':config,'protocol_file_sha256':digest_file(Path(protocol_path)),
                'stage':'before source admission, tier filtering and all method scoring'})
            inputs = read_json(Path(inputs_path))
            if not inputs.get('records'):
                report.update(status='blocked_no_reviewed_inputs',reason='No T02-admitted local records; approved shortlist identities do not admit source bytes')
            else:
                _, sources, _ = reviewed_sources(Path(inputs_path))
                if any(source.review['annotation_status']!='approved' for source in sources):
                    raise InputError('E3 requires annotation review approval before source parsing')
                inputs,proteins,provenance,excluded = load_proteins(Path(inputs_path))
                report.update(provenance=provenance,input_exclusions=excluded)
                evidence = inputs.get('evidence_kind','unclassified')
                if evidence not in ('admitted_research','fictional_software_fixture','unclassified'):
                    raise InputError('Unknown explicit evidence_kind')
                report['evidence_kind'] = evidence
                labels,label_exclusions = reviewed_labels(inputs,proteins)
                archetypes = {s['record_id']:s['archetype'] for s in provenance['sources']}
                labels = {i:{**l,'archetype':archetypes.get(proteins[i].construct_id,'unknown')} for i,l in labels.items()}
                report['label_exclusions'] = label_exclusions
                tiers = freeze_reference_tiers(config,proteins,labels,inputs.get('query_ids'),inputs.get('reference_ids'),
                    linkage_policy=_policy(inputs,'linkage_policy',inputs_path),
                    relatedness_policy=_policy(inputs,'relatedness_policy',inputs_path))
                write_json(output/'reference-tier-freeze.json',tiers)
                report['reference_tiers'] = tiers
                first = tiers['tiers']['tier1_t07_invariants']['panel']
                report['tiers'] = {name:{'status':tier['status'],'methods':{},'panel_sha256':tier['panel']['panel_sha256'] if tier.get('panel') else None}
                                   for name,tier in tiers['tiers'].items()}
                if first['status']!='frozen':
                    report.update(status='blocked_insufficient_panel',reason='No query has an independent positive and competing family')
                else:
                    report['resolved_device'] = 'cpu'
                    score_sources,work = {},{}
                    if 'alignment' in config['models']:
                        try:
                            from .pilot_alignment import alignment_scores
                            started=time.process_time()
                            with Measure() as measured:
                                alignment=alignment_scores(proteins,first,max_residues=config['max_residues'])
                            common={'cold_load_seconds':None,'inference_pooling_seconds':measured.data['wall_seconds'],
                                'total_wall_seconds':measured.data['wall_seconds'],'model_calls':0,
                                'alignment_calls':alignment['pairwise_score_calls']+alignment['reference_self_score_calls'],
                                'unique_encoded_sequences':0,'processed_residues':None,'cpu_seconds':time.process_time()-started,
                                'gpu_seconds':None,'peak_process_rss_bytes':measured.data['sampled_peak_process_rss_bytes'],
                                'peak_cuda_allocated_bytes':None,'cache_hits':0}
                            score_sources['alignment_raw']=alignment['raw_scores']
                            score_sources['alignment_reference_normalized']=alignment['normalized_scores']
                            work['alignment_raw']=work['alignment_reference_normalized']=common
                            report['alignment']={k:v for k,v in alignment.items() if k not in ('raw_scores','normalized_scores')}
                        except Exception as exc:
                            report['alignment']={'status':'blocked',**failure(exc)}
                    for name in (m for m in config['models'] if m!='alignment'):
                        adapter,cached=None,None
                        info={'status':'blocked','inference_executed':False}
                        report['models'][name]=info
                        load=Measure()
                        try:
                            resolved=prepare_device(device,config['seed'])
                            info['resolved_device']=resolved
                            if resolved.startswith('cuda'):
                                import torch
                                info['gpu']={'name':torch.cuda.get_device_name(resolved),
                                             'total_memory_bytes':torch.cuda.get_device_properties(resolved).total_memory}
                            report['resolved_device']=resolved
                            revision=config['model_revisions'].get(name,ESM2_REVISION if name=='esm2' else 'main')
                            load=Measure(resolved)
                            with load:
                                adapter=adapter_factory(name,resolved,revision,config['max_residues'],offline=not allow_download)
                            info['metadata']=adapter.metadata
                            write_json(output/('model-freeze-'+name+'.json'),adapter.metadata)
                            cached=_CachedAdapter(adapter)
                            cpu_start=time.process_time()
                            with Measure(resolved) as inference:
                                vectors,details=_vectors(cached,proteins,set(first['query_ids'])|set(first['reference_ids']))
                                score_sources[name]=_vector_scores(vectors,first)
                            successful=[event for event in cached.events if event['status']=='ok' and not event['cache_hit']]
                            resources=inference.data
                            work[name]={'cold_load_seconds':load.data['wall_seconds'],
                                'inference_pooling_seconds':resources['wall_seconds'],
                                'total_wall_seconds':load.data['wall_seconds']+resources['wall_seconds'],
                                'model_calls':len(successful),'alignment_calls':0,
                                'unique_encoded_sequences':len(successful),'processed_residues':sum(e['length_aa'] for e in successful),
                                'cpu_seconds':time.process_time()-cpu_start,'gpu_seconds':None,
                                'peak_process_rss_bytes':max(x for x in [resources['sampled_peak_process_rss_bytes'],load.data['sampled_peak_process_rss_bytes']] if x is not None) if any(x is not None for x in [resources['sampled_peak_process_rss_bytes'],load.data['sampled_peak_process_rss_bytes']]) else None,
                                'peak_cuda_allocated_bytes':resources.get('peak_cuda_allocated_bytes'),
                                'cache_hits':sum(item.get('cache_hit',False) for item in details)}
                            info.update(status='attempted',inputs=details,inference_resources=resources)
                        except Exception as exc:
                            info.update(status='blocked',**failure(exc))
                        finally:
                            info['loading_resources']=load.data
                            if cached is not None:
                                success=sum(e['status']=='ok' and not e['cache_hit'] for e in cached.events)
                                info.update(inference_attempted=bool(cached.events),inference_executed=bool(success),forward_requests=cached.events)
                                report['inference_executed'] |= bool(success)
                                cached.cache.clear()
                            if adapter is not None:
                                adapter.close()
                    expected=['alignment_raw','alignment_reference_normalized']+[m for m in config['models'] if m!='alignment']
                    for tier_name,tier in tiers['tiers'].items():
                        tier_report=report['tiers'][tier_name]
                        panel=tier.get('panel')
                        if panel is None:
                            continue
                        for method in expected:
                            if method not in score_sources:
                                tier_report['methods'][method]={'status':'blocked_method_unavailable','rows':[],'scored_queries':0}
                                continue
                            scored=score_reference_tier(_tier_scores(score_sources[method],panel),tier,labels)
                            tier_report['methods'][method]=scored
                            records.append(_cost_record(report,method,tier_name,panel,scored,work[method],config))
                        methods=tier_report['methods']
                        tier_report['paired_comparisons']=_compare_methods(methods,panel,config)
                        common_ids=set(panel['query_ids'])
                        for result in methods.values():
                            common_ids &= {r['query_id'] for r in result.get('rows',[])}
                        tier_report['common_scored_query_ids']=sorted(common_ids)
                        eligible={q for q in panel['query_ids'] if q not in panel['query_exclusions']}
                        any_scores=any(m.get('scored_queries',0) for m in methods.values())
                        tier_report['status']=('blocked_no_scores' if not any_scores else 'descriptive_only_budget_saturated_panel'
                            if eligible and eligible==set(panel['budget_saturated_query_ids']) else 'scored' if common_ids else 'partial_no_common_queries')
                    any_rows=any(m.get('scored_queries',0) for t in report['tiers'].values() for m in t['methods'].values())
                    all_comparable=all(t['status']=='scored' and bool(t.get('common_scored_query_ids')) for t in report['tiers'].values())
                    report.update(scientific_scores_produced=bool(any_rows and evidence=='admitted_research'),
                                  measurement_results_produced=bool(any_rows),status='ok' if all_comparable else 'partial' if any_rows else 'blocked_no_scores')
                    if evidence!='admitted_research' and any_rows:
                        report['status']='completed_software_fixture' if evidence=='fictional_software_fixture' else 'unclassified_measurements'
        except Exception as exc:
            report.update(status='blocked',**failure(exc))
    report['resources_total']=total.data
    report['completed_at_utc']=datetime.now(timezone.utc).isoformat()
    write_json(output/'cost-quality-records.json',{'records':records})
    write_json(output/'report.json',report)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--inputs',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--device',choices=('cpu','cuda','auto'),default='cpu')
    parser.add_argument('--allow-download',action='store_true',help='Model weights only; biological sources are never downloaded')
    parser.add_argument('--measurement-records',type=Path,nargs='*')
    args=parser.parse_args(argv)
    try:
        config=load_extension_config(args.protocol)
    except InputError as exc:
        parser.error(str(exc))
    extension=config['extension']
    if extension!='E4' and args.inputs is None:
        parser.error('--inputs is required for E1-E3, even when it is an explicitly empty pending configuration')
    log=configure_logging('INFO')
    try:
        if extension=='E1':
            from .annotation_quality import run_annotation_quality
            output=local_output(args.output)
            if output.exists() and any(output.iterdir()):
                raise InputError('Extension output must be a new empty ignored directory')
            report=run_metadata('T07_v1.1_E1',config)
            write_json(output/'protocol-freeze.json',{'protocol':config,'protocol_file_sha256':digest_file(args.protocol),
                'stage':'before source admission and annotation-quality audit',
                'frozen_at_utc':datetime.now(timezone.utc).isoformat()})
            result=run_annotation_quality(args.inputs,output/'annotation-quality',max_residues=config['max_residues'])
            report.update(status=result['status'],inference_executed=False,
                          scientific_scores_produced=result.get('scientific_scores_produced',False),
                          evidence_class='source_annotation_agreement_only',quality_report=result)
            write_json(output/'report.json',report)
        elif extension=='E2':
            from .multi_cds import run_multi_cds
            report=run_multi_cds(args.inputs,args.protocol,args.output,device=args.device,allow_download=args.allow_download)
        elif extension=='E3':
            report=run_reference_generalization(args.inputs,args.protocol,args.output,device=args.device,allow_download=args.allow_download)
        else:
            from .cost_quality import write_cost_quality_report
            records=[]
            for path in args.measurement_records or []:
                value=read_json(path)
                batch=value.get('records')
                if not isinstance(batch,list):
                    raise InputError('Measurement files must contain a records list')
                records.extend(batch)
            report=write_cost_quality_report(records,args.output,protocol=config)
    except InputError as exc:
        parser.error(str(exc))
    log.info('event=extension_end extension=%s status=%s',extension,report['status'])
    print('Extension report:',args.output/'report.json')
    return 0 if report['status'] in ('ok','completed','scored') else 1
