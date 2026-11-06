#!/usr/bin/env python3
"""Fair common-query comparisons of checksum-verified completed retrieval jobs.

Read-only source processing: no model inference, alignment or parameter tuning.
Clean sanity checks remain separate from combined-query headline results.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import datetime as dt
import io
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import UnitStore, atomic_bytes, canonical, digest, file_hash
from iri_constructs.session_benchmark import evaluate_query,evaluate_oracle_query,load_suite
from session_catalogue_analysis import SummaryAccumulator,latest_completion,source_arms
from session_independent_audit import verify_unit,validate_ranking


def require(condition,message):
    if not condition:raise ValueError(message)


def read(path):return json.loads(Path(path).read_text())


def primary_arm(mode):
    return {'model':'whole','alignment':'local_alignment','alignment_normalization':'normalized_alignment'}[mode]


def comparison_role(mode):
    return 'exploratory_normalization' if mode=='alignment_normalization' else 'primary'


def display_name(identity):
    if identity['mode']=='alignment':return 'Smith-Waterman raw'
    if identity['mode']=='alignment_normalization':return 'SW reference-self-score (exploratory)'
    config=identity.get('config',{})
    model=config.get('model_id',identity.get('method',{}).get('model_id','protein model')).split('/')[-1]
    layer=config.get('layer','last')
    return ' '.join(map(str,[model,config.get('device','unknown-device'),config.get('dtype','unknown-dtype'),f'layer={layer}']))


def new_bundle():return {'combined':SummaryAccumulator(),'clean':SummaryAccumulator()}


def add_bundle(bundle,evaluation):
    bundle['clean' if evaluation['component_count']==1 else 'combined'].add(evaluation)


def flatten_summaries(job,summaries,population):
    component_rows,query_rows,headlines=[],[],[]
    for arm,bundle in sorted(summaries.items()):
        for scope,accumulator in bundle.items():
            summary=accumulator.result()
            if not summary['query_count']:continue
            for kind,value in summary['kinds'].items():
                metadata={'job_id':job['job_id'],'label':job['label'],'mode':job['mode'],
                          'comparison_role':comparison_role(job['mode']),'arm':arm,
                          'retrieval_budget':summary['retrieval_budget'],'evaluation_population':population,
                          'scope':scope,'label_kind':kind}
                def component_row(row,dimension,key):
                    fields={k:v for k,v in row.items() if not k.startswith('by_') and k!='query_metrics'}
                    return {**metadata,'stratum_dimension':dimension,'stratum':key,**fields}
                component_rows.append(component_row(value,'overall','all'))
                query_rows.append({**metadata,'stratum_dimension':'overall','stratum':'all',**value['query_metrics']})
                for field,dimension in [('by_component_count','component_count'),('by_fraction','fraction'),
                                        ('by_position','position'),('by_target','target')]:
                    for key,row in value[field].items():
                        component_rows.append(component_row(row,dimension,key))
                        if field=='by_component_count':
                            query_rows.append({**metadata,'stratum_dimension':dimension,'stratum':key,**row['query_metrics']})
                if scope=='combined' and kind=='exact' and arm==primary_arm(job['mode']):
                    headlines.append({**metadata,**component_row(value,'overall','all'),**value['query_metrics']})
    return component_rows,query_rows,headlines


def _observation(job,query,component,ranking):
    row=component['exact'];accession=component['accession']
    rows=ranking[accession] if isinstance(ranking,dict) else ranking
    return {'job_id':job['job_id'],'label':job['label'],'mode':job['mode'],
            'comparison_role':comparison_role(job['mode']),'query_id':query['query_id'],
            'group_id':query['group_id'],'suite':query['suite'],'component_count':len(query['components']),
            'accession':accession,'component_index':component['component_index'],'position':component['position'],
            'fraction':component['fraction'],'fraction_bin':component['fraction_bin'],
            'start0':component['start0'],'end0':component['end0'],'length_aa':query.get('length_aa'),
            'components_in_order':list(query['components']),'rank':row['rank'],
            'top10_accessions':[r['accession'] for r in rows[:10]],
            'near_tie_count':row.get('near_tie_count'),'exact_tie_width':row.get('tie_worst_rank',0)-row.get('tie_best_rank',0)+1}


def order_movements(observations):
    """Within one unordered group, compare all observed orders descriptively."""
    groups=defaultdict(list)
    for row in observations:
        if row['component_count']>1:groups[(row['job_id'],row['arm'],row['group_id'],row['accession'])].append(row)
    output=[]
    for (job,arm,group,accession),rows in sorted(groups.items()):
        if len(rows)<2:continue
        best=min(rows,key=lambda r:(r['rank'],r['query_id']))
        worst=min(rows,key=lambda r:(-r['rank'],r['query_id']))
        positions={p:[r['rank'] for r in rows if r['position']==p] for p in ('N','internal','C')}
        means={p:math.fsum(values)/len(values) if values else None for p,values in positions.items()}
        output.append({'job_id':job,'label':rows[0]['label'],'arm':arm,
                       'retrieval_budget':rows[0]['retrieval_budget'],'group_id':group,'accession':accession,
                       'component_count':rows[0]['component_count'],'fraction':rows[0]['fraction'],
                       'observed_order_count':len(rows),'best_rank':best['rank'],'worst_rank':worst['rank'],
                       'rank_range':worst['rank']-best['rank'],'mean_rank':math.fsum(r['rank'] for r in rows)/len(rows),
                       'best_query_id':best['query_id'],'worst_query_id':worst['query_id'],
                       'N_mean_rank':means['N'],'internal_mean_rank':means['internal'],'C_mean_rank':means['C'],
                       'C_minus_N_mean_rank':means['C']-means['N'] if means['C'] is not None and means['N'] is not None else None,
                       'N_order_count':len(positions['N']),'internal_order_count':len(positions['internal']),
                       'C_order_count':len(positions['C'])})
    return output


def _distinct_group_cases(rows,limit):
    chosen=[];seen=set()
    for row in rows:
        if row['group_id'] in seen:continue
        chosen.append(row);seen.add(row['group_id'])
        if len(chosen)>=limit:break
    return chosen


def deterministic_cases(primary_observations,limit=5):
    """Outcome-described examples under fixed rules, without independence claims."""
    by_job=defaultdict(list)
    for row in primary_observations:
        if row['component_count']>1:by_job[row['job_id']].append(row)
    cases=[]
    for job,rows in sorted(by_job.items()):
        worst=sorted(rows,key=lambda r:(-r['rank'],r['query_id'],r['accession']))
        best=sorted(rows,key=lambda r:(r['rank'],r['fraction'],r['query_id'],r['accession']))
        for rule,ordered in [('worst_component_rank',worst),('best_rank_then_smallest_component_fraction',best)]:
            cases.extend({'selection_rule':rule,'case_type':'within_method','observation':row}
                         for row in _distinct_group_cases(ordered,limit))
    for left,right in itertools.combinations(sorted(by_job),2):
        a={(r['query_id'],r['accession']):r for r in by_job[left]}
        b={(r['query_id'],r['accession']):r for r in by_job[right]}
        require(set(a)==set(b),'Method cases must use identical component observations')
        differences=[{'group_id':a[key]['group_id'],'left':a[key],'right':b[key],
                      'right_minus_left_rank':b[key]['rank']-a[key]['rank']} for key in a]
        differences=[row for row in differences if row['right_minus_left_rank']!=0]
        differences.sort(key=lambda r:(-abs(r['right_minus_left_rank']),r['left']['query_id'],r['left']['accession']))
        cases.extend({'selection_rule':'largest_absolute_rank_disagreement','case_type':'between_methods',**row}
                     for row in _distinct_group_cases(differences,limit))
    return cases


def paired_method_comparisons(primary_observations):
    """Rank/recall differences only; no cross-representation cosine comparison."""
    by_job=defaultdict(dict)
    for row in primary_observations:
        if row['component_count']>1:by_job[row['job_id']][(row['query_id'],row['accession'])]=row
    results=[]
    for left,right in itertools.combinations(sorted(by_job),2):
        require(set(by_job[left])==set(by_job[right]),'Paired method query/component set differs')
        groupings=defaultdict(list)
        for key,a in by_job[left].items():
            b=by_job[right][key]
            groupings['combined'].append((a,b));groupings['count'+str(a['component_count'])].append((a,b))
        for scope,pairs in sorted(groupings.items()):
            differences=[b['rank']-a['rank'] for a,b in pairs]
            n=len(pairs)
            results.append({'left_job':left,'right_job':right,'left_label':pairs[0][0]['label'],
                            'right_label':pairs[0][1]['label'],'scope':scope,'component_observations':n,
                            'queries':len({a['query_id'] for a,b in pairs}),
                            'distinct_groups':len({a['group_id'] for a,b in pairs}),
                            'right_rank_better':sum(x<0 for x in differences),'same_rank':sum(x==0 for x in differences),
                            'right_rank_worse':sum(x>0 for x in differences),
                            'mean_right_minus_left_rank':math.fsum(differences)/n,
                            'mean_absolute_rank_difference':math.fsum(abs(x) for x in differences)/n,
                            **{f'right_minus_left_recall_at_{k}':sum(int(b['rank']<=k)-int(a['rank']<=k) for a,b in pairs)/n for k in (1,5,10)},
                            'contains_exploratory_method':any(row['comparison_role']!='primary' for pair in pairs[:1] for row in pair),
                            'inference':'descriptive paired differences; shared proteins/orders are correlated, no p-values'})
    return results


def mechanism_rows(job,query,evaluations):
    if job['mode']!='model' or len(query['components'])==1:return []
    contrasts=[('whole','isolated_mixture','joined_context_under_whole_pooling'),
               ('contextual_component','isolated_component','joined_context_with_known_component_pooling'),
               ('isolated_mixture','isolated_component','mixing_isolated_components'),
               ('whole','contextual_component','mixing_contextual_components')]
    result=[]
    for affected,reference,name in contrasts:
        left={c['accession']:c['exact'] for c in evaluations[affected]['components']}
        right={c['accession']:c['exact'] for c in evaluations[reference]['components']}
        for accession in query['components']:
            a,b=left[accession]['rank'],right[accession]['rank']
            result.append({'job_id':job['job_id'],'label':job['label'],'query_id':query['query_id'],
                           'group_id':query['group_id'],'component_count':len(query['components']),
                           'accession':accession,'contrast':name,'affected_arm':affected,'reference_arm':reference,
                           'affected_rank':a,'reference_rank':b,'rank_change':a-b,
                           'affected_budget':evaluations[affected]['retrieval_budget'],
                           'reference_budget':evaluations[reference]['retrieval_budget']})
    return result


def summarize_mechanism(rows):
    groups=defaultdict(list)
    for row in rows:
        groups[(row['job_id'],row['contrast'],'combined')].append(row)
        groups[(row['job_id'],row['contrast'],'count'+str(row['component_count']))].append(row)
    result=[]
    for (job,contrast,scope),values in sorted(groups.items()):
        n=len(values)
        result.append({'job_id':job,'label':values[0]['label'],'contrast':contrast,'scope':scope,
                       'component_observations':n,'distinct_groups':len({v['group_id'] for v in values}),
                       'queries':len({v['query_id'] for v in values}),
                       'ranks_improved':sum(v['rank_change']<0 for v in values),
                       'ranks_unchanged':sum(v['rank_change']==0 for v in values),
                       'ranks_worsened':sum(v['rank_change']>0 for v in values),
                       'mean_rank_change':math.fsum(v['rank_change'] for v in values)/n,
                       **{f'change_in_recall_at_{k}':sum(int(v['affected_rank']<=k)-int(v['reference_rank']<=k) for v in values)/n for k in (1,5,10)},
                       'affected_arm':values[0]['affected_arm'],'reference_arm':values[0]['reference_arm'],
                       'affected_budget':values[0]['affected_budget'],'reference_budget':values[0]['reference_budget'],
                       'interpretation':'Within-model representation contrast; ranks are nonlinear and contrasts are not an additive causal decomposition.'})
    return result


def analyse_jobs(jobs,equivalence,full_query_ids,query_metadata):
    """Consume each job's records iterator once, retaining only compact evidence."""
    require(jobs and len({j['job_id'] for j in jobs})==len(jobs),'Need distinct source jobs')
    common=set.intersection(*(set(j['query_ids']) for j in jobs))
    require(common,'Source jobs have no shared queries')
    require(common<=set(full_query_ids),'Common query set lies outside frozen suite')
    tables={'coverage':[],'component_metrics':[],'query_metrics':[],'headlines':[]}
    observations=[];primary=[];mechanism=[];control_checks=[]
    for job in jobs:
        own=set(job['query_ids']);require(len(own)==len(job['query_ids']),'Duplicate source query ID')
        require(own<=set(full_query_ids),'Source query outside frozen suite')
        tables['coverage'].append({'job_id':job['job_id'],'label':job['label'],'mode':job['mode'],
                                  'comparison_role':comparison_role(job['mode']),'full_suite_queries':len(full_query_ids),
                                  'job_completed_queries':len(own),'common_queries':len(common),
                                  'job_only_queries':len(own-common),'missing_from_full_suite':len(set(full_query_ids)-own),
                                  'completed_by_suite':dict(Counter(query_metadata[q]['suite'] for q in own)),
                                  'common_by_suite':dict(Counter(query_metadata[q]['suite'] for q in common)),
                                  'config':job.get('config',{})})
        full,shared={},{};seen=set()
        for record in job['records']:
            query=record['query'];qid=query['query_id']
            require(qid in own and qid not in seen,'Missing/unexpected/duplicated source record')
            seen.add(qid)
            arms=source_arms({'mode':job['mode']},record);evaluations={}
            for arm,(ranking,pointer) in arms.items():
                evaluator=evaluate_oracle_query if isinstance(ranking,dict) else evaluate_query
                evaluation=evaluator(ranking,query,equivalence)
                expected_key='accession105' if job['mode']=='alignment_normalization' else arm
                if 'evaluations' in record:
                    require(record['evaluations'][expected_key]==evaluation,'Saved/native reevaluation differs; investigate before comparison')
                evaluations[arm]=evaluation
                add_bundle(full.setdefault(arm,new_bundle()),evaluation)
                if qid in common:
                    add_bundle(shared.setdefault(arm,new_bundle()),evaluation)
                    for component in evaluation['components']:
                        require(component['exact']['scorable'],'Core comparison unexpectedly lacks target')
                        row={**_observation(job,query,component,ranking),'arm':arm,
                             'retrieval_budget':evaluation['retrieval_budget']}
                        observations.append(row)
                        if arm==primary_arm(job['mode']):primary.append(row)
            if qid in common:
                mechanism.extend(mechanism_rows(job,query,evaluations))
                if job['mode']=='model':
                    controls=record.get('controls',{})
                    for metric in ('whole_weighted_context_max_abs','mixture_reverse_max_abs','independent_clean_reference_max_abs'):
                        if metric in controls:
                            require(math.isfinite(controls[metric]) and 0<=controls[metric]<=1e-6,
                                    'Saved model control identity exceeds its frozen tolerance')
                            control_checks.append({'job_id':job['job_id'],'query_id':qid,'metric':metric,
                                                   'value':controls[metric],'tolerance':1e-6})
        require(seen==own,'Source records silently omitted a declared query')
        for population,bundles in [('common_queries',shared),('full_job_coverage',full)]:
            comp,queries,head=flatten_summaries(job,bundles,population)
            tables['component_metrics'].extend(comp);tables['query_metrics'].extend(queries);tables['headlines'].extend(head)
    tables['order_rank_movements']=order_movements(observations)
    tables['paired_method_comparisons']=paired_method_comparisons(primary)
    tables['mechanism_contrasts']=summarize_mechanism(mechanism)
    tables['mechanism_component_observations']=mechanism
    tables['control_identity_checks']=control_checks
    return {'common_query_ids':sorted(common),'tables':tables,'cases':deterministic_cases(primary),
            'case_selection_rules':{'limit_per_category_method_or_method_pair':5,'group_repetition':'at most one case per unordered group within each category',
                                    'worst':'descending exact component rank; then query ID, accession',
                                    'best':'ascending rank, fraction, query ID, accession',
                                    'discordant':'descending absolute rank difference, then query ID, accession; zero differences excluded'},
            'interpretation':['Clean sanity queries excluded from combined headline summaries.',
                              'Only common query IDs and the identical full105 catalogue support direct comparisons; per-job coverage summaries are separate.',
                              'Oracle arms use k per known component and remain labelled with their larger review budget.',
                              'Ranks/recall compared across models; cosine magnitudes are not compared across representations.',
                              'Repeated proteins and query permutations are correlated descriptive observations; no p-values or independence claim.',
                              'Model checkpoint differences include architecture/training and do not isolate a causal model-size effect.',
                              'Reference-self-score normalized alignment is exploratory and separate from primary raw alignment.']}


def csv_bytes(rows):
    if not rows:return b''
    fields=list(dict.fromkeys(key for row in rows for key in row))
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
    for row in rows:
        writer.writerow({key:json.dumps(value,sort_keys=True,separators=(',',':')) if isinstance(value,(list,dict)) else value
                         for key,value in row.items()})
    return stream.getvalue().encode()


def plot_common_recall(report,path):
    """Standard static scientific plot; called only with explicit --plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=[r for r in report['tables']['component_metrics'] if r['evaluation_population']=='common_queries'
          and r['scope']=='combined' and r['label_kind']=='exact' and r['stratum_dimension']=='component_count'
          and r['arm']==primary_arm(r['mode'])]
    grouped=defaultdict(list)
    for row in rows:grouped[row['job_id']].append(row)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    for _,values in sorted(grouped.items()):
        values.sort(key=lambda r:int(r['stratum']))
        style='--' if values[0]['comparison_role']!='primary' else '-'
        ax.plot([int(v['stratum']) for v in values],[v['recall_at_5'] for v in values],style,marker='o',label=values[0]['label'])
    ax.set(xlabel='Number of complete components in query',ylabel='Component recall@5',ylim=(0,1.02),xticks=[2,3,4],
           title='Shared-query retrieval against the same 105-accession catalogue')
    ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
    stream=io.BytesIO();fig.savefig(stream,format='png',dpi=180);plt.close(fig)
    atomic_bytes(path,stream.getvalue())


def run(args):
    started=time.monotonic();session=args.session.resolve()
    deadline=dt.datetime.fromisoformat(read(session/'launch.json')['hard_stop_utc'])
    def check_deadline():require(dt.datetime.now(dt.timezone.utc)<deadline,'Session hard stop reached')
    check_deadline();suite=load_suite(ROOT,session/'inputs');accessions=sorted(suite['sequences'])
    require(len(accessions)==105,'Core comparison requires the fixed105 catalogue')
    sequences={a:digest(suite['sequences'][a].encode()) for a in accessions}
    queries={q['query_id']:q for q in suite['queries']};jobs=[];sources=[]
    if args.label:require(len(args.label)==len(args.source_job),'Supply exactly one label per source job')
    for index,path in enumerate(args.source_job):
        source=path.resolve();identity=read(source/'identity.json');source_hash=digest(canonical(identity))
        require(identity['mode'] in ('model','alignment','alignment_normalization'),'Unsupported source mode')
        require(identity['catalogue_accessions']==accessions and identity['catalogue_sequence_sha256']==sequences,
                'Fair comparison requires identical catalogue accessions and source sequences')
        require(identity['challenge_sha256']==suite['verification']['challenge_manifest_sha256'],'Source challenge differs')
        qids=identity['query_ids'];require(qids and len(qids)==len(set(qids)) and set(qids)<=set(queries),'Invalid source queries')
        completion_path=latest_completion(source);completion=read(completion_path)
        expected={'query-'+qid for qid in qids}|{'summary'}
        if identity['mode']=='model':expected|={'catalogue-'+a for a in accessions}
        if identity['mode']=='alignment_normalization':expected.add('reference-selfscores')
        require(completion['identity_hash']==source_hash and completion['expected_count']==len(expected)
                and set(completion['units'])==expected,'Source completion missing/inconsistent units')
        paths={key:verify_unit(source,key,entry,source_hash) for key,entry in completion['units'].items()}
        def records(paths=paths,qids=qids,identity=identity):
            for qid in qids:
                check_deadline();record=read(paths['query-'+qid]/'record.json')
                require(record['query']=={k:v for k,v in queries[qid].items() if k!='sequence'},'Source query metadata changed')
                for _,(ranking,_) in source_arms(identity,record).items():
                    for rows in ranking.values() if isinstance(ranking,dict) else [ranking]:validate_ranking(rows,accessions)
                yield record
        label=args.label[index] if args.label else display_name(identity)
        jobs.append({'job_id':source.name,'label':label,'mode':identity['mode'],'config':identity.get('config',{}),
                     'query_ids':qids,'records':records()})
        sources.append({'job_id':source.name,'label':label,'path':str(source),'identity_hash':source_hash,
                        'completion_path':str(completion_path),'completion_sha256':file_hash(completion_path),
                        'mode':identity['mode'],'config':identity.get('config',{}),'query_ids':qids})
    analysis_identity={'version':1,'mode':'common_model_comparison','sources':sources,
                       'challenge_sha256':suite['verification']['challenge_manifest_sha256'],
                       'catalogue_accessions':accessions,'catalogue_sequence_sha256':sequences,
                       'plots_requested':args.plots,'code_hashes':{str(p.relative_to(ROOT)):file_hash(p) for p in
                            [Path(__file__),ROOT/'scripts/session_catalogue_analysis.py',ROOT/'scripts/session_independent_audit.py',
                             ROOT/'src/iri_constructs/session_benchmark.py',ROOT/'src/iri_constructs/session_store.py']}}
    jobid=args.name+'-'+digest(canonical(analysis_identity))[:12]
    store=UnitStore(session/'jobs'/jobid,analysis_identity,session/'job_ledger.jsonl')
    store.event(event='job_started',state='running',job_id=jobid,pid=os.getpid())
    try:
        report=analyse_jobs(jobs,suite['equivalence'],set(queries),queries)
        report['source_provenance']=sources
        old=store.read('comparison-data')
        if old is None:store.publish('comparison-data',report)
        else:require(old==report,'Recomputed comparison differs from preserved results')
        # Export regenerable views of checksum-committed data; never overwrite
        # a differing existing export silently.
        exports=store.root/'exports';exports.mkdir(exist_ok=True)
        artifacts={}
        for name,rows in report['tables'].items():
            target=exports/(name+'.csv');data=csv_bytes(rows)
            if target.exists():require(target.read_bytes()==data,'Existing CSV export differs')
            else:atomic_bytes(target,data)
            artifacts[target.name]={'sha256':file_hash(target),'bytes':target.stat().st_size}
        target=exports/'cases.json';data=canonical({'rules':report['case_selection_rules'],'cases':report['cases']})+b'\n'
        if target.exists():require(target.read_bytes()==data,'Existing case export differs')
        else:atomic_bytes(target,data)
        artifacts[target.name]={'sha256':file_hash(target),'bytes':target.stat().st_size}
        if args.plots:
            target=exports/'component_recall5_by_count.png'
            if not target.exists():plot_common_recall(report,target)
            from PIL import Image
            with Image.open(target) as picture:picture.verify()
            artifacts[target.name]={'sha256':file_hash(target),'bytes':target.stat().st_size}
        index={'path':str(exports),'artifacts':artifacts,'regeneration':'CSV/JSON views derive from immutable comparison-data unit; optional PNG uses plot_common_recall.'}
        old=store.read('exports')
        if old is None:store.publish('exports',index)
        else:require(old==index,'Saved export hashes differ')
        for name,evidence in artifacts.items():require(file_hash(exports/name)==evidence['sha256'],'Export changed before completion')
        complete=store.complete(['comparison-data','exports'],{'common_queries':len(report['common_query_ids']),
                                                               'source_jobs':len(jobs),'elapsed_seconds':time.monotonic()-started})
        print(json.dumps({'status':'complete','job':str(store.root),'common_queries':len(report['common_query_ids']),
                          'units':complete['expected_count'],'elapsed_seconds':time.monotonic()-started}),flush=True)
    except BaseException as exc:
        store.event(event='job_failed',state='failed',job_id=jobid,exception=f'{type(exc).__name__}: {exc}')
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--session',type=Path,required=True);p.add_argument('--source-job',type=Path,action='append',required=True)
    p.add_argument('--label',action='append');p.add_argument('--name',default='fair-common-model-comparison')
    p.add_argument('--plots',action='store_true');run(p.parse_args())
