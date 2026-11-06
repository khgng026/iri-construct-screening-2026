#!/usr/bin/env python3
"""Render a frozen descriptive failure profile from existing comparison exports.

No inference, ranking, fitting, hypothesis tests, or new order/case search.
The external design hash and completed export chain must match before outcomes
are read. Run only in a lead-assigned lightweight reporting slot.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import datetime as dt
import hashlib
import io
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from session_independent_audit import verify_unit

MAX_INPUT_BYTES=32*2**20
FILTER={'mode':'model','arm':'whole','label_kind':'exact','evaluation_population':'common_queries',
        'scope':'combined','retrieval_budget':'k total per original query'}
COUNTS=('original_components','scorable_components','unscorable_components','parent_queries',
        'distinct_groups','scorable_targets','hits_at_5')
DENOMINATORS=COUNTS[:-1]
SOURCE_RULES={'limit_per_category_method_or_method_pair':5,
    'group_repetition':'at most one case per unordered group within each category',
    'worst':'descending exact component rank; then query ID, accession',
    'best':'ascending rank, fraction, query ID, accession',
    'discordant':'descending absolute rank difference, then query ID, accession; zero differences excluded'}


def require(condition,message):
    if not condition:raise ValueError(message)


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def sha(path):
    path=Path(path)
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<=MAX_INPUT_BYTES,
            'Hash input must be a bounded regular file')
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(2**20),b''):h.update(chunk)
    return h.hexdigest()


def read(path):
    path=Path(path);require(path.is_file() and not path.is_symlink(),'Expected regular source file')
    require(path.stat().st_size<=MAX_INPUT_BYTES,'Input exceeds bounded reporting size')
    return json.loads(path.read_text())


def confined(path,root):
    path=Path(path).absolute();root=Path(root).resolve()
    require(path.is_relative_to(root),'Source path escapes pinned root')
    current=root
    for part in path.relative_to(root).parts:
        require(part not in ('.','..'),'Unsafe source path')
        current/=part;require(not current.is_symlink(),'Symlink source is not eligible')
    return path


def bind_source(source,design):
    """Validate metadata and hashes; no outcome CSV/JSON parsing in this step."""
    source=source.resolve(strict=True)
    require(source.name==design['source_job'],'Source job differs from frozen reporting design')
    identity_path=confined(source/'identity.json',source)
    require(sha(identity_path)==design['identity_file_sha256'],'Source identity file changed')
    identity=read(identity_path);identity_hash=hashlib.sha256(canonical(identity)).hexdigest()
    require(identity['mode']=='common_model_comparison' and len(identity['catalogue_accessions'])==105,
            'Source is not the original fixed105 comparison')
    models={row['job_id'] for row in identity['sources'] if row['mode']=='model'}
    require(models=={row['job_id'] for row in design['models']} and len(models)==5,'Five-model roster differs')
    require(identity['code_hashes']['scripts/session_compare_models.py']==design['producer_code_sha256_at_design'],
            'Producer implementation differs from frozen design')
    completion_path=confined(source/design['completion_file'],source)
    require(sha(completion_path)==design['completion_sha256'],'Pinned completion bytes changed')
    completion=read(completion_path)
    require(completion['identity_hash']==identity_hash and completion['expected_count']==2
            and set(completion['units'])=={'comparison-data','exports'},'Comparison completion is incomplete')
    for key,reference in completion['units'].items():
        unit=confined(source/reference['path'],source)
        require(sha(unit/'complete.json')==reference['descriptor_sha256'],'Completion descriptor changed')
        descriptor=read(unit/'complete.json')
        require(descriptor['key']==key and descriptor['identity_hash']==identity_hash,'Descriptor identity changed')
    reference=completion['units']['exports']
    require(reference['descriptor_sha256']==design['exports_descriptor_sha256'],'Pinned exports descriptor changed')
    unit=verify_unit(source,'exports',reference,identity_hash)
    index=read(unit/'record.json')
    require(Path(index['path'])==source/'exports','Exports path differs from completed source')
    require(set(design['allowed_exports'])=={'component_metrics.csv','cases.json'},'Reporting export allowlist changed')
    paths={};evidence={}
    for name in design['allowed_exports']:
        path=confined(source/'exports'/name,source);expected=index['artifacts'][name]
        require(path.stat().st_size==expected['bytes']<=MAX_INPUT_BYTES and sha(path)==expected['sha256'],
                'Completed export bytes/hash changed: '+name)
        paths[name]=path;evidence[name]=dict(expected)
    return paths,{'source_job':str(source),'identity_file_sha256':sha(identity_path),
        'identity_hash':identity_hash,'completion_file':str(completion_path),'completion_sha256':sha(completion_path),
        'exports_descriptor_sha256':reference['descriptor_sha256'],'exports_record_sha256':sha(unit/'record.json'),
        'exports':evidence,'comparison_data_body_read':False,'raw_model_records_read':False}


def count(value,field):
    require(isinstance(value,(str,int)) and not isinstance(value,bool),'Invalid integer '+field)
    parsed=int(value);require(str(parsed)==str(value) and parsed>=0,'Invalid nonnegative integer '+field)
    return parsed


def parse_metric(row,index):
    result={key:row[key] for key in ('job_id','label','stratum_dimension','stratum')}
    result.update({key:count(row[key],key) for key in COUNTS})
    n=result['scorable_components'];hits=result['hits_at_5'];raw=row['recall_at_5']
    value=None if raw in ('',None) else float(raw)
    require(result['original_components']==n+result['unscorable_components'] and hits<=n,'Invalid component denominator/hits')
    require((n==0 and value is None) or (n>0 and value is not None and math.isfinite(value)
            and 0<=value<=1 and math.isclose(value,hits/n,rel_tol=0,abs_tol=1e-12)),
            'Saved recall disagrees with its original numerator/denominator')
    require(result['parent_queries']<=result['original_components'] and result['distinct_groups']<=result['parent_queries']
            and result['scorable_targets']<=n,'Invalid stratum support counts')
    return {**result,'recall_at_5':value,'source_csv_row':index,'source_row_present':True,**FILTER}


def extract_profile(rows,design):
    """Select existing summaries; validate partitions without rescoring queries."""
    models={m['job_id']:m for m in design['models']};selected={}
    categories={'overall':['all'],'fraction':design['fraction_bins'],'position':design['positions']}
    for index,row in enumerate(rows,2):
        if any(row.get(key)!=value for key,value in FILTER.items()):continue
        require(row.get('job_id') in models,'Unexpected whole-model source in common summary')
        dimension=row.get('stratum_dimension')
        if dimension not in categories:continue
        require(row['stratum'] in categories[dimension],'Unknown frozen plotting category')
        key=(row['job_id'],dimension,row['stratum'])
        require(key not in selected,'Duplicate source summary row')
        selected[key]=parse_metric(row,index)
    output=[]
    for job,model in models.items():
        overall=selected.get((job,'overall','all'));require(overall is not None,'Missing overall model summary')
        require(overall['original_components']==overall['scorable_components']==design['population']['component_observations_per_model']
                and overall['unscorable_components']==0
                and overall['parent_queries']==design['population']['joined_queries_per_model'],
                'Common joined panel or exact-reference scorable denominator changed')
        for dimension,values in categories.items():
            partition=[]
            for stratum in values:
                row=selected.get((job,dimension,stratum))
                if row is None:
                    require(dimension!='overall','Missing overall summary')
                    row={'job_id':job,'label':overall['label'],'stratum_dimension':dimension,'stratum':stratum,
                         **{key:0 for key in COUNTS},'recall_at_5':None,'source_csv_row':None,
                         'source_row_present':False,**FILTER}
                require(all(row[key]<=overall[key] for key in DENOMINATORS),'Stratum support exceeds overall support')
                partition.append(row);output.append({**row,'short_label':model['short_label']})
            if dimension!='overall':
                for field in ('original_components','scorable_components','unscorable_components','hits_at_5'):
                    require(sum(row[field] for row in partition)==overall[field],
                            dimension+' partition differs from original overall '+field)
    for dimension,values in categories.items():
        for stratum in values:
            support=[tuple(row[key] for key in (*DENOMINATORS,'source_row_present')) for row in output
                     if row['stratum_dimension']==dimension and row['stratum']==stratum]
            require(len(support)==len(models) and len(set(support))==1,'Models do not share identical stratum denominators')
    return output


def validate_observation(row,models):
    require(row['job_id'] in models and row['mode']=='model' and row['arm']=='whole'
            and row['retrieval_budget']==FILTER['retrieval_budget'],'Case is outside whole-model budget')
    require(row['component_count'] in (2,3,4) and row['query_id'] and row['group_id'] and row['accession'],
            'Case is outside joined-query population')
    require(isinstance(row['rank'],int) and not isinstance(row['rank'],bool) and 1<=row['rank']<=105,'Invalid saved case rank')


def select_cases(document,design):
    require(document['rules']==SOURCE_RULES,'Existing deterministic case rules changed')
    models={row['job_id'] for row in design['models']};counts=Counter();groups={};selected=[]
    require(len(document['cases'])<=1000,'Case input exceeds bounded source expectation')
    for index,case in enumerate(document['cases']):
        rule=case['selection_rule']
        if rule=='worst_component_rank':
            row=case['observation']
            if row['job_id'] not in models:continue
            require(case['case_type']=='within_method','Case type differs from frozen rule')
            validate_observation(row,models);key=(rule,row['job_id']);limit=2;group=row['group_id']
        elif rule=='largest_absolute_rank_disagreement':
            left,right=case['left'],case['right']
            if left['job_id'] not in models or right['job_id'] not in models:continue
            require(case['case_type']=='between_methods' and left['job_id']<right['job_id'],'Invalid saved method pair')
            validate_observation(left,models);validate_observation(right,models)
            require(all(left[k]==right[k] for k in ('query_id','accession','group_id','component_count'))
                    and case['group_id']==left['group_id'],'Discordant case has unmatched component observations')
            require(case['right_minus_left_rank']==right['rank']-left['rank']!=0,'Saved rank disagreement changed')
            key=(rule,left['job_id'],right['job_id']);limit=1;group=case['group_id']
        else:continue
        seen=groups.setdefault(key,set())
        require(group not in seen,'Source repeats an unordered group within an existing case category');seen.add(group)
        if counts[key]>=limit:continue
        selected.append({'source_case_index':index,'case':case});counts[key]+=1
    require(len(selected)<=20,'Frozen case budget exceeded')
    categories=[('worst_component_rank',job) for job in sorted(models)]
    categories += [('largest_absolute_rank_disagreement',a,b) for a,b in itertools.combinations(sorted(models),2)]
    return {'source_rules':document['rules'],'reporting_rules':design['case_rules'],'selected_cases':selected,
            'category_counts':[{'category':list(key),'selected':counts[key]} for key in categories],
            'interpretation':'Examples retain original deterministic order; worst rank need not exceed 5, and rank disagreement need not be a top5 hit/miss difference.'}


def csv_bytes(rows):
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return stream.getvalue().encode()


def plot_profile(profile,design,directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42}):
        fig,axes=plt.subplots(1,2,figsize=(10.8,5.4),sharey=True,gridspec_kw={'width_ratios':[1.5,1]})
        fig.subplots_adjust(left=.075,right=.98,top=.77,bottom=.30,wspace=.12)
        for ax,dimension,labels,title in zip(axes,['fraction','position'],
                [design['fraction_labels'],['N-terminal','Internal','C-terminal']],
                ['A  Component fraction','B  Position: unadjusted descriptives']):
            keys=design['fraction_bins'] if dimension=='fraction' else design['positions']
            n=[]
            for model in design['models']:
                rows=[next(r for r in profile if r['job_id']==model['job_id'] and r['stratum_dimension']==dimension
                           and r['stratum']==key) for key in keys]
                n=[r['scorable_components'] for r in rows]
                ax.plot(range(len(keys)),[float('nan') if r['recall_at_5'] is None else r['recall_at_5'] for r in rows],
                        color=model['color'],marker=model['marker'],markersize=5,linewidth=1.5,label=model['short_label'])
            ax.set_xticks(range(len(keys)),[f'{label}\nn={number:,}' for label,number in zip(labels,n)],fontsize=8)
            ax.set_ylim(*design['figure']['y_range']);ax.set_yticks([0,.25,.5,.75,1]);ax.set_title(title,loc='left',fontsize=10,pad=10)
            ax.grid(axis='y',alpha=.22);ax.spines[['top','right']].set_visible(False)
        axes[0].set_ylabel('Exact component recall@5')
        fig.suptitle('Whole-query retrieval against the original 105-accession catalogue',y=.98,fontsize=12)
        fig.text(.5,.91,'Same 1,130 joined queries and 3,520 component observations per model',ha='center',fontsize=9)
        handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.88),ncol=5,frameon=False,fontsize=8.5)
        fig.text(.075,.15,'Position and fraction are descriptive: target identity, component count and component fraction covary.\n'
                 'Repeated proteins and permutations are correlated; no causal position effect, fitted trend or uncertainty interval is inferred.\n'
                 'n counts component observations per category. Connecting lines only guide the eye; clean/oracle arms are excluded.',
                 fontsize=8,va='top',linespacing=1.5)
        fig.savefig(directory/'failure_profile.png',dpi=200,metadata={'Software':'session_failure_profile.py'})
        fig.savefig(directory/'failure_profile.pdf',metadata={'Creator':'session_failure_profile.py','CreationDate':None,'ModDate':None})
        plt.close(fig)
    for name in ('failure_profile.png','failure_profile.pdf'):
        with (directory/name).open('rb') as stream:os.fsync(stream.fileno())
    return matplotlib.__version__


def write_new(path,data):
    with Path(path).open('xb') as stream:stream.write(data);stream.flush();os.fsync(stream.fileno())


def run(args):
    started=time.monotonic();require(0<args.timeout<=300,'Reporting timeout must be 1–300 seconds')
    def check():require(time.monotonic()-started<args.timeout,'Bounded reporting time exceeded')
    require(sha(args.design)==args.design_sha256,'Externally reviewed design SHA256 must match')
    design=read(args.design)
    require(design['version']=='descriptive-failure-profile-v1' and design['frozen_before_outcome_row_reads'],
            'Reporting design was not frozen before outcome reads')
    source=Path(args.source_job);paths,provenance=bind_source(source,design);check()
    with paths['component_metrics.csv'].open(newline='') as stream:
        profile=extract_profile(csv.DictReader(stream),design)
    cases=select_cases(read(paths['cases.json']),design);check()
    identity={'version':1,'mode':'descriptive_failure_profile','design_sha256':args.design_sha256,
              'source':provenance,'script_sha256':sha(Path(__file__)),
              'verifier_sha256':sha(ROOT/'scripts/session_independent_audit.py')}
    token=hashlib.sha256(canonical(identity)).hexdigest()
    output=args.output_root.resolve()/('failure-profile-v1-'+token[:12])
    require(not output.exists(),'Preserve existing output; do not overwrite a completed or partial profile')
    output.mkdir(parents=True)
    write_new(output/'identity.json',canonical(identity)+b'\n')
    write_new(output/'profile_metrics.csv',csv_bytes(profile))
    write_new(output/'cases.json',canonical(cases)+b'\n')
    version=plot_profile(profile,design,output);check()
    # Recheck the exact consumed source bytes before marking this view complete.
    for name,path in paths.items():require(sha(path)==provenance['exports'][name]['sha256'],'Source export changed during reporting')
    report={'design':design,'source':provenance,'profile':profile,'selected_case_count':len(cases['selected_cases']),
            'case_category_counts':cases['category_counts'],'matplotlib_version':version,
            'inference_or_ranking_performed':False,'new_order_analysis_performed':False,
            'elapsed_seconds':time.monotonic()-started}
    write_new(output/'report.json',canonical(report)+b'\n')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(output.iterdir()) if p.is_file()}
    completion={'status':'complete','identity_hash':token,'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
                'files':files,'profile_rows':len(profile),'selected_cases':len(cases['selected_cases'])}
    write_new(output/'completion.json',canonical(completion)+b'\n')
    print(json.dumps({'status':'complete','output':str(output),'completion_sha256':sha(output/'completion.json'),
                      'rows':len(profile),'cases':len(cases['selected_cases'])}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-job',type=Path,required=True)
    parser.add_argument('--design',type=Path,required=True)
    parser.add_argument('--design-sha256',required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=180)
    run(parser.parse_args())
