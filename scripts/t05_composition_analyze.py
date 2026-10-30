#!/usr/bin/env python3
"""Descriptive summaries and publication-exportable figures, without new selection."""
import argparse
from collections import defaultdict, Counter
import csv
import json
from pathlib import Path
import statistics
import sys
from time import perf_counter
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.t05_composition import load_reviewed, verify_freeze
from iri_constructs.t05_data import sha256
from t05_composition_experiment import read,write,csv_write,BASE,PLAN,MANIFEST,FREEZE,safe_id


def mean(x):return float(np.mean(x))


def averaged_ranks(values):
    return np.array([sum(t<v for t in values)+(sum(t==v for t in values)+1)/2 for v in values])


def spearman(x,y):
    a,b=averaged_ranks(x),averaged_ranks(y)
    if np.std(a)==0 or np.std(b)==0:return None
    return float(np.corrcoef(a,b)[0,1])


def summarize(rows):
    result={'component_instances':len(rows),'distinct_proteins':len({r['accession'] for r in rows}),
            'distinct_pairs':len({r['pair_id'] for r in rows if r['pair_id']}),
            'ordered_queries':len({r['query_id'] for r in rows}),
            'mean_fraction':mean([r['fraction_in_parent'] for r in rows]),
            'mean_rank':mean([r['rank'] for r in rows]),'median_rank':float(np.median([r['rank'] for r in rows])),
            'mean_rank_change':mean([r['rank_change'] for r in rows]),
            'mrr':mean([r['reciprocal_rank'] for r in rows]),
            'mean_cosine':mean([r['cosine_similarity'] for r in rows]),
            'mean_cosine_change':mean([r['cosine_change'] for r in rows]),
            'improved_rank':sum(r['rank_change']<0 for r in rows),
            'unchanged_rank':sum(r['rank_change']==0 for r in rows),
            'worsened_rank':sum(r['rank_change']>0 for r in rows),
            'target_exact_tie_instances':sum(r['exact_score_tie_count']>1 for r in rows),
            'target_near_tie_instances':sum(r['near_tie_count_1e6']>1 for r in rows)}
    for k in (1,5,10):
        result[f'recall_at_{k}']=mean([r[f'hit_at_{k}'] for r in rows])
        result[f'sequence_equivalent_recall_at_{k}']=mean([r[f'sequence_equivalent_hit_at_{k}'] for r in rows])
    result['fraction_vs_rank_spearman']=spearman([r['fraction_in_parent'] for r in rows],[r['rank'] for r in rows])
    return result


def group_summary(rows,keys):
    groups=defaultdict(list)
    for r in rows:groups[tuple(r[k] for k in keys)].append(r)
    return [{**dict(zip(keys,key)),**summarize(values)} for key,values in sorted(groups.items())]


def relatedness(sequences,plan,out):
    from Bio.Align import PairwiseAligner,substitution_matrices
    aligner=PairwiseAligner(mode='global'); aligner.substitution_matrix=substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
    ids=sorted(sequences);rows=[];started=perf_counter()
    for i,a in enumerate(ids):
        for b in ids[i+1:]:
            alignment=aligner.align(sequences[a],sequences[b])[0]
            matches=sum(x==y and x!='-' for x,y in zip(alignment[0],alignment[1],strict=True))
            rows.append({'a':a,'b':b,'identical_residues':matches,'alignment_columns':alignment.shape[1],
                         'global_identity':matches/alignment.shape[1],'exact_duplicate':sequences[a]==sequences[b]})
    summary={'method':'Global alignment: BLOSUM62; gap open -10, extension -0.5, including termini; first optimum; identity denominator includes gap columns',
             'manifest_sha256':sha256(MANIFEST.read_bytes()),'alignment_seconds':perf_counter()-started,'pairs':len(rows),
             'catalogue_neighbour_counts':{str(c):sum(any(a in (r['a'],r['b']) and r['global_identity']>=c for r in rows) for a in ids) for c in [.3,.5,.9,1.]},'targets':[]}
    for a in plan['targets']:
        others=[r for r in rows if a in (r['a'],r['b'])];best=max(others,key=lambda r:r['global_identity'])
        summary['targets'].append({'accession':a,'gene':plan['target_genes'][a],
            'nearest_other':best['b'] if best['a']==a else best['a'],'identity':best['global_identity'],
            'neighbours_at_least_50pct':sum(r['global_identity']>=.5 for r in others)})
    csv_write(out/'catalogue_pairwise_identity.csv',rows);write(out/'relatedness.json',summary)
    return summary,rows


def analyze(identifier):
    run=BASE/safe_id(identifier);data=read(run/'results.json')
    if data['status']!='passed':raise ValueError('Analyze only a completed run; inspect failure record separately')
    verify_freeze(ROOT,read(FREEZE))
    if data['freeze_sha256']!=sha256(FREEZE.read_bytes()):raise ValueError('Run does not match frozen sources')
    plan,manifest=read(PLAN),read(MANIFEST);seq=load_reviewed(ROOT,manifest)
    out=run/'analysis';out.mkdir(exist_ok=True)
    rel,identity_rows=relatedness(seq,plan,out)
    whole=[r for r in data['component_rows'] if r['arm']=='whole']
    clean=[r for r in data['component_rows'] if r['arm']=='clean']
    oracle=[r for r in data['component_rows'] if r['arm'].startswith('oracle')]
    summary={'catalogue_size':data['catalogue_size'],'counts':data['counts'],
             'whole':summarize(whole),'clean':summarize(clean),'oracle':summarize(oracle),
             'by_fraction':group_summary(whole,['fraction_bin']),
             'by_position':group_summary(whole,['position_in_parent']),
             'by_fraction_and_position':group_summary(whole,['fraction_bin','position_in_parent']),
             'by_target':group_summary(whole,['gene','accession']),
             'clean_failures':data['clean_query_failures'],'clean_rows':clean,'relatedness':rel,
             'analysis_code_sha256':sha256(Path(__file__).read_bytes())}
    for kind in ['by_fraction','by_position','by_fraction_and_position','by_target']:
        csv_write(out/(kind+'.csv'),summary[kind])
    summary['equal_target_weighted']={k:mean([t[k] for t in summary['by_target']]) for k in ['mrr','recall_at_1','recall_at_5','recall_at_10','mean_rank']}
    qwhole=[e for e in data['evaluations'] if e['arm']=='whole']
    query_rows=[{'query_id':e['query_id'],'pair_id':e['pair_id'],**e['retrieval']} for e in qwhole]
    csv_write(out/'query_retrieval.csv',query_rows)
    summary['query_metrics']={k:mean([q[k] for q in query_rows]) for k in query_rows[0] if k not in ['query_id','pair_id']}
    bypair=defaultdict(dict)
    for r in whole:bypair[(r['pair_id'],r['accession'])][r['position_in_parent']]=r
    reversals=[]
    for (pair,a),arms in sorted(bypair.items()):
        n,c=arms['N'],arms['C']
        reversals.append({'pair_id':pair,'accession':a,'gene':n['gene'],'partner':n['partner'],
            'fraction':n['fraction_in_parent'],'rank_N':n['rank'],'rank_C':c['rank'],
            'C_minus_N_rank':c['rank']-n['rank'],'C_minus_N_cosine':c['cosine_similarity']-n['cosine_similarity']})
    csv_write(out/'position_reversals.csv',reversals)
    deltas=[r['C_minus_N_rank'] for r in reversals]
    summary['position_reversal']={'protein_pair_comparisons':len(deltas),'distinct_pairs':len({r['pair_id'] for r in reversals}),
        'C_better':sum(d<0 for d in deltas),'N_better':sum(d>0 for d in deltas),'same_rank':sum(d==0 for d in deltas),
        'mean_C_minus_N_rank':mean(deltas),'median_C_minus_N_rank':float(np.median(deltas)),
        'mean_abs_rank_change':mean([abs(d) for d in deltas]),'max_abs_rank_change':max(map(abs,deltas)),
        'largest_changes':sorted(reversals,key=lambda r:-abs(r['C_minus_N_rank']))[:6]}
    summary['unexpected_examples']={
        'small_fraction_rank1':sorted([r for r in whole if r['rank']==1],key=lambda r:r['fraction_in_parent'])[:5],
        'large_fraction_low_rank':sorted([r for r in whole if r['fraction_in_parent']>=.6],key=lambda r:-r['rank'])[:5],
        'worst_ranks':sorted(whole,key=lambda r:-r['rank'])[:5]}
    # Describe incorrect top hits relative to the TWO valid components, not just one.
    relations={frozenset((r['a'],r['b'])):r['global_identity'] for r in identity_rows}
    miss_rows=[]
    for e in qwhole:
        top=e['ranking'][0]['accession']
        if top not in e['expected']:
            miss_rows.append({'query_id':e['query_id'],'top_accession':top,
                'max_global_identity_to_present_component':max(relations.get(frozenset((top,a)),1.) for a in e['expected']),
                'best_present_rank':min(r['rank'] for r in e['ranking'] if r['accession'] in e['expected'])})
    csv_write(out/'noncomponent_top_hits.csv',miss_rows)
    summary['noncomponent_top_hits']={'count':len(miss_rows),'at_least_50pct_identity_to_a_component':sum(r['max_global_identity_to_present_component']>=.5 for r in miss_rows),
                                     'most_frequent':Counter(r['top_accession'] for r in miss_rows).most_common(8)}
    summary['ties']={'evaluations_with_any_exact_score_tie':sum(any(r['exact_score_tie_count']>1 for r in e['ranking']) for e in data['evaluations']),
        'duplicate_sequence_groups':manifest['duplicate_sequence_groups'],
        'target_exact_ties':sum(r['exact_score_tie_count']>1 for r in data['component_rows']),
        'target_near_ties':sum(r['near_tie_count_1e6']>1 for r in data['component_rows'])}
    # Explicitly post hoc sensitivity: primary frozen accession ranking stays intact.
    sequence_hash={a:sha256(s.encode()) for a,s in seq.items()}
    collapsed=[]
    for e in qwhole:
        seen={}
        for r in e['ranking']:
            h=sequence_hash[r['accession']]
            if h not in seen:seen[h]=len(seen)+1
            if r['accession'] in e['expected']:
                collapsed.append({'query_id':e['query_id'],'accession':r['accession'],
                    'primary_rank':r['rank'],'unique_sequence_rank':seen[h]})
    sensitivity={'label':'Post hoc duplicate-collapse sensitivity; no re-encoding or primary catalogue change',
        'unique_sequence_catalogue_size':len(set(seq.values())),
        'instances_with_lower_numerical_rank':sum(r['unique_sequence_rank']<r['primary_rank'] for r in collapsed),
        'max_rank_shift':max(r['primary_rank']-r['unique_sequence_rank'] for r in collapsed),
        'hit_counts':{str(k):sum(r['unique_sequence_rank']<=k for r in collapsed) for k in (1,5,10)}}
    write(out/'duplicate-collapse-sensitivity.json',{**sensitivity,'rows':collapsed})
    summary['duplicate_collapse_sensitivity']=sensitivity
    write(out/'summary.json',summary)
    plots(out,whole,summary,reversals,plan)
    print(json.dumps({k:summary[k] for k in ['whole','query_metrics','position_reversal','noncomponent_top_hits']},indent=2))


def plots(out,whole,summary,reversals,plan):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    def save(fig,name):
        for suffix in ['png','pdf','svg']:fig.savefig(out/(name+'.'+suffix),bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for position,color,marker in [('N','#2479ad','o'),('C','#d67b20','^')]:
        rows=[r for r in whole if r['position_in_parent']==position]
        axes[0].scatter([r['fraction_in_parent'] for r in rows],[r['rank'] for r in rows],s=19,alpha=.55,c=color,marker=marker,label=position+'-terminal component')
    axes[0].set(xlabel='Component fraction of total residues',ylabel='Component rank among 105 references (log scale)',yscale='log',ylim=(.8,130))
    axes[0].legend(fontsize=8);axes[0].set_title('Length fraction is associated with retrieval')
    bins=summary['by_fraction'];x=np.arange(len(bins))
    for k,color in [(1,'#283c63'),(5,'#24836f'),(10,'#cc7b23')]:
        axes[1].plot(x,[r[f'recall_at_{k}'] for r in bins],'o-',label=f'Recall@{k}',color=color)
    axes[1].set_xticks(x,[r['fraction_bin']+'\nn='+str(r['component_instances']) for r in bins],fontsize=8)
    axes[1].set(xlabel='Prespecified fraction bins; n = component instances',ylabel='Fraction of components retrieved',ylim=(-.02,1.04))
    axes[1].legend(fontsize=8);axes[1].set_title('460 correlated instances from 115 pairs')
    fig.suptitle('ESM-2 8M • unchanged proteins • direct artificial concatenation',fontsize=13)
    save(fig,'fraction_retrieval')
    targets=plan['targets'];a=np.full((len(targets),len(targets)),np.nan)
    for i,t in enumerate(targets):
        a[i,i]=next(r['rank'] for r in summary['clean_rows'] if r['accession']==t)
        for j,p in enumerate(targets):
            rows=[r for r in whole if r['accession']==t and r['partner']==p]
            if rows:a[i,j]=mean([r['rank'] for r in rows])
    fig,ax=plt.subplots(figsize=(10,8),layout='constrained');cmap=plt.get_cmap('viridis_r').copy();cmap.set_bad('#dadada')
    im=ax.imshow(a,norm=LogNorm(vmin=1,vmax=105),cmap=cmap)
    labels=[plan['target_genes'][t] for t in targets]
    ax.set_xticks(range(len(targets)),labels,rotation=60,ha='right');ax.set_yticks(range(len(targets)),labels)
    ax.set(xlabel='Partner protein (targets ordered by length)',ylabel='Component being retrieved',title='Component identity matters: mean rank across both orders\nDiagonal = clean control; grey = excluded length combination')
    for i in range(len(targets)):
        for j in range(len(targets)):
            if np.isfinite(a[i,j]):
                rgb=cmap(im.norm(a[i,j]))[:3]
                luminance=sum(w*c for w,c in zip((.299,.587,.114),rgb))
                ax.text(j,i,f'{a[i,j]:g}',ha='center',va='center',fontsize=7,color='black' if luminance>.52 else 'white')
    fig.colorbar(im,ax=ax,label='Mean accession rank / 105 (log colour scale)',shrink=.75)
    save(fig,'target_partner_heatmap')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    axes[0].scatter([r['rank_N'] for r in reversals],[r['rank_C'] for r in reversals],alpha=.5,s=23,color='#4664a5')
    axes[0].plot([1,105],[1,105],'k--',lw=1)
    axes[0].set(xscale='log',yscale='log',xlabel='Rank with component first (N)',ylabel='Rank with component second (C)',title='Same proteins and length fraction; reversed order')
    targetrows={r['accession']:r for r in summary['by_target']};x=np.arange(len(targets))
    axes[1].bar(x,[targetrows[t]['recall_at_5'] for t in targets],color='#24836f')
    axes[1].set_xticks(x,labels,rotation=65,ha='right',fontsize=8)
    axes[1].set(ylim=(0,1.05),ylabel='Component recall@5 across its partners',title='Protein-specific outcomes; unequal eligible pair counts')
    fig.suptitle('Order and identity • descriptive results, not independent biological observations',fontsize=12)
    save(fig,'position_and_identity')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run_id');analyze(p.parse_args().run_id)
