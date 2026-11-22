"""Recompute paired metrics from new saved scores, with explicit review budgets."""
import argparse
import csv
from collections import defaultdict
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'sources/nishit-working/src'))
from pilot.align import best_target_rank
from budget_compare import atomic_json, digest, fixed_ranking, rank_lists, coverage
from run_comparison import configuration, seq_key, valid_array, implementation_signature


def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        if rows:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def summarize(rows,keys):
    cells=defaultdict(list)
    for row in rows:cells[tuple(row[k] for k in keys)].append(row)
    out=[]
    for cell,rs in sorted(cells.items()):
        payloads=sorted({r['payload'] for r in rs})
        macro=np.mean([np.mean([r['hit'] for r in rs if r['payload']==a]) for a in payloads])
        out.append({**dict(zip(keys,cell)),'n_items':len(rs),'n_payloads':len(payloads),
                    'n_unique_queries':len({r['sequence_sha256'] for r in rs}),
                    'n_unscorable':sum(r['n_targets']==0 for r in rs),
                    'hits':sum(r['hit'] for r in rs),'recall':sum(r['hit'] for r in rs)/len(rs),
                    'payload_macro_recall':float(macro),'mrr_within_review_budget':float(np.mean([r['rr'] for r in rs])),
                    'returned_review_slots':sum(r['review_slots'] for r in rs),
                    'returned_unique_hits':sum(r['unique_returns'] for r in rs),
                    'non_target_labelled_returns':sum(r['non_target_returns'] for r in rs),
                    'abstentions':sum(r['abstained'] for r in rs),'mean_sequence_coverage':float(np.mean([r['sequence_coverage'] for r in rs])),
                    'mean_segment_count':float(np.mean([r['n_segments'] for r in rs])),
                    'bhs_fallbacks':sum(r['fallback'] for r in rs),
                    'worst_tie_vs_actual_list_disagreements':sum(r['tie_disagreement'] for r in rs)})
    return out


def main(models,split):
    cfg=configuration();freeze=cfg['freeze_sha256'];ref_ids=[r['accession'] for r in cfg['references']]
    signature=implementation_signature()
    def load_checked(path,shape):
        value=valid_array(path,shape,cfg,signature)
        if value is None:raise ValueError('Incomplete result: '+str(path))
        return value
    embdir=ROOT/'results/embeddings'/freeze
    refmatrix=np.stack([load_checked(embdir/(r['sequence_sha256']+'.npy'),(1280,)) for r in cfg['references']]) if 'esm2_650m' in models else None
    records=[];details=[]
    for item in cfg['items']:
        if split!='all' and item['split']!=split:continue
        g=json.loads((ROOT/'results/geometry'/(digest(item['construct_id'])[:24]+'.json')).read_text())
        if (g['freeze_sha256']!=freeze or g['implementation_signature']!=signature or
            g['input_sha256']!=item['sequence_sha256']):raise ValueError('Geometry provenance mismatch')
        oracle=g['spans']['oracle'][0];os,oe=oracle['start'],oracle['end']
        for c in g['conditions']:
            allowed=[i for i in range(len(ref_ids)) if i not in set(c['mask'])]
            targets=np.array(c['targets'],dtype=int); mask=np.array(c['mask'],dtype=int)
            for arm,spans in {**g['spans'],'bhs':c['bhs']}.items():
                target_cov=max(max(0,min(s['end'],oe)-max(s['start'],os))/(oe-os) for s in spans)
                target_purity=max(max(0,min(s['end'],oe)-max(s['start'],os))/(s['end']-s['start']) for s in spans)
                cov=coverage([(s['start'],s['end']) for s in spans],len(item['sequence']))
                for model in models:
                    if model=='alignment':
                        scores=np.stack([load_checked(ROOT/'results/alignment_scores'/freeze/(s['sequence_sha256']+'.npy'),(len(ref_ids),)) for s in spans])
                    elif model=='esm2_650m':
                        vectors=np.stack([load_checked(embdir/(s['sequence_sha256']+'.npy'),(1280,)) for s in spans]); scores=vectors@refmatrix.T
                    else:raise ValueError(model)
                    if not np.isfinite(scores).all():raise ValueError('Nonfinite incomplete panel')
                    per_ranks=[best_target_rank(row,targets,mask) for row in scores]
                    per_rank=min([r for r in per_ranks if r is not None],default=None)
                    lists=rank_lists(scores,ref_ids,allowed)
                    for budget in (cfg['secondary_budget'],cfg['primary_budget']):
                        policies=['max_score','round_robin','per_segment_diagnostic'] if len(spans)>1 else ['max_score','round_robin','per_segment_diagnostic']
                        for policy in policies:
                            if policy=='per_segment_diagnostic':
                                returned=list(dict.fromkeys(i for ranking in lists for i in ranking[:budget]))
                                rank=per_rank
                                hit=int(rank is not None and rank<=budget)
                                actual_hit=int(bool(set(returned)&set(c['targets'])))
                                slots=sum(min(budget,len(lst)) for lst in lists)
                                rr=1/rank if rank and rank<=budget else 0.0
                                tie_disagreement=int(actual_hit!=hit)
                            else:
                                returned=fixed_ranking(scores,ref_ids,allowed,budget,policy)
                                rank=next((j+1 for j,i in enumerate(returned) if i in c['targets']),None)
                                hit=int(rank is not None);slots=len(returned);rr=1/rank if rank else 0.0
                                original_rank=best_target_rank(np.max(scores,axis=0),targets,mask)
                                tie_disagreement=int((original_rank is not None and original_rank<=budget)!=bool(hit)) if policy=='max_score' else 0
                            base={'item_id':item['construct_id'],'sequence_sha256':item['sequence_sha256'],
                                  'payload':item['payload_accession'],'split':item['split'],'stratum':item['stratum'],
                                  'threshold':c['threshold'],'model':model,'arm':arm,'aggregation':policy,'budget':budget,
                                  'n_targets':len(c['targets']),'n_references':len(ref_ids),'n_allowed':len(allowed),
                                  'best_identity_left':c['best_identity_left'],'hit':hit,'rr':rr,'target_rank_in_output':rank,
                                  'rank_semantics':'minimum pessimistic rank across separate lists' if policy=='per_segment_diagnostic' else 'first target rank in actual fixed output',
                                  'review_slots':slots,'unique_returns':len(returned),
                                  'non_target_returns':sum(i not in c['targets'] for i in returned),
                                  'abstained':int(not returned),'n_segments':len(spans),'sequence_coverage':cov,
                                  'max_target_span_coverage':target_cov,'max_segment_target_purity':target_purity,
                                  'fallback':int(arm=='bhs' and c['bhs_fallback']),
                                  'tie_disagreement':tie_disagreement}
                            records.append(base)
                            details.append({**base,'returned_accessions':[ref_ids[i] for i in returned],
                                            'per_segment_worst_tie_ranks':per_ranks,
                                            'segment_coordinates':[[s['start'],s['end']] for s in spans]})
    out=ROOT/'results'/('comparison_'+split+'_'+'_'.join(models));out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'per_item.csv',records)
    atomic_json(out/'returned_lists.json',details)
    keys=['split','stratum','threshold','model','arm','aggregation','budget']
    write_csv(out/'summary.csv',summarize(records,keys))
    write_csv(out/'per_payload.csv',summarize(records,keys+['payload']))
    # Paired policy changes; never turn correlated items into independent trials.
    lookup={(r['item_id'],r['threshold'],r['model'],r['arm'],r['budget'],r['aggregation']):r for r in records}
    pairs=[]
    for r in records:
        if r['aggregation']!='max_score':continue
        b=lookup[(r['item_id'],r['threshold'],r['model'],r['arm'],r['budget'],'round_robin')]
        diag=lookup[(r['item_id'],r['threshold'],r['model'],r['arm'],r['budget'],'per_segment_diagnostic')]
        pairs.append({k:r[k] for k in ['item_id','payload','split','stratum','threshold','model','arm','budget']}|
                     {'max_hit':r['hit'],'round_robin_hit':b['hit'],'separate_lists_hit':diag['hit'],
                      'gain':int(not r['hit'] and b['hit']),'loss':int(r['hit'] and not b['hit']),
                      'absent_reference':int(r['n_targets']==0),'boundary_target_coverage':r['max_target_span_coverage'],
                      'candidate_local_hit_lost_by_max':int(diag['hit'] and not r['hit']),
                      'round_robin_recovers_local_hit':int(diag['hit'] and not r['hit'] and b['hit'])})
    write_csv(out/'paired_aggregation.csv',pairs)
    atomic_json(out/'manifest.json',{'freeze_sha256':freeze,'models':models,'split':split,'n_metric_rows':len(records),
                'n_items':len({r['item_id'] for r in records}),
                'tie_policy':'fixed lists descending score then accession; original per-segment diagnostic pessimistic tied target rank',
                'scope':cfg['claim_scope'],'complete_recovery':cfg['complete_recovery_limit'],
                'false_positive_limit':cfg['negative_label_limit'],
                'all_absent_targets_kept':True,'calibrated_abstention_threshold':None})
    print(out)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--models',nargs='+',default=['alignment','esm2_650m']);p.add_argument('--split',default='all',choices=['all','development','evaluation']);a=p.parse_args();main(a.models,a.split)
