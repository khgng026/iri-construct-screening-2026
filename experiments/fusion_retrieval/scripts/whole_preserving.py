"""Frozen total-budget mergers; target labels never enter this module."""
import numpy as np
from budget_compare import rank_lists


def merge(scores, whole_scores, reference_ids, allowed, budget, policy):
    if type(budget) is not int or budget < 1:
        raise ValueError('Positive integer budget required')
    scores=np.asarray(scores)
    if scores.ndim==1 and scores.size==0:
        scores=np.empty((0,len(reference_ids)))
    whole_scores=np.asarray(whole_scores)
    if whole_scores.shape!=(len(reference_ids),):
        raise ValueError('Incomplete whole-query scores')
    whole=rank_lists(whole_scores[None,:],reference_ids,allowed)[0]
    lists=rank_lists(scores,reference_ids,allowed)
    if policy not in {'max_score','round_robin','whole_half'}:
        raise ValueError('Unknown frozen policy')
    if not lists:
        return [(i,-1) for i in whole[:budget]]
    if policy=='max_score':
        pooled=scores.max(axis=0)
        ranking=sorted(allowed,key=lambda i:(-float(pooled[i]),reference_ids[i]))[:budget]
        # Earliest span is the deterministic insertion source; all tied supports
        # are recorded separately by the evaluation exporter.
        return [(i,int(np.flatnonzero(scores[:,i]==pooled[i])[0])) for i in ranking]
    out=[];seen=set()
    if policy=='whole_half':
        for i in whole[:budget//2]:
            out.append((i,-1));seen.add(i)
    positions=[0]*len(lists)
    while len(out)<min(budget,len(allowed)):
        progress=False
        for j,ranking in enumerate(lists):
            while positions[j]<len(ranking) and ranking[positions[j]] in seen:
                positions[j]+=1
            if positions[j]<len(ranking):
                i=ranking[positions[j]];positions[j]+=1
                out.append((i,j));seen.add(i);progress=True
                if len(out)==budget:break
        if not progress:break
    return out
