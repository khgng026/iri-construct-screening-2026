"""Within-model representation metrics and independent-family retrieval."""
from __future__ import annotations

import math
from .runtime import InputError, require_text


def mean_pool(rows: list[list[float]]) -> list[float]:
    if not rows or not rows[0] or any(len(r) != len(rows[0]) for r in rows):
        raise InputError('Empty or inconsistent residue representation')
    result = [math.fsum(column) / len(rows) for column in zip(*rows)]
    if not all(math.isfinite(x) for x in result):
        raise InputError('Representation contains nonfinite values')
    return result


def cosine(a: list[float], b: list[float]) -> float:
    if not a or len(a) != len(b) or not all(math.isfinite(x) for x in a + b):
        raise InputError('Cannot compare incompatible representations')
    norm = math.sqrt(math.fsum(x*x for x in a) * math.fsum(x*x for x in b))
    if norm <= 1e-12:
        raise InputError('Cannot compare zero-norm representations')
    return max(-1.0, min(1.0, math.fsum(x*y for x, y in zip(a, b)) / norm))


def reviewed_labels(config: dict, proteins: dict) -> tuple[dict, list[dict]]:
    entries = config.get('families', [])
    if not isinstance(entries, list):
        raise InputError('families must be a list')
    labels, excluded, seen = {}, [], set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise InputError('Family entries must be objects')
        iid = require_text(entry, 'input_id')
        if iid in seen:
            raise InputError('Duplicate family annotation for input_id')
        seen.add(iid)
        if (entry.get('annotation_status') not in ('known', 'unknown', 'ambiguous', 'pending_review') or
                entry.get('review_status') not in ('pending', 'approved', 'uncertain', 'rejected')):
            raise InputError('Family annotation/review status is invalid')
        if iid not in proteins:
            excluded.append({'input_id': iid, 'reason': 'No validated protein unit for family annotation'})
            continue
        if entry.get('annotation_status') != 'known' or entry.get('review_status') != 'approved':
            excluded.append({'input_id': iid, 'reason': 'Family annotation is unknown, ambiguous or unreviewed'})
            continue
        for key in ('family', 'evidence', 'independence_group'):
            require_text(entry, key)
        labels[iid] = entry
    return labels, excluded


def candidates(query: str, labels: dict, proteins: dict) -> list[str]:
    """Preserve prototype leave-construct-out policy and exclude exact copies/linked groups."""
    q = proteins[query]
    return sorted(r for r in labels if r != query
                  and proteins[r].construct_id != q.construct_id
                  and proteins[r].sequence_sha256 != q.sequence_sha256
                  and labels[r]['independence_group'] != labels[query]['independence_group'])


def retrieve(vectors: dict, labels: dict, proteins: dict, query_ids=None,
             reference_ids=None) -> dict:
    """Use a fixed panel; a missing embedding causes abstention, never a smaller reference pool."""
    queries = sorted(labels if query_ids is None else query_ids)
    panel = set(labels if reference_ids is None else reference_ids)
    predictions, abstentions = [], []
    for q in queries:
        refs = [r for r in candidates(q, labels, proteins) if r in panel]
        same = [r for r in refs if labels[r]['family'] == labels[q]['family']]
        other = [r for r in refs if labels[r]['family'] != labels[q]['family']]
        reason = None
        if not same or not other:
            reason = 'No independent same-family reference and competing family'
        elif q not in vectors or any(r not in vectors for r in refs):
            reason = 'Missing query/reference representation in the fixed panel'
        if reason:
            abstentions.append({'input_id': q, 'reason': reason})
            continue
        ranked = sorted(((cosine(vectors[q], vectors[r]), r) for r in refs), key=lambda x: (-x[0], x[1]))
        value, neighbor = ranked[0]
        predictions.append({'input_id': q, 'family': labels[q]['family'],
                            'retrieved_family': labels[neighbor]['family'], 'neighbor': neighbor,
                            'cosine': value, 'margin_same_minus_other':
                            max(cosine(vectors[q], vectors[r]) for r in same) -
                            max(cosine(vectors[q], vectors[r]) for r in other),
                            'reference_ids': refs,
                            'tie_at_top': sum(abs(s-value) < 1e-12 for s, _ in ranked) > 1})
    f1 = []
    for family in sorted({p['family'] for p in predictions} | {p['retrieved_family'] for p in predictions}):
        tp = sum(p['family'] == family and p['retrieved_family'] == family for p in predictions)
        fp = sum(p['family'] != family and p['retrieved_family'] == family for p in predictions)
        fn = sum(p['family'] == family and p['retrieved_family'] != family for p in predictions)
        f1.append(2*tp/(2*tp+fp+fn))
    return {'predictions': predictions, 'abstentions': abstentions, 'query_ids': queries,
            'reference_ids': sorted(panel), 'scored': len(predictions), 'requested': len(queries),
            'coverage': len(predictions)/len(queries) if queries else 0.0,
            'annotation_family_accuracy': sum(p['family'] == p['retrieved_family'] for p in predictions)/len(predictions) if predictions else None,
            'macro_f1': math.fsum(f1)/len(f1) if f1 else None}
