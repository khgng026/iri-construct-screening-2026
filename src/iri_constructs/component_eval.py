"""Component-attributed retrieval evaluation for multi-unit parents (pure Python).

A parent construct contributes several query units (validated proteins). Every
unit is scored against the same reference panel. These functions turn the
per-unit score rows into one parent-level hit list under a fixed total budget
and say *which unit produced each hit*, so that a family is only counted as
characterized when the component that belongs to it retrieved it.

No model, sequence or file access happens here.
"""
from __future__ import annotations

import math

from .runtime import InputError


def _check(unit_scores: dict, k: int) -> list[str]:
    if type(k) is not int or k < 1:
        raise InputError('Hit budget must be a positive integer')
    if not unit_scores or any(not row for row in unit_scores.values()):
        raise InputError('Each unit needs a nonempty score row')
    references = set(next(iter(unit_scores.values())))
    for row in unit_scores.values():
        if set(row) != references:
            raise InputError('All units of a parent must be scored against the same reference panel')
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in row.values()):
            raise InputError('Scores must be finite numbers')
    return sorted(references)


def ranked(row: dict) -> list[str]:
    """References by descending score; ties broken by reference ID."""
    return sorted(row, key=lambda r: (-row[r], r))


def merge_single(unit_scores: dict, unit: str, k: int) -> list[dict]:
    _check(unit_scores, k)
    row = unit_scores[unit]
    return [{'reference_id': r, 'origin_unit': unit, 'score': row[r], 'origin_rank': i + 1}
            for i, r in enumerate(ranked(row)[:k])]


def merge_max(unit_scores: dict, k: int) -> list[dict]:
    """Historical rule: best score per reference over all units, then one top-k list.

    The origin of a hit is the unit holding that best score (ties: unit ID).
    """
    references = _check(unit_scores, k)
    best = {}
    for r in references:
        unit = min(unit_scores, key=lambda u: (-unit_scores[u][r], u))
        best[r] = (unit_scores[unit][r], unit)
    order = sorted(references, key=lambda r: (-best[r][0], r))[:k]
    ranks = {u: {r: i + 1 for i, r in enumerate(ranked(row))} for u, row in unit_scores.items()}
    return [{'reference_id': r, 'origin_unit': best[r][1], 'score': best[r][0],
             'origin_rank': ranks[best[r][1]][r]} for r in order]


def merge_round_robin(unit_scores: dict, unit_order: list[str], k: int) -> list[dict]:
    """Component-balanced rule: units take turns contributing their next-best unseen reference.

    ``unit_order`` is fixed by the caller before scoring. The same total budget
    ``k`` applies, so with more units than ``k`` later units contribute nothing.
    Scores are never compared across units.
    """
    _check(unit_scores, k)
    if sorted(unit_order) != sorted(unit_scores):
        raise InputError('unit_order must list every scored unit exactly once')
    orders = {u: ranked(unit_scores[u]) for u in unit_order}
    position = {u: 0 for u in unit_order}
    hits, taken = [], set()
    while len(hits) < k:
        progressed = False
        for u in unit_order:
            while position[u] < len(orders[u]) and orders[u][position[u]] in taken:
                position[u] += 1
            if position[u] < len(orders[u]):
                r = orders[u][position[u]]
                taken.add(r)
                hits.append({'reference_id': r, 'origin_unit': u, 'score': unit_scores[u][r],
                             'origin_rank': position[u] + 1})
                position[u] += 1
                progressed = True
                if len(hits) == k:
                    break
        if not progressed:
            break
    return hits


def attributed_recall(hits: list[dict], unit_families: dict, reference_family: dict) -> dict:
    """Score supported components by hits that *they* produced.

    ``unit_families`` maps each supported unit to the set of families it is
    documented to contain. Units absent from it are unlabeled: their hits use
    budget and are counted as unscorable, never as right or wrong.
    """
    recovered, composition = set(), {'attributed_correct': 0, 'attributed_wrong': 0, 'from_unlabeled_unit': 0}
    for hit in hits:
        families = unit_families.get(hit['origin_unit'])
        if families is None:
            composition['from_unlabeled_unit'] += 1
        elif reference_family[hit['reference_id']] in families:
            composition['attributed_correct'] += 1
            recovered.add(hit['origin_unit'])
        else:
            composition['attributed_wrong'] += 1
    supported = sorted(unit_families)
    return {'supported_components': len(supported), 'recovered_components': len(recovered),
            'recovered_units': sorted(recovered),
            'attributed_recall': len(recovered) / len(supported) if supported else None,
            'hit_composition': composition, 'returned_hits': len(hits)}


def family_recall(hits: list[dict], unit_families: dict, reference_family: dict) -> dict:
    """Historical-style endpoint: an expected family counts whichever unit retrieved it.

    ``unattributed_credit`` lists expected families present among the hits only
    through a unit that does not belong to that family.
    """
    expected = set().union(*unit_families.values()) if unit_families else set()
    present = {reference_family[h['reference_id']] for h in hits} & expected
    through_owner = {reference_family[h['reference_id']] for h in hits
                     if reference_family[h['reference_id']] in unit_families.get(h['origin_unit'], ())}
    return {'expected_families': len(expected), 'recovered_families': len(present),
            'family_recall': len(present) / len(expected) if expected else None,
            'unattributed_credit': sorted(present - through_owner)}


def unit_retrieval(row: dict, families: set, reference_family: dict, k: int) -> dict:
    """One unit queried alone: where does its first same-family reference rank?"""
    order = ranked(row)
    first = next((i + 1 for i, r in enumerate(order) if reference_family[r] in families), None)
    return {'first_same_family_rank': first, 'hit_at_1': first == 1, 'hit_at_k': first is not None and first <= k,
            'reciprocal_rank': 1 / first if first else 0.0,
            'same_family_in_top_k': sum(reference_family[r] in families for r in order[:k])}
