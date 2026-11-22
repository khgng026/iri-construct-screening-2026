"""Fixed review-budget adapters; no labels are accepted by ranking functions.

Nishit Kalani's upstream pilot code supplies spans, embedding and target-rank
rules (MIT; see sources/nishit-working/LICENSE). This module adds a total budget
and audit records without changing those upstream functions or benchmark files.
"""
import hashlib
import json
from pathlib import Path
import numpy as np


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def rank_lists(scores, reference_ids, allowed):
    scores = np.asarray(scores)
    if scores.ndim != 2 or scores.shape[1] != len(reference_ids):
        raise ValueError('Expected complete segment by reference matrix')
    if len(reference_ids) != len(set(reference_ids)):
        raise ValueError('Duplicate reference IDs')
    if len(allowed) != len(set(allowed)) or any(i < 0 or i >= len(reference_ids) for i in allowed):
        raise ValueError('Invalid eligible reference indices')
    if not np.isfinite(scores[:, allowed]).all():
        raise ValueError('Nonfinite scores require abstention, never dropping a reference')
    return [sorted(allowed, key=lambda i: (-float(row[i]), reference_ids[i])) for row in scores]


def fixed_ranking(scores, reference_ids, allowed, budget, policy='max_score'):
    if type(budget) is not int or budget < 1:
        raise ValueError('Positive integer total budget required')
    lists = rank_lists(scores, reference_ids, allowed)
    if not lists:
        return []
    if policy == 'max_score':
        pooled = np.max(scores, axis=0)
        return sorted(allowed, key=lambda i: (-float(pooled[i]), reference_ids[i]))[:budget]
    if policy != 'round_robin':
        raise ValueError('Unknown frozen aggregation policy')
    # Coordinate-ordered segment lists, one novel hit from each on every round.
    # Duplicate references do not consume budget; no target or boundary labels used.
    positions, out, seen = [0] * len(lists), [], set()
    while len(out) < min(budget, len(allowed)):
        changed = False
        for j, ranking in enumerate(lists):
            while positions[j] < len(ranking) and ranking[positions[j]] in seen:
                positions[j] += 1
            if positions[j] < len(ranking):
                i = ranking[positions[j]]
                positions[j] += 1
                out.append(i)
                seen.add(i)
                changed = True
                if len(out) == budget:
                    break
        if not changed:
            break
    return out


def coverage(spans, length):
    covered = set()
    for start, end in spans:
        if not 0 <= start < end <= length:
            raise ValueError('Invalid 0-based half-open span')
        covered.update(range(start, end))
    return len(covered) / length if length else 0.0
