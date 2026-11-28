"""Grouping and paired analysis for the Q2 reporter replication.

``group`` clusters documented components into lineage groups by within-family protein
identity (single linkage at >= 90%). Lineage groups are the primary equal-weight unit;
the family label is the stricter grouping used for sensitivity. Output is sequence-free.

``analyze`` consumes the frozen grouping and the verifier report and computes, per arm:
pooled and equal-lineage-group recall, the gene-caller minus minimal gain with a
whole-group bootstrap 95% interval, the reporter/marker role split, grouping and budget
sensitivity, the coverage funnel, and resource use.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

REPORTER_FAMILIES = {'gfp_like', 'firefly_luciferase', 'renilla_luciferase', 'beta_galactosidase'}
IDENTITY_THRESHOLD = 0.90


def read_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def write_json(path: Path, value, *, exist_ok=False):
    with path.open('w' if exist_ok else 'x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def _identity(aligner, a: str, b: str) -> float:
    if a == b:
        return 1.0
    alignment = aligner.align(a, b)[0]
    try:
        counts = alignment.counts()
        identities = counts.identities
    except AttributeError:  # pragma: no cover - older Biopython
        top, bottom = str(alignment).split('\n')[0], str(alignment).split('\n')[2]
        identities = sum(x == y and x != '-' for x, y in zip(top, bottom))
    return identities / max(len(a), len(b))


class _Union:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def group(answer_key: Path, annotation_units: Path, out: Path):
    from Bio.Align import PairwiseAligner, substitution_matrices
    aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-10.0, extend_gap_score=-0.5)

    key = read_json(answer_key)
    families, proteins = {}, {}
    for construct in key['constructs']:
        for component in construct['components']:
            families[component['id']] = component['family']
    for construct_id, units in read_json(annotation_units).items():
        for unit in units:
            proteins[unit['id']] = unit['protein']

    assignments = {}
    by_family: dict[str, list[str]] = {}
    for cid, family in families.items():
        by_family.setdefault(family, []).append(cid)

    lineage = {}
    for family, ids in sorted(by_family.items()):
        ids = sorted(ids)
        union = _Union(ids)
        for i, a in enumerate(ids):
            for b in ids[i + 1:]:
                if _identity(aligner, proteins[a], proteins[b]) >= IDENTITY_THRESHOLD:
                    union.union(a, b)
        roots = {}
        for cid in ids:
            root = union.find(cid)
            roots.setdefault(root, []).append(cid)
        for k, (root, members) in enumerate(sorted(roots.items(), key=lambda kv: min(kv[1])), 1):
            for cid in members:
                lineage[cid] = f'{family}#{k}'

    for cid, family in families.items():
        assignments[cid] = {'family': family, 'role': 'reporter' if family in REPORTER_FAMILIES else 'marker',
                            'lineage_group': lineage[cid], 'family_group': family}

    n_lineage = len({v['lineage_group'] for v in assignments.values()})
    write_json(out, {'study': 'q2-reporter-replication', 'identity_threshold': IDENTITY_THRESHOLD,
                     'components': len(assignments), 'lineage_groups': n_lineage,
                     'family_groups': len(by_family), 'assignments': assignments}, exist_ok=True)
    print(f'grouped {len(assignments)} components into {n_lineage} lineage groups '
          f'({len(by_family)} families)')


# --------------------------------------------------------------------------- analysis


def _recovered_by_cap(run_report):
    """method -> cap -> set(recovered component ids), plus resource + funnel sums."""
    out = {}
    for method, caps in run_report['runs'].items():
        out[method] = {}
        for cap, payload in caps.items():
            recovered, resources = set(), {
                'queries': 0, 'input_residues': 0, 'candidate_count': 0,
                'characterize_seconds': 0.0, 'preprocess_seconds': 0.0,
                'components_without_covering_query': 0, 'covered_components_not_recovered': 0}
            for row in payload['per_construct']:
                recovered.update(row['recovered_component_ids'])
                for k in ('queries', 'input_residues', 'candidate_count', 'characterize_seconds',
                          'preprocess_seconds', 'components_without_covering_query',
                          'covered_components_not_recovered'):
                    resources[k] += row[k]
            out[method][int(cap)] = {'recovered': recovered, 'resources': resources,
                                     'eligible': payload['overall']['eligible_components']}
    return out


def _subset(assignments, role):
    if role == 'all':
        return set(assignments)
    return {cid for cid, a in assignments.items() if a['role'] == role}


def _group_map(assignments, ids, key):
    groups: dict[str, list[str]] = {}
    for cid in ids:
        groups.setdefault(assignments[cid][key], []).append(cid)
    return groups


def _pooled(recovered, ids):
    return len(recovered & ids) / len(ids) if ids else 0.0


def _equal_group(recovered, groups):
    if not groups:
        return 0.0
    return sum(len(recovered & set(m)) / len(m) for m in groups.values()) / len(groups)


def _bootstrap_gain(recovered_a, recovered_b, groups, reps, seed):
    """Equal-group and pooled gain (a - b) resampling whole groups with replacement."""
    keys = list(groups)
    rng = random.Random(seed)
    eq, pooled = [], []
    members = {k: set(groups[k]) for k in keys}
    for _ in range(reps):
        sample = [keys[rng.randrange(len(keys))] for _ in keys]
        eq_a = sum(len(recovered_a & members[k]) / len(members[k]) for k in sample) / len(sample)
        eq_b = sum(len(recovered_b & members[k]) / len(members[k]) for k in sample) / len(sample)
        eq.append(eq_a - eq_b)
        tot = sum(len(members[k]) for k in sample)
        pa = sum(len(recovered_a & members[k]) for k in sample) / tot
        pb = sum(len(recovered_b & members[k]) for k in sample) / tot
        pooled.append(pa - pb)
    eq.sort(); pooled.sort()
    lo, hi = int(0.025 * reps), int(0.975 * reps) - 1
    return {'equal_group_interval': [eq[lo], eq[hi]], 'pooled_interval': [pooled[lo], pooled[hi]]}


def _effect(recovered_a, recovered_b, assignments, ids, group_key, reps, seed):
    groups = _group_map(assignments, ids, group_key)
    eq = _equal_group(recovered_a, groups) - _equal_group(recovered_b, groups)
    pooled = _pooled(recovered_a, ids) - _pooled(recovered_b, ids)
    ci = _bootstrap_gain(recovered_a, recovered_b, groups, reps, seed)
    return {'equal_group_gain': eq, 'pooled_gain': pooled, 'groups': len(groups), 'components': len(ids),
            **ci}


def analyze(report: Path, groups: Path, out: Path, baseline='minimal', candidate='gene_caller',
            reps=4000, seed=20261124):
    run = read_json(report)
    grouping = read_json(groups)
    assignments = grouping['assignments']
    data = _recovered_by_cap(run)
    protocol = run['protocol']
    primary_cap = protocol['primary_query_cap']
    target = protocol['minimum_interesting_gain']

    arms = sorted(data)
    recall_table = {}
    for method in arms:
        recall_table[method] = {}
        for cap, payload in data[method].items():
            rec = payload['recovered']
            row = {}
            for role in ('all', 'reporter', 'marker'):
                ids = _subset(assignments, role)
                row[role] = {'pooled': _pooled(rec, ids),
                             'equal_lineage_group': _equal_group(rec, _group_map(assignments, ids, 'lineage_group')),
                             'recovered': len(rec & ids), 'eligible': len(ids)}
            row['resources'] = payload['resources']
            recall_table[method][cap] = row

    rec_c = data[candidate]
    rec_b = data[baseline]
    effects = {}
    for role in ('all', 'reporter', 'marker'):
        ids = _subset(assignments, role)
        effects[role] = {}
        for cap in sorted(rec_b):
            e = _effect(rec_c[cap]['recovered'], rec_b[cap]['recovered'], assignments, ids,
                        'lineage_group', reps, seed)
            if cap == primary_cap:
                e['meets_target'] = e['equal_group_interval'][0] >= target
            effects[role][cap] = e

    # Sensitivity at the primary cap, all components.
    ids_all = _subset(assignments, 'all')
    a, b = rec_c[primary_cap]['recovered'], rec_b[primary_cap]['recovered']
    lineage_groups = _group_map(assignments, ids_all, 'lineage_group')
    dominant = max(lineage_groups, key=lambda k: len(lineage_groups[k]))
    ids_no_dom = {cid for cid in ids_all if assignments[cid]['lineage_group'] != dominant}
    sensitivity = {
        'family_group_weighting': _effect(a, b, assignments, ids_all, 'family_group', reps, seed),
        'omit_dominant_lineage': {'dropped_group': dominant, 'dropped_components': len(lineage_groups[dominant]),
                                  **_effect(a, b, assignments, ids_no_dom, 'lineage_group', reps, seed + 1)},
        'reporter_only_lineage': _effect(a, b, assignments, _subset(assignments, 'reporter'),
                                         'lineage_group', reps, seed + 2),
    }

    summary = {
        'study': 'q2-reporter-replication', 'baseline': baseline, 'candidate': candidate,
        'primary_cap': primary_cap, 'target': target,
        'denominators': {'components': grouping['components'], 'lineage_groups': grouping['lineage_groups'],
                         'family_groups': grouping['family_groups'],
                         'reporter_components': len(_subset(assignments, 'reporter')),
                         'marker_components': len(_subset(assignments, 'marker')),
                         'constructs': len(run['runs'][baseline][str(primary_cap)]['per_construct'])},
        'recall': recall_table, 'effects': effects, 'sensitivity': sensitivity,
        'bootstrap': {'replicates': reps, 'seed': seed, 'interval': 0.95},
    }
    write_json(out, summary, exist_ok=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)

    g = sub.add_parser('group')
    g.add_argument('--answer-key', type=Path, required=True)
    g.add_argument('--annotation-units', type=Path, required=True)
    g.add_argument('--out', type=Path, required=True)

    a = sub.add_parser('analyze')
    a.add_argument('--report', type=Path, required=True)
    a.add_argument('--groups', type=Path, required=True)
    a.add_argument('--out', type=Path, required=True)

    args = parser.parse_args()
    if args.stage == 'group':
        group(args.answer_key, args.annotation_units, args.out)
    else:
        summary = analyze(args.report, args.groups, args.out)
        print(json.dumps({k: summary[k] for k in ('denominators',)}, indent=2))


if __name__ == '__main__':
    main()
