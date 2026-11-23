"""Analyze crossed cohort traces against one evaluator key, without rerunning inference.

Inputs are the traces written by ``cohort_study run`` for both scorers, the evaluator
answer key, the private construct map (groups, families, tiers) and the group
manifests. The verifier scores each scorer's traces; this module adds the paired
group analysis, uncertainty intervals, failure breakdown and strata.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .verify import _covered, _location, _require, verify

SCORERS = ('esm2', 'blosum62')
BASELINE = 'minimal'
CANDIDATES = ('start_ranked', 'gene_caller')
REFERENCE = 'annotation_reference'
MAIN_ARMS = (BASELINE, *CANDIDATES, REFERENCE)
BOOTSTRAP = {'replicates': 4000, 'seed': 20261123, 'interval': 0.95}
# Decision thresholds fixed in cohort_protocol.json; repeated here for the verdict text.
TARGET = 0.05


def read(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def load_runs(run_dir: Path, arms):
    protocols = {s: read(run_dir / s / 'protocol.json') for s in SCORERS}
    traces = {s: {a: read(run_dir / s / f'{a}.json') for a in arms} for s in SCORERS}
    return protocols, traces


def check_matched_inputs(protocols, traces, arms):
    """Both scorers must see the same protocol settings and the same selected inputs for every arm."""
    common = [{k: v for k, v in protocols[s].items() if k not in ('characterizer', 'characterization_id')}
              for s in SCORERS]
    _require(common[0] == common[1], 'scorers must share the evaluation protocol')
    fingerprints = {}
    for arm in arms:
        left, right = (traces[s][arm] for s in SCORERS)
        for s in SCORERS:
            _require(traces[s][arm]['characterizer'] == s and traces[s][arm]['method'] == arm, 'wrong trace')
        _require(left['selected_inputs_sha256'] == right['selected_inputs_sha256']
                 and left['input_sha256'] == right['input_sha256'], f'{arm}: selected inputs differ between scorers')
        _require(left['candidate_pools'] == right['candidate_pools'], f'{arm}: candidate pools differ')
        fingerprints[arm] = left['selected_inputs_sha256']
    return fingerprints


def component_table(protocols, key, traces, construct_map, arms):
    """One row per eligible documented component with its outcome in every scorer, arm and cap."""
    reports = {s: verify(protocols[s], key, traces[s][BASELINE],
                         [traces[s][a] for a in arms if a not in (BASELINE, REFERENCE)],
                         traces[s].get(REFERENCE)) for s in SCORERS}
    overlap = protocols['esm2']['minimum_component_overlap']
    info = {row['id']: row for row in construct_map}
    rows = []
    for construct in key['constructs']:
        cid, length = construct['id'], construct['length']
        meta = info[cid]
        components = {c['id']: c for c in meta['components']}
        for component in construct['components']:
            locus = _location(component, length, component['id'])
            extra = components[component['id']]
            row = {'construct': cid, 'component': component['id'], 'family': component['family'],
                   'archetype': construct['archetype'], 'lineage_group': meta['lineage_group'],
                   'conservative_group': meta['conservative_group'],
                   'similarity_tier': extra['similarity']['tier'], 'fusion_class': extra['fusion_class'],
                   'protein_sha256': extra['protein_sha256'], 'outcome': {}}
            for s in SCORERS:
                for arm in arms:
                    trace = traces[s][arm]
                    pool = trace['candidate_pools'][cid]
                    in_pool = any(_covered(locus, (u['strand'], sorted(map(tuple, u['spans'])))) >= overlap
                                  for u in pool)
                    for tier in trace['tiers']:
                        cap = tier['query_cap']
                        trow = next(r for r in tier['constructs'] if r['id'] == cid)
                        scored = next(r for r in reports[s]['runs'][arm][cap]['per_construct']
                                      if r['construct_id'] == cid)
                        covering = [q['id'] for q in trow['queries']
                                    if _covered(locus, (q['strand'], sorted(map(tuple, q['spans'])))) >= overlap]
                        ranks = [trow['query_family_rankings'][q].index(component['family']) + 1
                                 for q in covering if component['family'] in trow['query_family_rankings'][q]]
                        recovered = component['id'] in scored['recovered_component_ids']
                        if recovered:
                            status = 'recovered'
                        elif covering:
                            status = 'covered_label_missed'
                        elif in_pool:
                            status = 'compatible_candidate_below_cap'
                        else:
                            status = 'no_compatible_candidate'
                        row['outcome'][(s, arm, cap)] = {'recovered': recovered, 'status': status,
                                                         'best_rank_in_covering_query': min(ranks) if ranks else None}
            rows.append(row)
    return rows, reports


# ------------------------------------------------------------- estimators

def _group_arrays(rows, group_field, a, b):
    groups = sorted({r[group_field] for r in rows})
    index = {g: i for i, g in enumerate(groups)}
    n, sa, sb = np.zeros(len(groups)), np.zeros(len(groups)), np.zeros(len(groups))
    for r in rows:
        i = index[r[group_field]]
        n[i] += 1
        sa[i] += r['outcome'][a]['recovered']
        sb[i] += r['outcome'][b]['recovered']
    return groups, n, sa, sb


def _interval(values):
    alpha = (1 - BOOTSTRAP['interval']) / 2
    return [float(np.quantile(values, alpha)), float(np.quantile(values, 1 - alpha))]


def paired_effect(rows, a, b, group_field='lineage_group', family_bootstrap=True):
    """Recall of a minus recall of b on the same components, pooled and with equal group weights."""
    if not rows:
        return None
    groups, n, sa, sb = _group_arrays(rows, group_field, a, b)
    d = sa - sb
    pooled = float(d.sum() / n.sum())
    equal = float(np.mean(d / n))
    out = {'components': int(n.sum()), 'groups': len(groups),
           'largest_group_share': float(n.max() / n.sum()),
           'recall_a_pooled': float(sa.sum() / n.sum()), 'recall_b_pooled': float(sb.sum() / n.sum()),
           'recall_a_equal_group': float(np.mean(sa / n)), 'recall_b_equal_group': float(np.mean(sb / n)),
           'pooled_gain': pooled, 'equal_group_gain': equal,
           'groups_improved': int((d > 0).sum()), 'groups_declined': int((d < 0).sum()),
           'groups_tied': int((d == 0).sum()),
           'components_only_a': sum(r['outcome'][a]['recovered'] and not r['outcome'][b]['recovered'] for r in rows),
           'components_only_b': sum(r['outcome'][b]['recovered'] and not r['outcome'][a]['recovered'] for r in rows)}
    if len(groups) > 1:
        omit = (d.sum() - d) / (n.sum() - n)
        omit_equal = (equal * len(groups) - d / n) / (len(groups) - 1)
        out['leave_one_group_out_pooled_range'] = [float(omit.min()), float(omit.max())]
        out['leave_one_group_out_equal_group_range'] = [float(omit_equal.min()), float(omit_equal.max())]
        largest = int(np.argmax(n))
        out['pooled_gain_without_largest_group'] = float(omit[largest])
        rng = np.random.default_rng(BOOTSTRAP['seed'])
        draws = rng.integers(0, len(groups), size=(BOOTSTRAP['replicates'], len(groups)))
        out['group_bootstrap'] = {
            'pooled_gain_interval': _interval(d[draws].sum(1) / n[draws].sum(1)),
            'equal_group_gain_interval': _interval((d / n)[draws].mean(1))}
    if family_bootstrap and len(groups) > 1:
        out['two_way_bootstrap_pooled_gain_interval'] = _two_way(rows, group_field, a, b)
    return out


def _two_way(rows, group_field, a, b):
    """Pigeonhole bootstrap: resample groups and families independently and weight each cell by both counts."""
    groups = sorted({r[group_field] for r in rows})
    families = sorted({r['family'] for r in rows})
    gi, fi = {g: i for i, g in enumerate(groups)}, {f: i for i, f in enumerate(families)}
    cells = defaultdict(lambda: [0, 0])
    for r in rows:
        cell = cells[(gi[r[group_field]], fi[r['family']])]
        cell[0] += 1
        cell[1] += r['outcome'][a]['recovered'] - r['outcome'][b]['recovered']
    keys = list(cells)
    g_idx = np.array([k[0] for k in keys])
    f_idx = np.array([k[1] for k in keys])
    n = np.array([cells[k][0] for k in keys], dtype=float)
    d = np.array([cells[k][1] for k in keys], dtype=float)
    rng = np.random.default_rng(BOOTSTRAP['seed'] + 1)
    values = []
    for _ in range(BOOTSTRAP['replicates']):
        wg = np.bincount(rng.integers(0, len(groups), len(groups)), minlength=len(groups))
        wf = np.bincount(rng.integers(0, len(families), len(families)), minlength=len(families))
        w = wg[g_idx] * wf[f_idx]
        if w @ n:
            values.append(float((w @ d) / (w @ n)))
    return _interval(values)


def verdict(effect, target=TARGET):
    """Apply the locked decision rule to the equal-group-weighted paired gain."""
    if effect is None or 'group_bootstrap' not in effect:
        return 'descriptive only: too few groups for an interval'
    low, high = effect['group_bootstrap']['equal_group_gain_interval']
    point = effect['equal_group_gain']
    consistent = effect['leave_one_group_out_equal_group_range'][0] >= target
    if low >= target and consistent:
        return 'practically meaningful benefit (interval lower bound at least the target)'
    if low >= target:
        return 'interval above target but sensitive to omitting one group'
    if high < target:
        return 'gain of the target size unsupported (interval upper bound below the target)'
    if point > 0 and low > 0:
        return 'positive signal; interval spans the target, so a practical improvement is not established'
    return 'inconclusive (interval spans zero and the target)'


# --------------------------------------------------------------- summaries

def work_table(traces, key, arms, reference_count):
    constructs = len(key['constructs'])
    out = {}
    for s in SCORERS:
        for arm in arms:
            trace = traces[s][arm]
            for tier in trace['tiers']:
                cap = tier['query_cap']
                rows = tier['constructs']
                queries = sum(len(r['queries']) for r in rows)
                out[f'{s}|{arm}|{cap}'] = {
                    'allowed_inputs': cap * constructs, 'actual_inputs': queries,
                    'input_residues': sum(q['input_residues'] for r in rows for q in r['queries']),
                    'candidate_pool': sum(r['candidate_count'] for r in rows),
                    'reference_comparisons': sum(r['work']['reference_comparisons'] for r in rows),
                    'alignment_cells': sum(r['work']['alignment_cells'] for r in rows) if s == 'blosum62' else None,
                    'encoder_calls': (queries + reference_count) if s == 'esm2' else None,
                    'characterize_seconds': sum(r['characterize_seconds'] for r in rows)}
            out[f'{s}|{arm}|run'] = {'unique_scorer_inputs': trace['unique_scorer_inputs'],
                                     'unique_scoring_seconds': trace['unique_scoring_seconds'],
                                     'wall_seconds': trace['wall_seconds']}
    return out


def failure_breakdown(rows, arms, cap):
    out = {}
    for s in SCORERS:
        for arm in arms:
            statuses = Counter(r['outcome'][(s, arm, cap)]['status'] for r in rows)
            ranks = Counter(r['outcome'][(s, arm, cap)]['best_rank_in_covering_query'] for r in rows
                            if r['outcome'][(s, arm, cap)]['status'] == 'covered_label_missed')
            out[f'{s}|{arm}'] = {'statuses': dict(statuses),
                                 'covered_missed_best_rank': {str(k): v for k, v in sorted(ranks.items(),
                                                                                           key=lambda x: (x[0] is None, x[0] or 0))}}
    # Where do the baseline's misses go in the other arms?
    transitions = {}
    for s in SCORERS:
        for arm in arms:
            if arm == BASELINE:
                continue
            transitions[f'{s}|minimal->{arm}'] = dict(Counter(
                f"{r['outcome'][(s, BASELINE, cap)]['status']}->{r['outcome'][(s, arm, cap)]['status']}"
                for r in rows if not r['outcome'][(s, BASELINE, cap)]['recovered']))
    return {'by_arm': out, 'baseline_miss_transitions': transitions}


def recall_curves(rows, arms, caps, group_field='lineage_group'):
    out = {}
    for s in SCORERS:
        for arm in arms:
            for cap in caps:
                k = (s, arm, cap)
                groups = defaultdict(list)
                for r in rows:
                    groups[r[group_field]].append(r['outcome'][k]['recovered'])
                out[f'{s}|{arm}|{cap}'] = {
                    'recovered': sum(r['outcome'][k]['recovered'] for r in rows), 'eligible': len(rows),
                    'pooled_recall': sum(r['outcome'][k]['recovered'] for r in rows) / len(rows) if rows else None,
                    'equal_group_recall': float(np.mean([np.mean(v) for v in groups.values()])) if rows else None}
    return out


def effects(rows, cap, arms, group_field='lineage_group'):
    e = {}
    for s in SCORERS:
        for arm in arms:
            if arm != BASELINE:
                e[f'{s}: {arm} - minimal'] = paired_effect(rows, (s, arm, cap), (s, BASELINE, cap), group_field)
    for arm in arms:
        e[f'esm2 - blosum62: {arm}'] = paired_effect(rows, ('esm2', arm, cap), ('blosum62', arm, cap), group_field)
    if 'gene_caller' in arms:
        # Difference in preprocessing gain between scorers, as a paired contrast per component.
        e['difference in gain (gene_caller - minimal, esm2 - blosum62)'] = _interaction(rows, cap, group_field)
        e['esm2: annotation_reference - gene_caller'] = paired_effect(
            rows, ('esm2', REFERENCE, cap), ('esm2', 'gene_caller', cap), group_field)
        e['blosum62: annotation_reference - gene_caller'] = paired_effect(
            rows, ('blosum62', REFERENCE, cap), ('blosum62', 'gene_caller', cap), group_field)
    return e


def _interaction(rows, cap, group_field):
    if not rows:
        return None
    groups = sorted({r[group_field] for r in rows})
    index = {g: i for i, g in enumerate(groups)}
    n, d = np.zeros(len(groups)), np.zeros(len(groups))
    for r in rows:
        o = r['outcome']
        i = index[r[group_field]]
        n[i] += 1
        d[i] += ((o[('esm2', 'gene_caller', cap)]['recovered'] - o[('esm2', BASELINE, cap)]['recovered'])
                 - (o[('blosum62', 'gene_caller', cap)]['recovered'] - o[('blosum62', BASELINE, cap)]['recovered']))
    out = {'components': int(n.sum()), 'groups': len(groups), 'pooled': float(d.sum() / n.sum()),
           'equal_group': float(np.mean(d / n))}
    if len(groups) > 1:
        rng = np.random.default_rng(BOOTSTRAP['seed'] + 2)
        draws = rng.integers(0, len(groups), size=(BOOTSTRAP['replicates'], len(groups)))
        out['equal_group_interval'] = _interval((d / n)[draws].mean(1))
        out['pooled_interval'] = _interval(d[draws].sum(1) / n[draws].sum(1))
    return out


def analyze(run_dir: Path, key_path: Path, map_path: Path, arms=MAIN_ARMS):
    protocols, traces = load_runs(run_dir, arms)
    fingerprints = check_matched_inputs(protocols, traces, arms)
    key, construct_map = read(key_path), read(map_path)
    rows, reports = component_table(protocols, key, traces, construct_map, arms)
    protocol = protocols['esm2']
    caps, cap = protocol['query_caps'], protocol['primary_query_cap']
    strata = {}
    for field in ('archetype', 'similarity_tier'):
        for value in sorted({r[field] for r in rows}):
            subset = [r for r in rows if r[field] == value]
            strata[f'{field}={value}'] = {
                'components': len(subset), 'lineage_groups': len({r['lineage_group'] for r in subset}),
                'curves': recall_curves(subset, arms, [cap]),
                'effects': {k: v for k, v in effects(subset, cap, arms).items()
                            if k.endswith('gene_caller - minimal') or k.startswith('esm2 - blosum62')}}
    report = {
        'schema': 'paired-cohort-analysis-1',
        'bootstrap': BOOTSTRAP, 'primary_query_cap': cap, 'target': protocol['minimum_interesting_gain'],
        'matched_inputs_sha256': fingerprints,
        'denominator': {'components': len(rows), 'constructs': len(key['constructs']),
                        'lineage_groups': len({r['lineage_group'] for r in rows}),
                        'conservative_groups': len({r['conservative_group'] for r in rows}),
                        'families': len({r['family'] for r in rows})},
        'curves': recall_curves(rows, arms, caps),
        'curves_conservative': recall_curves(rows, arms, caps, 'conservative_group'),
        'effects_by_cap': {c: effects(rows, c, arms) for c in caps},
        'effects_conservative_primary_cap': effects(rows, cap, arms, 'conservative_group'),
        'work': work_table(traces, key, arms, traces['esm2'][arms[0]]['reference_count']),
        'failures_primary_cap': failure_breakdown(rows, arms, cap),
        'strata_primary_cap': strata,
    }
    primary = report['effects_by_cap'][cap]
    report['verdicts'] = {k: verdict(primary[k]) for k in ('esm2: gene_caller - minimal',
                                                           'blosum62: gene_caller - minimal',
                                                           'esm2: start_ranked - minimal',
                                                           'blosum62: start_ranked - minimal')
                          if k in primary}
    return report, rows, reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--answer-key', type=Path, required=True)
    parser.add_argument('--construct-map', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report, _, reports = analyze(args.run, args.answer_key, args.construct_map)
    report['input_sha256'] = {
        'answer_key': hashlib.sha256(args.answer_key.read_bytes()).hexdigest(),
        'construct_map': hashlib.sha256(args.construct_map.read_bytes()).hexdigest(),
        'traces': {f'{s}/{a}': hashlib.sha256((args.run / s / f'{a}.json').read_bytes()).hexdigest()
                   for s in SCORERS for a in MAIN_ARMS}}
    with args.out.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True, default=str)
        handle.write('\n')
    with args.out.with_suffix('.verifier.json').open('x', encoding='utf-8') as handle:
        json.dump(reports, handle, indent=2, sort_keys=True, default=str)
        handle.write('\n')


if __name__ == '__main__':
    main()
