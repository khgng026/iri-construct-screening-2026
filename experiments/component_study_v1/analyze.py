"""Evaluate the frozen conditions and endpoints of SPEC.md on the private score matrices.

    python experiments/component_study_v1/analyze.py

Writes source-free results (IDs, counts, rates) for Git. Pure Python.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
from statistics import mean, median

from iri_constructs.component_eval import (attributed_recall, family_recall, merge_max, merge_round_robin,
                                           merge_single, unit_retrieval)

ROOT = Path(__file__).resolve().parents[2]
K = 5
CONDITIONS = ('single_longest', 'all_max', 'all_round_robin', 'single_longest_supported')


def rate(numerator, denominator):
    return {'numerator': numerator, 'denominator': denominator,
            'value': round(numerator / denominator, 4) if denominator else None}


def avg(values):
    values = list(values)
    return round(mean(values), 4) if values else None


def clusters(parents, sequence_of):
    """Connected components of parents that share any eligible unit sequence."""
    owner, group = {}, {}
    def find(x):
        while group[x] != x:
            group[x] = group[group[x]]
            x = group[x]
        return x
    for p in parents:
        group[p['construct_id']] = p['construct_id']
    for p in parents:
        for unit in p['eligible_unit_ids']:
            s = sequence_of[unit]
            if s in owner:
                group[find(p['construct_id'])] = find(owner[s])
            else:
                owner[s] = p['construct_id']
    return len({find(p['construct_id']) for p in parents})


def evaluate(freeze, scores):
    units = {u['unit_id']: u for u in freeze['units']}
    sequence_of = {u['unit_id']: u['sequence_sha256'] for u in freeze['units']}
    reference_family = freeze['reference_family']
    out = {}

    # ---- per-unit retrieval over distinct supported sequences (identical in every condition) ----
    distinct = {}
    for u in freeze['units']:
        if u['supported']:
            distinct.setdefault(u['sequence_sha256'], u)
    rows = []
    for s, u in sorted(distinct.items()):
        r = unit_retrieval(scores[s], set(u['supported_families']), reference_family, K)
        per_family = {f: unit_retrieval(scores[s], {f}, reference_family, K)['hit_at_k'] for f in u['supported_families']}
        rows.append({'sequence_sha256': s, 'label_class': u['label_class'], 'families': u['supported_families'],
                     'length_aa': u['length_aa'], **r, 'families_hit_at_k': per_family})
    def unit_summary(subset):
        return {'distinct_sequences': len(subset), 'hit_at_1': rate(sum(r['hit_at_1'] for r in subset), len(subset)),
                'hit_at_5': rate(sum(r['hit_at_k'] for r in subset), len(subset)),
                'mean_reciprocal_rank': avg(r['reciprocal_rank'] for r in subset)}
    by_family = defaultdict(list)
    for r in rows:
        if r['label_class'] == 'full':
            by_family[r['families'][0]].append(r)
    multi = [r for r in rows if r['label_class'] == 'multi_family']
    out['per_unit'] = {
        'all_supported': unit_summary(rows),
        'by_label_class': {c: unit_summary([r for r in rows if r['label_class'] == c]) for c in ('full', 'embedded', 'multi_family')},
        'full_by_family': {f: {**unit_summary(v), 'first_same_family_ranks': sorted(r['first_same_family_rank'] for r in v)}
                           for f, v in sorted(by_family.items())},
        'multi_family_units': {'units': len(multi), 'family_slots': sum(len(r['families']) for r in multi),
                               'family_slots_hit_at_5': sum(sum(r['families_hit_at_k'].values()) for r in multi),
                               'units_with_every_family_hit_at_5': sum(all(r['families_hit_at_k'].values()) for r in multi)},
        'misses_at_5': [{'sequence_sha256': r['sequence_sha256'][:16], 'label_class': r['label_class'], 'families': r['families'],
                         'length_aa': r['length_aa'], 'first_same_family_rank': r['first_same_family_rank']}
                        for r in rows if not r['hit_at_k']]}

    # ---- parent-level conditions ----
    parent_rows = []
    for p in freeze['parents']:
        if p['analysis_set'] == 'no_supported_component':
            continue
        unit_scores = {u: scores[sequence_of[u]] for u in p['eligible_unit_ids']}
        families = {u: set(units[u]['supported_families']) for u in p['supported_unit_ids']}
        hits = {'single_longest': merge_single(unit_scores, p['baseline_unit_id'], K),
                'all_max': merge_max(unit_scores, K),
                'all_round_robin': merge_round_robin(unit_scores, p['round_robin_order'], K),
                'single_longest_supported': merge_single(unit_scores, p['baseline_supported_unit_id'], K)}
        row = {'construct_id': p['construct_id'], 'analysis_set': p['analysis_set'],
               'eligible_units': len(p['eligible_unit_ids']), 'supported_components': len(families),
               'baseline_is_supported': p['baseline_unit_id'] in families,
               'signature': '|'.join(sorted(sequence_of[u] for u in p['eligible_unit_ids'])),
               'ceiling_components': sum(unit_retrieval(unit_scores[u], f, reference_family, K)['hit_at_k'] for u, f in families.items()),
               'query_residues_all': sum(units[u]['length_aa'] for u in p['eligible_unit_ids']),
               'query_residues_single': units[p['baseline_unit_id']]['length_aa'], 'conditions': {}}
        # Exploratory only (not in SPEC.md): the same rules at other budgets, to see how far the
        # max-merge loss depends on the panel holding five references per family.
        row['budget_sensitivity'] = {}
        for k in (3, 10):
            alternatives = {'single_longest_supported': merge_single(unit_scores, p['baseline_supported_unit_id'], k),
                            'all_max': merge_max(unit_scores, k),
                            'all_round_robin': merge_round_robin(unit_scores, p['round_robin_order'], k)}
            row['budget_sensitivity'][str(k)] = {
                **{name: attributed_recall(h, families, reference_family)['recovered_components'] for name, h in alternatives.items()},
                'ceiling': sum(unit_retrieval(unit_scores[u], f, reference_family, k)['hit_at_k'] for u, f in families.items())}
        for name, hit_list in hits.items():
            a = attributed_recall(hit_list, families, reference_family)
            f = family_recall(hit_list, families, reference_family)
            origins = Counter(h['origin_unit'] for h in hit_list)
            row['conditions'][name] = {
                'recovered_components': a['recovered_components'], 'attributed_recall': a['attributed_recall'],
                'hit_composition': a['hit_composition'], 'family_recall': f['family_recall'],
                'recovered_families': f['recovered_families'], 'expected_families': f['expected_families'],
                'unattributed_family_credit': len(f['unattributed_credit']),
                'distinct_origin_units': len(origins), 'one_unit_owns_all_hits': max(origins.values()) == len(hit_list)}
        parent_rows.append(row)

    def summarize(subset):
        if not subset:
            return None
        summary = {'parents': len(subset), 'supported_components': sum(r['supported_components'] for r in subset),
                   'per_unit_ceiling': rate(sum(r['ceiling_components'] for r in subset), sum(r['supported_components'] for r in subset)),
                   'mean_ceiling_recall': avg(r['ceiling_components'] / r['supported_components'] for r in subset), 'conditions': {}}
        for name in CONDITIONS:
            cells = [r['conditions'][name] for r in subset]
            composition = Counter()
            for c in cells:
                composition.update(c['hit_composition'])
            summary['conditions'][name] = {
                'mean_attributed_recall': avg(c['attributed_recall'] for c in cells),
                'components_recovered': rate(sum(c['recovered_components'] for c in cells), sum(r['supported_components'] for r in subset)),
                'parents_with_every_component_recovered': rate(sum(c['attributed_recall'] == 1 for c in cells), len(cells)),
                'mean_family_recall_origin_agnostic': avg(c['family_recall'] for c in cells),
                'families_recovered_origin_agnostic': rate(sum(c['recovered_families'] for c in cells), sum(c['expected_families'] for c in cells)),
                'family_recoveries_without_attribution': sum(c['unattributed_family_credit'] for c in cells),
                'hit_composition': dict(sorted(composition.items())),
                'parents_where_one_unit_owns_all_hits': sum(c['one_unit_owns_all_hits'] for c in cells),
                'mean_distinct_origin_units': avg(c['distinct_origin_units'] for c in cells)}
        def paired(a, b):
            diffs = [r['conditions'][a]['attributed_recall'] - r['conditions'][b]['attributed_recall'] for r in subset]
            return {'better': sum(d > 0 for d in diffs), 'equal': sum(d == 0 for d in diffs), 'worse': sum(d < 0 for d in diffs),
                    'mean_difference': avg(diffs)}
        summary['paired_attributed_recall'] = {'all_max_minus_single_longest': paired('all_max', 'single_longest'),
                                               'all_round_robin_minus_single_longest': paired('all_round_robin', 'single_longest'),
                                               'all_round_robin_minus_all_max': paired('all_round_robin', 'all_max'),
                                               'all_max_minus_single_longest_supported': paired('all_max', 'single_longest_supported'),
                                               'all_round_robin_minus_single_longest_supported': paired('all_round_robin', 'single_longest_supported')}
        summary['query_cost'] = {
            'single_longest': {'model_calls': len(subset), 'residues': sum(r['query_residues_single'] for r in subset)},
            'all_units': {'model_calls': sum(r['eligible_units'] for r in subset), 'residues': sum(r['query_residues_all'] for r in subset)}}
        summary['query_cost']['all_over_single'] = {
            'model_calls': round(summary['query_cost']['all_units']['model_calls'] / len(subset), 3),
            'residues': round(summary['query_cost']['all_units']['residues'] / summary['query_cost']['single_longest']['residues'], 3)}
        return summary

    def dedupe(subset):
        seen, kept = set(), []
        for r in subset:
            if r['signature'] not in seen:
                seen.add(r['signature'])
                kept.append(r)
        return kept

    multi_rows = [r for r in parent_rows if r['analysis_set'] == 'multi_unit']
    single_rows = [r for r in parent_rows if r['analysis_set'] == 'single_unit']
    by_id = {p['construct_id']: p for p in freeze['parents']}
    out['multi_unit'] = {
        'all_parents': summarize(multi_rows), 'distinct_unit_sets': summarize(dedupe(multi_rows)),
        'dependence_clusters': clusters([by_id[r['construct_id']] for r in multi_rows], sequence_of),
        'strata_distinct_unit_sets': {
            'two_or_more_supported_components': summarize([r for r in dedupe(multi_rows) if r['supported_components'] >= 2]),
            'one_supported_component': summarize([r for r in dedupe(multi_rows) if r['supported_components'] == 1]),
            'baseline_unit_is_supported': summarize([r for r in dedupe(multi_rows) if r['baseline_is_supported']]),
            'baseline_unit_is_unlabeled': summarize([r for r in dedupe(multi_rows) if not r['baseline_is_supported']])}}
    distinct_multi = dedupe(multi_rows)
    out['exploratory_budget_sensitivity'] = {
        'status': 'Exploratory, added after the registered K=5 results were seen; not a registered endpoint.',
        'scope': 'distinct multi-unit unit sets', 'supported_components': sum(r['supported_components'] for r in distinct_multi),
        'components_recovered': {k: {name: sum(r['budget_sensitivity'][k][name] for r in distinct_multi)
                                     for name in ('single_longest_supported', 'all_max', 'all_round_robin', 'ceiling')} for k in ('3', '10')}}
    out['single_unit_control'] = {
        'parents': len(single_rows), 'distinct_unit_sets': len(dedupe(single_rows)),
        'all_conditions_identical': all(len({json.dumps(r['conditions'][c], sort_keys=True) for c in CONDITIONS}) == 1 for r in single_rows),
        'components_recovered': rate(sum(r['conditions']['single_longest']['recovered_components'] for r in single_rows), len(single_rows))}
    out['parent_rows'] = [{k: v for k, v in r.items() if k != 'signature'} for r in parent_rows]
    return out


def coverage(freeze):
    parents, units = freeze['parents'], freeze['units']
    withheld = Counter(u['reason_code'] for u in units if u['extraction_state'] == 'withheld')
    fractions = [p['exported_cds_nt_fraction'] for p in parents]
    return {
        'parents': len(parents), 'source_features': sum(p['source_features'] for p in parents),
        'noncoding_features': sum(p['noncoding_features'] for p in parents), 'source_cds': len(units),
        'exported_cds': rate(sum(u['extraction_state'] == 'exported' for u in units), len(units)),
        'withheld_by_reason': dict(sorted(withheld.items())),
        'exported_but_over_model_limit': sum(u['extraction_state'] == 'exported' and not u['eligible'] for u in units),
        'eligible_units': rate(sum(u['eligible'] for u in units), len(units)),
        'label_class_of_exported': freeze['counts']['label_class'],
        'supported_units': rate(sum(u['supported'] for u in units), len(units)),
        'supported_of_eligible': rate(sum(u['supported'] for u in units), sum(u['eligible'] for u in units)),
        'labeled_but_family_not_scorable': sum(bool(u['families']) and u['eligible'] and not u['supported'] for u in units),
        'parents_without_any_cds': sum(p['source_cds'] == 0 for p in parents),
        'parents_without_exported_cds': sum(p['exported_cds'] == 0 for p in parents),
        'parents_with_withheld_cds': sum(bool(p['withheld_by_reason']) for p in parents),
        'analysis_sets': freeze['counts']['analysis_sets'],
        'construct_length_inside_exported_cds': {'mean_fraction': avg(fractions), 'median_fraction': round(median(fractions), 4),
                                                 'min': min(fractions), 'max': max(fractions)},
        'note': 'Everything outside exported CDS segments (promoters, origins, terminators, withheld CDSs, unannotated sequence) '
                'is not seen by the protein model in any condition.'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--scores', type=Path, default=ROOT / 'results' / 'component_study_v1')
    parser.add_argument('--record', type=Path, default=ROOT / 'docs' / 'results' / 'component-study-v1')
    args = parser.parse_args(argv)
    freeze = json.loads((args.record / 'panel.freeze.json').read_text())
    results = {'schema_version': 'component-study-results-1', 'panel_sha256': freeze['panel_sha256'],
               'spec_sha256': freeze['spec_sha256'], 'budget_k': K, 'coverage': coverage(freeze), 'methods': {}}
    for method in ('esm2', 'alignment'):
        private = json.loads((args.scores / f'scores.{method}.json').read_text())
        if private['panel_sha256'] != freeze['panel_sha256']:
            raise SystemExit('Scores do not belong to the committed freeze')
        record = json.loads((args.record / f'scoring.{method}.json').read_text())
        results['methods'][method] = {'score_matrix_sha256': record['score_matrix_sha256'], **evaluate(freeze, private['scores'])}
    (args.record / 'results.json').write_text(json.dumps(results, indent=1, sort_keys=True) + '\n')
    for method, r in results['methods'].items():
        print('==', method)
        print(' per-unit', json.dumps(r['per_unit']['all_supported']))
        for scope in ('all_parents', 'distinct_unit_sets'):
            s = r['multi_unit'][scope]
            print(f"  {scope}: parents {s['parents']} components {s['supported_components']} ceiling {s['per_unit_ceiling']['value']}")
            for name, c in s['conditions'].items():
                print(f"    {name:26s} attr {c['mean_attributed_recall']:.3f} comps {c['components_recovered']['numerator']}/{c['components_recovered']['denominator']}"
                      f" fam(agn) {c['mean_family_recall_origin_agnostic']:.3f} unattributed {c['family_recoveries_without_attribution']} comp {c['hit_composition']}")
        print('  clusters', r['multi_unit']['dependence_clusters'], 'control', json.dumps(r['single_unit_control']))


if __name__ == '__main__':
    main()
