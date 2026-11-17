"""Registered analyses 2 and 3 of SPEC_EXTENSION.md (pure Python; no model).

    python experiments/component_study_v1/analyze_extension.py
"""
from __future__ import annotations

import json
from pathlib import Path

from analyze import CONDITIONS, K, evaluate
from iri_constructs.component_eval import unit_retrieval
from iri_constructs.pilot_scoring import paired_cluster_bootstrap

ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / 'docs' / 'results' / 'component-study-v1'
EXT = V1 / 'extension-1'
PRIVATE = {'v1': ROOT / 'results' / 'component_study_v1', 'ext': ROOT / 'results' / 'component_study_v1_ext1'}
SEED, REPLICATES = 20261118, 2000


def signatures(freeze):
    sequence_of = {u['unit_id']: u['sequence_sha256'] for u in freeze['units']}
    return {p['construct_id']: '|'.join(sorted(sequence_of[u] for u in p['eligible_unit_ids'])) for p in freeze['parents']}, sequence_of


def cluster_ids(parents, sequence_of):
    group, owner = {p['construct_id']: p['construct_id'] for p in parents}, {}
    def find(x):
        while group[x] != x:
            group[x] = group[group[x]]
            x = group[x]
        return x
    for p in parents:
        for unit in p['eligible_unit_ids']:
            s = sequence_of[unit]
            if s in owner:
                group[find(p['construct_id'])] = find(owner[s])
            else:
                owner[s] = p['construct_id']
    return {p['construct_id']: 'cluster-' + find(p['construct_id']) for p in parents}


def per_unit(freeze, scores, exclude):
    family = freeze['reference_family']
    distinct = {}
    for u in freeze['units']:
        if u['supported'] and u['sequence_sha256'] not in exclude:
            distinct.setdefault(u['sequence_sha256'], u)
    out = {}
    for label_class in ('full', 'embedded', 'multi_family'):
        rows = [(u, unit_retrieval(scores[s], set(u['supported_families']), family, K)) for s, u in distinct.items() if u['label_class'] == label_class]
        slots = [unit_retrieval(scores[u['sequence_sha256']], {f}, family, K)['hit_at_k'] for u, _ in rows for f in u['supported_families']]
        out[label_class] = {'distinct_sequences': len(rows), 'hit_at_1': sum(r['hit_at_1'] for _, r in rows),
                            'hit_at_5': sum(r['hit_at_k'] for _, r in rows), 'family_slots': len(slots), 'family_slots_hit_at_5': sum(slots)}
    return out


def main():
    v1, ext = json.loads((V1 / 'panel.freeze.json').read_text()), json.loads((EXT / 'panel.freeze.json').read_text())
    v1_signatures, _ = signatures(v1)
    ext_signatures, sequence_of = signatures(ext)
    seen_sets = set(v1_signatures.values())
    v1_sequences = {u['sequence_sha256'] for u in v1['units'] if u['eligible']}
    repeated = [p['construct_id'] for p in ext['parents'] if p['analysis_set'] == 'multi_unit' and ext_signatures[p['construct_id']] in seen_sets]
    new_only = {**ext, 'parents': [p for p in ext['parents'] if ext_signatures[p['construct_id']] not in seen_sets]}
    results = {'schema_version': 'component-study-extension-1', 'spec': 'experiments/component_study_v1/SPEC_EXTENSION.md',
               'panel_sha256': ext['panel_sha256'], 'reference_panel_reused_from': ext['reference_panel_reused_from'],
               'multi_unit_parents_in_extension': sum(p['analysis_set'] == 'multi_unit' for p in ext['parents']),
               'multi_unit_parents_removed_because_unit_set_seen_in_v1': len(repeated), 'methods': {}}
    for method in ('esm2', 'alignment'):
        scores = json.loads((PRIVATE['ext'] / f'scores.{method}.json').read_text())['scores']
        v1_scores = json.loads((PRIVATE['v1'] / f'scores.{method}.json').read_text())['scores']
        evaluated = evaluate(new_only, scores)
        rows, seen = [], set()
        for r in evaluated['parent_rows']:                      # one row per distinct unit set, multi-unit only
            signature = ext_signatures[r['construct_id']]
            if r['analysis_set'] == 'multi_unit' and signature not in seen:
                seen.add(signature)
                rows.append(r)
        by_id = {p['construct_id']: p for p in new_only['parents']}
        clusters = cluster_ids([by_id[r['construct_id']] for r in rows], sequence_of)
        values = {c: {r['construct_id']: r['conditions'][c]['attributed_recall'] for r in rows} for c in CONDITIONS}
        intervals = {}
        for left, right in (('single_longest_supported', 'all_round_robin'), ('single_longest_supported', 'all_max'),
                            ('all_max', 'all_round_robin'), ('single_longest', 'all_round_robin')):
            b = paired_cluster_bootstrap(values[left], values[right], clusters, seed=SEED, replicates=REPLICATES, min_clusters=8)
            intervals[f'{right}_minus_{left}'] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in b.items()
                                                  if k in ('status', 'independent_clusters', 'paired_cluster_mean_difference', 'ci_lower', 'ci_upper', 'unit')}
        sizes = sorted((sum(1 for c in clusters.values() if c == g) for g in set(clusters.values())), reverse=True)
        new_units, v1_units = per_unit(ext, scores, v1_sequences), per_unit(v1, v1_scores, set())
        fusion = {}
        for label_class in ('embedded', 'multi_family'):
            fusion[label_class] = {part: {k: source[label_class][k] for k in ('distinct_sequences', 'hit_at_5', 'family_slots', 'family_slots_hit_at_5')}
                                   for part, source in (('v1', v1_units), ('extension_new', new_units))}
        slots = sum(fusion[c][p]['family_slots'] for c in fusion for p in fusion[c])
        results['methods'][method] = {
            'replication_on_unseen_unit_sets': {**evaluated['multi_unit']['distinct_unit_sets'],
                                                'dependence_clusters': len(set(clusters.values())), 'largest_cluster_sizes': sizes[:5]},
            'paired_cluster_bootstrap': {'seed': SEED, 'replicates': REPLICATES, 'intervals': intervals,
                'note': 'Each dependence cluster has equal weight, however many parents it holds; 95% percentile interval. '
                        'Population: records matching the fixed NCBI query. One cluster holds most parents, so the '
                        'equal-weight cluster mean differs from the per-parent mean.'},
            'per_unit_new_sequences': new_units,
            'fusion_classes_v1_plus_extension': {**fusion, 'family_slots_total': slots,
                'family_slots_hit_at_5_total': sum(fusion[c][p]['family_slots_hit_at_5'] for c in fusion for p in fusion[c])},
            'single_unit_control_all_conditions_identical': evaluated['single_unit_control']['all_conditions_identical']}
    (EXT / 'results.extension.json').write_text(json.dumps(results, indent=1, sort_keys=True) + '\n')
    for method, r in results['methods'].items():
        rep = r['replication_on_unseen_unit_sets']
        print('==', method, 'parents', rep['parents'], 'components', rep['supported_components'], 'ceiling', rep['per_unit_ceiling'], 'clusters', rep['dependence_clusters'], rep['largest_cluster_sizes'])
        print('  ', {c: rep['conditions'][c]['components_recovered']['numerator'] for c in CONDITIONS}, {c: rep['conditions'][c]['mean_attributed_recall'] for c in CONDITIONS})
        print('   boot', json.dumps(r['paired_cluster_bootstrap']['intervals']))
        print('   new per-unit', json.dumps(r['per_unit_new_sequences']))
        print('   fusion', json.dumps(r['fusion_classes_v1_plus_extension']))
    print('removed (unit set seen in v1):', len(repeated))


if __name__ == '__main__':
    main()
