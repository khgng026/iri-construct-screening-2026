"""Consolidated exploratory comparison of all scorers on the locked cohort v2 holdout.

Reads the verifier traces for BLOSUM62 (locked), mean-pooled and residue-level ESM-2 at
8M and 150M (all on identical selected inputs), builds the per-component outcome table once,
and emits a single reproducible JSON: recall per scorer/arm/cap, recall by similarity tier
and archetype, the oracle union of the ESM configurations with alignment, and the
complementarity counts. No inference; the traces are the inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from . import cohort_analysis as ca

SCORERS = ['blosum62', 'esm2_8m_mean', 'esm2_8m_residue', 'esm2_150m_mean', 'esm2_150m_residue']
ESM = [s for s in SCORERS if s.startswith('esm2')]
ARMS = ('minimal', 'start_ranked', 'gene_caller', 'annotation_reference')


def _equal_group(rows, pred):
    groups = defaultdict(list)
    for r in rows:
        groups[r['lineage_group']].append(bool(pred(r)))
    return (sum(pred(r) for r in rows), len(rows),
            float(np.mean([np.mean(v) for v in groups.values()])) if rows else None)


def compare(run_dir: Path, key_path: Path, map_path: Path):
    ca.use_scorers(SCORERS)
    protocols = {s: ca.read(run_dir / s / 'protocol.json') for s in SCORERS}
    traces = {s: {a: ca.read(run_dir / s / f'{a}.json') for a in ARMS} for s in SCORERS}
    # All scorers must share the selected inputs for every arm (they differ only in scoring).
    for arm in ARMS:
        fingerprints = {traces[s][arm]['selected_inputs_sha256'] for s in SCORERS}
        pools = {json.dumps(traces[s][arm]['candidate_pools'], sort_keys=True) for s in SCORERS}
        ca._require(len(fingerprints) == 1 and len(pools) == 1, f'{arm}: scorers used different selected inputs')
    key, construct_map = ca.read(key_path), ca.read(map_path)
    rows, _ = ca.component_table(protocols, key, traces, construct_map, ARMS)
    caps = protocols['blosum62']['query_caps']

    def recovered(r, s, arm, cap):
        return r['outcome'][(s, arm, cap)]['recovered']

    report = {'scorers': SCORERS, 'components': len(rows),
              'lineage_groups': len({r['lineage_group'] for r in rows}),
              'recall': {}, 'by_tier': {}, 'by_archetype': {}, 'complementarity': {}, 'oracle': {}}
    for s in SCORERS:
        for arm in ARMS:
            for cap in caps:
                total, n, eq = _equal_group(rows, lambda r: recovered(r, s, arm, cap))
                report['recall'][f'{s}|{arm}|{cap}'] = {'recovered': total, 'eligible': n,
                                                        'pooled': total / n if n else None, 'equal_group': eq}
    cap = protocols['blosum62']['primary_query_cap']
    for field, bucket in (('similarity_tier', 'by_tier'), ('archetype', 'by_archetype')):
        for value in sorted({r[field] for r in rows}):
            sub = [r for r in rows if r[field] == value]
            report[bucket][value] = {'components': len(sub), 'lineage_groups': len({r['lineage_group'] for r in sub}),
                                     'gene_caller_recovered': {s: sum(recovered(r, s, 'gene_caller', cap) for r in sub)
                                                               for s in SCORERS}}
    # Oracle union of each ESM config with BLOSUM62, and complementarity, at gene_caller / primary cap.
    for s in ESM:
        only_esm = sum(recovered(r, s, 'gene_caller', cap) and not recovered(r, 'blosum62', 'gene_caller', cap)
                       for r in rows)
        only_bl = sum(recovered(r, 'blosum62', 'gene_caller', cap) and not recovered(r, s, 'gene_caller', cap)
                      for r in rows)
        union_total, _, union_eq = _equal_group(
            rows, lambda r: recovered(r, s, 'gene_caller', cap) or recovered(r, 'blosum62', 'gene_caller', cap))
        report['complementarity'][s] = {'only_this_esm': only_esm, 'only_blosum62': only_bl}
        report['oracle'][f'{s} OR blosum62'] = {'recovered': union_total, 'equal_group': union_eq}
    # Union of every ESM config with BLOSUM62.
    all_total, _, all_eq = _equal_group(rows, lambda r: recovered(r, 'blosum62', 'gene_caller', cap)
                                        or any(recovered(r, s, 'gene_caller', cap) for s in ESM))
    report['oracle']['any_esm OR blosum62'] = {'recovered': all_total, 'equal_group': all_eq}
    report['blosum62_gene_caller_recovered'] = sum(recovered(r, 'blosum62', 'gene_caller', cap) for r in rows)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--answer-key', type=Path, required=True)
    parser.add_argument('--construct-map', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.run, args.answer_key, args.construct_map)
    report['input_sha256'] = {'answer_key': hashlib.sha256(args.answer_key.read_bytes()).hexdigest(),
                              'construct_map': hashlib.sha256(args.construct_map.read_bytes()).hexdigest()}
    with args.out.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True, default=str)
        handle.write('\n')
    print(f"{report['components']} components, {report['lineage_groups']} groups")
    print('gene_caller recall (equal-group) at primary cap:')
    cap = 2
    for s in SCORERS:
        print(f"  {s:20s} {report['recall'][f'{s}|gene_caller|{cap}']['recovered']}")
    print('insufficient_coverage tier:', report['by_tier'].get('insufficient_coverage', {}).get('gene_caller_recovered'))
    print('complementarity vs blosum62:', report['complementarity'])


if __name__ == '__main__':
    main()
