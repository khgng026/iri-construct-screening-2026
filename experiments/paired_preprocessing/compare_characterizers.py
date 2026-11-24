"""Compare two scorers on identical reporter inputs and an identical answer key."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from .reporter_pilot import read, write
from .verify import _require, verify

BASELINE = 'minimal'
CANDIDATES = ('start_aware', 'start_ranked', 'gene_caller')
REFERENCE = 'annotation_reference'
METHODS = (BASELINE, *CANDIDATES, REFERENCE)


def _recovered(tier):
    return {(row['construct_id'], component) for row in tier['per_construct']
            for component in row['recovered_component_ids']}


def compare(protocols, runs, answer_key):
    """Verify input equivalence and report paired recall at each cap."""
    names = ('esm2', 'blosum62')
    common = [{k: v for k, v in protocols[name].items() if k not in ('characterizer', 'characterization_id')}
              for name in names]
    _require(common[0] == common[1], 'scorers must share the evaluation protocol and reference panel')
    for name in names:
        _require(protocols[name].get('characterizer', 'esm2') == name, 'protocol uses the wrong scorer')
        _require(runs[name].keys() == set(METHODS), 'expected all reporter comparison arms')

    matched_inputs = {}
    for method in METHODS:
        left, right = (runs[name][method] for name in names)
        for name, run in zip(names, (left, right), strict=True):
            _require(run['method'] == method and run['characterizer'] == name, 'trace uses the wrong arm or scorer')
            _require(run['reference_panel_sha256'] == common[0]['reference_panel_sha256'], 'trace uses the wrong reference panel')
        digest = left['selected_inputs_sha256']
        _require(isinstance(digest, str) and re.fullmatch(r'[0-9a-f]{64}', digest) is not None
                 and digest == right['selected_inputs_sha256'] and left['input_sha256'] == right['input_sha256'],
                 f'{method}: selected inputs differ between scorers')
        matched_inputs[method] = digest

    verified = {name: verify(protocols[name], answer_key, runs[name][BASELINE],
                             [runs[name][method] for method in CANDIDATES], runs[name][REFERENCE]) for name in names}
    matched = {}
    for cap in common[0]['query_caps']:
        at_budget = {}
        for method in METHODS:
            esm, alignment = (verified[name]['runs'][method][cap] for name in names)
            esm_hits, alignment_hits = _recovered(esm), _recovered(alignment)
            esm_total, alignment_total = esm['overall'], alignment['overall']
            for field in ('queries', 'input_residues', 'eligible_components', 'components_without_covering_query'):
                _require(esm_total[field] == alignment_total[field], f'{method}: incompatible {field} counts')
            at_budget[method] = {
                'eligible_components': esm_total['eligible_components'],
                'queries_per_scorer': esm_total['queries'],
                'input_residues_per_scorer': esm_total['input_residues'],
                'recovered': {'esm2': len(esm_hits), 'blosum62': len(alignment_hits)},
                'scorer_only_recoveries': {'esm2': len(esm_hits - alignment_hits),
                                           'blosum62': len(alignment_hits - esm_hits)},
            }
            if method in CANDIDATES:
                gains = {name: verified[name]['candidate_minus_baseline'][method][cap]['recall_gain'] for name in names}
                at_budget[method]['recall_gain_over_minimal'] = gains
                at_budget[method]['difference_in_gain'] = gains['esm2'] - gains['blosum62']
        matched[cap] = at_budget
    return {'primary_query_cap': common[0]['primary_query_cap'], 'matched_inputs_sha256': matched_inputs,
            'by_query_cap': matched, 'scorer_reports': verified}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--esm2-run', type=Path, required=True)
    parser.add_argument('--blosum62-run', type=Path, required=True)
    parser.add_argument('--answer-key', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    paths = {'esm2': args.esm2_run, 'blosum62': args.blosum62_run}
    try:
        protocols = {name: read(path / 'protocol.json') for name, path in paths.items()}
        runs = {name: {method: read(path / f'{method}.json') for method in METHODS} for name, path in paths.items()}
        report = compare(protocols, runs, read(args.answer_key))
        write(args.out, report)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.exit(2, f'Cannot compare characterizers: {exc}\n')


if __name__ == '__main__':
    main()
