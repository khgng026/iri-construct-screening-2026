"""Fusion subregion analysis: intact candidates, one windowing arm and documented regions.

Endpoints per documented subregion, for each scorer, arm and query cap:

* ``region_supported``: the verifier's endpoint on the subregion key. A correct family
  call comes from a query on the same strand covering at least the protocol fraction of
  the documented subregion's bases.
* ``family_recovered``: a correct family call from any query that overlaps the subregion.
* For misses, whether no query reached the region, a query reached it but no single query
  covered it (the split cut across the region), or a covering query reached the scorer
  without a correct call.

Also reports query counts, input residues and how much of each supporting input lies
outside the documented region.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from .verify import SCHEMA, _covered, _location, verify

SCORERS = ('esm2', 'blosum62')
ARMS = ('gene_caller', 'gene_caller_windows', 'annotation_reference', 'documented_regions')


def read(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def anchor_region_key(freeze_path: Path, construct_map_path: Path, out: Path):
    """Integration-only key from anchor regions of development fusion components (not evidence)."""
    freeze = read(freeze_path)
    ids = {row['accession']: row['id'] for row in read(construct_map_path)}
    constructs = []
    for record in freeze['records']:
        if record['accession'] not in ids:
            continue
        components = [{'id': f'{c["id"]}-anchor', 'family': c['family'], 'strand': c['strand'],
                       'spans': c['anchor_region']['spans']}
                      for c in record['components'] if c['fusion_class'] and c['anchor_region']['usable']]
        constructs.append({'id': ids[record['accession']], 'archetype': record['archetype'],
                           'length': record['length'], 'components': components})
    with out.open('x', encoding='utf-8') as handle:
        json.dump({'schema': SCHEMA, 'constructs': sorted(constructs, key=lambda c: c['id']),
                   'note': 'anchor-region boundaries; integration check only'}, handle, indent=2, sort_keys=True)


def _score(protocol, key, trace):
    """Score one trace with the verifier; input identity is checked separately by the caller."""
    fixed = {k: v for k, v in protocol.items() if k != 'sequence_input_sha256'}
    return verify(fixed, key, trace, [dict(trace, method=f'{trace["method"]}#copy')])['runs'][trace['method']]


def analyze(run_dir: Path, key_path: Path, arms=ARMS, scorers=SCORERS):
    key = read(key_path)
    truth = {c['id']: c for c in key['constructs']}
    out = {'regions': sum(len(c['components']) for c in key['constructs']), 'by_arm': {}}
    for scorer in scorers:
        protocol = read(run_dir / scorer / 'protocol.json')
        overlap = protocol['minimum_component_overlap']
        traces = {a: read(run_dir / scorer / f'{a}.json') for a in arms if (run_dir / scorer / f'{a}.json').exists()}
        for arm, trace in traces.items():
            if arm.startswith('gene_caller') and trace['input_sha256'] != protocol['sequence_input_sha256']:
                raise ValueError(f'{arm}: sequence input differs from the protocol')
        for arm, trace in traces.items():
            scored = _score(protocol, key, trace)
            for tier in trace['tiers']:
                cap = tier['query_cap']
                rows = {r['id']: r for r in tier['constructs']}
                recovered = {(r['construct_id'], c) for r in scored[cap]['per_construct']
                             for c in r['recovered_component_ids']}
                statuses, outside, residues, queries = Counter(), [], 0, 0
                family_hits = 0
                for cid, construct in truth.items():
                    row = rows[cid]
                    residues += sum(q['input_residues'] for q in row['queries'])
                    queries += len(row['queries'])
                    for component in construct['components']:
                        locus = _location(component, construct['length'], component['id'])
                        region_bases = sum(b - a for a, b in locus[1])
                        loci = {q['id']: (q['strand'], sorted(map(tuple, q['spans']))) for q in row['queries']}
                        fractions = {qid: _covered(locus, ql) for qid, ql in loci.items()}
                        touching = {qid for qid, f in fractions.items() if f > 0}
                        covering = {qid for qid, f in fractions.items() if f >= overlap}
                        callers = {p['query_id'] for p in row['predictions'] if p['family'] == component['family']}
                        family = bool(callers & touching)
                        family_hits += family
                        if (cid, component['id']) in recovered:
                            status = 'region_supported'
                            for qid in callers & covering:
                                query_bases = sum(b - a for a, b in loci[qid][1])
                                outside.append(1 - region_bases * fractions[qid] / query_bases)
                        elif family:
                            status = 'family_from_partial_region'
                        elif covering:
                            status = 'scorer_missed_with_covering_input'
                        elif touching:
                            status = 'split_cut_across_region'
                        else:
                            status = 'region_not_reached'
                        statuses[status] += 1
                out['by_arm'][f'{scorer}|{arm}|{cap}'] = {
                    'region_supported': statuses['region_supported'], 'family_recovered': family_hits,
                    'statuses': dict(statuses), 'queries': queries, 'input_residues': residues,
                    'mean_fraction_of_supporting_input_outside_region':
                        sum(outside) / len(outside) if outside else None}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    key = sub.add_parser('anchor-key')
    key.add_argument('--freeze', type=Path, required=True)
    key.add_argument('--construct-map', type=Path, required=True)
    key.add_argument('--out', type=Path, required=True)
    run = sub.add_parser('analyze')
    run.add_argument('--run', type=Path, required=True)
    run.add_argument('--key', type=Path, required=True)
    run.add_argument('--out', type=Path, required=True)
    run.add_argument('--scorers', nargs='+', default=list(SCORERS))
    args = parser.parse_args()
    if args.stage == 'anchor-key':
        anchor_region_key(args.freeze, args.construct_map, args.out)
    else:
        result = analyze(args.run, args.key, scorers=args.scorers)
        with args.out.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write('\n')


if __name__ == '__main__':
    main()
