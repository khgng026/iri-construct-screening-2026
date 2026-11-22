"""Score paired construct runs against held-out component annotations."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

SCHEMA = 'paired-preprocessing-1'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _int(value, minimum):
    return type(value) is int and value >= minimum


def _number(value, low, high=math.inf):
    return type(value) in (float, int) and math.isfinite(value) and low <= value <= high


def _index(rows, key, context):
    result = {}
    for row in rows:
        value = row[key]
        _require(type(value) is str and value.strip(), f'{context}: missing {key}')
        _require(value not in result, f'{context}: duplicate {key} {value!r}')
        result[value] = row
    return result


def _location(row, length, context):
    strand = row['strand']
    _require(type(strand) is int and strand in (-1, 1), f'{context}: expected strand -1 or 1')
    spans = row['spans']
    _require(isinstance(spans, list) and spans, f'{context}: missing spans')
    for start, end in spans:
        _require(_int(start, 0) and _int(end, 1) and start < end <= length,
                 f'{context}: span outside construct or empty')
    spans = sorted(map(tuple, spans))
    _require(all(a[1] <= b[0] for a, b in zip(spans, spans[1:])), f'{context}: overlapping spans')
    return strand, spans


def _covered(component, query):
    component_strand, spans = component
    query_strand, query_spans = query
    if component_strand != query_strand:
        return 0
    shared = sum(max(0, min(b, d) - max(a, c)) for a, b in spans for c, d in query_spans)
    return shared / sum(end - start for start, end in spans)


def _match(edges):
    """Assign each prediction and documented component at most once."""
    owner = {}

    def assign(call, seen):
        for component in edges[call]:
            if component in seen:
                continue
            seen.add(component)
            previous = owner.get(component)
            if previous is None or assign(previous, seen):
                owner[component] = call
                return True
        return False

    for call in range(len(edges)):
        assign(call, set())
    return owner


def _score_construct(truth, trace, cap, call_limit, overlap):
    construct_id = truth['id']
    length = truth['length']
    _require(_int(length, 1), f'{construct_id}: invalid construct length')
    components = _index(truth['components'], 'id', construct_id)
    loci = {cid: _location(row, length, cid) for cid, row in components.items()}
    for cid, component in components.items():
        _require(isinstance(component['family'], str) and component['family'].strip(),
                 f'{cid}: missing family label')

    queries = _index(trace['queries'], 'id', construct_id)
    _require(len(queries) <= cap, f'{construct_id}: query count exceeds cap {cap}')
    _require(_int(trace['candidate_count'], len(queries)), f'{construct_id}: invalid candidate count')
    query_loci = {qid: _location(row, length, qid) for qid, row in queries.items()}
    covered_by_query = {
        qid: {cid for cid, locus in loci.items() if _covered(locus, qlocus) >= overlap}
        for qid, qlocus in query_loci.items()
    }
    covered = set().union(*covered_by_query.values())
    _require(all(_int(row['input_residues'], 1) for row in queries.values()),
             f'{construct_id}: invalid model input length')
    for field in ('preprocess_seconds', 'characterize_seconds'):
        _require(_number(trace[field], 0), f'{construct_id}: invalid {field}')

    calls = trace['predictions']
    _require(len(calls) <= call_limit, f'{construct_id}: prediction count exceeds limit {call_limit}')
    _require(trace['status'] in ('ok', 'failed'), f'{construct_id}: invalid status')
    _require(trace['status'] != 'failed' or not calls, f'{construct_id}: failed trace has predictions')
    edges, seen = [], set()
    for call in calls:
        qid, label = call['query_id'], call['family']
        _require(qid in queries, f'{construct_id}: prediction from unqueried input {qid!r}')
        _require(isinstance(label, str) and label.strip(), f'{construct_id}: empty predicted family')
        _require((qid, label) not in seen, f'{construct_id}: duplicate prediction for {qid!r}')
        seen.add((qid, label))
        edges.append(sorted(cid for cid in covered_by_query[qid] if components[cid]['family'] == label))
    recovered = _match(edges)
    return {
        'construct_id': construct_id,
        'archetype': truth['archetype'],
        'status': trace['status'],
        'eligible_components': len(components),
        'recovered_components': len(recovered),
        'recovered_component_ids': sorted(recovered),
        'covered_component_ids': sorted(covered),
        'components_without_covering_query': len(components) - len(covered),
        'covered_components_not_recovered': len(covered) - len(recovered),
        'candidate_count': trace['candidate_count'],
        'queries': len(queries),
        'input_residues': sum(query['input_residues'] for query in queries.values()),
        'returned_predictions': len(calls),
        'preprocess_seconds': trace['preprocess_seconds'],
        'characterize_seconds': trace['characterize_seconds'],
    }


def _summarize(rows):
    out = {field: sum(row[field] for row in rows) for field in (
        'eligible_components', 'recovered_components', 'candidate_count', 'queries', 'input_residues',
        'returned_predictions', 'preprocess_seconds', 'characterize_seconds',
        'components_without_covering_query', 'covered_components_not_recovered')}
    out['constructs'] = len(rows)
    out['failed_constructs'] = sum(row['status'] == 'failed' for row in rows)
    n = out['eligible_components']
    out['recall'] = out['recovered_components'] / n if n else None
    return out


def _score_run(run, truth, caps, call_limit, overlap):
    by_cap = {}
    for tier in run['tiers']:
        cap = tier['query_cap']
        _require(_int(cap, 1) and cap in caps and cap not in by_cap, f'{run["method"]}: invalid or repeated cap')
        traces = _index(tier['constructs'], 'id', run['method'])
        _require(traces.keys() == truth.keys(), f'{run["method"]} at cap {cap}: construct set differs from key')
        rows = [_score_construct(truth[cid], traces[cid], cap, call_limit, overlap) for cid in sorted(truth)]
        groups = defaultdict(list)
        for row in rows:
            groups[row['archetype']].append(row)
        by_cap[cap] = {
            'overall': _summarize(rows),
            'by_archetype': {key: _summarize(group) for key, group in sorted(groups.items())},
            'per_construct': rows,
        }
    _require(by_cap.keys() == set(caps), f'{run["method"]}: missing budget tiers')
    return by_cap


def verify(protocol, answer_key, baseline, candidates, annotation_reference=None):
    """Return recall curves and within-cap deltas for the same documented components."""
    runs = [baseline, *candidates]
    if annotation_reference is not None:
        runs.append(annotation_reference)
    for document in (protocol, answer_key, *runs):
        _require(document['schema'] == SCHEMA, f'expected schema {SCHEMA}')
    _require(candidates, 'at least one preprocessing run is required')
    methods = _index(runs, 'method', 'runs')
    fixed_workflow = protocol['characterization_id']
    _require(isinstance(fixed_workflow, str) and fixed_workflow.strip(), 'missing characterization ID')
    _require(all(run['characterization_id'] == fixed_workflow for run in runs),
             'all runs must use the protocol characterization workflow')
    if 'sequence_input_sha256' in protocol:
        _require(all(run.get('input_sha256') == protocol['sequence_input_sha256'] for run in [baseline, *candidates]),
                 'sequence-only runs must use the same frozen input bytes')
    caps = protocol['query_caps']
    _require(caps and all(_int(cap, 1) for cap in caps) and len(set(caps)) == len(caps), 'invalid query caps')
    cap = protocol['primary_query_cap']
    _require(_int(cap, 1) and cap in caps, 'primary cap must be one of the query caps')
    call_limit, overlap = protocol['prediction_limit'], protocol['minimum_component_overlap']
    _require(_int(call_limit, 1), 'invalid prediction limit')
    _require(_number(overlap, 0, 1) and overlap > 0, 'minimum component overlap must be in (0, 1]')
    target = protocol['minimum_interesting_gain']
    _require(_number(target, 0, 1), 'minimum interesting gain must be a fraction in [0, 1]')
    truth = _index(answer_key['constructs'], 'id', 'answer key')
    _require(truth and any(row['components'] for row in truth.values()), 'no eligible documented components')
    _require(all(isinstance(row['archetype'], str) and row['archetype'].strip() for row in truth.values()),
             'each construct needs an archetype')

    scored = {method: _score_run(run, truth, caps, call_limit, overlap) for method, run in methods.items()}
    baseline_id = baseline['method']
    deltas = {}
    for run in candidates:
        name = run['method']
        deltas[name] = {}
        for tier in sorted(caps):
            base = scored[baseline_id][tier]['overall']
            observed = scored[name][tier]['overall']
            gain = (observed['recovered_components'] - base['recovered_components']) / base['eligible_components']
            deltas[name][tier] = {'recall_gain': gain, 'primary_budget': tier == cap}
            if tier == cap:
                deltas[name][tier]['meets_effect_target_by_point_estimate'] = gain >= target
    return {
        'schema': SCHEMA,
        'protocol': protocol,
        'baseline': baseline_id,
        'annotation_reference': annotation_reference['method'] if annotation_reference is not None else None,
        'runs': scored,
        'candidate_minus_baseline': deltas,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--answer-key', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, action='append', required=True)
    parser.add_argument('--annotation-reference', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    digests = {}

    def read(role, path):
        data = path.read_bytes()
        digests[role] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    try:
        report = verify(
            read('protocol', args.protocol), read('answer_key', args.answer_key), read('baseline', args.baseline),
            [read(f'candidate_{i}', path) for i, path in enumerate(args.candidate, 1)],
            read('annotation_reference', args.annotation_reference) if args.annotation_reference else None)
        report['input_sha256'] = digests
        with args.out.open('x', encoding='utf-8') as output:
            output.write(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + '\n')
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.exit(2, f'Cannot score runs: {error}\n')


if __name__ == '__main__':
    main()
