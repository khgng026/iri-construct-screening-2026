"""Verified benign challenge loading and explicit retrieval denominators.

No network, model inference, input alteration, or writing occurs in this module.
Sequence-bearing values returned by ``load_suite`` belong only in ignored storage.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path

from .t05_composition import equivalence_map, load_reviewed, verify_freeze

KS = (1, 5, 10)


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _check(condition, message):
    if not condition:
        raise ValueError(message)


def _read(path):
    return json.loads(Path(path).read_text())


def _verified_path(root, relative, expected):
    path = (root / relative).resolve()
    _check(path.is_relative_to(root.resolve()), 'Frozen path escapes its root')
    _check(_digest(path.read_bytes()) == expected, f'Frozen hash mismatch: {relative}')


def _stable(seed, values, salt):
    return _digest((seed + '|' + salt + '|' + '|'.join(values)).encode())


def _selection_bin(value):
    return 'lt_0.10' if value < .1 else ('0.10_to_lt_0.20' if value < .2 else 'ge_0.20')


def fraction_bin(value):
    _check(math.isfinite(value) and 0 < value <= 1, 'Invalid component fraction')
    for upper, name in ((.25, 'lt_0.25'), (.4, '0.25_to_lt_0.40'),
                        (.6, '0.40_to_lt_0.60'), (.75, '0.60_to_lt_0.75')):
        if value < upper:
            return name
    return 'ge_0.75'


def _validate_groups(manifest, lengths):
    """Independently enumerate every combination, selection, and allowed order."""
    seed, cap, targets = (manifest['selection_seed'], manifest['shared_max_residues'],
                          manifest['targets'])
    expected_orders = {}
    group_counts = {}
    for size, budget in ((3, 90), (4, 45)):
        design = manifest['group_selection'][str(size)]
        bins = {name: [] for name in ('lt_0.10', '0.10_to_lt_0.20', 'ge_0.20')}
        excluded = []
        for values in itertools.combinations(targets, size):
            total = sum(lengths[a] for a in values)
            group = {'group_id': '__'.join(values), 'components': list(values),
                     'length_aa': total, 'minimum_fraction': min(lengths[a] for a in values) / total}
            if total > cap:
                excluded.append({**group, 'reason': 'shared_length_cap'})
            else:
                bins[_selection_bin(group['minimum_fraction'])].append(group)
        for name in bins:
            bins[name].sort(key=lambda g: _stable(seed, g['components'], f'count{size}:{name}'))
        original_counts = {name: len(groups) for name, groups in bins.items()}
        selected = []
        while len(selected) < budget and any(bins.values()):
            for name in sorted(bins):
                if bins[name] and len(selected) < budget:
                    selected.append(bins[name].pop(0))
        rebuilt = {'group_budget': budget, 'selected_groups': selected,
                   'eligible_bin_counts': original_counts,
                   'selected_bin_counts': {name: sum(_selection_bin(g['minimum_fraction']) == name
                                                    for g in selected) for name in bins},
                   'eligible_omitted_group_ids': [g['group_id'] for groups in bins.values() for g in groups],
                   'excluded_groups': excluded,
                   'order_rule': 'all_permutations' if size == 3 else 'cyclic_rotations_and_reversals'}
        _check(design == rebuilt, f'Group selection differs from deterministic design: {size}')
        for group in selected:
            values = tuple(group['components'])
            if size == 3:
                orders = set(itertools.permutations(values))
            else:
                rotations = [values[i:] + values[:i] for i in range(size)]
                orders = set(rotations + [x[::-1] for x in rotations])
            expected_orders[group['group_id']] = orders
        group_counts[str(size)] = {'considered': math.comb(len(targets), size),
                                  'selected': len(selected), 'excluded': len(excluded),
                                  'eligible_omitted': sum(len(v) for v in bins.values())}
    observed = defaultdict(set)
    for query in manifest['queries']:
        if query['component_count'] >= 3:
            observed[query['group_id']].add(tuple(query['components']))
    _check(dict(observed) == expected_orders, 'Missing, duplicated, or unexpected higher-order query')
    return group_counts


def _validate_subsets(manifest, sequences):
    groups = defaultdict(list)
    for accession in sorted(sequences):
        groups[_digest(sequences[accession].encode())].append(accession)
    targets = {_digest(sequences[a].encode()) for a in manifest['targets']}
    decoys = sorted(set(groups) - targets)
    expected = []
    for replicate in range(3):
        ordered = sorted(decoys, key=lambda h: _stable(manifest['selection_seed'], [h], f'decoys:{replicate}'))
        for size in (25, 50, len(groups)):
            chosen = sorted(targets) + ordered[:size-len(targets)]
            expected.append({'subset_id': f'unique{size}_rep{replicate}', 'replicate': replicate,
                             'unique_sequence_count': size,
                             'entries': [{'sequence_sha256': h, 'representative_accession': min(groups[h]),
                                          'equivalent_accessions': groups[h]} for h in chosen]})
    _check(manifest['catalogue_subsets'] == expected, 'Catalogue subset design mismatch')
    return len(groups)


def load_suite(repo_root, packet_dir):
    """Validate original bytes, reviews and all design entries; return local inputs.

    ``packet_dir`` is the original prepared packet, or a complete exact copy.
    This does not mutate or refresh sources, nor treat historical review as new
    institutional approval. Any mismatch raises before returning scored inputs.
    """
    root, packet = Path(repo_root).resolve(), Path(packet_dir).resolve()
    manifest = _read(packet / 'challenge_manifest.json')
    checks = _read(packet / 'preparation_checks.json')
    _check(_digest((packet / 'challenge_manifest.json').read_bytes()) == checks['challenge_manifest_sha256'],
           'Challenge manifest checksum mismatch')
    _check(manifest['source_files_sha256'] == checks['source_files_sha256'], 'Source hash lists disagree')
    for relative, expected in checks['source_files_sha256'].items():
        _verified_path(root, relative, expected)
    for relative, expected in checks['prepared_document_sha256'].items():
        _verified_path(packet, relative, expected)
    _verified_path(packet, 'build_challenge_manifest.py', checks['builder_sha256'])
    freeze = _read(root / 'config/t05_expanded_freeze.json')
    verify_freeze(root, freeze)
    reviewed = _read(root / 'data/manifests/t05_expanded_catalogue.json')
    sequences = load_reviewed(root, reviewed)
    records = {record['accession']: record for record in reviewed['records']}
    _check(manifest['catalogue_accessions'] == sorted(sequences), 'Catalogue differs from reviewed sources')
    original = _read(root / 'config/t05_expanded_queries.json')
    _check(manifest['targets'] == sorted(original['targets']), 'Targets changed')
    _check(manifest['original_pair_exclusions'] == original['excluded_pairs'], 'Pair exclusions changed')
    target_hashes = {_digest(sequences[a].encode()) for a in manifest['targets']}
    _check(len(target_hashes) == len(manifest['targets']), 'Targets are not unique sequences')
    lengths = {a: len(s) for a, s in sequences.items()}
    # Recheck all 120 original unordered pairs and the admitted/excluded sets.
    eligible_pairs, excluded_pairs = set(), set()
    for pair in itertools.combinations(manifest['targets'], 2):
        (eligible_pairs if sum(lengths[a] for a in pair) <= manifest['shared_max_residues']
         else excluded_pairs).add(pair)
    _check({tuple(p['components']) for p in original['pairs']} == eligible_pairs, 'Original pair eligibility mismatch')
    _check({tuple(p['components']) for p in original['excluded_pairs']} == excluded_pairs, 'Original pair exclusion mismatch')
    original_queries = {q['query_id']: q for q in original['queries']}
    queries, ids, suites = [], set(), Counter()
    query_hashes = {}
    for q in manifest['queries']:
        qid, components, spans = q['query_id'], q['components'], q['spans']
        _check(qid not in ids, 'Duplicate query ID')
        ids.add(qid)
        _check(len(components) == q['component_count'] == len(spans) and 1 <= len(components) <= 4,
               'Component-count mismatch')
        _check(len(set(components)) == len(components) and set(components) <= set(manifest['targets']),
               'Invalid query components')
        _check(q['artificial'] == (len(components) > 1), 'Artificial-input designation mismatch')
        expected_suite = {1: 'clean', 2: 'pairs', 3: 'triples', 4: 'quadruples'}[len(components)]
        _check(q['suite'] == expected_suite, 'Suite designation mismatch')
        sequence = ''.join(sequences[a] for a in components)
        _check(len(sequence) == q['length_aa'] <= manifest['shared_max_residues'], 'Query length mismatch')
        pos = 0
        for a, span in zip(components, spans):
            expected = {'accession': a, 'start0': pos, 'end0': pos + lengths[a],
                        'fraction': lengths[a] / len(sequence), 'sequence_sha256': _digest(sequences[a].encode())}
            _check(span == expected, 'Component boundary, fraction or checksum mismatch')
            _check(sequence[span['start0']:span['end0']] == sequences[a], 'Component construction mismatch')
            pos = span['end0']
        if len(components) <= 2:
            oq = original_queries.get(qid)
            _check(oq is not None and oq['components'] == components, 'Original query changed')
            _check(q['group_id'] == (oq.get('pair_id') or qid), 'Original group changed')
        else:
            _check(q['group_id'] == '__'.join(sorted(components)), 'Group ID mismatch')
            _check(qid == f'join{len(components)}_' + '__'.join(components), 'Higher-order query ID mismatch')
        query_hashes[qid] = _digest(sequence.encode())
        queries.append({**q, 'sequence': sequence, 'joined_sequence_sha256': query_hashes[qid]})
        suites[q['suite']] += 1
    _check(set(original_queries) == {q['query_id'] for q in queries if q['component_count'] <= 2},
           'Original clean/pair coverage changed')
    _check(dict(suites) == manifest['suite_counts'] == checks['suite_counts'], 'Suite counts disagree')
    groups = _validate_groups(manifest, lengths)
    unique_sequences = _validate_subsets(manifest, sequences)
    repeat = list(original['repeat_query_ids'])
    for suite_name in ('triples', 'quadruples'):
        candidates = [q for q in queries if q['suite'] == suite_name]
        repeat.extend(q['query_id'] for q in sorted(candidates,
                      key=lambda q: _stable(manifest['selection_seed'], [q['query_id']], 'repeat'))[:4])
    _check(repeat == manifest['repeat_query_ids'] and len(set(repeat)) == len(repeat), 'Repeat subset mismatch')
    observed_checks = {'reviewed_records_checked': len(sequences), 'unique_sequences': unique_sequences,
                       'target_count': len(manifest['targets']), 'suite_counts': dict(suites),
                       'query_count': len(queries), 'max_query_length': max(q['length_aa'] for q in queries),
                       'combined_component_observations': sum(q['component_count'] for q in queries if q['artificial']),
                       'catalogue_subset_count': len(manifest['catalogue_subsets']), 'repeat_query_count': len(repeat),
                       'existing_freeze_hashes_verified': len(freeze['file_sha256'])}
    for key, value in observed_checks.items():
        _check(checks[key] == value, f'Preparation count mismatch: {key}')
    verification = {'status': 'passed', **observed_checks,
                    'all_component_observations': sum(q['component_count'] for q in queries),
                    'challenge_manifest_sha256': checks['challenge_manifest_sha256'],
                    'source_files_sha256': manifest['source_files_sha256'],
                    'group_verification': groups, 'query_sequence_sha256': query_hashes,
                    'review_scope': reviewed['scope'], 'attribution': reviewed['attribution'],
                    'review_note': 'All existing record-specific approvals and original bytes revalidated; not new institutional sign-off.'}
    return {'manifest': manifest, 'sequences': sequences,
            'catalogue': {a: {'sequence': sequences[a], 'metadata': records[a]} for a in sorted(sequences)},
            'queries': queries, 'equivalence': equivalence_map(sequences), 'verification': verification}


def rank_scores(scores, accessions, score_name='score'):
    """Rank descending finite scores, then accession ascending for exact ties."""
    _check(len(scores) == len(accessions), 'Scores and catalogue sizes differ')
    _check(len(set(accessions)) == len(accessions), 'Duplicate ranking accession')
    values = [float(value) for value in scores]
    _check(all(math.isfinite(v) for v in values), 'Nonfinite retrieval score')
    counts = Counter(values)
    return [{'rank': rank, 'accession': accessions[i], 'score': values[i], score_name: values[i],
             'exact_score_tie_count': counts[values[i]]}
            for rank, i in enumerate(sorted(range(len(values)), key=lambda i: (-values[i], accessions[i])), 1)]


def cosine_ranking(vector, reference_vectors, accessions, *, dtype=None):
    """NumPy cosine with input dtype preserved unless explicitly overridden."""
    import numpy as np
    query, refs = np.asarray(vector, dtype=dtype), np.asarray(reference_vectors, dtype=dtype)
    _check(query.ndim == 1 and refs.shape == (len(accessions), len(query)), 'Incompatible cosine shapes')
    _check(len(accessions) > 0 and np.isfinite(query).all() and np.isfinite(refs).all(), 'Invalid cosine vectors')
    qnorm, norms = np.linalg.norm(query), np.linalg.norm(refs, axis=1)
    _check(qnorm > 0 and (norms > 0).all(), 'Zero cosine vector')
    return rank_scores(refs @ query / (norms * qnorm), accessions, 'cosine_similarity')


def filter_ranking(ranking, included_accessions, *, tolerance=1e-6):
    """Apply a predeclared common catalogue mask and recompute ranks/tie counts."""
    included = set(included_accessions)
    available = {row['accession'] for row in ranking}
    _check(included <= available, 'Catalogue mask contains unavailable accessions')
    rows = [row for row in ranking if row['accession'] in included]
    counts = Counter(_score(row) for row in rows)
    filtered = [{**row, 'rank': rank, 'exact_score_tie_count': counts[_score(row)]}
                for rank, row in enumerate(rows, 1)]
    for row in filtered:
        score = _score(row)
        tied = [r['rank'] for r in filtered if _score(r) == score]
        near = [r['rank'] for r in filtered if abs(_score(r)-score) <= tolerance]
        row.update(tie_best_rank=min(tied), tie_worst_rank=max(tied),
                   near_tie_count=len(near), near_tie_best_rank=min(near), near_tie_worst_rank=max(near),
                   near_tie_tolerance=tolerance)
        # Older alignment records carry this fixed-tolerance spelling.
        if 'near_tie_count_1e6' in row:
            row['near_tie_count_1e6'] = sum(abs(_score(r)-score) <= 1e-6 for r in filtered)
    return filtered


def catalogue_masks(suite):
    """Unique representatives only; identical full endpoints share one mask."""
    masks = {'accession105': sorted(suite['sequences'])}
    seen = {}
    aliases = {}
    for subset in suite['manifest']['catalogue_subsets']:
        accessions = sorted(entry['representative_accession'] for entry in subset['entries'])
        key = tuple(accessions)
        if key in seen:
            aliases[subset['subset_id']] = seen[key]
        else:
            masks[subset['subset_id']] = accessions
            seen[key] = subset['subset_id']
    return {'masks': masks, 'aliases': aliases}


def common_reference_removal_mask(query, equivalence, catalogue_accessions):
    removed = set().union(*(set(equivalence[a]) for a in query['components']))
    return sorted(set(catalogue_accessions) - removed)


def _score(row):
    for name in ('score', 'cosine_similarity', 'alignment_score'):
        if name in row:
            return float(row[name])
    raise ValueError('Ranking row has no supported score')


def _match_metrics(ranking, accepted, tolerance):
    surviving = set(accepted) & {row['accession'] for row in ranking}
    if not surviving:
        return {'scorable': False, 'accepted_surviving': [], 'rank': None,
                'reciprocal_rank': None, **{f'hit_at_{k}': None for k in KS}}
    row = next(row for row in ranking if row['accession'] in surviving)
    score = _score(row)
    tied = [r['rank'] for r in ranking if _score(r) == score]
    near = [r['rank'] for r in ranking if abs(_score(r)-score) <= tolerance]
    return {'scorable': True, 'accepted_surviving': sorted(surviving),
            'best_accession': row['accession'], 'rank': row['rank'], 'score': score,
            'reciprocal_rank': 1 / row['rank'], 'tie_best_rank': min(tied), 'tie_worst_rank': max(tied),
            'near_tie_count': len(near), 'near_tie_best_rank': min(near), 'near_tie_worst_rank': max(near),
            **{f'hit_at_{k}': int(row['rank'] <= k) for k in KS}}


def _query_metrics(components, kind):
    scored = [component[kind] for component in components if component[kind]['scorable']]
    n, original = len(scored), len(components)
    result = {'component_count': original, 'scorable_components': n, 'unscorable_components': original-n,
              'all_components_scorable': n == original}
    for k in KS:
        hits = sum(row[f'hit_at_{k}'] for row in scored)
        result.update({f'component_hits_at_{k}': hits,
                       f'component_recall_at_{k}': hits/n if n else None,
                       f'any_component_hit_at_{k}': int(hits > 0) if n == original else None,
                       f'any_scorable_component_hit_at_{k}': int(hits > 0) if n else None,
                       f'all_components_hit_at_{k}': int(hits == original) if n == original else None})
    return result


def evaluate_query(ranking, query, equivalence, acceptable_matches=None, *, tolerance=1e-6):
    """Evaluate all components, retaining absent-label denominators explicitly.

    A supplied acceptable-match map is a separately frozen family/domain label
    policy, never inferred from scores. Missing/empty labels remain unscorable.
    Exact accession and duplicate equivalence are always reported separately.
    """
    _check([r['rank'] for r in ranking] == list(range(1, len(ranking)+1)), 'Noncontiguous ranking')
    accessions = [r['accession'] for r in ranking]
    _check(len(set(accessions)) == len(accessions), 'Duplicate ranking accession')
    _check(all(math.isfinite(_score(r)) for r in ranking), 'Nonfinite ranking')
    _check(ranking == sorted(ranking, key=lambda r: (-_score(r), r['accession'])), 'Unsorted ranking')
    expected = query['components']
    _check(expected and len(set(expected)) == len(expected), 'Components must be distinct')
    _check(len(query['spans']) == len(expected), 'Missing component spans')
    components = []
    for i, (accession, span) in enumerate(zip(expected, query['spans'])):
        _check(span['accession'] == accession, 'Component/span accession mismatch')
        position = 'single' if len(expected) == 1 else ('N' if i == 0 else ('C' if i == len(expected)-1 else 'internal'))
        result = {'accession': accession, 'component_index': i, 'position': position,
                  'start0': span['start0'], 'end0': span['end0'], 'fraction': span['fraction'],
                  'fraction_bin': fraction_bin(span['fraction']),
                  'exact': _match_metrics(ranking, [accession], tolerance),
                  'duplicate_equivalent': _match_metrics(ranking, equivalence.get(accession, [accession]), tolerance)}
        if acceptable_matches is not None:
            result['labelled'] = _match_metrics(ranking, acceptable_matches.get(accession, []), tolerance)
        components.append(result)
    kinds = ('exact', 'duplicate_equivalent') + (('labelled',) if acceptable_matches is not None else ())
    return {'query_id': query['query_id'], 'suite': query.get('suite'), 'group_id': query.get('group_id', query['query_id']),
            'component_count': len(expected), 'catalogue_count': len(ranking), 'components': components,
            'query_metrics': {kind: _query_metrics(components, kind) for kind in kinds},
            'exact_recall_ceiling': {str(k): min(1, k/len(expected)) for k in KS},
            'retrieval_budget': 'k total per original query',
            'returned_hits_ceiling': {str(k): min(k, len(ranking)) for k in KS},
            'near_tie_tolerance': tolerance}


def evaluate_oracle_query(rankings_by_accession, query, equivalence, acceptable_matches=None, *, tolerance=1e-6):
    """Combine boundary-known component retrievals with their larger budget.

    Each component is evaluated using its own list. A hit obtained only from a
    different component's list is not borrowed. This is a supplied-boundary
    oracle, and its per-component budget must not be called k total per query.
    """
    expected = set(query['components'])
    _check(set(rankings_by_accession) == expected, 'Oracle rankings must cover exactly the query components')
    evaluations = [evaluate_query(rankings_by_accession[a], query, equivalence,
                                 acceptable_matches, tolerance=tolerance) for a in query['components']]
    masks = [{r['accession'] for r in rankings_by_accession[a]} for a in query['components']]
    _check(all(mask == masks[0] for mask in masks), 'Oracle components must use a common catalogue mask')
    result = dict(evaluations[0])
    components = [e['components'][i] for i, e in enumerate(evaluations)]
    result.update(components=components,
                  query_metrics={kind: _query_metrics(components, kind) for kind in result['query_metrics']},
                  retrieval_budget='k per known component',
                  exact_recall_ceiling={str(k): 1.0 for k in KS},
                  returned_hits_ceiling={str(k): sum(min(k, len(rows)) for rows in rankings_by_accession.values()) for k in KS})
    return result


def summarize(evaluations):
    """Aggregate one fixed model/method/arm/catalogue; reject duplicate queries.

    This function intentionally does not silently mix configurations. Caller must
    partition them before calling. Component strata retain parent query IDs.
    """
    evaluations = list(evaluations)
    ids = [e['query_id'] for e in evaluations]
    _check(len(set(ids)) == len(ids), 'Duplicate query evaluations in summary')
    if not evaluations:
        return {'query_count': 0, 'component_count': 0, 'kinds': {}}
    kinds = tuple(evaluations[0]['query_metrics'])
    _check(all(tuple(e['query_metrics']) == kinds for e in evaluations), 'Mixed label policies')
    _check(len({e['retrieval_budget'] for e in evaluations}) == 1, 'Mixed retrieval budgets')
    output = {'query_count': len(evaluations), 'distinct_group_count': len({e['group_id'] for e in evaluations}),
              'component_count': sum(e['component_count'] for e in evaluations),
              'retrieval_budget': evaluations[0]['retrieval_budget'], 'kinds': {}}
    for kind in kinds:
        entries = [(e, c) for e in evaluations for c in e['components']]
        def component_summary(rows):
            scored = [c[kind] for _, c in rows if c[kind]['scorable']]
            n = len(scored)
            targets = defaultdict(list)
            for _, c in rows:
                if c[kind]['scorable']:
                    targets[c['accession']].append(c[kind])
            return {'original_components': len(rows), 'scorable_components': n,
                    'unscorable_components': len(rows)-n,
                    'parent_queries': len({e['query_id'] for e, _ in rows}),
                    'distinct_groups': len({e['group_id'] for e, _ in rows}),
                    'scorable_targets': len(targets),
                    'mean_rank': sum(r['rank'] for r in scored)/n if n else None,
                    'mrr': sum(r['reciprocal_rank'] for r in scored)/n if n else None,
                    **{f'hits_at_{k}': sum(r[f'hit_at_{k}'] for r in scored) for k in KS},
                    **{f'recall_at_{k}': sum(r[f'hit_at_{k}'] for r in scored)/n if n else None for k in KS},
                    **{f'per_target_macro_recall_at_{k}':
                       sum(sum(r[f'hit_at_{k}'] for r in rs)/len(rs) for rs in targets.values())/len(targets)
                       if targets else None for k in KS}}
        def query_summary(rows):
            distinct = {e['query_id']: e for e, _ in rows}
            complete = [e for e in distinct.values() if e['query_metrics'][kind]['all_components_scorable']]
            summary = {'original_queries': len(distinct), 'all_components_scorable_queries': len(complete),
                       'excluded_unscorable_queries': len(distinct)-len(complete)}
            for k in KS:
                for metric in ('any_component_hit', 'all_components_hit'):
                    count = sum(e['query_metrics'][kind][f'{metric}_at_{k}'] for e in complete)
                    summary[f'{metric}_count_at_{k}'] = count
                    summary[f'{metric}_at_{k}'] = count/len(complete) if complete else None
            return summary
        result = component_summary(entries)
        result['query_metrics'] = query_summary(entries)
        for field, getter in (('by_component_count', lambda e,c: str(e['component_count'])),
                              ('by_fraction', lambda e,c: c['fraction_bin']),
                              ('by_position', lambda e,c: c['position']),
                              ('by_target', lambda e,c: c['accession'])):
            grouped = defaultdict(list)
            for e, c in entries:
                grouped[getter(e,c)].append((e,c))
            result[field] = {key: component_summary(rows) for key, rows in sorted(grouped.items())}
            if field == 'by_component_count':
                for key, rows in grouped.items():
                    result[field][key]['query_metrics'] = query_summary(rows)
        output['kinds'][kind] = result
    return output
