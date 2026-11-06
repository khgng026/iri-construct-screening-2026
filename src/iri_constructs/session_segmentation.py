"""Label-free candidate segmentation and fixed-budget ranking aggregation.

Prediction functions take only sequence/length, catalogue and fixed parameters.
Constructed target spans enter ``interval_metrics`` only after prediction.
The thresholds below generate candidates; they are not calibrated detection.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
import hashlib
import math
from time import perf_counter

from .session_benchmark import rank_scores


@dataclass(frozen=True)
class SegmentationConfig:
    minimum_span: int = 32
    minimum_reference_coverage: float = .5
    minimum_identity_fraction: float = .25
    minimum_score: float = 60.
    maximum_segments: int = 6
    empty_prediction_fallback: str = 'whole_query'

    def __post_init__(self):
        if not isinstance(self.minimum_span, int) or self.minimum_span < 1:
            raise ValueError('Minimum span must be a positive integer')
        if not isinstance(self.maximum_segments, int) or self.maximum_segments < 1:
            raise ValueError('Segment cap must be a positive integer')
        if not (0 <= self.minimum_reference_coverage <= 1 and 0 <= self.minimum_identity_fraction <= 1):
            raise ValueError('Coverage and identity thresholds must be in [0,1]')
        if not math.isfinite(self.minimum_score) or self.minimum_score < 0:
            raise ValueError('Score threshold must be finite and nonnegative')
        if self.empty_prediction_fallback != 'whole_query':
            raise ValueError('Frozen empty-prediction fallback is whole_query')


def fixed_windows(length, width, stride=None):
    """Cover both termini using fixed overlapping windows, without labels."""
    stride = width // 2 if stride is None else stride
    if not all(isinstance(v, int) and v > 0 for v in (length, width, stride)) or stride > width:
        raise ValueError('Length/width/stride must be positive integers, stride <= width')
    if length <= width:
        starts = [0]
    else:
        starts = list(range(0, length-width+1, stride))
        if starts[-1] != length-width:
            starts.append(length-width)
    spans = [{'start0': start, 'end0': min(length, start+width)} for start in starts]
    return {'segments': spans, 'segment_count': len(spans), 'query_length': length,
            'total_queried_residues': sum(s['end0']-s['start0'] for s in spans),
            'width': width, 'stride': stride, 'cap_bound': False,
            'method': 'fixed_overlapping_windows', 'boundaries_supplied': False}


def alignment_guided_segments(sequence, catalogue, config=None, *, ranker=None):
    """Iteratively remove strongest accepted local span from disjoint intervals.

    Every remaining interval is aligned exhaustively against the same catalogue.
    Its best *eligible* hit is considered, rather than assuming its highest raw
    hit passes the fixed candidate filters. The globally strongest candidate is
    accepted; score ties use accession then start/end coordinates. No expected
    component count, accession, boundary or family label enters this function.

    ``ranker`` is an injectable pure sequence/catalogue interface for tests.
    Default: the pinned exhaustive Smith-Waterman comparator. Returned proposal
    scores belong to the source fragment, not an independently encoded segment.
    Re-encode returned spans for end-to-end embedding retrieval and fusion.
    """
    from .session_alignment import rank_alignment
    if not isinstance(sequence, str) or not sequence or not catalogue:
        raise ValueError('A nonempty sequence and catalogue are required')
    config = config or SegmentationConfig()
    ranker = ranker or rank_alignment
    started = perf_counter()
    active = [(0, len(sequence))]
    cache, attempts, accepted, discarded = {}, [], [], []
    queried = 0
    alignments = 0

    def interval_candidate(interval):
        nonlocal queried, alignments
        if interval in cache:
            return cache[interval]
        start, end = interval
        result = ranker(sequence[start:end], catalogue)
        ranking = result['ranking'] if isinstance(result, dict) else result
        if {row['accession'] for row in ranking} != set(catalogue) or len(ranking) != len(catalogue):
            raise ValueError('Candidate ranking must cover the full fixed catalogue')
        eligible = []
        for row in ranking:
            score = float(row['score'])
            if not math.isfinite(score):
                raise ValueError('Nonfinite candidate score')
            span = row.get('query_span')
            if span is None:
                continue
            q0, q1 = span
            if not (isinstance(q0, int) and isinstance(q1, int) and 0 <= q0 < q1 <= end-start):
                raise ValueError('Invalid local candidate coordinates')
            coverage, identity = float(row['reference_coverage']), float(row['identity_fraction'])
            if not (math.isfinite(coverage) and math.isfinite(identity) and 0 <= coverage <= 1 and 0 <= identity <= 1):
                raise ValueError('Invalid candidate coverage or identity')
            if (q1-q0 >= config.minimum_span and score >= config.minimum_score and
                    coverage >= config.minimum_reference_coverage and identity >= config.minimum_identity_fraction):
                eligible.append({'start0': start+q0, 'end0': start+q1,
                                 'source_interval': [start, end], 'seed_hit': dict(row)})
        eligible.sort(key=lambda hit: (-hit['seed_hit']['score'], hit['seed_hit']['accession'],
                                      hit['start0'], hit['end0']))
        candidate = eligible[0] if eligible else None
        attempts.append({'start0': start, 'end0': end, 'catalogue_count': len(ranking),
                         'eligible_reference_count': len(eligible),
                         'best_eligible_accession': candidate['seed_hit']['accession'] if candidate else None,
                         'best_eligible_score': candidate['seed_hit']['score'] if candidate else None})
        queried += end-start
        alignments += len(ranking)
        cache[interval] = candidate
        return candidate

    while active and len(accepted) < config.maximum_segments:
        candidates = []
        for interval in list(active):
            if interval[1]-interval[0] < config.minimum_span:
                active.remove(interval)
                discarded.append({'start0': interval[0], 'end0': interval[1], 'reason': 'below_minimum_span'})
                continue
            candidate = interval_candidate(interval)
            if candidate is None:
                active.remove(interval)
                discarded.append({'start0': interval[0], 'end0': interval[1], 'reason': 'no_eligible_reference'})
            else:
                candidates.append(candidate)
        if not candidates:
            break
        candidates.sort(key=lambda hit: (-hit['seed_hit']['score'], hit['seed_hit']['accession'],
                                        hit['start0'], hit['end0']))
        best = candidates[0]
        accepted.append({**best, 'selection_order': len(accepted)})
        left, right = best['source_interval']
        active.remove((left, right))
        for interval in ((left, best['start0']), (best['end0'], right)):
            if interval[1] > interval[0]:
                if interval[1]-interval[0] >= config.minimum_span:
                    active.append(interval)
                else:
                    discarded.append({'start0': interval[0], 'end0': interval[1], 'reason': 'below_minimum_span'})
    # This conservative cap indicator means an admissibly long interval remains,
    # not proof that it contains another accepted hit; no unbudgeted lookahead.
    cap_bound = len(accepted) == config.maximum_segments and bool(active)
    spans = sorted(accepted, key=lambda hit: (hit['start0'], hit['end0']))
    if any(a['end0'] > b['start0'] for a, b in zip(spans, spans[1:])):
        raise ValueError('Iterative segmentation produced overlapping spans')
    fallback = not spans
    if fallback:
        spans = [{'start0': 0, 'end0': len(sequence), 'source_interval': [0, len(sequence)],
                  'seed_hit': None, 'selection_order': None, 'fallback': True}]
    return {'method': 'iterative_alignment_candidates', 'segments': spans, 'segment_count': len(spans),
            'accepted_alignment_segment_count': len(accepted), 'whole_query_fallback_used': fallback,
            'query_length': len(sequence), 'config': asdict(config), 'boundaries_supplied': False,
            'threshold_status': 'prespecified candidate filters; not calibrated detection',
            'total_queried_residues': queried,
            'total_predicted_segment_residues': sum(s['end0']-s['start0'] for s in spans),
            'alignment_comparisons': alignments, 'intervals_aligned': len(attempts),
            'interval_attempts': attempts, 'discarded_intervals': discarded,
            'unprocessed_intervals': [list(interval) for interval in sorted(active)],
            'cap_bound': cap_bound, 'cap_bound_definition': 'segment cap reached with an interval >= minimum_span remaining',
            'seconds': perf_counter()-started}


def _validated_rankings(rankings):
    rankings = list(rankings)
    if not rankings:
        return [], []
    expected = {row['accession'] for row in rankings[0]}
    if not expected:
        raise ValueError('Empty segment catalogue')
    for rows in rankings:
        if len(rows) != len(expected) or {r['accession'] for r in rows} != expected:
            raise ValueError('Segment rankings must use one complete common catalogue')
        if [r['rank'] for r in rows] != list(range(1, len(rows)+1)):
            raise ValueError('Segment ranks must be contiguous')
    return rankings, sorted(expected)


def fuse_rankings(rankings, method='rrf', *, rrf_constant=60, catalogue_accessions=None):
    """Label-free complete ranking; a fused top k uses k hits per input query.

    RRF is the sum of 1/(60+rank) over all segments, with accession tie breaks.
    Max-score must be applied only to same-model cosine rankings, never scores
    pooled across representations. No expected answer selects a best segment.
    """
    rankings, accessions = _validated_rankings(rankings)
    if method not in ('rrf', 'max_score') or rrf_constant != 60:
        raise ValueError('Frozen fusion choices are RRF constant 60 or max_score')
    if not rankings:
        # A no-candidate query is explicitly retained. Alphabetical ties are not
        # evidence-backed predictions and must not be scored as retrieval hits.
        return {'ranking': [], 'segment_count': 0, 'candidate_status': 'no_candidates',
                'catalogue_accessions': sorted(catalogue_accessions or []), 'method': method,
                'retrieval_budget': 'k total per original query', 'rrf_constant': rrf_constant}
    if catalogue_accessions is not None and set(catalogue_accessions) != set(accessions):
        raise ValueError('Fusion catalogue mismatch')
    values = {accession: [] for accession in accessions}
    for rows in rankings:
        for row in rows:
            if method == 'rrf':
                value = 1/(rrf_constant+row['rank'])
            else:
                if 'cosine_similarity' not in row:
                    raise ValueError('max_score requires explicit same-model cosine scores')
                value = float(row['cosine_similarity'])
                if not math.isfinite(value):
                    raise ValueError('Nonfinite segment cosine')
            values[row['accession']].append(value)
    scores = [math.fsum(values[a]) if method == 'rrf' else max(values[a]) for a in accessions]
    return {'ranking': rank_scores(scores, accessions, 'fusion_score'), 'segment_count': len(rankings),
            'candidate_status': 'candidates_available', 'method': method,
            'rrf_constant': rrf_constant, 'retrieval_budget': 'k total per original query'}


def per_segment_union(rankings, k):
    """Explicit larger-budget candidate list, with duplicates/review volume."""
    if not isinstance(k, int) or k < 1:
        raise ValueError('Returned hit budget must be a positive integer')
    rankings, _ = _validated_rankings(rankings)
    occurrences = defaultdict(list)
    for segment, rows in enumerate(rankings):
        for row in rows[:k]:
            occurrences[row['accession']].append({'segment_index': segment, 'rank': row['rank']})
    volume = sum(min(k, len(rows)) for rows in rankings)
    return {'accessions': sorted(occurrences), 'occurrences': dict(occurrences),
            'per_segment_k': k, 'segment_count': len(rankings),
            'returned_hit_occurrences': volume, 'unique_returned_hits': len(occurrences),
            'duplicate_returned_hits': volume-len(occurrences),
            'retrieval_budget': 'k per predicted segment'}


def interval_metrics(predicted, targets):
    """Evaluation only: overlap diagnostics using supplied ground-truth spans."""
    def span(value):
        a, b = value['start0'], value['end0']
        if not isinstance(a, int) or not isinstance(b, int) or not 0 <= a < b:
            raise ValueError('Invalid evaluation interval')
        return a, b
    predicted_spans = [span(value) for value in predicted]
    rows = []
    for target in targets:
        start, end = span(target)
        overlaps = []
        intersections = []
        for a, b in predicted_spans:
            overlap = max(0, min(b, end)-max(a, start))
            overlaps.append(overlap / (max(b, end)-min(a, start)))
            if overlap:
                intersections.append((max(a,start), min(b,end)))
        covered = 0
        right = start
        for a, b in sorted(intersections):
            covered += max(0, b-max(a, right))
            right = max(right, b)
        rows.append({'accession': target.get('accession'), 'start0': start, 'end0': end,
                     'best_interval_iou': max(overlaps, default=0),
                     'covered_fraction_by_union': covered/(end-start),
                     'exact_boundary_match': (start,end) in predicted_spans})
    return {'components': rows, 'predicted_segment_count': len(predicted_spans),
            'target_component_count': len(rows), 'evaluation_only': True}


def grouped_partition(queries, seed='iri-segmentation-development-v1'):
    """SHA-256 20% group split; all permutations share an assignment."""
    groups = {}
    assignments = {}
    component_groups = {}
    for query in queries:
        group = query['group_id']
        component_set = tuple(sorted(query['components']))
        if group in groups and groups[group]['components'] != component_set:
            raise ValueError('One group ID refers to different component sets')
        if component_set in component_groups and component_groups[component_set] != group:
            raise ValueError('One unordered component set has multiple group IDs')
        component_groups[component_set] = group
        digest = hashlib.sha256((seed+'|'+group).encode()).hexdigest()
        partition = 'dev' if int(digest,16) < (1 << 256)//5 else 'test'
        groups[group] = {'components': component_set, 'partition': partition, 'sha256': digest}
        assignments[query['query_id']] = partition
    if len(assignments) != len(queries):
        raise ValueError('Duplicate query IDs')
    return {'seed': seed, 'rule': 'int(SHA256(seed|group_id),16) < floor(2^256/5) assigns dev',
            'groups': groups, 'queries': assignments}


def select_configuration(candidate_results, design):
    """Dev-only prespecified group-macro selection; rejects missing/extra IDs.

    Each candidate supplies rows with query_id, component_recall_at_5,
    all_components_hit_at_5 and total_queried_residues. The primary score is
    the average of within-group mean component recall@5, then similarly all-hit
    @5, then lower queried residues per query, then configuration ID ascending.
    """
    expected = set(design['selection_query_ids'])
    if not expected:
        raise ValueError('Development selection panel is empty')
    lookup = {qid: group for group, item in design['group_queries'].items() for qid in item}
    declared = {c['configuration_id'] for c in design['candidate_configurations']}
    if {c['configuration_id'] for c in candidate_results} != declared or len(candidate_results) != len(declared):
        raise ValueError('Must compare every prespecified configuration exactly once')
    scored = []
    for candidate in candidate_results:
        rows = candidate['rows']
        if len(rows) != len(expected) or {r['query_id'] for r in rows} != expected:
            raise ValueError('Selection rows must cover exactly the frozen combined dev panel')
        groups = defaultdict(list)
        for row in rows:
            if design['partition']['queries'][row['query_id']] != 'dev':
                raise ValueError('Test query cannot influence parameter selection')
            recall, allhit, work = (row['component_recall_at_5'], row['all_components_hit_at_5'],
                                    row['total_queried_residues'])
            if not (math.isfinite(recall) and 0 <= recall <= 1 and allhit in (0,1) and
                    math.isfinite(work) and work >= 0):
                raise ValueError('Invalid development summary metric')
            groups[lookup[row['query_id']]].append(row)
        macro_recall = math.fsum(math.fsum(r['component_recall_at_5'] for r in rs)/len(rs)
                                 for rs in groups.values())/len(groups)
        macro_allhit = math.fsum(math.fsum(r['all_components_hit_at_5'] for r in rs)/len(rs)
                                 for rs in groups.values())/len(groups)
        work = math.fsum(r['total_queried_residues'] for r in rows)/len(rows)
        scored.append({'configuration_id': candidate['configuration_id'],
                       'group_macro_component_recall_at_5': macro_recall,
                       'group_macro_all_components_hit_at_5': macro_allhit,
                       'mean_total_queried_residues': work, 'dev_query_count': len(rows),
                       'dev_group_count': len(groups)})
    scored.sort(key=lambda r: (-r['group_macro_component_recall_at_5'],
                              -r['group_macro_all_components_hit_at_5'],
                              r['mean_total_queried_residues'], r['configuration_id']))
    segmenters = {c['configuration_id']: c.get('segmenter', 'unspecified') for c in design['candidate_configurations']}
    selected_by_segmenter = {}
    for row in scored:
        selected_by_segmenter.setdefault(segmenters[row['configuration_id']], row['configuration_id'])
    return {'selected_configuration_id': scored[0]['configuration_id'],
            'selected_by_segmenter': selected_by_segmenter, 'development_scores': scored,
            'selection_scope': 'combined dev queries only; test remains untouched',
            'must_freeze_before_test_scoring': True}
