"""Frozen reviewed-panel scoring shared by pilot experiments; standard library only.

No sequences are serialized and no model scores are manufactured. Family retrieval
is annotation agreement, distinct from identifying biochemical function or proving
the identity of every biological component in a construct.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import random
import statistics

from .metrics import candidates, reviewed_labels
from .runtime import InputError, read_json, require_text

ARCHETYPES = {'single_cds', 'multi_cds', 'tagged_fusion', 'backbone_heavy', 'unknown'}
MODELS = {'alignment', 'esm2', 'esmc', 'esm3'}


def _hash(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode('utf-8')).hexdigest()


def _query_clusters(query_ids, fingerprints):
    """Conservatively join query observations with documented shared provenance.

    Distinct protein-family groups do not make two CDSs from one construct
    independent. Exact query copies and indirect parent/group links likewise
    belong to one connected cluster, fixed before any method is scored.
    """
    components, tokens_by_query = [], {}
    for qid in sorted(query_ids):
        fp = fingerprints[qid]
        tokens = {'parent:' + fp['construct_id'],
                  'group:' + fp['independence_group'],
                  'sequence:' + fp['sequence_sha256']}
        tokens_by_query[qid] = tokens.copy()
        for component in list(components):
            if tokens & component:
                tokens |= component
                components.remove(component)
        components.append(tokens)
    return {qid: 'query-cluster-' + _hash(sorted(next(
                c for c in components if c & tokens_by_query[qid])))[:20]
            for qid in sorted(query_ids)}


def query_independence_clusters(panel):
    """Derive cluster membership from frozen evidence, including legacy panels."""
    validate_panel(panel)
    return _query_clusters(panel['query_ids'], panel['input_fingerprints'])


def load_scoring_config(path: str | Path) -> dict:
    config = read_json(Path(path))
    validate_scoring_config(config)
    return config


def validate_scoring_config(config: dict) -> None:
    if not isinstance(config, dict) or config.get('schema_version') != 'pilot-scoring-1':
        raise InputError('Expected pilot-scoring-1 configuration')
    for key in ('protocol_id', 'protocol_version'):
        require_text(config, key)
    if config.get('experiment') not in ('A', 'B', 'C', 'D'):
        raise InputError('Pilot experiment must be A, B, C or D')
    if config.get('design_status') != 'design_frozen_pending_admission':
        raise InputError('Protocol design must be frozen before panel admission/scoring')
    limits = {'top_k': (1, 100), 'retrieval_budget': (1, 100), 'seed': (0, 2**32-1),
              'min_independent_clusters': (2, 10000), 'bootstrap_replicates': (100, 10000),
              'max_residues': (1, 1022)}
    for key, (low, high) in limits.items():
        value = config.get(key)
        if type(value) is not int or not low <= value <= high:
            raise InputError(f'Invalid pilot configuration field: {key}')
    if config['top_k'] != config['retrieval_budget']:
        raise InputError('top_k and total retrieval_budget must match')
    models = config.get('models')
    if not isinstance(models, list) or any(m not in MODELS for m in models) or len(models) != len(set(models)):
        raise InputError('models must be a list of unique supported methods')
    if config['experiment'] != 'A' and not models:
        raise InputError('Inference/retrieval experiments require an explicit method list')
    if config['experiment'] == 'B' and models != ['esm2']:
        raise InputError('Prioritized fusion comparison uses one fixed ESM-2 adapter')
    if config['experiment'] == 'A' and models:
        raise InputError('Annotation-preservation experiment A does not load protein models')
    if config['experiment'] == 'D' and ('alignment' in models or len(models) < 2):
        raise InputError('Model-cost experiment D requires at least two protein model adapters')
    if config['experiment'] == 'C':
        if 'alignment' not in models or len(models) < 2:
            raise InputError('Experiment C requires alignment and at least one protein embedding method')
        alignment = config.get('alignment', {})
        if (alignment.get('mode') != 'local' or alignment.get('substitution_matrix') != 'BLOSUM62'
                or alignment.get('open_gap_score') != -10.0 or alignment.get('extend_gap_score') != -0.5
                or alignment.get('primary') != 'raw_score'
                or alignment.get('sensitivity') != 'reference_self_normalized'):
            raise InputError('Alignment settings must match the preregistered local BLOSUM62 comparison')
    if not isinstance(config.get('model_revisions', {}), dict):
        raise InputError('model_revisions must be an object')


def freeze_panel(config: dict, proteins: dict, labels: dict, query_ids=None,
                 reference_ids=None, component_families=None) -> dict:
    """Bind exact identities and annotations before scores are produced.

    Query/reference roles may overlap. Existing leave-construct-out, exact-copy
    and independence-group exclusions still apply. The caller must enforce T02
    source approval/checksums before providing proteins; this function validates
    the separate known/approved annotation labels rather than approving sources.
    """
    validate_scoring_config(config)
    if any(key != entry.get('input_id') for key, entry in labels.items()):
        raise InputError('Annotation map keys must match their reviewed input_id')
    reviewed, excluded = reviewed_labels({'families': list(labels.values())}, proteins)
    if excluded or set(reviewed) != set(labels):
        raise InputError('Frozen scoring panel contains unreviewed or missing protein annotations')
    for ids in (query_ids, reference_ids):
        if ids is not None and (not isinstance(ids, list) or any(not isinstance(i, str) for i in ids)):
            raise InputError('Panel roles must be explicit lists of reviewed input IDs')
    queries = list(sorted(labels) if query_ids is None else query_ids)
    references = list(sorted(labels) if reference_ids is None else reference_ids)
    for ids in (queries, references):
        if len(ids) != len(set(ids)) or any(i not in reviewed for i in ids):
            raise InputError('Panel IDs must be unique and have reviewed known annotations')
    queries, references = sorted(queries), sorted(references)
    expected = {}
    for qid in queries:
        families = [labels[qid]['family']] if component_families is None else component_families.get(qid)
        if not isinstance(families, list) or not families or any(not isinstance(f, str) or not f.strip() for f in families):
            raise InputError('Every query requires a frozen nonempty documented-family expectation')
        expected[qid] = sorted(set(families))
    if component_families is not None and set(component_families) != set(queries):
        raise InputError('Documented-family expectations must match the exact query panel')
    fingerprints = {}
    for iid in sorted(set(queries) | set(references)):
        archetype = labels[iid].get('archetype', 'unknown')
        if archetype not in ARCHETYPES:
            raise InputError('Invalid annotation archetype in frozen panel')
        fingerprints[iid] = {'construct_id': proteins[iid].construct_id,
            'sequence_sha256': proteins[iid].sequence_sha256, 'family': labels[iid]['family'],
            'source_sha256': proteins[iid].source_sha256, 'feature_id': proteins[iid].feature_id,
            'independence_group': labels[iid]['independence_group'], 'archetype': archetype}
    allowed, exclusions = {}, {}
    for qid in queries:
        allowed[qid] = [r for r in candidates(qid, labels, proteins) if r in references]
        observed = {labels[r]['family'] for r in allowed[qid]}
        if not set(expected[qid]) <= observed:
            exclusions[qid] = 'No independent reference for every documented query family'
        elif not observed - set(expected[qid]):
            exclusions[qid] = 'No independent competing-family reference'
    status = 'frozen' if any(q not in exclusions for q in queries) else 'blocked_insufficient_panel'
    panel = {'schema_version': 'pilot-panel-1', 'status': status,
        'protocol_id': config['protocol_id'], 'protocol_version': config['protocol_version'],
        'protocol_sha256': _hash(config), 'query_ids': queries, 'reference_ids': references,
        'eligible_reference_ids': allowed, 'query_exclusions': exclusions,
        'budget_saturated_query_ids': [q for q in queries if q not in exclusions and len(allowed[q]) <= config['top_k']],
        'budget_saturation_note': 'Returning every eligible reference is descriptive only, excluded from method-advantage intervals.',
        'expected_families_by_query': expected, 'input_fingerprints': fingerprints,
        'top_k': config['top_k'], 'retrieval_budget': config['retrieval_budget'],
        'tie_policy': 'descending method score, then ascending stable reference ID',
        'claim_scope': 'Reviewed annotation-family recovery only; no biological-function inference'}
    panel['query_independence_clusters'] = _query_clusters(queries, fingerprints)
    panel['query_cluster_rules'] = ['shared_source_parent', 'reviewed_independence_group',
                                    'exact_query_sequence_copy', 'transitive_links']
    panel['panel_sha256'] = _hash(panel)
    return panel


def validate_panel(panel: dict) -> None:
    if not isinstance(panel, dict) or panel.get('schema_version') != 'pilot-panel-1':
        raise InputError('Expected frozen pilot-panel-1 object')
    digest = panel.get('panel_sha256')
    value = {key: val for key, val in panel.items() if key != 'panel_sha256'}
    if digest != _hash(value):
        raise InputError('Frozen panel changed after admission; refuse scoring')


def rank_scores(scores_by_query: dict, panel: dict, top_k=None) -> dict:
    """Rank a complete fixed eligible score panel; missing/nonfinite rows abstain.

    Scores for excluded references are an error, not silently accepted leakage.
    All method/reference scores must exist before ranking so a failed model does
    not benefit from dropping a difficult reference. This function has no model.
    """
    validate_panel(panel)
    if top_k is not None and top_k != panel['top_k']:
        raise InputError('Scoring cannot change the frozen return budget')
    if not isinstance(scores_by_query, dict) or set(scores_by_query) - set(panel['query_ids']):
        raise InputError('Scores contain queries outside the frozen panel')
    rankings = {}
    for qid in panel['query_ids']:
        if qid in panel['query_exclusions']:
            continue
        scores = scores_by_query.get(qid)
        if scores is None:
            continue
        if not isinstance(scores, dict):
            raise InputError('Reference scores must be an object keyed by frozen reference ID')
        allowed = panel['eligible_reference_ids'][qid]
        if set(scores) - set(allowed):
            raise InputError('Reference scores contain leakage or IDs outside the frozen eligible panel')
        if set(scores) != set(allowed):
            continue
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in scores.values()):
            continue
        rankings[qid] = sorted(scores, key=lambda rid: (-scores[rid], rid))[:panel['top_k']]
    return rankings


def _aggregate(rows: list[dict]) -> dict:
    n = len(rows)
    return {'scored_queries': n, 'family_hit_at_k': statistics.mean(r['family_hit_at_k'] for r in rows) if rows else None,
        'documented_family_recall_at_k': statistics.mean(r['documented_family_recall_at_k'] for r in rows) if rows else None,
        'all_documented_families_recovered_rate': statistics.mean(r['all_documented_families_recovered'] for r in rows) if rows else None,
        'returned_hits': sum(r['returned_hits'] for r in rows),
        'documented_family_observations': sum(r['documented_family_count'] for r in rows),
        'recovered_family_observations': sum(r['recovered_family_count'] for r in rows)}


def score_rankings(rankings: dict, panel: dict, labels: dict, *, component_families=None) -> dict:
    """Score reviewed family relationships, never biochemical function/presence.

    Per-construct summaries average protein-target queries and state total hits;
    callers comparing whole versus decomposed constructs must merge unit scores
    before truncating one parent-level ranking to the frozen total budget.
    """
    validate_panel(panel)
    clusters = query_independence_clusters(panel)
    if not isinstance(rankings, dict) or set(rankings) - set(panel['query_ids']):
        raise InputError('Rankings contain queries outside the frozen panel')
    for iid, fingerprint in panel['input_fingerprints'].items():
        label = labels.get(iid, {})
        if (label.get('family') != fingerprint['family'] or
                label.get('independence_group') != fingerprint['independence_group'] or
                label.get('annotation_status') != 'known' or label.get('review_status') != 'approved'):
            raise InputError('Scoring annotations differ from the frozen reviewed panel')
    if component_families is not None:
        if {k: sorted(set(v)) for k,v in component_families.items()} != panel['expected_families_by_query']:
            raise InputError('Documented-family expectations changed after panel freeze')
    rows, abstentions = [], []
    for qid in panel['query_ids']:
        hits = rankings.get(qid)
        if qid in panel['query_exclusions']:
            if hits:
                raise InputError('An ineligible query was assigned a scored ranking')
            abstentions.append({'query_id': qid, 'reason': panel['query_exclusions'][qid]})
            continue
        if hits is None or hits == []:
            abstentions.append({'query_id': qid, 'reason': 'Missing/incompatible method input or incomplete score panel'})
            continue
        if not isinstance(hits, list) or len(hits) != len(set(hits)):
            raise InputError('Ranking must contain unique reference IDs')
        if len(hits) > panel['retrieval_budget']:
            raise InputError('Ranking exceeds the frozen total return budget')
        allowed = panel['eligible_reference_ids'][qid]
        if any(r not in allowed for r in hits):
            raise InputError('Ranking contains same-parent, duplicate, linked-group or out-of-panel leakage')
        if len(hits) != min(panel['top_k'], len(allowed)):
            abstentions.append({'query_id': qid, 'reason': 'Method returned fewer hits than the fixed available budget'})
            continue
        expected = set(panel['expected_families_by_query'][qid])
        retrieved = {labels[r]['family'] for r in hits}
        matched = expected & retrieved
        fp = panel['input_fingerprints'][qid]
        rows.append({'query_id': qid, 'construct_id': fp['construct_id'], 'archetype': fp['archetype'],
            'independence_group': fp['independence_group'], 'expected_families': sorted(expected),
            'independence_cluster': clusters[qid],
            'budget_saturated': qid in panel.get('budget_saturated_query_ids', []),
            'retrieved_families': sorted(retrieved), 'reference_ids': hits,
            'family_hit_at_k': int(fp['family'] in retrieved),
            'documented_family_count': len(expected), 'recovered_family_count': len(matched),
            'documented_family_recall_at_k': len(matched)/len(expected),
            'all_documented_families_recovered': int(expected <= retrieved), 'returned_hits': len(hits)})
    overall = _aggregate(rows)
    overall['micro_documented_family_recall_at_k'] = (overall['recovered_family_observations']/overall['documented_family_observations']
                                                    if overall['documented_family_observations'] else None)
    constructs = {cid: _aggregate([r for r in rows if r['construct_id'] == cid])
                  for cid in sorted({fp['construct_id'] for fp in panel['input_fingerprints'].values()})
                  if any(r['construct_id'] == cid for r in rows)}
    per_archetype = {a: _aggregate([r for r in rows if r['archetype'] == a])
                     for a in sorted({panel['input_fingerprints'][q]['archetype'] for q in panel['query_ids']})}
    return {'status': 'scored' if rows else 'blocked_no_scorable_queries', 'panel_sha256': panel['panel_sha256'],
        'requested_queries': len(panel['query_ids']), 'eligible_queries': len(panel['query_ids'])-len(panel['query_exclusions']),
        'scored_queries': len(rows), 'coverage': len(rows)/len(panel['query_ids']) if panel['query_ids'] else 0.0,
        'overall': overall, 'per_construct': constructs, 'per_archetype': per_archetype,
        'macro_construct_documented_family_recall_at_k': statistics.mean(x['documented_family_recall_at_k'] for x in constructs.values()) if constructs else None,
        'rows': rows, 'abstentions': abstentions, 'claim_scope': panel['claim_scope']}


def paired_cluster_bootstrap(left: dict, right: dict, clusters: dict, *, seed=0,
                             replicates=2000, min_clusters=8) -> dict:
    """Percentile CI for right-minus-left mean on matched independent clusters.

    Each independent group has equal weight, regardless of how many correlated
    protein observations it contains. Missing method pairs are reported rather
    than counted as zero. The estimate is exploratory for a small reviewed pilot.
    """
    if type(seed) is not int or not 0 <= seed < 2**32 or type(replicates) is not int or not 100 <= replicates <= 10000:
        raise InputError('Invalid bootstrap seed or replicate count')
    if type(min_clusters) is not int or min_clusters < 2:
        raise InputError('Bootstrap requires at least two independent clusters')
    paired = sorted(set(left) & set(right))
    group_deltas = {}
    for qid in paired:
        a, b = left[qid], right[qid]
        if (isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int,float)) or not isinstance(b, (int,float))
                or not math.isfinite(a) or not math.isfinite(b) or not 0 <= a <= 1 or not 0 <= b <= 1):
            raise InputError('Bootstrap values must be finite per-query proportions in 0..1')
        group = clusters.get(qid)
        if not isinstance(group, str) or not group:
            raise InputError('Every paired query needs a documented independence group')
        group_deltas.setdefault(group, []).append(b-a)
    means = [statistics.mean(group_deltas[g]) for g in sorted(group_deltas)]
    result = {'status': 'estimated' if len(means) >= min_clusters else 'insufficient_independent_clusters',
        'paired_query_ids': paired, 'excluded_unpaired_query_ids': sorted(set(left) ^ set(right)),
        'independent_clusters': len(means), 'min_independent_clusters': min_clusters,
        'paired_cluster_mean_difference': statistics.mean(means) if means else None,
        'confidence_level': 0.95, 'ci_lower': None, 'ci_upper': None,
        'seed': seed, 'bootstrap_replicates': replicates, 'direction': 'right minus left',
        'unit': 'Equal-weight independent-group mean of paired query proportions'}
    if len(means) < min_clusters:
        return result
    rng = random.Random(seed)
    estimates = sorted(statistics.mean(rng.choices(means, k=len(means))) for _ in range(replicates))
    def percentile(p):
        position = (len(estimates)-1)*p
        low = math.floor(position)
        high = math.ceil(position)
        return estimates[low] + (estimates[high]-estimates[low])*(position-low)
    result.update(ci_lower=percentile(0.025), ci_upper=percentile(0.975))
    return result
