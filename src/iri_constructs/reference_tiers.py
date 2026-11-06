"""Pre-score independent-reference tiers extending the approved T07 panel guards.

Every tier retains exact-copy, same-parent and reviewed independence-group
exclusions. Additional vector/linkage or stricter relatedness exclusions require
predeclared, human-reviewed evidence. No identity threshold or source approval is
inferred here. Callers must enforce T02 source gates before supplying proteins.
"""
from __future__ import annotations

import copy
from datetime import date
import hashlib
import json

from .pilot_scoring import freeze_panel, rank_scores, score_rankings, validate_panel
from .runtime import InputError, require_text


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode('utf-8')).hexdigest()


def _base_config(config):
    if config.get('schema_version') == 'pilot-scoring-1':
        return config
    from .extension_config import base_scoring_config
    return base_scoring_config(config, experiment='C')


def _reviewed_pairs(policy, scope, proteins, queries, references):
    """Validate a policy declaration, not the scientific truth of its evidence."""
    if policy is None:
        return None, 'No independently reviewed predeclared exclusion policy supplied'
    try:
        if not isinstance(policy, dict) or policy.get('schema_version') != 'reviewed-relatedness-policy-1':
            raise InputError('Expected reviewed-relatedness-policy-1 metadata')
        if policy.get('policy_kind') != 'reviewed_pair_exclusions':
            raise InputError('Unsupported relatedness policy; no automatic identity threshold is admitted')
        if policy.get('policy_scope') != scope:
            raise InputError('Reviewed exclusion policy has the wrong declared scope')
        if policy.get('review_status') != 'approved' or policy.get('declared_before_results') is not True:
            raise InputError('Exclusion policy requires human approval and pre-score declaration')
        for key in ('policy_id', 'reviewer', 'reviewed_on', 'evidence'):
            require_text(policy, key)
        date.fromisoformat(policy['reviewed_on'])
        entries = policy.get('excluded_pairs')
        if not isinstance(entries, list):
            raise InputError('Reviewed exclusions must be an explicit pair list')
        pairs, seen = {}, set()
        for item in entries:
            if not isinstance(item, dict):
                raise InputError('Every reviewed exclusion must be an object')
            qid, rid = require_text(item, 'query_id'), require_text(item, 'reference_id')
            if qid not in queries or rid not in references or (qid,rid) in seen:
                raise InputError('Reviewed exclusion IDs must be unique members of the frozen roles')
            seen.add((qid,rid))
            for key in ('reason', 'evidence'):
                require_text(item, key)
            if (item.get('query_sequence_sha256') != proteins[qid].sequence_sha256 or
                    item.get('reference_sequence_sha256') != proteins[rid].sequence_sha256):
                raise InputError('Reviewed relatedness evidence does not match the frozen sequence fingerprints')
            pairs[(qid,rid)] = {'reason': item['reason'], 'evidence': item['evidence']}
        return {'pairs': pairs, 'policy_sha256': _hash(policy), 'policy_id': policy['policy_id'],
                'policy_scope': scope}, None
    except (InputError, ValueError, TypeError):
        # Keep a bad stricter tier blocked while independent valid tiers can run.
        # Do not include arbitrary values or source sequence text in an error.
        return None, 'Unsupported, incomplete, unapproved or fingerprint-mismatched relatedness policy'


def _base_removals(panel):
    removals = {}
    for qid in panel['query_ids']:
        query = panel['input_fingerprints'][qid]
        removed = {}
        for rid in panel['reference_ids']:
            ref = panel['input_fingerprints'][rid]
            reasons = []
            if qid == rid:
                reasons.append('query_self_reference')
            if query['sequence_sha256'] == ref['sequence_sha256']:
                reasons.append('exact_sequence_copy')
            if query['construct_id'] == ref['construct_id']:
                reasons.append('same_source_parent')
            if query['independence_group'] == ref['independence_group']:
                reasons.append('reviewed_independence_group')
            if reasons:
                removed[rid] = [{'rule': r, 'evidence': 'Frozen T07 source/label provenance'} for r in reasons]
        removals[qid] = removed
    return removals


def _make_tier(base, name, removals, policy_status, policy_sha256=None):
    panel = copy.deepcopy(base)
    allowed, exclusions, counts = {}, {}, {}
    for qid in base['query_ids']:
        expected = set(base['expected_families_by_query'][qid])
        allowed[qid] = [r for r in base['reference_ids'] if r not in removals[qid]]
        positive = [r for r in allowed[qid] if base['input_fingerprints'][r]['family'] in expected]
        competing = [r for r in allowed[qid] if base['input_fingerprints'][r]['family'] not in expected]
        observed = {base['input_fingerprints'][r]['family'] for r in positive}
        if not expected <= observed:
            exclusions[qid] = 'No independent reference for every frozen documented family at this tier'
        elif not competing:
            exclusions[qid] = 'No independent competing-family reference at this tier'
        counts[qid] = {'catalog_references': len(base['reference_ids']), 'eligible_references': len(allowed[qid]),
            'removed_references': len(removals[qid]), 'positive_references': len(positive),
            'competing_references': len(competing),
            'expected_documented_families': len(expected), 'supported_documented_families': len(observed),
            'eligible_reference_parents': len({base['input_fingerprints'][r]['construct_id'] for r in allowed[qid]}),
            'eligible_reference_independence_groups': len({base['input_fingerprints'][r]['independence_group'] for r in allowed[qid]})}
    panel.update(eligible_reference_ids=allowed, query_exclusions=exclusions,
        status='frozen' if any(q not in exclusions for q in panel['query_ids']) else 'blocked_insufficient_panel',
        reference_tier=name, relatedness_policy_sha256=policy_sha256,
        budget_saturated_query_ids=[q for q in panel['query_ids'] if q not in exclusions and len(allowed[q]) <= panel['top_k']])
    panel.pop('panel_sha256', None)
    panel['panel_sha256'] = _hash(panel)
    return {'tier': name, 'status': panel['status'], 'panel': panel, 'removals': removals,
        'reference_counts_by_query': counts, 'additional_policy_status': policy_status,
        'relatedness_policy_sha256': policy_sha256, 'requested_queries': len(panel['query_ids']),
        'eligible_queries': len(panel['query_ids'])-len(exclusions),
        'eligible_query_ids': [q for q in panel['query_ids'] if q not in exclusions],
        'query_parent_count': len({panel['input_fingerprints'][q]['construct_id'] for q in panel['query_ids']}),
        'query_independence_group_count': len({panel['input_fingerprints'][q]['independence_group'] for q in panel['query_ids']})}


def _add_pairs(removals, policy, rule):
    result = copy.deepcopy(removals)
    for (qid,rid), evidence in sorted(policy['pairs'].items()):
        result[qid].setdefault(rid, []).append({'rule': rule, 'policy_id': policy['policy_id'], **evidence})
    return result


def freeze_reference_tiers(config, proteins, labels, query_ids=None, reference_ids=None,
                           *, linkage_policy=None, relatedness_policy=None, component_families=None):
    """Freeze all tier identities/rules/removals before observing method results.

    Tier 1 is the full original T07 guard, not a relaxation that permits known
    same-parent/linked references. Tier 2 adds optional reviewed vector/linkage
    pair exclusions, otherwise transparently repeats tier 1. Tier 3 additionally
    needs a reviewed sequence-relatedness pair policy; absent evidence blocks it.
    """
    base = freeze_panel(_base_config(config), proteins, labels, query_ids, reference_ids, component_families)
    queries, references = base['query_ids'], base['reference_ids']
    first_removals = _base_removals(base)
    tier1 = _make_tier(base, 'tier1_t07_invariants', first_removals, 'Original T07 guards mandatory')
    linkage, linkage_error = _reviewed_pairs(linkage_policy, 'vector_linkage', proteins, queries, references)
    if linkage_policy is None:
        second_removals, second_status, second_hash = first_removals, 'No additional policy; redundant with tier 1', None
    elif linkage is None:
        second_removals, second_status, second_hash = None, linkage_error, None
    else:
        second_removals = _add_pairs(first_removals, linkage, 'reviewed_vector_linkage')
        second_status, second_hash = 'Reviewed additional vector/linkage exclusions', linkage['policy_sha256']
    if second_removals is None:
        tier2 = {'tier': 'tier2_additional_linkage', 'status': 'blocked_relatedness_policy', 'panel': None,
                 'reason': second_status}
    else:
        tier2 = _make_tier(base, 'tier2_additional_linkage', second_removals, second_status, second_hash)
    strict, strict_error = _reviewed_pairs(relatedness_policy, 'sequence_relatedness', proteins, queries, references)
    if second_removals is None or strict is None:
        tier3 = {'tier': 'tier3_reviewed_relatedness', 'status': 'blocked_relatedness_policy', 'panel': None,
                 'reason': 'Prior linkage policy is invalid' if second_removals is None else strict_error}
    else:
        strict_removals = _add_pairs(second_removals, strict, 'reviewed_stricter_relatedness')
        tier3 = _make_tier(base, 'tier3_reviewed_relatedness', strict_removals,
                          'Reviewed stricter pair-relatedness exclusions', strict['policy_sha256'])
    tiers = {t['tier']: t for t in (tier1,tier2,tier3)}
    previous = set(tier1['eligible_query_ids'])
    for name, tier in tiers.items():
        current = set(tier.get('eligible_query_ids', []))
        tier['query_attrition_from_prior_tier'] = sorted(previous-current)
        tier['query_attrition_from_tier1'] = sorted(set(tier1['eligible_query_ids'])-current)
        if tier.get('panel') is not None:
            previous = current
    result = {'schema_version': 'reference-tiers-1', 'stage': 'before all method results',
        'source_admission': 'Caller must enforce unchanged T02 gates; policies never approve sources',
        'base_protocol_sha256': base['protocol_sha256'], 'query_ids': queries, 'reference_ids': references,
        'mandatory_rules': ['query_self_reference','exact_sequence_copy','same_source_parent','reviewed_independence_group'],
        'tiers': tiers, 'limitations': ['Tier/candidate changes affect catalogue composition and score denominators.',
            'Reviewed relatedness does not establish unseen-protein or pretraining independence.',
            'No scored inference is authorized by a metadata-only tier freeze.']}
    result['tier_freeze_sha256'] = _hash(result)
    return result


def assert_same_tier_panel(method_panels):
    """Verify exact method query/reference/eligibility and budget parity."""
    hashes = []
    for panel in method_panels.values():
        validate_panel(panel)
        hashes.append(panel['panel_sha256'])
    if not hashes or len(set(hashes)) != 1:
        raise InputError('Methods do not share the exact frozen tier panel and return budget')
    return hashes[0]


def score_reference_tier(scores_by_query, tier_or_panel, labels):
    """Use the existing frozen scorer and expose bounded wrong-family hit rates."""
    panel = tier_or_panel.get('panel') if 'panel' in tier_or_panel else tier_or_panel
    if panel is None:
        return {'status': 'blocked_relatedness_policy', 'scored_queries': 0,
                'false_family_hit_fraction': None, 'queries_with_false_family_match_rate': None,
                'false_family_hit_count': 0, 'returned_hit_denominator': 0}
    result = score_rankings(rank_scores(scores_by_query, panel), panel, labels)
    false_hits = false_queries = returned = 0
    for row in result['rows']:
        expected = set(panel['expected_families_by_query'][row['query_id']])
        wrong = sum(labels[r]['family'] not in expected for r in row['reference_ids'])
        row.update(false_family_hit_count=wrong, returned_hit_denominator=len(row['reference_ids']))
        false_hits += wrong
        false_queries += bool(wrong)
        returned += len(row['reference_ids'])
    result.update(false_family_hit_count=false_hits, returned_hit_denominator=returned,
        false_family_hit_fraction=false_hits/returned if returned else None,
        queries_with_false_family_match_rate=false_queries/len(result['rows']) if result['rows'] else None,
        false_family_claim_scope='Wrong reviewed-family labels among bounded returned references, not a biological or operational screening false-positive rate')
    return result
