"""T05 ESM-2 isolated/full-fusion/contextual-target comparison."""
from __future__ import annotations

import math
from .metrics import cosine, mean_pool
from .runtime import InputError, Measure, failure, require_text


def validate_case(case: dict, proteins: dict) -> tuple:
    if not isinstance(case, dict):
        raise InputError('Fusion case must be an object')
    require_text(case, 'case_id')
    require_text(case, 'coordinate_evidence')
    if case.get('review_status') != 'approved':
        raise InputError('Fusion coordinates have not been reviewed')
    try:
        fusion = proteins[case['fusion_id']]
        target = proteins[case['target_id']]
    except (KeyError, TypeError) as exc:
        raise InputError('Fusion case references a missing validated protein input') from exc
    start, end = case.get('start0'), case.get('end0')
    if type(start) is not int or type(end) is not int or not 0 <= start < end <= len(fusion.sequence):
        raise InputError('Fusion target coordinates must be 0-based half-open amino-acid bounds')
    if fusion.sequence[start:end] != target.sequence:
        raise InputError('Reviewed target span does not exactly match the isolated target sequence')
    return fusion, target, start, end


def compare_case(case: dict, proteins: dict, adapter, *, pooled_sink: dict | None = None) -> dict:
    fusion, target, start, end = validate_case(case, proteins)
    result = {'case_id': case['case_id'], 'fusion_id': fusion.input_id, 'target_id': target.input_id,
              'start0': start, 'end0': end, 'coordinate_system': '0-based half-open amino acids',
              'coordinate_evidence': case['coordinate_evidence'], 'conditions': {}}
    raw = {}
    # One shared model instance, no concurrent inference. Context is sliced after full fusion inference.
    for condition, sequence in [('isolated_target', target.sequence), ('complete_fusion', fusion.sequence)]:
        with Measure(adapter.device) as measured:
            try:
                raw[condition] = adapter.residues(sequence)
                result['conditions'][condition] = {'status': 'ok', 'length_aa': len(sequence)}
            except Exception as exc:
                result['conditions'][condition] = {'status': 'failed', **failure(exc)}
        result['conditions'][condition]['resources'] = measured.data
    if 'complete_fusion' in raw:
        raw['contextual_target'] = raw['complete_fusion'][start:end]
        result['conditions']['contextual_target'] = {'status': 'ok', 'length_aa': end-start,
            'inference_reused_from': 'complete_fusion'}
    else:
        result['conditions']['contextual_target'] = {'status': 'unavailable', 'reason': 'Full fusion inference failed'}
    pooled = {key: mean_pool(rows) for key, rows in raw.items()}
    # Callers may reuse pooled values for a separately frozen retrieval panel.
    # Vectors remain in memory and are never inserted in the sequence-free report.
    if pooled_sink is not None:
        pooled_sink.update(pooled)
    result['pairwise_cosine'] = {a + '__' + b: cosine(pooled[a], pooled[b])
        for a, b in [('isolated_target', 'complete_fusion'), ('isolated_target', 'contextual_target'),
                     ('complete_fusion', 'contextual_target')] if a in pooled and b in pooled}
    if 'isolated_target' in raw and 'contextual_target' in raw:
        isolated, contextual = raw['isolated_target'], raw['contextual_target']
        if len(isolated) != len(contextual) or any(len(a) != len(b) for a, b in zip(isolated, contextual)):
            raise InputError('Isolated and contextual target dimensions differ')
        squared = math.fsum((x-y)**2 for a,b in zip(isolated, contextual) for x,y in zip(a,b))
        result['context_difference'] = {
            'residue_representation_rms': math.sqrt(squared / sum(map(len, isolated))),
            'mean_residue_cosine': math.fsum(cosine(a,b) for a,b in zip(isolated, contextual))/len(isolated),
            'pooled_cosine_distance': 1-cosine(pooled['isolated_target'], pooled['contextual_target'])}
    result['status'] = 'ok' if len(raw) == 3 else 'partial'
    return result
