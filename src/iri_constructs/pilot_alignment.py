"""Fixed local protein alignment scores on the same admitted panel as embeddings.

PR #1 motivates the raw/length-normalized comparison, not implementation imports.
No boundary inference, sequence reconstruction or annotation assignment occurs.
"""
from __future__ import annotations

import math
import hashlib

from .bench_inputs import validate_sequence
from .runtime import InputError, Measure


def alignment_scores(proteins: dict, panel: dict, *, max_residues: int = 1022) -> dict:
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner(mode='local')
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.epsilon = 1e-6
    raw, normalized, failures = {}, {}, []
    self_scores, calls, self_calls, work = {}, 0, 0, 0
    with Measure() as measured:
        for qid in panel['query_ids']:
            references = panel['eligible_reference_ids'].get(qid, [])
            if not references:
                continue
            try:
                query = proteins[qid].sequence
                validate_sequence(query, 'protein')
                if len(query) > max_residues:
                    raise InputError('Alignment query exceeds the same experimental cap as embeddings')
                qraw, qnorm = {}, {}
                # The entire frozen candidate panel must be available; never rank a partial panel.
                for rid in references:
                    reference = proteins[rid].sequence
                    validate_sequence(reference, 'protein')
                    if len(reference) > max_residues:
                        raise InputError('Alignment reference exceeds the common experimental cap')
                    score = float(aligner.score(query, reference))
                    calls += 1
                    work += len(query) * len(reference)
                    if rid not in self_scores:
                        self_scores[rid] = float(aligner.score(reference, reference))
                        self_calls += 1
                    if not math.isfinite(score) or self_scores[rid] <= 0:
                        raise InputError('Invalid local alignment score or reference self score')
                    qraw[rid], qnorm[rid] = score, score / self_scores[rid]
                raw[qid], normalized[qid] = qraw, qnorm
            except InputError as exc:
                failures.append({'input_id': qid, 'reason': str(exc)})
    return {
        'status': 'ok' if raw and not failures else 'partial' if raw else 'no_scorable_queries', 'raw_scores': raw,
        'normalized_scores': normalized, 'failures': failures, 'resources': measured.data,
        'configuration': {'algorithm': 'Biopython PairwiseAligner local affine-gap scoring',
                          'selected_algorithm': aligner.algorithm,
                          'substitution_matrix': 'BLOSUM62', 'open_gap_score': -10.0,
                          'extend_gap_score': -0.5, 'max_residues': max_residues,
                          'epsilon': 1e-6, 'gap_cost': 'open + (length - 1) * extend',
                          'matrix_sha256': hashlib.sha256(aligner.substitution_matrix.tobytes()).hexdigest(),
                          'matrix_alphabet': aligner.substitution_matrix.alphabet,
                          'normalization': 'Raw primary; reference-self score normalization predeclared sensitivity'},
        'pairwise_score_calls': calls, 'reference_self_score_calls': self_calls,
        'query_reference_length_products': work,
        'work_note': 'Length products are a work proxy, not measured operation counts or runtime.',
    }
