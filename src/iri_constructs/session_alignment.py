"""Exhaustive local-alignment comparator for the benign IRI digital benchmark.

Prediction receives only the query and catalogue; construction labels/boundaries
are deliberately absent from this API. Coordinates are zero-based, half-open.
This is retrieval, not a calibrated biological-function or hazard classifier.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import hashlib
import math
from typing import Any

CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class AlignmentConfig:
    matrix: str = "BLOSUM62"
    gap_open: float = 10.0
    gap_extend: float = 0.5
    epsilon: float = 1e-6

    def __post_init__(self):
        if self.matrix != "BLOSUM62":
            raise ValueError("The frozen alignment arm requires BLOSUM62")
        if (self.gap_open, self.gap_extend) != (10.0, 0.5):
            raise ValueError("The frozen alignment arm requires open=10, extend=0.5")
        if self.epsilon != 1e-6:
            raise ValueError("The frozen alignment arm requires epsilon=1e-6")


def _validate_sequence(sequence: str) -> None:
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("Alignment requires a nonempty sequence string")
    if not set(sequence) <= CANONICAL_AA:
        raise ValueError("Alignment input contains noncanonical residue symbols")


def make_aligner(config: AlignmentConfig | None = None):
    from Bio.Align import PairwiseAligner, substitution_matrices
    config = config or AlignmentConfig()
    aligner = PairwiseAligner(mode="local")
    aligner.substitution_matrix = substitution_matrices.load(config.matrix)
    aligner.open_gap_score = -config.gap_open
    aligner.extend_gap_score = -config.gap_extend
    aligner.epsilon = config.epsilon
    return aligner


def alignment_metadata(config: AlignmentConfig | None = None) -> dict[str, Any]:
    import Bio
    config = config or AlignmentConfig()
    aligner = make_aligner(config)
    matrix = aligner.substitution_matrix
    return {
        "implementation": "Bio.Align.PairwiseAligner",
        "implementation_version": Bio.__version__,
        "algorithm": aligner.algorithm,
        "mode": "local",
        **asdict(config),
        "matrix_sha256": hashlib.sha256(matrix.tobytes()).hexdigest(),
        "matrix_alphabet": matrix.alphabet,
        "gap_cost": "gap_open + (gap_length - 1) * gap_extend",
        "ranking": "descending raw local score, accession ascending for exact ties",
        "coordinate_convention": "0-based half-open; query is alignment row 0",
        "alignment_tie_coordinates": "first alignment yielded by pinned Biopython; alternative optimal coordinates are not enumerated",
        "coverage_convention": "consumed residues in local span / full sequence length; includes residues opposite gaps",
        "identity_convention": "identical residue pairs / all alignment columns including gaps",
        "zero_score": "no positive local alignment; null spans and zero coverage",
        "normalization": "none; raw scores favour longer high-scoring matches",
    }


def align_pair(query: str, reference: str, *, aligner=None) -> dict[str, Any]:
    """Return the highest raw local score and one optimal alignment's geometry."""
    _validate_sequence(query)
    _validate_sequence(reference)
    aligner = aligner or make_aligner()
    alignments = aligner.align(query, reference)
    score = float(alignments.score)
    if not math.isfinite(score):
        raise ValueError("Nonfinite local alignment score")
    result = {
        "score": score,
        "query_length": len(query),
        "reference_length": len(reference),
        "query_span": None,
        "reference_span": None,
        "coordinates": [[], []],
        "aligned_length": 0,
        "aligned_residue_pairs": 0,
        "identities": 0,
        "identity_fraction": 0.0,
        "query_coverage": 0.0,
        "reference_coverage": 0.0,
        "query_paired_coverage": 0.0,
        "reference_paired_coverage": 0.0,
        "gap_columns": 0,
        "gap_runs": 0,
    }
    if score == 0.0:
        return result
    alignment = next(iter(alignments))
    coords = alignment.coordinates.tolist()
    columns = pairs = identities = gap_columns = gap_runs = 0
    reconstructed_score = 0.0
    matrix = aligner.substitution_matrix
    for j in range(len(coords[0]) - 1):
        q0, q1 = coords[0][j:j + 2]
        r0, r1 = coords[1][j:j + 2]
        dq, dr = q1 - q0, r1 - r0
        if min(dq, dr) < 0 or max(dq, dr) == 0:
            raise ValueError("Unexpected alignment coordinate block")
        columns += max(dq, dr)
        if dq and dr:
            if dq != dr:
                raise ValueError("Unequal ungapped block lengths")
            pairs += dq
            for a, b in zip(query[q0:q1], reference[r0:r1]):
                identities += int(a == b)
                reconstructed_score += float(matrix[a, b])
        else:
            gap_columns += max(dq, dr)
            gap_runs += 1
            reconstructed_score += float(aligner.open_gap_score) + (max(dq, dr) - 1) * float(aligner.extend_gap_score)
    if abs(reconstructed_score - score) > 1e-6:
        raise ValueError("Alignment coordinates fail independent score reconstruction")
    qspan, rspan = [coords[0][0], coords[0][-1]], [coords[1][0], coords[1][-1]]
    result.update({
        "query_span": qspan,
        "reference_span": rspan,
        "coordinates": coords,
        "aligned_length": columns,
        "aligned_residue_pairs": pairs,
        "identities": identities,
        "identity_fraction": identities / columns,
        "query_coverage": (qspan[1] - qspan[0]) / len(query),
        "reference_coverage": (rspan[1] - rspan[0]) / len(reference),
        "query_paired_coverage": pairs / len(query),
        "reference_paired_coverage": pairs / len(reference),
        "gap_columns": gap_columns,
        "gap_runs": gap_runs,
        "reconstructed_score": reconstructed_score,
    })
    return result


def _catalogue_sequences(catalogue):
    if not isinstance(catalogue, Mapping) or not catalogue:
        raise ValueError("Catalogue must be a nonempty accession mapping")
    for accession, item in catalogue.items():
        if not isinstance(accession, str) or not accession:
            raise ValueError("Catalogue accession must be a nonempty string")
        sequence = item if isinstance(item, str) else item["sequence"]
        yield accession, sequence


def rank_alignment(sequence: str, catalogue: Mapping, *, config: AlignmentConfig | None = None) -> dict[str, Any]:
    """Exhaustively compare one unsplit query with every supplied reference."""
    _validate_sequence(sequence)
    aligner = make_aligner(config)
    rows = [{"accession": accession, **align_pair(sequence, reference, aligner=aligner)}
            for accession, reference in _catalogue_sequences(catalogue)]
    rows.sort(key=lambda row: (-row["score"], row["accession"]))
    # Keep exact ties distinct from near ties; neither changes primary ordering.
    for rank, row in enumerate(rows, 1):
        exact = [i for i, other in enumerate(rows, 1) if other["score"] == row["score"]]
        row.update({
            "rank": rank,
            "exact_score_tie_count": len(exact),
            "tie_best_rank": min(exact),
            "tie_worst_rank": max(exact),
            "near_tie_count_1e6": sum(abs(other["score"] - row["score"]) <= 1e-6 for other in rows),
        })
    return {"ranking": rows, "method": alignment_metadata(config)}


def independent_affine_score(query: str, reference: str) -> float:
    """Slow independent Python DP for audit samples, not benchmark execution.

    Uses separate local-score and horizontal/vertical affine-gap recurrences.
    Space is O(reference length); no Biopython alignment routine is called.
    """
    from Bio.Align import substitution_matrices
    _validate_sequence(query)
    _validate_sequence(reference)
    matrix = substitution_matrices.load("BLOSUM62")
    width = len(reference) + 1
    previous = [0.0] * width
    vertical = [float("-inf")] * width
    best = 0.0
    for qa in query:
        current = [0.0] * width
        horizontal = float("-inf")
        for j, rb in enumerate(reference, 1):
            vertical[j] = max(previous[j] - 10.0, vertical[j] - 0.5)
            horizontal = max(current[j - 1] - 10.0, horizontal - 0.5)
            current[j] = max(0.0, previous[j - 1] + float(matrix[qa, rb]), vertical[j], horizontal)
            best = max(best, current[j])
        previous = current
    return best
