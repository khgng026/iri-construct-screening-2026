"""Sequence alphabets shared by preprocessing and the historical experiment code.

Kept free of manifest, approval and model imports so the installed toolkit can
validate a local record without any research-admission machinery.
"""
from __future__ import annotations

from .runtime import InputError

AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
DNA = frozenset('ACGTRYSWKMBDHVN')


def validate_sequence(sequence: str, molecule: str) -> None:
    alphabet = AA if molecule == 'protein' else DNA if molecule == 'dna' else None
    if alphabet is None or not sequence or set(sequence) - alphabet:
        raise InputError('Empty sequence or unsupported alphabet; protein inference requires canonical amino acids')
