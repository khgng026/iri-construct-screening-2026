"""Numerical fixtures test logic only; no biological or model measurements."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.fusion import compare_case, validate_case
from iri_constructs.metrics import cosine, mean_pool
from iri_constructs.runtime import InputError


class NumericalAdapter:
    device = 'cpu'
    def __init__(self):
        self.calls = []
    def residues(self, sequence):
        self.calls.append(len(sequence))
        return [[float(i + 1), float(len(sequence))] for i in range(len(sequence))]


class FusionTests(unittest.TestCase):
    def setUp(self):
        self.proteins = {'f': ProteinInput('f', 'parent-f', 'A' * 6, 'fake'),
                         't': ProteinInput('t', 'parent-t', 'A' * 2, 'fake')}
        self.case = dict(case_id='numerical-test', fusion_id='f', target_id='t', start0=2, end0=4,
                         coordinate_evidence='Fictional test indices only', review_status='approved')

    def test_three_conditions_use_two_forward_passes(self):
        adapter = NumericalAdapter()
        result = compare_case(self.case, self.proteins, adapter)
        self.assertEqual(adapter.calls, [2, 6])
        self.assertEqual(result['status'], 'ok')
        self.assertAlmostEqual(result['pairwise_cosine']['isolated_target__contextual_target'],
                               cosine([1.5, 2.0], [3.5, 6.0]))
        self.assertGreater(result['context_difference']['residue_representation_rms'], 0)

    def test_bad_coordinates_and_review_fail_before_inference(self):
        for change in ({'start0': True}, {'start0': -1}, {'end0': 7}, {'end0': 2},
                       {'review_status': 'pending'}, {'coordinate_evidence': ''}):
            with self.subTest(change=change), self.assertRaises(InputError):
                validate_case({**self.case, **change}, self.proteins)

    def test_exact_target_required(self):
        self.proteins['t'] = ProteinInput('t', 'p', 'C' * 2, 'fake')
        with self.assertRaises(InputError):
            validate_case(self.case, self.proteins)

    def test_full_fusion_failure_preserves_isolated_coverage(self):
        class Limited(NumericalAdapter):
            def residues(self, sequence):
                if len(sequence) > 3:
                    raise InputError('Test residue cap')
                return super().residues(sequence)
        result = compare_case(self.case, self.proteins, Limited())
        self.assertEqual(result['status'], 'partial')
        self.assertEqual(result['conditions']['isolated_target']['status'], 'ok')
        self.assertEqual(result['conditions']['contextual_target']['status'], 'unavailable')
        self.assertEqual(result['pairwise_cosine'], {})

    def test_invalid_vectors(self):
        for a,b in [([0.0], [1.0]), ([1.0], [1.0, 2.0]), ([float('nan')], [1.0])]:
            with self.assertRaises(InputError):
                cosine(a,b)
        with self.assertRaises(InputError):
            mean_pool([])
