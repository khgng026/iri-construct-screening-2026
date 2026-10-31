"""Catalogue expansion must preserve saved baseline geometry and precision."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_expanded_model import validate_vector, validate_base_reproduction
from iri_constructs.session_store import digest
from iri_constructs.session_benchmark import rank_scores


class ExpandedModelTests(unittest.TestCase):
    def test_vector_requires_exact_source_hash_length_precision_dimension(self):
        sequence = '0123'  # Nonbiological fixture.
        record = {'sequence_sha256': digest(sequence.encode()), 'length': len(sequence)}
        arrays = {'mean': np.array([1., 2., 3.], dtype=np.float32)}
        self.assertEqual(validate_vector(record, arrays, sequence, 3).shape, (3,))
        with self.assertRaises(ValueError): validate_vector(record, arrays, '0132', 3)
        with self.assertRaises(ValueError): validate_vector(record, arrays, sequence, 4)
        with self.assertRaises(ValueError): validate_vector(record, {'mean': arrays['mean'].astype(np.float64)}, sequence, 3)
        with self.assertRaises(ValueError): validate_vector(record, {'mean': np.array([1., np.nan, 3.], dtype=np.float32)}, sequence, 3)

    def test_source_baseline_rerank_preserves_ties_after_alias_removal(self):
        parent = rank_scores([1., 1., .8, .7], ['A', 'A_alias', 'B', 'C'], 'cosine_similarity')
        expanded = rank_scores([1., .8, .7], ['A', 'B', 'C'], 'cosine_similarity')
        self.assertEqual(validate_base_reproduction(parent, expanded, ['A', 'B', 'C']), 0.)
        broken = deepcopy(expanded); broken[0]['score'] += 1e-5
        with self.assertRaises(ValueError): validate_base_reproduction(parent, broken, ['A', 'B', 'C'])
        with self.assertRaises(ValueError): validate_base_reproduction(parent, expanded[::-1], ['A', 'B', 'C'])


if __name__ == '__main__':
    unittest.main()
