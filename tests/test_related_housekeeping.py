"""Check that local embedding alignment responds to residue order."""
import importlib.util
import unittest

from experiments.benign_protein_pilot.related_study import LOCAL_EXTEND, LOCAL_OPEN, _residue_local_score


@unittest.skipUnless(importlib.util.find_spec('Bio') and importlib.util.find_spec('numpy'), 'requires parser dependencies')
class ResidueLocalScoreTests(unittest.TestCase):
    def test_distinguishes_vectors_with_identical_means_and_different_order(self):
        import numpy as np
        from Bio.Align import PairwiseAligner

        query = np.eye(6, dtype=np.float32)
        reversed_order = query[::-1]
        aligner = PairwiseAligner(mode='local', open_gap_score=LOCAL_OPEN, extend_gap_score=LOCAL_EXTEND)
        np.testing.assert_array_equal(query.mean(axis=0), reversed_order.mean(axis=0))
        aligned = _residue_local_score(query, query, aligner)
        scrambled = _residue_local_score(query, reversed_order, aligner)
        self.assertGreater(aligned, scrambled)
        self.assertAlmostEqual(aligned, 5 / 6)


if __name__ == '__main__':
    unittest.main()
