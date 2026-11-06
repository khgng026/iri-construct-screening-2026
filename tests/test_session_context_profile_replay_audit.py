from pathlib import Path
import math
import sys
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_context_profile_replay_audit import independent_residue_changes, selected_queries, guard_admission


class IndependentReplayArithmeticTests(unittest.TestCase):
    def test_independent_residue_dot_norm_and_change(self):
        left = np.array([[1., 0.], [1., 1.], [2., 0.]], dtype=np.float32)
        right = np.array([[1., 0.], [0., 1.], [1., 0.]], dtype=np.float32)
        cosine, relative = independent_residue_changes(left, right)
        np.testing.assert_allclose(cosine, [1., 1 / math.sqrt(2), 1.], rtol=0, atol=1e-15)
        np.testing.assert_allclose(relative, [0., 1., 1.], rtol=0, atol=1e-15)
        with self.assertRaisesRegex(ValueError, 'zero-norm'):
            independent_residue_changes(left, np.zeros_like(right))

    def test_three_query_selection_is_fixed_by_lexicographic_id_within_count(self):
        queries = [{'query_id': name, 'component_count': count}
                   for name, count in [('z2', 2), ('z4', 4), ('a3', 3), ('a4', 4), ('a2', 2), ('z3', 3)]]
        self.assertEqual([q['query_id'] for q in selected_queries(queries)], ['a2', 'a3', 'a4'])
        with self.assertRaises(ValueError):
            selected_queries([q for q in queries if q['component_count'] != 4])

    def test_fresh_mode_refuses_missing_recorded_guard(self):
        with patch.dict('os.environ', {}, clear=True):
            with self.assertRaisesRegex(ValueError, 'prepared guard'):
                guard_admission()


if __name__ == '__main__': unittest.main()
