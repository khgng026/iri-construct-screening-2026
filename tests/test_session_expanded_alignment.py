from copy import deepcopy
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_expanded_alignment import rerank, representative_query, validate_base_reproduction, validate_raw_rows


class ExpandedAlignmentTests(unittest.TestCase):
    def test_reduced_baseline_recomputes_ties_and_preserves_geometry(self):
        rows = [{'accession': a, 'score': score, 'query_span': [0, 2], 'reference_span': [0, 2]}
                for a, score in [('A', 10.), ('A_alias', 10.), ('B', 8.)]]
        source = rerank(rows, ['A', 'A_alias', 'B'])
        reduced = rerank([r for r in rows if r['accession'] != 'A_alias'], ['A', 'B'])
        validate_base_reproduction(source, reduced, ['A', 'B'])
        self.assertEqual(reduced[0]['exact_score_tie_count'], 1)
        altered = deepcopy(reduced); altered[0]['query_span'] = [1, 3]
        with self.assertRaisesRegex(ValueError, 'geometry'):
            validate_base_reproduction(source, altered, ['A', 'B'])
        with self.assertRaisesRegex(ValueError, 'repeated'):
            rerank([rows[0], rows[0]], ['A', 'B'])

    def test_representative_mapping_checks_exact_component_string(self):
        query = {'query_id': 'toy', 'sequence': '012345', 'components': ['alias', 'B'], 'component_count': 2,
                 'spans': [{'accession': 'alias', 'start0': 0, 'end0': 3}, {'accession': 'B', 'start0': 3, 'end0': 6}]}
        mapped = representative_query(query, {'alias': 'A', 'B': 'B'}, {'A': '012', 'B': '345'})
        self.assertEqual(mapped['components'], ['A', 'B'])
        self.assertNotIn('sequence', mapped)
        self.assertEqual(query['components'], ['alias', 'B'])
        with self.assertRaisesRegex(ValueError, 'sequence differs'):
            representative_query(query, {'alias': 'A', 'B': 'B'}, {'A': '021', 'B': '345'})
        with self.assertRaisesRegex(ValueError, 'collapse'):
            representative_query(query, {'alias': 'A', 'B': 'A'}, {'A': '012'})

    def test_cached_raw_score_is_bounded_and_geometry_in_range(self):
        row = {'accession': 'A', 'score': 8., 'query_length': 4, 'reference_length': 2,
               'identity_fraction': 1., 'query_coverage': .5, 'reference_coverage': 1.,
               'query_span': [1, 3], 'reference_span': [0, 2]}
        valid = rerank([row], ['A'])
        validate_raw_rows(valid, {'A': '01'}, 4, {'A': 8.})
        broken = deepcopy(valid); broken[0]['score'] = 8.5
        with self.assertRaisesRegex(ValueError, 'upper bound'):
            validate_raw_rows(broken, {'A': '01'}, 4, {'A': 8.})
        broken = deepcopy(valid); broken[0]['query_span'] = [1, 5]
        with self.assertRaisesRegex(ValueError, 'outside'):
            validate_raw_rows(broken, {'A': '01'}, 4, {'A': 8.})


if __name__ == '__main__': unittest.main()
