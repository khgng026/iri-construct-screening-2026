import importlib.util
import unittest

from experiments.paired_preprocessing.distant_reporter_study import (_check_structure_mapping, _domain_locus,
                                                                   _score_targets, _signature)


class DistantReporterStudyTest(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython structure mapping dependency')
    def test_structure_mapping_checks_the_documented_in_frame_variant(self):
        target = {'missing_structure_reference_residues': 1, 'substitutions_to_structure_reference': 1}
        _check_structure_mapping('ACEFGAHIK', 'ACDEFGPHIK', target)
        with self.assertRaises(ValueError):
            _check_structure_mapping('ACEFGAAIK', 'ACDEFGPHIK', target)

    def test_archived_domain_offsets_project_from_the_coding_start(self):
        self.assertEqual(_domain_locus([[100, 130]], 1, 1, 4), [[103, 112]])
        self.assertEqual(_domain_locus([[100, 130]], -1, 1, 4), [[118, 127]])
        with self.assertRaises(ValueError):
            _domain_locus([[100, 107], [110, 130]], 1, 1, 4)

    def test_one_call_cannot_recover_two_targets(self):
        query = {'strand': 1, 'spans': [[0, 60]], 'protein': 'candidate'}
        targets = [dict(query, id=key, family='reporter') for key in ('first', 'second')]
        outcomes = _score_targets(targets, [query], [query], {'method': {_signature('candidate'): 'reporter'}})
        self.assertEqual(sum(case['correct']['method'] for case in outcomes.values()), 1)
        self.assertTrue(all(case['selected_covering'] for case in outcomes.values()))

    def test_coverage_can_be_lost_at_the_query_cap(self):
        query = {'strand': 1, 'spans': [[0, 60]], 'protein': 'candidate'}
        target = dict(query, id='first', family='reporter')
        case, = _score_targets([target], [query], [], {'method': {}}).values()
        self.assertTrue(case['pool_covering'])
        self.assertFalse(case['selected_covering'])
        self.assertFalse(case['selected_exact_protein'])
        self.assertFalse(case['correct']['method'])
