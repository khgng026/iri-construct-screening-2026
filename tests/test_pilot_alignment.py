"""Fictional generated protein strings verify alignment contracts, not biological findings."""
import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.bench_inputs import ProteinInput


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython alignment tests')
class PilotAlignmentTests(unittest.TestCase):
    def setUp(self):
        self.proteins = {key: ProteinInput(key,key,sequence,'fictional-test') for key,sequence in
                         [('q','A'*5),('same','A'*4+'C'),('other','W'*5)]}
        self.panel = {'query_ids':['q'], 'eligible_reference_ids':{'q':['same','other']}}

    def test_raw_and_self_normalized_use_identical_panel(self):
        from iri_constructs.pilot_alignment import alignment_scores
        result = alignment_scores(self.proteins,self.panel)
        self.assertEqual(result['status'],'ok')
        self.assertGreater(result['raw_scores']['q']['same'],result['raw_scores']['q']['other'])
        self.assertEqual(set(result['raw_scores']['q']),set(result['normalized_scores']['q']))
        self.assertEqual(result['pairwise_score_calls'],2)
        self.assertEqual(result['reference_self_score_calls'],2)

    def test_cap_failure_abstains_on_entire_query_panel(self):
        from iri_constructs.pilot_alignment import alignment_scores
        result = alignment_scores(self.proteins,self.panel,max_residues=3)
        self.assertEqual(result['raw_scores'],{})
        self.assertEqual(result['status'],'no_scorable_queries')
        self.assertEqual(len(result['failures']),1)

    def test_empty_panel_has_no_scores_or_calls(self):
        from iri_constructs.pilot_alignment import alignment_scores
        result = alignment_scores(self.proteins,{'query_ids':[],'eligible_reference_ids':{}})
        self.assertEqual(result['pairwise_score_calls'],0)
        self.assertEqual(result['raw_scores'],{})
