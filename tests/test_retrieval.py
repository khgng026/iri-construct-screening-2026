"""Hand-computable vectors exercise retrieval; no pretrained measurements."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.metrics import candidates, retrieve, reviewed_labels
from iri_constructs.runtime import InputError


class RetrievalTests(unittest.TestCase):
    def setUp(self):
        self.proteins = {k: ProteinInput(k, k, residue * 3, 'fictional')
                         for k,residue in zip('abcd', 'ACDE')}
        self.labels = {k: dict(input_id=k, family='one' if k in 'ab' else 'two',
                              evidence='Fictional test', independence_group=k,
                              annotation_status='known', review_status='approved') for k in 'abcd'}
        self.vectors = dict(a=[1.,0.], b=[1.,.1], c=[0.,1.], d=[.1,1.])

    def test_perfect_distinct_independent_family_retrieval(self):
        result = retrieve(self.vectors, self.labels, self.proteins)
        self.assertEqual(result['scored'], 4)
        self.assertEqual(result['macro_f1'], 1.0)
        self.assertTrue(all(p['margin_same_minus_other'] > 0 for p in result['predictions']))
        self.assertTrue(all(p['input_id'] != p['neighbor'] for p in result['predictions']))

    def test_no_independent_positive_abstains(self):
        self.labels['b']['independence_group'] = 'a'
        result = retrieve(self.vectors, self.labels, self.proteins)
        self.assertEqual(result['scored'], 2)
        self.assertEqual({a['input_id'] for a in result['abstentions']}, {'a','b'})

    def test_exact_duplicate_and_same_construct_excluded(self):
        self.proteins['b'] = ProteinInput('b', 'b', self.proteins['a'].sequence, 'fake')
        self.assertNotIn('b', candidates('a', self.labels, self.proteins))
        self.proteins['c'] = ProteinInput('c', 'a', 'D' * 3, 'fake')
        self.assertNotIn('c', candidates('a', self.labels, self.proteins))

    def test_failed_reference_does_not_shrink_fixed_panel(self):
        del self.vectors['d']
        result = retrieve(self.vectors, self.labels, self.proteins)
        self.assertEqual(result['scored'], 0)
        self.assertEqual(result['requested'], 4)
        self.assertIsNone(result['macro_f1'])

    def test_explicit_common_panel_ids(self):
        result = retrieve(self.vectors, self.labels, self.proteins, ['a','c'], ['a','b','c'])
        self.assertEqual(result['query_ids'], ['a','c'])
        self.assertEqual(result['reference_ids'], ['a','b','c'])
        self.assertEqual(result['scored'], 1)

    def test_unreviewed_ambiguous_and_duplicate_labels(self):
        entry = {**self.labels['a'], 'annotation_status': 'ambiguous'}
        labels, excluded = reviewed_labels({'families': [entry]}, self.proteins)
        self.assertEqual(labels, {})
        self.assertEqual(len(excluded), 1)
        with self.assertRaises(InputError):
            reviewed_labels({'families': [self.labels['a'], self.labels['a']]}, self.proteins)
