"""Regression tests for the Q2 reporter-replication curation and analysis helpers.

These cover the dependency-free logic: family mapping from deposited product strings,
and the pooled / equal-group / bootstrap math used to interpret the paired runs. Network
fetches, toolkit validation, and Biopython clustering are exercised by the committed
run artifacts, not here.
"""
import unittest

from experiments.q2_reporter_replication import analyze, curate


class FamilyMappingTests(unittest.TestCase):
    def test_reporter_products_map_to_expected_families(self):
        cases = {
            'green fluorescent protein GFP': ('gfp_like', 'reporter'),
            'enhanced GFP': ('gfp_like', 'reporter'),
            'mCherry': ('gfp_like', 'reporter'),
            'firefly luciferase': ('firefly_luciferase', 'reporter'),
            'Renilla luciferase': ('renilla_luciferase', 'reporter'),
            'beta-galactosidase': ('beta_galactosidase', 'reporter'),
        }
        for product, expected in cases.items():
            self.assertEqual(curate.map_family(product, None), expected, product)

    def test_marker_and_regulator_products(self):
        self.assertEqual(curate.map_family('beta-lactamase', None), ('beta_lactamase', 'marker'))
        self.assertEqual(curate.map_family('neomycin phosphotransferase', None),
                         ('aminoglycoside_phosphotransferase', 'marker'))
        self.assertEqual(curate.map_family('arabinose operon regulatory protein', None),
                         ('araC', 'regulator'))

    def test_gene_fallback_and_unknown(self):
        self.assertEqual(curate.map_family(None, 'bla'), ('beta_lactamase', 'marker'))
        self.assertEqual(curate.map_family('hypothetical protein', None), (None, None))


class AnalysisMathTests(unittest.TestCase):
    def setUp(self):
        # Two lineage groups: a big near-duplicate group (4 components) and a singleton.
        self.assignments = {
            'a1': {'family': 'f', 'role': 'reporter', 'lineage_group': 'g1', 'family_group': 'f'},
            'a2': {'family': 'f', 'role': 'reporter', 'lineage_group': 'g1', 'family_group': 'f'},
            'a3': {'family': 'f', 'role': 'reporter', 'lineage_group': 'g1', 'family_group': 'f'},
            'a4': {'family': 'f', 'role': 'reporter', 'lineage_group': 'g1', 'family_group': 'f'},
            'b1': {'family': 'h', 'role': 'marker', 'lineage_group': 'g2', 'family_group': 'h'},
        }

    def test_pooled_versus_equal_group_weighting(self):
        ids = set(self.assignments)
        groups = analyze._group_map(self.assignments, ids, 'lineage_group')
        recovered = {'a1', 'a2', 'a3', 'a4'}  # the whole big group, not the singleton
        # Pooled: 4/5 = 0.8; equal-group: (1.0 + 0.0)/2 = 0.5
        self.assertAlmostEqual(analyze._pooled(recovered, ids), 0.8)
        self.assertAlmostEqual(analyze._equal_group(recovered, groups), 0.5)

    def test_effect_and_bootstrap_interval_orientation(self):
        ids = set(self.assignments)
        a = {'a1', 'a2', 'a3', 'a4', 'b1'}  # candidate recovers everything
        b = {'a1', 'a2', 'a3', 'a4'}        # baseline misses the singleton group
        effect = analyze._effect(a, b, self.assignments, ids, 'lineage_group', reps=200, seed=1)
        self.assertAlmostEqual(effect['equal_group_gain'], 0.5)   # singleton group 0 -> 1
        self.assertAlmostEqual(effect['pooled_gain'], 0.2)        # 4/5 -> 5/5
        lo, hi = effect['equal_group_interval']
        self.assertLessEqual(lo, effect['equal_group_gain'])
        self.assertGreaterEqual(hi, 0.0)


if __name__ == '__main__':
    unittest.main()
