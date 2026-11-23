"""Regression tests for the cohort v2 freeze, runner helpers and paired analysis (fictional data)."""
import importlib.util
import random
import unittest

from experiments.paired_preprocessing import cohort_study as cs

# The analysis module needs NumPy, which the dependency-free CI jobs do not install.
HAS_NUMPY = importlib.util.find_spec('numpy') is not None
if HAS_NUMPY:
    from experiments.paired_preprocessing import cohort_analysis as ca
from experiments.paired_preprocessing.verify import _covered


class Ref:
    def __init__(self, authors='', title='', journal=''):
        self.authors, self.title, self.journal = authors, title, journal


def _random_dna(n, seed):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


class GroupingTests(unittest.TestCase):
    def test_kmers_ignore_strand_and_rotation(self):
        dna = _random_dna(600, 1)
        reverse = dna.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        rotated = dna[200:] + dna[:200]
        base = cs.sampled_kmers(dna, 25, 8)
        self.assertTrue(base)
        self.assertEqual(base, cs.sampled_kmers(reverse, 25, 8))
        self.assertEqual(base, cs.sampled_kmers(rotated, 25, 8))

    def test_near_copy_requires_mutual_containment(self):
        a = _random_dna(2000, 2)
        sets = {'A': cs.sampled_kmers(a, 25, 8), 'B': cs.sampled_kmers(a[:1900], 25, 8),
                'C': cs.sampled_kmers(a + _random_dna(4000, 3), 25, 8)}
        pairs = cs.near_copy_pairs(sets, 0.9)
        self.assertIn(('A', 'B'), pairs)
        # C contains A, but A covers only a third of C: lineage by containment alone is not a near copy.
        self.assertNotIn(('A', 'C'), pairs)

    def test_lineage_keys_separate_submission_publication_and_catalog(self):
        keys = cs.lineage_keys([Ref('Doe,J.', 'Direct Submission', 'Submitted (01-JAN-2000) Lab'),
                                Ref('Doe,J.', 'A vector paper', 'Gene 1 (1), 1-2 (2000)'),
                                Ref('Doe,J.', 'Vector catalog', 'Unpublished')])
        self.assertEqual({k[0] for k in keys}, {'submission', 'publication_title', 'unpublished_title'})

    def test_union_find_root_is_deterministic(self):
        uf = cs.UnionFind(['b', 'a', 'c'])
        uf.union('c', 'b')
        uf.union('b', 'a')
        self.assertEqual(uf.groups(), {'a': ['a', 'b', 'c']})


class UnitTests(unittest.TestCase):
    def test_windows_keep_forward_and_reverse_coordinates(self):
        protein = 'M' + 'A' * 699
        for strand, spans in ((1, [[100, 2200]]), (-1, [[100, 2200]])):
            windows = cs._windows([{'id': 'q1', 'strand': strand, 'spans': spans, 'protein': protein}])
            self.assertEqual([len(w['protein']) for w in windows], [300] * 4)
            self.assertEqual([w['id'] for w in windows], ['q1w0', 'q1w150', 'q1w300', 'q1w400'])
            first = windows[0]['spans'][0]
            self.assertEqual(first, [100, 1000] if strand == 1 else [1300, 2200])
            last = windows[-1]['spans'][0]
            self.assertEqual(last, [1300, 2200] if strand == 1 else [100, 1000])

    def test_short_candidates_are_not_windowed(self):
        unit = {'id': 'q1', 'strand': 1, 'spans': [[0, 900]], 'protein': 'A' * 300}
        self.assertEqual(cs._windows([unit]), [unit])

    def test_chunk_unit_maps_joined_spans(self):
        protein = 'A' * 1100
        unit = {'id': 'u', 'strand': 1, 'spans': [[0, 1500], [2000, 3800]], 'protein': protein}
        chunks = cs.chunk_unit(unit)
        self.assertEqual([len(c['protein']) for c in chunks], [1022, 78])
        self.assertEqual(chunks[0]['spans'], [[0, 1500], [2000, 3566]])
        reverse = dict(unit, strand=-1)
        first = cs.chunk_unit(reverse)[0]
        # Translation starts at the 3' end of the last segment on the reverse strand.
        self.assertEqual(first['spans'], [[234, 1500], [2000, 3800]])

    def test_documented_unit_over_limit_reaches_the_cap(self):
        units = [{'id': 'u', 'strand': 1, 'spans': [[0, 3180]], 'protein': 'A' * 1060}]
        pool = [u for u in cs.candidate_units('annotation_reference', units, None)
                if cs.MIN_RESIDUES <= len(u['protein']) <= cs.MAX_RESIDUES]
        self.assertEqual(len(pool[0]['protein']), 1022)
        self.assertGreaterEqual(_covered((1, [(0, 3180)]), (1, [tuple(pool[0]['spans'][0])])), 0.8)


def _rows(spec):
    """spec: list of (group, family, a_recovered, b_recovered)."""
    a, b = ('s', 'A', 2), ('s', 'B', 2)
    return [{'lineage_group': g, 'family': f, 'outcome': {a: {'recovered': x}, b: {'recovered': y}}}
            for g, f, x, y in spec]


@unittest.skipUnless(HAS_NUMPY, 'cohort analysis tests require NumPy')
class EffectTests(unittest.TestCase):
    def test_pooled_equal_group_and_leave_one_out(self):
        rows = _rows([('g1', 'f1', 1, 0), ('g1', 'f1', 1, 0), ('g1', 'f2', 0, 0), ('g1', 'f2', 0, 0),
                      ('g2', 'f1', 0, 1), ('g3', 'f2', 1, 1)])
        e = ca.paired_effect(rows, ('s', 'A', 2), ('s', 'B', 2))
        self.assertAlmostEqual(e['pooled_gain'], 1 / 6)
        self.assertAlmostEqual(e['equal_group_gain'], (0.5 - 1 + 0) / 3)
        self.assertEqual((e['groups_improved'], e['groups_declined'], e['groups_tied']), (1, 1, 1))
        self.assertAlmostEqual(e['pooled_gain_without_largest_group'], -1 / 2)
        low, high = e['group_bootstrap']['pooled_gain_interval']
        self.assertLessEqual(low, high)
        self.assertEqual(len(e['two_way_bootstrap_pooled_gain_interval']), 2)

    def test_verdict_rules(self):
        base = {'equal_group_gain': 0.1, 'leave_one_group_out_equal_group_range': [0.06, 0.12]}
        self.assertTrue(ca.verdict({**base, 'group_bootstrap': {'equal_group_gain_interval': [0.06, 0.14]}})
                        .startswith('practically meaningful'))
        self.assertTrue(ca.verdict({**base, 'group_bootstrap': {'equal_group_gain_interval': [0.01, 0.2]}})
                        .startswith('positive signal'))
        self.assertTrue(ca.verdict({**base, 'equal_group_gain': 0.01,
                                    'group_bootstrap': {'equal_group_gain_interval': [-0.01, 0.04]}})
                        .startswith('gain of the target size unsupported'))
        self.assertTrue(ca.verdict({**base, 'group_bootstrap': {'equal_group_gain_interval': [-0.05, 0.2]}})
                        .startswith('inconclusive'))
        self.assertTrue(ca.verdict(None).startswith('descriptive'))


if __name__ == '__main__':
    unittest.main()


class ScorerModuleTests(unittest.TestCase):
    def test_scorer_and_partition_tables_are_consistent(self):
        from experiments.paired_preprocessing import cohort_scorers as cse
        for model, _ in cse.SCORERS.values():
            self.assertIn(model, cse.MODELS)
        for spec in cse.PARTITIONS.values():
            self.assertTrue(set(spec['units']).issubset(spec['arms']))

    def test_precomputed_characterizer_returns_stored_rankings(self):
        from experiments.paired_preprocessing import cohort_scorers as cse
        protein = 'MAAAA'
        rankings = {'model': {'model': 'esm2_8m'},
                    'mean_rankings': {cse.sha(protein): ['famA', 'famB']},
                    'residue_rankings': {cse.sha(protein): ['famB']}}
        with cse._precomputed(rankings, 'mean') as (rank, meta):
            self.assertEqual(rank(protein), ['famA', 'famB'])
            self.assertEqual(meta['representation'], 'mean')
        with cse._precomputed(rankings, 'residue') as (rank, _):
            self.assertEqual(rank(protein), ['famB'])

    @unittest.skipUnless(HAS_NUMPY, 'representation helpers require NumPy')
    def test_normalized_rows_are_unit_length(self):
        import numpy as np
        from experiments.paired_preprocessing import cohort_scorers as cse
        states = np.array([[3.0, 4.0], [0.0, 2.0]], dtype=np.float32)
        rows = cse._normalized_rows(states)
        self.assertTrue(np.allclose(np.linalg.norm(rows, axis=1), 1.0))
        self.assertEqual(len(cse._mean_vector(states)), 2)


@unittest.skipUnless(HAS_NUMPY, 'ablation analysis requires NumPy')
class AblationTests(unittest.TestCase):
    def _doc(self):
        ref_family = {'r1': 'A', 'r2': 'A', 'r3': 'B'}
        def per(idA1, idA2, idB, bl):
            return {'r1': {'identity': idA1, 'blosum': bl[0], 'mean': bl[0], 'residue': bl[0]},
                    'r2': {'identity': idA2, 'blosum': bl[1], 'mean': bl[1], 'residue': bl[1]},
                    'r3': {'identity': idB, 'blosum': bl[2], 'mean': bl[2], 'residue': bl[2]}}
        comps = [{'construct': 'c1', 'component': 'c1-A', 'family': 'A', 'lineage_group': 'g1',
                  'protein_sha256': 'p1', 'fusion_class': False}]
        # p1: close famA ref r1 wins; once r1 is stripped, famB ref r3 outscores the distant famA r2.
        scores = {'p1': per(0.9, 0.3, 0.25, (0.9, 0.2, 0.3))}
        return {'model': {'model': 'esm2_8m'}, 'reference_family': ref_family, 'components': comps,
                'scores_by_protein': scores, 'encode_seconds': 0}

    def test_ablation_strips_close_references_and_reports_recall(self):
        import json
        import tempfile
        from experiments.paired_preprocessing import cohort_ablation as ab
        from pathlib import Path
        d = Path(tempfile.mkdtemp())
        (d / 'scores.json').write_text(json.dumps(self._doc()))
        ab.analyze([str(d / 'scores.json')], d / 'out.json')
        out = json.loads((d / 'out.json').read_text())
        # Above every retained-reference identity: the close famA ref makes BLOSUM correct.
        self.assertEqual(out['by_threshold']['1.01']['recall']['blosum62']['recovered'], 1)
        # Once references >= 0.80 identity are removed, the panel's best is famB: BLOSUM misses.
        self.assertEqual(out['by_threshold']['0.80']['recall']['blosum62']['recovered'], 0)
        # Below the distant famA reference's identity (0.30) nothing same-family survives.
        self.assertEqual(out['by_threshold']['0.25']['scorable_components'], 0)

    def test_top_family_takes_best_reference_per_family(self):
        from experiments.paired_preprocessing import cohort_ablation as ab
        per_ref = {'r1': {'s': 0.2}, 'r2': {'s': 0.9}, 'r3': {'s': 0.5}}
        self.assertEqual(ab._top_family(per_ref, {'r1': 'A', 'r2': 'A', 'r3': 'B'}, 's', ['r1', 'r2', 'r3']), 'A')
        self.assertEqual(ab._top_family(per_ref, {'r1': 'A', 'r2': 'A', 'r3': 'B'}, 's', ['r1', 'r3']), 'B')
        self.assertIsNone(ab._top_family(per_ref, {'r1': 'A'}, 's', []))
