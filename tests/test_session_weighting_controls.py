"""Fictional vectors and outcomes only; no real parent arrays or models."""
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'scripts'), str(ROOT/'src')]
import numpy as np
import session_weighting_controls as control


def rows(scores):
    counts = control.Counter(scores.values())
    return [{'accession': a, 'rank': i+1, 'score': scores[a], 'exact_score_tie_count': counts[scores[a]]}
            for i, a in enumerate(sorted(scores, key=lambda a: (-scores[a], a)))]


def fixture(clean=False):
    refs = {'A': np.array([10, 0], np.float32), 'B': np.array([0, 1], np.float32),
            'C': np.array([1, 1], np.float32)}
    components = ['A'] if clean else ['A', 'B']
    lengths = {'A': 1} if clean else {'A': 1, 'B': 3}; total = sum(lengths.values())
    spans = []; offset = 0
    for a in components:
        spans.append({'accession': a, 'start0': offset, 'end0': offset+lengths[a], 'fraction': lengths[a]/total})
        offset += lengths[a]
    query = {'query_id': 'q', 'components': components, 'component_count': len(components), 'spans': spans,
             'length_aa': total, 'group_id': '__'.join(components), 'suite': 'clean' if clean else 'pairs'}
    isolated = {a: refs[a] for a in components}; contextual = {a: refs[a].copy() for a in components}
    arrays = {'whole': control.aggregate(contextual, lengths)[0], 'isolated_mixture': control.aggregate(isolated, lengths)[0],
              **{'isolated_'+a: isolated[a] for a in components}, **{'contextual_'+a: contextual[a] for a in components}}
    matrix = np.stack([refs[a] for a in sorted(refs)])
    rankings = {arm: control.cosine_ranking(arrays[arm], matrix, sorted(refs), dtype=np.float64)
                for arm in ('whole', 'isolated_mixture')}
    rankings.update(isolated_component={}, contextual_component={})
    return query, {'query': deepcopy(query), 'rankings': rankings}, arrays, refs


class ArithmeticTests(unittest.TestCase):
    def test_equal_mean_preserves_raw_norms_and_differs_from_length_weighting(self):
        vectors = {'A': np.array([10, 0], np.float32), 'B': np.array([0, 1], np.float32)}
        equal, raw = control.aggregate(vectors)
        np.testing.assert_array_equal(equal, np.array([5, .5], np.float32))
        np.testing.assert_array_equal(raw, [5, .5])
        weighted, _ = control.aggregate(vectors, {'A': 1, 'B': 3})
        np.testing.assert_array_equal(weighted, [2.5, .75])
        self.assertFalse(np.array_equal(equal, [.5, .5]))

    def test_float64_reduction_precedes_single_float32_cast(self):
        vectors = {'C': np.array([-1e8], np.float32), 'A': np.array([1e8], np.float32), 'B': np.array([1], np.float32)}
        stored, exact = control.aggregate(vectors)
        self.assertEqual(exact[0], 1/3)
        self.assertEqual(stored.dtype, np.float32)
        self.assertEqual(stored[0], np.float32(1/3))
        np.testing.assert_array_equal(control.aggregate(dict(reversed(list(vectors.items()))))[0], stored)

    def test_invalid_vectors_lengths_and_nonfinite_rejected(self):
        for vectors, lengths in [({'A': np.array([1.], np.float64)}, None),
                                 ({'A': np.array([np.nan], np.float32)}, None),
                                 ({'A': np.array([1.], np.float32)}, {'A': 0}),
                                 ({'A': np.array([1.], np.float32)}, {'B': 1})]:
            with self.assertRaises(ValueError): control.aggregate(vectors, lengths)

    def test_clean_controls_reproduce_saved_single_component_baselines(self):
        q, record, arrays, refs = fixture(True)
        result, vectors = control.derive(q, record, arrays, refs)
        self.assertEqual(result['checks']['status'], 'passed')
        for arm, baseline in control.MATCHED.items():
            np.testing.assert_array_equal(vectors[arm], arrays[baseline])
            self.assertEqual(result['evaluations'][arm], result['evaluations'][baseline])
            self.assertEqual(result['checks']['clean_cosine'][arm]['status'], 'passed')

    def test_clean_small_vector_error_can_still_fail_separate_score_gate(self):
        q, record, arrays, refs = fixture(True)
        arrays['contextual_A'][1] = np.float32(1e-7)
        result, vectors = control.derive(q, record, arrays, refs)
        self.assertLessEqual(result['checks']['clean_vector_max_abs']['equal_contextual'], control.VECTOR_TOL)
        self.assertLessEqual(result['checks']['weighted_vector_max_abs']['whole'], control.VECTOR_TOL)
        self.assertEqual(result['checks']['baseline_cosine']['whole']['status'], 'passed')
        self.assertGreater(result['checks']['clean_cosine']['equal_contextual']['maximum_absolute_score_error'], control.SCORE_TOL)
        self.assertEqual(result['checks']['clean_cosine']['equal_contextual']['status'], 'failed')
        self.assertEqual(result['checks']['status'], 'failed')
        self.assertEqual(vectors['equal_contextual'][1], np.float32(1e-7))
        self.assertEqual(record['rankings']['whole'], result['baseline_top10']['whole'])

    def test_weighted_reconstruction_and_new_controls_use_same_catalogue(self):
        q, record, arrays, refs = fixture()
        result, vectors = control.derive(q, record, arrays, refs)
        self.assertEqual(result['checks']['status'], 'passed')
        self.assertEqual(result['checks']['weighted_vector_max_abs'], {'whole': 0., 'isolated_mixture': 0.})
        for arm in control.MATCHED:
            np.testing.assert_array_equal(vectors[arm], [5, .5])
            self.assertEqual({r['accession'] for r in result['new_rankings'][arm]}, set(refs))
            self.assertEqual(len(result['new_rankings'][arm]), len(refs))

    def test_reconstruction_failure_is_retained_not_forced_to_pass(self):
        q, record, arrays, refs = fixture(); arrays['whole'] = arrays['whole'] + .01
        record['rankings']['whole'] = control.cosine_ranking(arrays['whole'], np.stack(list(refs.values())), list(refs), dtype=np.float64)
        result, _ = control.derive(q, record, arrays, refs)
        self.assertEqual(result['checks']['status'], 'failed')
        self.assertGreater(result['checks']['weighted_vector_max_abs']['whole'], control.VECTOR_TOL)

    def test_source_isolated_vector_substitution_rejected(self):
        q, record, arrays, refs = fixture(); arrays['isolated_A'] = arrays['isolated_A'] + .001
        with self.assertRaisesRegex(ValueError, 'catalogue'):
            control.derive(q, record, arrays, refs)

    def test_baseline_rank_inversions_retained_with_fixed_near_tie_policy(self):
        saved = rows({'A': .5, 'B': .5-1e-13})
        near = rows({'A': .5-1e-13, 'B': .5})
        result = control.baseline_reproduction(saved, near, ['A', 'B'])
        self.assertEqual(result['status'], 'passed'); self.assertEqual(len(result['inverted_pairs']), 1)
        self.assertTrue(result['saved_baseline_ranks_preserved'])
        self.assertEqual(saved[0]['accession'], 'A')
        far = control.baseline_reproduction(saved, rows({'A': .4, 'B': .6}), ['A', 'B'])
        self.assertEqual(far['status'], 'failed')


class MetricTests(unittest.TestCase):
    def test_exact_accession_total_budget_all_any_and_original_geometry(self):
        q, _, _, _ = fixture()
        result = control.evaluate(q, rows({'A': .9, 'B': .7, 'C': .8}))
        self.assertEqual(result['query']['1'], {'any': 1, 'all': 0})
        self.assertEqual(result['query']['5'], {'any': 1, 'all': 1})
        self.assertEqual([c['fraction'] for c in result['components']], [.25, .75])
        self.assertEqual([c['ordinal_position'] for c in result['components']], [1, 2])
        self.assertEqual([c['position'] for c in result['components']], ['N', 'C'])
        self.assertEqual(result['components'][1]['rank'], 3)

    def test_gain_loss_pairs_retain_both_directions_and_query_membership(self):
        q, _, _, _ = fixture()
        before = control.evaluate(q, rows({'A': .9, 'B': .7, 'C': .8}))
        after = control.evaluate(q, rows({'A': .7, 'B': .9, 'C': .8}))
        result = control.paired(before, after)
        self.assertEqual([c['outcomes']['1'] for c in result['components']], ['loss', 'gain'])
        self.assertEqual(result['query']['1'], {'any': 'both_hit', 'all': 'neither'})
        totals = control.Totals(['A','B'], True); totals.add(q, result); summary = totals.result()
        self.assertEqual(summary['by_k']['1']['component'], {'denominator': 2, 'both_hit': 0, 'gain': 1,
            'loss': 1, 'neither': 0, 'net_gain': 0, 'net_gain_rate': 0.})
        self.assertEqual(summary['by_k']['1']['any']['denominator'], 1)

    def test_empty_strata_null_and_component_denominators_are_preserved(self):
        q, _, _, _ = fixture(); totals = control.Totals(['A','B','unobserved'])
        totals.add(q, control.evaluate(q, rows({'A': .9, 'B': .7, 'C': .8})))
        result = totals.result()
        self.assertEqual(result['queries'], 1); self.assertEqual(result['components'], 2)
        self.assertEqual(result['by_k']['1']['component'], {'hits': 1, 'denominator': 2, 'recall': .5})
        self.assertIsNone(result['component_strata']['by_fraction']['lt_0.25']['1']['recall'])
        self.assertIsNone(result['component_strata']['by_target']['unobserved']['5']['recall'])
        self.assertIsNone(result['query_by_component_count']['4']['5']['all']['recall'])
        self.assertEqual(result['by_k']['1']['target_macro'], {'support': 2, 'recall': .5})

    def test_missing_exact_accession_and_altered_fraction_rejected(self):
        q, _, _, _ = fixture()
        with self.assertRaisesRegex(ValueError, 'absent'):
            control.evaluate(q, rows({'A': .9, 'C': .8}))
        q['spans'][0]['fraction'] = .1
        with self.assertRaisesRegex(ValueError, 'fraction'):
            control.evaluate(q, rows({'A': .9, 'B': .7, 'C': .8}))

    def test_semantic_resumption_checks_record_and_vector_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            store = control.UnitStore(Path(folder)/'job', {'mode': 'toy'})
            arrays = {'equal_isolated': np.array([1.,2.], np.float32)}
            control.persist(store, 'q', {'value': 1}, arrays)
            control.persist(store, 'q', {'value': 1}, arrays)
            with self.assertRaisesRegex(ValueError, 'Resumed'):
                control.persist(store, 'q', {'value': 2}, arrays)
            with self.assertRaisesRegex(ValueError, 'Resumed'):
                control.persist(store, 'q', {'value': 1}, {'equal_isolated': np.array([1.,3.], np.float32)})
            control.persist(store, 'summary', {'passed': True})
            with self.assertRaisesRegex(ValueError, 'Resumed'):
                control.persist(store, 'summary', {'passed': False})


if __name__ == '__main__': unittest.main()
