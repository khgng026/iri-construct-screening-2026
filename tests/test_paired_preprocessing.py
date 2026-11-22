"""Regression tests for scoring without passing truth IDs to experimental arms."""
import copy
import unittest

from experiments.paired_preprocessing.verify import SCHEMA, verify


def _component(name, spans, strand=1, family='blue'):
    return {'id': name, 'spans': spans, 'strand': strand, 'family': family}


def _query(name, spans, strand=1):
    return {'id': name, 'spans': spans, 'strand': strand, 'input_residues': 10}


def _trace(name, queries, families):
    return {'id': name, 'queries': queries, 'candidate_count': len(queries),
            'predictions': [{'query_id': qid, 'family': f} for qid, f in families],
            'status': 'ok', 'preprocess_seconds': 0.1, 'characterize_seconds': 0.2}


def _run(name, one, two):
    return {'schema': SCHEMA, 'method': name, 'characterization_id': 'fixed',
            'tiers': [{'query_cap': 1, 'constructs': [one]}, {'query_cap': 2, 'constructs': [two]}]}


class VerificationTests(unittest.TestCase):
    def setUp(self):
        self.protocol = {'schema': SCHEMA, 'characterization_id': 'fixed', 'query_caps': [1, 2],
                         'primary_query_cap': 1, 'prediction_limit': 4, 'minimum_component_overlap': 0.8,
                         'minimum_interesting_gain': 0.05}
        self.truth = {'schema': SCHEMA, 'constructs': [
            {'id': 'c', 'archetype': 'multi_cds', 'length': 100,
             'components': [_component('a', [[0, 30]]), _component('b', [[40, 70]])]},
        ]}
        # A broad blue call could cover either label instance. A narrower blue call covers only a.
        broad = _query('wide', [[0, 80]])
        narrow = _query('narrow', [[0, 30]])
        self.one = _trace('c', [broad], [('wide', 'blue'), ('wide', 'unlabelled')])
        self.two = _trace('c', [broad, narrow], [('wide', 'blue'), ('narrow', 'blue')])

    def test_additional_queries_do_not_count_as_a_matched_budget_gain(self):
        base = _run('baseline', self.one, self.two)
        alternative = _run('candidate', self.one, self.two)
        report = verify(self.protocol, self.truth, base, [alternative])
        self.assertEqual(report['runs']['baseline'][1]['overall']['recall'], 0.5)
        self.assertEqual(report['runs']['candidate'][2]['overall']['recall'], 1.0)
        self.assertEqual(report['candidate_minus_baseline']['candidate'][1]['recall_gain'], 0)
        self.assertFalse(report['candidate_minus_baseline']['candidate'][1]['meets_effect_target_by_point_estimate'])
        self.assertEqual(report['runs']['baseline'][1]['overall']['returned_predictions'], 2)

    def test_one_to_one_assignment_uses_narrow_call_for_repeated_label(self):
        base = _run('baseline', self.one, self.one)
        alternative = _run('candidate', self.one, self.two)
        self.protocol['primary_query_cap'] = 2
        report = verify(self.protocol, self.truth, base, [alternative])
        scored = report['runs']['candidate'][2]['per_construct'][0]
        self.assertEqual(scored['recovered_component_ids'], ['a', 'b'])
        self.assertEqual(report['candidate_minus_baseline']['candidate'][2]['recall_gain'], 0.5)
        self.assertTrue(report['candidate_minus_baseline']['candidate'][2]['meets_effect_target_by_point_estimate'])

    def test_wrong_strand_and_bounding_interval_gap_get_no_credit(self):
        self.truth['constructs'][0]['components'] = [_component('join', [[95, 100], [0, 20]], -1)]
        gap = _query('in_gap', [[30, 70]], -1)
        wrong = _query('wrong_strand', [[0, 20], [95, 100]], 1)
        right = _query('right', [[0, 20], [95, 100]], -1)
        one = _trace('c', [gap], [('in_gap', 'blue')])
        two = _trace('c', [gap, wrong], [('in_gap', 'blue'), ('wrong_strand', 'blue')])
        candidate = _trace('c', [right], [('right', 'blue')])
        report = verify(self.protocol, self.truth, _run('baseline', one, two),
                        [_run('candidate', candidate, candidate)])
        self.assertEqual(report['runs']['baseline'][2]['overall']['recall'], 0)
        self.assertEqual(report['runs']['candidate'][1]['overall']['recall'], 1)
        self.assertEqual(report['runs']['baseline'][2]['overall']['components_without_covering_query'], 1)
        self.assertEqual(report['runs']['candidate'][1]['overall']['components_without_covering_query'], 0)

    def test_matching_locus_without_correct_label_is_tracked_separately(self):
        no_match = _trace('c', [_query('q', [[40, 70]])], [('q', 'unlabelled')])
        report = verify(self.protocol, self.truth, _run('base', no_match, no_match),
                        [_run('candidate', self.one, self.two)])
        scores = report['runs']['base'][1]['overall']
        self.assertEqual(scores['components_without_covering_query'], 1)
        self.assertEqual(scores['covered_components_not_recovered'], 1)
        self.assertEqual(scores['recovered_components'], 0)

    def test_failed_construct_is_a_miss_with_same_denominator(self):
        fail = _trace('c', [], [])
        fail['status'] = 'failed'
        report = verify(self.protocol, self.truth, _run('base', fail, self.one),
                        [_run('candidate', self.one, self.two)])
        self.assertEqual(report['runs']['base'][1]['overall']['failed_constructs'], 1)
        self.assertEqual(report['runs']['base'][1]['overall']['eligible_components'], 2)
        self.assertEqual(report['runs']['base'][1]['overall']['recall'], 0)

    def test_enforces_caps_same_panel_and_output_limit(self):
        base = _run('baseline', self.one, self.two)
        alternative = _run('candidate', self.one, self.two)
        for edit, error in [
            (lambda run: run.update(characterization_id='different'), 'characterization workflow'),
            (lambda run: run['tiers'].pop(), 'missing budget tiers'),
            (lambda run: run['tiers'][0].update(constructs=[]), 'construct set differs'),
            (lambda run: run['tiers'][0].update(constructs=[self.two]), 'exceeds cap'),
            (lambda run: run['tiers'][0]['constructs'][0].update(predictions=self.one['predictions'] * 3),
             'prediction count exceeds'),
            (lambda run: run['tiers'][0]['constructs'][0]['predictions'][0].update(query_id='unknown'),
             'unqueried input'),
        ]:
            with self.subTest(error=error):
                bad = copy.deepcopy(alternative)
                edit(bad)
                with self.assertRaisesRegex(ValueError, error):
                    verify(self.protocol, self.truth, base, [bad])


if __name__ == '__main__':
    unittest.main()
