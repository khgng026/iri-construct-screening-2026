"""Pure numerical/metadata fixtures test pilot scoring; no research measurements."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.pilot_scoring import (freeze_panel, load_scoring_config,
    paired_cluster_bootstrap, rank_scores, score_rankings, validate_scoring_config)
from iri_constructs.runtime import InputError


def configuration():
    return {'schema_version': 'pilot-scoring-1', 'protocol_id': 'fictional-test',
        'protocol_version': '1.0', 'experiment': 'D',
        'design_status': 'design_frozen_pending_admission', 'top_k': 2,
        'retrieval_budget': 2, 'seed': 0, 'min_independent_clusters': 8,
        'bootstrap_replicates': 2000, 'max_residues': 1022, 'models': ['esm2', 'esmc']}


class PilotScoringTests(unittest.TestCase):
    def setUp(self):
        self.config = configuration()
        self.proteins = {key: ProteinInput(key, key, aa*3, 'fictional')
                         for key,aa in zip('abcd', 'ACDE')}
        self.labels = {key: {'input_id': key, 'family': 'one' if key in 'ab' else 'two',
            'evidence': 'Fictional unit-test evidence only', 'independence_group': key,
            'annotation_status': 'known', 'review_status': 'approved',
            'archetype': 'tagged_fusion' if key in 'ab' else 'single_cds'} for key in 'abcd'}
        self.panel = freeze_panel(self.config, self.proteins, self.labels)

    def test_configs_load_before_research_inputs_exist(self):
        configs = sorted((ROOT / 'config' / 'pilot').glob('experiment_*.json'))
        self.assertEqual(len(configs), 4)
        self.assertEqual({load_scoring_config(p)['experiment'] for p in configs}, set('ABCD'))

    def test_checked_in_design_freeze_matches_protocol_and_configs(self):
        freeze=json.loads((ROOT/'config'/'pilot'/'design_freeze.json').read_text())
        self.assertEqual(freeze['schema_version'],'pilot-design-freeze-1')
        self.assertEqual(freeze['status'],'design_frozen_pending_admission')
        for path,expected in freeze['files'].items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT/path).read_text(encoding='utf-8').encode('utf-8')).hexdigest(),expected)

    def test_invalid_or_unfrozen_configuration_rejected(self):
        for change in ({'top_k': True}, {'retrieval_budget': 3}, {'models': ['esm2', 'esm2']},
                       {'design_status': 'tuned_after_results'}, {'max_residues': 5000}):
            with self.subTest(change=change), self.assertRaises(InputError):
                validate_scoring_config({**self.config, **change})

    def test_roles_and_method_comparisons_must_be_explicit(self):
        for roles in ('abcd', {'a':True}):
            with self.assertRaises(InputError):
                freeze_panel(self.config,self.proteins,self.labels,roles)
        with self.assertRaises(InputError):
            validate_scoring_config({**self.config,'models':['esm2']})

    def test_budget_saturation_is_visible_without_invented_advantage(self):
        config = {**self.config,'top_k':5,'retrieval_budget':5}
        panel=freeze_panel(config,self.proteins,self.labels)
        self.assertEqual(panel['budget_saturated_query_ids'],list('abcd'))

    def test_ranking_and_family_scoring_are_deterministic(self):
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['b','c','d'])
        rankings = rank_scores({'a': {'b': 1.0, 'c': .5, 'd': .5}}, panel)
        self.assertEqual(rankings, {'a': ['b','c']})
        result = score_rankings(rankings, panel, self.labels)
        self.assertEqual(result['overall']['family_hit_at_k'], 1.0)
        self.assertEqual(result['overall']['documented_family_recall_at_k'], 1.0)
        self.assertEqual(result['overall']['returned_hits'], 2)
        self.assertEqual(result['per_archetype']['tagged_fusion']['scored_queries'], 1)

    def test_empty_or_insufficient_panels_do_not_become_success(self):
        empty = freeze_panel(self.config, {}, {}, [], [])
        result = score_rankings({}, empty, {})
        self.assertEqual(empty['status'], 'blocked_insufficient_panel')
        self.assertEqual(result['status'], 'blocked_no_scorable_queries')
        self.assertIsNone(result['overall']['family_hit_at_k'])
        self.assertEqual(result['coverage'], 0.0)
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['c','d'])
        self.assertEqual(panel['status'], 'blocked_insufficient_panel')
        self.assertEqual(len(score_rankings({}, panel, self.labels)['abstentions']), 1)

    def test_reference_failures_do_not_shrink_the_score_panel(self):
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['b','c','d'])
        for scores in ({'b': 1.0,'c': .5}, {'b': 1.0,'c': .5,'d': float('nan')}):
            self.assertEqual(rank_scores({'a': scores}, panel), {})
        result = score_rankings({}, panel, self.labels)
        self.assertEqual(result['scored_queries'], 0)
        self.assertIsNone(result['macro_construct_documented_family_recall_at_k'])

    def test_same_parent_duplicate_and_group_leakage_rejected(self):
        self.proteins['b'] = ProteinInput('b', 'a', self.proteins['b'].sequence, 'fictional')
        self.proteins['c'] = ProteinInput('c', 'c', self.proteins['a'].sequence, 'fictional')
        self.labels['d']['independence_group'] = self.labels['a']['independence_group']
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['b','c','d'])
        self.assertEqual(panel['eligible_reference_ids']['a'], [])
        with self.assertRaises(InputError):
            score_rankings({'a': ['b']}, panel, self.labels)

    def test_injected_reference_duplicate_hit_and_budget_change_rejected(self):
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['b','c','d'])
        with self.assertRaises(InputError):
            rank_scores({'a': {'a': 100, 'b': 1, 'c': .5,'d': .4}}, panel)
        with self.assertRaises(InputError):
            rank_scores({}, panel, top_k=1)
        for rankings in ({'a':['b','b']}, {'a':['b','c','d']}, {'a':['a','c']}):
            with self.subTest(rankings=rankings), self.assertRaises(InputError):
                score_rankings(rankings, panel, self.labels)

    def test_short_return_budget_is_explicit_abstention(self):
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['b','c','d'])
        result = score_rankings({'a':['b']}, panel, self.labels)
        self.assertEqual(result['scored_queries'], 0)
        self.assertIn('fewer hits', result['abstentions'][0]['reason'])

    def test_changed_annotations_or_panel_fail(self):
        changed = copy.deepcopy(self.labels)
        changed['a']['review_status'] = 'pending'
        with self.assertRaises(InputError):
            freeze_panel(self.config, self.proteins, changed)
        changed['a']['review_status'] = 'approved'
        changed['a']['family'] = 'posthoc-family'
        with self.assertRaises(InputError):
            score_rankings({}, self.panel, changed)
        altered = copy.deepcopy(self.panel)
        altered['reference_ids'].append('invented')
        with self.assertRaises(InputError):
            rank_scores({}, altered)

    def test_duplicate_role_ids_are_rejected(self):
        with self.assertRaises(InputError):
            freeze_panel(self.config, self.proteins, self.labels, ['a','a'], ['b','c'])
        changed = copy.deepcopy(self.labels)
        changed['a']['input_id'] = 'b'
        with self.assertRaisesRegex(InputError, 'map keys'):
            freeze_panel(self.config, self.proteins, changed)

    def test_frozen_expected_families_cannot_be_changed_after_ranking(self):
        panel = freeze_panel(self.config, self.proteins, self.labels, ['a'], ['b','c','d'], {'a':['one']})
        with self.assertRaises(InputError):
            score_rankings({}, panel, self.labels, component_families={'a':['two']})

    def test_overlapping_roles_exclude_query_self_reference(self):
        self.assertNotIn('a', self.panel['eligible_reference_ids']['a'])
        self.assertIn('b', self.panel['eligible_reference_ids']['a'])


class PairedUncertaintyTests(unittest.TestCase):
    def test_cluster_bootstrap_is_seeded_and_uses_group_weight(self):
        left = {str(i): .25 for i in range(10)}
        right = {str(i): .75 for i in range(10)}
        clusters = {str(i): 'group-'+str(i) for i in range(10)}
        result = paired_cluster_bootstrap(left, right, clusters)
        self.assertEqual(result, paired_cluster_bootstrap(left, right, clusters))
        self.assertEqual(result['ci_lower'], .5)
        self.assertEqual(result['ci_upper'], .5)
        self.assertEqual(result['independent_clusters'], 10)

    def test_missing_or_small_paired_sets_have_null_intervals(self):
        result = paired_cluster_bootstrap({'a':0.0,'b':0.0}, {'a':1.0}, {'a':'linked'})
        self.assertEqual(result['excluded_unpaired_query_ids'], ['b'])
        self.assertEqual(result['status'], 'insufficient_independent_clusters')
        self.assertIsNone(result['ci_lower'])
        self.assertEqual(result['paired_cluster_mean_difference'], 1.0)
        self.assertIsNone(paired_cluster_bootstrap({}, {}, {})['paired_cluster_mean_difference'])

    def test_invalid_proportions_and_missing_cluster_fail(self):
        for value in (float('nan'), True, 2.0):
            with self.subTest(value=value), self.assertRaises(InputError):
                paired_cluster_bootstrap({'a':0.0}, {'a':value}, {'a':'g'})
        with self.assertRaises(InputError):
            paired_cluster_bootstrap({'a':0.0}, {'a':1.0}, {})
