"""Metadata/denominator tests only; no saved retrieval outcomes or scoring."""
from pathlib import Path
from copy import deepcopy
import sys
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'src')]
import session_expanded_family as family


def label(target, survivors, common=True):
    return {'accession': target, 'accepted_survivors': survivors,
            'native_scorable': bool(survivors), 'common_scorable': common and bool(survivors)}


def ranking(*accessions):
    return [{'accession': a, 'rank': i+1, 'score': 1-i*.01} for i, a in enumerate(accessions)]


class FirstFamilyTests(unittest.TestCase):
    def test_exact_first_sentence_and_subfamily_evidence(self):
        self.assertEqual(family.tsv_family('SIMILARITY: Belongs to the example family. Narrow subfamily. {ECO:0000305}.'),
                         'Belongs to the example family')
        self.assertIsNone(family.tsv_family('SIMILARITY: Contains a domain. Belongs to the example family.'))
        self.assertIsNone(family.tsv_family('SIMILARITY: belongs to the example family.'))
        self.assertIsNone(family.tsv_family('SIMILARITY: Belongs to the example superfamily.'))
        self.assertIsNone(family.tsv_family(''))

    def test_multiple_first_families_block_freeze(self):
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            family.tsv_family('SIMILARITY: Belongs to the first family. SIMILARITY: Belongs to the second family.')
        with self.assertRaisesRegex(ValueError, 'transport'):
            family.tsv_family('Unmarked transport text')


class RemovalAndDenominatorTests(unittest.TestCase):
    def test_removal_is_sequence_based_and_query_wide(self):
        query = {'components': ['A', 'B']}
        mask = family.query_mask(query, ['a', 'same_a', 'b', 'x', 'y'], {'A': 'a', 'B': 'b'},
                                 {'a': 'h1', 'same_a': 'h1', 'b': 'h2', 'x': 'h3', 'y': 'h4'})
        self.assertEqual(mask['removed_representatives'], ['a', 'b', 'same_a'])
        self.assertEqual(mask['remaining_count'], 2)
        self.assertEqual(mask['included_representatives_sha256'], family.digest(family.canonical(['x', 'y'])))

    def test_common_stratum_keeps_ineligible_original_component(self):
        query = {'query_id': 'q', 'components': ['A', 'B'], 'group_id': 'A__B'}
        labels = {'components': [label('A', ['x']), label('B', ['y'], common=False)]}
        native = family.evaluate_filtered(query, ranking('x', 'y', 'z'), labels, 'native')
        common = family.evaluate_filtered(query, ranking('x', 'y', 'z'), labels, 'common11')
        self.assertTrue(native['all_components_scorable'])
        self.assertEqual(native['query_hits']['5']['all'], 1)
        self.assertFalse(common['all_components_scorable'])
        self.assertIsNone(common['query_hits']['5']['all'])
        self.assertEqual(len(common['components']), 2)
        self.assertFalse(common['components'][1]['scorable'])
        self.assertEqual(common['components'][1]['accepted_survivors'], ['y'])
        acc = family.Accumulator(); acc.add(query, common)
        result = acc.result(family.denominator([labels], 'common11'))
        self.assertEqual(result['original_components'], 2)
        self.assertEqual(result['scorable_components'], 1)
        self.assertEqual(result['excluded_unscorable_queries'], 1)
        self.assertIsNone(result['by_k']['5']['all_components_hit'])

    def test_empty_accepted_labels_are_unscorable_not_misses(self):
        query = {'components': ['A'], 'group_id': 'A'}; labels = {'components': [label('A', [])]}
        result = family.evaluate_filtered(query, ranking('x', 'y'), labels, 'native')
        self.assertIsNone(result['components'][0]['rank'])
        self.assertIsNone(result['components'][0]['hit']['5'])
        acc = family.Accumulator(); acc.add(query, result)
        summary = acc.result(family.denominator([labels], 'native'))
        self.assertIsNone(summary['by_k']['5']['component_recall'])
        self.assertIsNone(summary['by_k']['5']['all_components_hit'])

    def test_available_but_low_rank_is_a_scorable_miss(self):
        query = {'components': ['A'], 'group_id': 'A'}; labels = {'components': [label('A', ['a'])]}
        result = family.evaluate_filtered(query, ranking('v', 'w', 'x', 'y', 'z', 'a'), labels, 'native')
        self.assertTrue(result['components'][0]['scorable'])
        self.assertEqual(result['components'][0]['hit']['5'], 0)
        self.assertEqual(result['components'][0]['hit']['10'], 1)

    def test_missing_survivor_or_changed_component_order_rejected(self):
        with self.assertRaisesRegex(ValueError, 'missing'):
            family.evaluate_filtered({'components': ['A']}, ranking('x'), {'components': [label('A', ['a'])]}, 'native')
        with self.assertRaisesRegex(ValueError, 'order'):
            family.evaluate_filtered({'components': ['A', 'B']}, ranking('x', 'y'),
                                     {'components': [label('B', ['y']), label('A', ['x'])]}, 'native')

    def test_frozen_denominator_is_an_acceptance_gate(self):
        query = {'components': ['A'], 'group_id': 'A'}; labels = {'components': [label('A', ['a'])]}
        acc = family.Accumulator(); acc.add(query, family.evaluate_filtered(query, ranking('a'), labels, 'native'))
        wrong = family.denominator([labels], 'native'); wrong['scorable_components'] = 2
        with self.assertRaisesRegex(ValueError, 'frozen metadata'):
            acc.result(wrong)

    def test_macro_target_metric_does_not_weight_repeated_target_counts(self):
        acc = family.Accumulator(); labels = []
        for i, target in enumerate(['A', 'A', 'A', 'B']):
            row = {'components': [label(target, ['hit' if target == 'A' else 'miss'])]}; labels.append(row)
            query = {'components': [target], 'group_id': 'g' + str(i)}
            acc.add(query, family.evaluate_filtered(query, ranking('hit', 'v', 'w', 'x', 'y', 'miss'), row, 'native'))
        result = acc.result(family.denominator(labels, 'native'))
        self.assertEqual(result['by_k']['5']['component_recall'], .75)
        self.assertEqual(result['by_k']['5']['per_target_macro_recall'], .5)


class CompletionCoverageTests(unittest.TestCase):
    def test_full_unique_query_panel_required(self):
        identity = {'mode': 'expanded_model_retrieval_v1', 'query_ids': ['q'+str(i) for i in range(1146)]}
        self.assertEqual(len(family.expected_units(identity)), 1147)
        identity['query_ids'][-1] = identity['query_ids'][0]
        with self.assertRaises(ValueError): family.expected_units(identity)

    def test_alignment_requires_selfscores_and_every_query_shard(self):
        identity = {'mode': 'expanded_alignment_retrieval_v1', 'query_ids': ['q'+str(i) for i in range(1146)],
                    'new_reference_shards': [['a'], ['b']], 'reference_sequence_sha256': {'a': 'ha', 'b': 'hb'}}
        keys = family.expected_units(identity, ['a', 'b'], {'a': 'ha', 'b': 'hb'})
        self.assertEqual(len(keys), 1146*3+2)
        self.assertIn('raw-q1145-001', keys)
        self.assertIn('reference-selfscores', keys)

    def test_partial_method_roster_rejected_before_outcome_access(self):
        args = SimpleNamespace(model_job=[Path('one')] * 4, alignment_job=Path('alignment'))
        with self.assertRaisesRegex(ValueError, 'all five'):
            family.open_sources(args, Path('session'), {}, {}, {}, {})


class SourceAdmissionTests(unittest.TestCase):
    def setUp(self):
        self.order = ['ref' + str(i) for i in range(899)]
        self.hashes = {'base': 'base-sha', **{a: 'sha-' + a for a in self.order}}
        self.identity = {'mode': 'expanded_alignment_retrieval_v1',
                         'query_ids': ['q' + str(i) for i in range(1146)],
                         'new_reference_shards': [self.order[i:i+32] for i in range(0, 899, 32)],
                         'reference_sequence_sha256': self.hashes.copy()}

    def test_exact_899_partition_and_canonical_hashes_accept(self):
        family.validate_alignment_partition(self.identity, self.order, self.hashes)
        self.assertEqual(len(family.expected_units(self.identity, self.order, self.hashes)), 1146 * 30 + 2)

    def test_empty_missing_duplicate_reordered_and_oversize_shards_reject(self):
        variants = [[], [[]], self.identity['new_reference_shards'][:-1],
                    [[*self.order[:-1], self.order[0]]],
                    [self.order[32:64], self.order[:32], *self.identity['new_reference_shards'][2:]],
                    [self.order[:33], self.order[33:]]]
        for shards in variants:
            with self.subTest(shards=shards[:1]):
                changed = {**self.identity, 'new_reference_shards': shards}
                with self.assertRaises(ValueError):
                    family.expected_units(changed, self.order, self.hashes)
        missing = dict(self.identity); missing.pop('new_reference_shards')
        with self.assertRaises((KeyError, ValueError)):
            family.expected_units(missing, self.order, self.hashes)

    def test_duplicate_with_valid_shard_sizes_rejects(self):
        changed = deepcopy(self.identity)
        changed['new_reference_shards'][-1][-1] = self.order[0]
        with self.assertRaisesRegex(ValueError, 'order/coverage'):
            family.validate_alignment_partition(changed, self.order, self.hashes)

    def test_changed_missing_extra_reference_hash_rejects(self):
        variants = [{**self.hashes, 'base': 'changed'},
                    {k: v for k, v in self.hashes.items() if k != 'base'},
                    {**self.hashes, 'outside-prefix': 'extra'}]
        for hashes in variants:
            with self.subTest(size=len(hashes)):
                with self.assertRaisesRegex(ValueError, 'sequence hashes'):
                    family.validate_alignment_partition({**self.identity, 'reference_sequence_sha256': hashes},
                                                        self.order, self.hashes)

    def test_source_dependencies_bound_even_when_wrapper_unchanged(self):
        for mode, key in [('expanded_model_retrieval_v1', 'expanded_model_code_hashes'),
                          ('expanded_alignment_retrieval_v1', 'expanded_alignment_code_hashes')]:
            approved = {'wrapper.py': 'original-wrapper', 'scientific_dependency.py': 'original-dependency'}
            execution = {key: approved}
            family.validate_source_code({'mode': mode, 'code': approved.copy()}, execution)
            variants = [{**approved, 'scientific_dependency.py': 'changed'},
                        {'wrapper.py': 'original-wrapper'}, {**approved, 'extra.py': 'extra'}]
            for code in variants:
                with self.subTest(mode=mode, code=code):
                    with self.assertRaisesRegex(ValueError, 'dependency dictionary'):
                        family.validate_source_code({'mode': mode, 'code': code}, execution)


if __name__ == '__main__': unittest.main()
