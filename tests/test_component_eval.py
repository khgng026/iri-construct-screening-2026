"""Aggregation and attribution rules on hand-built score tables (no model, no biology)."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.component_eval import (attributed_recall, family_recall, merge_max, merge_round_robin,
                                           merge_single, unit_retrieval)
from iri_constructs.runtime import InputError

# Unit "g" scores every reference high; unit "o" is correct but on a lower scale.
REF_FAMILY = {'g1': 'G', 'g2': 'G', 'g3': 'G', 'o1': 'O', 'o2': 'O', 'x1': 'X'}
SCORES = {'g': {'g1': .99, 'g2': .98, 'g3': .97, 'o1': .90, 'o2': .89, 'x1': .88},
          'o': {'g1': .60, 'g2': .59, 'g3': .58, 'o1': .80, 'o2': .79, 'x1': .50}}
FAMILIES = {'g': {'G'}, 'o': {'O'}}


class ComponentEvalTests(unittest.TestCase):
    def test_max_merge_lets_one_unit_take_the_whole_budget(self):
        hits = merge_max(SCORES, 3)
        self.assertEqual([h['reference_id'] for h in hits], ['g1', 'g2', 'g3'])
        self.assertEqual({h['origin_unit'] for h in hits}, {'g'})
        result = attributed_recall(hits, FAMILIES, REF_FAMILY)
        self.assertEqual((result['recovered_components'], result['supported_components']), (1, 2))

    def test_round_robin_shares_the_same_budget_and_keeps_origin(self):
        hits = merge_round_robin(SCORES, ['g', 'o'], 3)
        self.assertEqual([(h['reference_id'], h['origin_unit'], h['origin_rank']) for h in hits],
                         [('g1', 'g', 1), ('o1', 'o', 1), ('g2', 'g', 2)])
        self.assertEqual(attributed_recall(hits, FAMILIES, REF_FAMILY)['attributed_recall'], 1.0)
        # Deduplication: a reference already taken is skipped, not counted twice.
        same = {'a': {'r1': 1.0, 'r2': 0.5}, 'b': {'r1': 0.9, 'r2': 0.1}}
        self.assertEqual([h['reference_id'] for h in merge_round_robin(same, ['a', 'b'], 2)], ['r1', 'r2'])
        # More units than budget: later units contribute nothing.
        self.assertEqual([h['origin_unit'] for h in merge_round_robin(SCORES, ['o', 'g'], 1)], ['o'])

    def test_origin_agnostic_endpoint_gives_credit_that_attribution_withholds(self):
        hits = merge_single(SCORES, 'g', 5)            # the single "g" unit also returns O-family references
        self.assertEqual(family_recall(hits, FAMILIES, REF_FAMILY),
                         {'expected_families': 2, 'recovered_families': 2, 'family_recall': 1.0, 'unattributed_credit': ['O']})
        attributed = attributed_recall(hits, FAMILIES, REF_FAMILY)
        self.assertEqual((attributed['recovered_units'], attributed['attributed_recall']), (['g'], 0.5))
        self.assertEqual(attributed['hit_composition'], {'attributed_correct': 3, 'attributed_wrong': 2, 'from_unlabeled_unit': 0})

    def test_unlabeled_units_use_budget_but_are_never_scored(self):
        hits = merge_max(SCORES, 3)
        result = attributed_recall(hits, {'o': {'O'}}, REF_FAMILY)    # "g" has no documented family
        self.assertEqual(result['hit_composition']['from_unlabeled_unit'], 3)
        self.assertEqual((result['recovered_components'], result['attributed_recall']), (0, 0.0))

    def test_unit_retrieval_and_ties(self):
        self.assertEqual(unit_retrieval(SCORES['o'], {'O'}, REF_FAMILY, 2),
                         {'first_same_family_rank': 1, 'hit_at_1': True, 'hit_at_k': True, 'reciprocal_rank': 1.0, 'same_family_in_top_k': 2})
        self.assertEqual(unit_retrieval(SCORES['g'], {'O'}, REF_FAMILY, 3)['first_same_family_rank'], 4)
        self.assertEqual(unit_retrieval(SCORES['g'], {'Z'}, REF_FAMILY, 3)['reciprocal_rank'], 0.0)
        tied = {'u': {'b': 1.0, 'a': 1.0}, 'v': {'b': 1.0, 'a': 1.0}}
        self.assertEqual([(h['reference_id'], h['origin_unit']) for h in merge_max(tied, 2)], [('a', 'u'), ('b', 'u')])

    def test_invalid_inputs_are_rejected(self):
        for scores, k in [({}, 1), (SCORES, 0), ({'a': {'r': 1.0}, 'b': {'s': 1.0}}, 1), ({'a': {'r': float('nan')}}, 1)]:
            with self.assertRaises(InputError):
                merge_max(scores, k)
        with self.assertRaises(InputError):
            merge_round_robin(SCORES, ['g'], 2)


if __name__ == '__main__':
    unittest.main()
