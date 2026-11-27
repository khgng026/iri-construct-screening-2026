import unittest

from experiments.benign_protein_pilot.fold_study import _admitted_folds, _fold_scores, comparisons


class FoldStudyTest(unittest.TestCase):
    def test_requires_two_competing_superfamilies(self):
        single = [{'family': family} for family in ('a.118.1.1', 'a.118.1.2')]
        competitor = [{'family': family} for family in ('a.118.8.1', 'a.118.8.7')]
        self.assertEqual(_admitted_folds(single), {})
        splits = _admitted_folds(single + competitor)
        for split in set(splits.values()):
            self.assertEqual({family.rsplit('.', 1)[0] for family, value in splits.items() if value != split},
                             {'a.118.1', 'a.118.8'})

    def test_superfamily_weight_is_equal_despite_unequal_domain_and_family_counts(self):
        rows = [{'fold': 'f', 'label': label, 'fine_family': family, 'correct': {'method': correct}}
                for label, family, correct in (('a', 'a1', 1), ('a', 'a1', 1), ('a', 'a2', 1), ('b', 'b1', 0))]
        self.assertEqual(_fold_scores(rows, ('method',)), {'f': {'method': .5}})

    def test_expanded_reference_library_keeps_the_whole_query_family_held_out(self):
        families = ('a.118.1.1', 'a.118.1.1', 'a.118.1.2', 'a.118.8.1', 'a.118.8.7')
        rows = [{'family': family} for family in families]
        splits = _admitted_folds(rows)
        entries = {str(index): dict(row, split=splits[row['family']]) for index, row in enumerate(rows)}
        sequences = {key: 'A' for key in entries}
        _, original, original_refs = next(comparisons(entries, sequences, 'three-way'))
        _, expanded, expanded_refs = next(comparisons(entries, sequences, 'leave-one-family-out'))
        self.assertEqual(set(original), {'0', '1', '3'})
        self.assertEqual(set(expanded), {'0', '1'})
        self.assertLess(original_refs.keys(), expanded_refs.keys())
