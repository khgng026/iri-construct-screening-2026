import unittest

from experiments.benign_protein_pilot.fold_study import _admitted_folds, _fold_scores


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
