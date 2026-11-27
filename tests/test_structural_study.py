import unittest

from experiments.benign_protein_pilot.reaction_study import _evidence
from experiments.benign_protein_pilot.structural_study import _identity_tier, _top_label, family_folds, superfamily


class StructuralStudyTest(unittest.TestCase):
    def test_split_leaves_out_whole_families_and_preserves_each_label(self):
        rows = [{'family': family} for family in (
            'c.1.2.1', 'c.1.2.1', 'c.1.2.3', 'c.1.2.4', 'c.1.2.6', 'c.1.10.1', 'c.1.10.2', 'c.1.11.1',
        )]
        folds = family_folds(rows)
        self.assertEqual(folds, family_folds(reversed(rows)))
        self.assertNotIn('c.1.11.1', folds)
        for fold in set(folds.values()):
            held_out = {family for family, value in folds.items() if value == fold}
            references = folds.keys() - held_out
            self.assertEqual({superfamily(family) for family in references}, {'c.1.2', 'c.1.10'})

    def test_repeated_references_agree_on_one_label_but_cross_label_ties_abstain(self):
        labels = {'r1': 'a', 'r2': 'a', 'r3': 'b'}
        self.assertEqual(_top_label({'r1': 2, 'r2': 2, 'r3': 1}, labels), 'a')
        self.assertIsNone(_top_label({'r1': 2, 'r2': 2, 'r3': 2}, labels))
        self.assertIsNone(_top_label({}, labels))

    def test_low_identity_boundary_is_strict(self):
        self.assertEqual(_identity_tier(0.24999), 'under_25')
        self.assertEqual(_identity_tier(0.25), '25_to_under_30')
        self.assertEqual(_identity_tier(0.30), '30_to_under_40')

    def test_reaction_evidence_requires_experimental_pubmed_support_for_the_exact_label(self):
        labelled_reaction = {'commentType': 'CATALYTIC ACTIVITY', 'reaction': {
            'ecNumber': '3.2.1.1', 'evidences': [
                {'evidenceCode': 'ECO:0000250', 'source': 'PubMed', 'id': '123'},
                {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '456'},
            ],
        }}
        self.assertEqual(_evidence({'comments': [labelled_reaction]}, '3.2.1.1'), ['456'])
        self.assertEqual(_evidence({'comments': [labelled_reaction]}, '3.2.1.2'), [])
        second_reaction = {'commentType': 'CATALYTIC ACTIVITY', 'reaction': {'ecNumber': '3.2.1.2'}}
        self.assertEqual(_evidence({'comments': [labelled_reaction, second_reaction]}, '3.2.1.1'), [])
