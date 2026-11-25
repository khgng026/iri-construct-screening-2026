"""Checks for candidate accounting in the archived preprocessing validation."""

from importlib.util import find_spec
from types import SimpleNamespace
import unittest

from experiments.preprocessing_exploration.hybrid import Candidate, distinct_supplements


class FallbackTest(unittest.TestCase):
    def test_deduplicates_fallbacks_without_removing_primary_or_other_frames(self):
        primary = Candidate("gene", 0, 1, 0, 300, 100)
        redundant = Candidate("fallback", 0, 1, 30, 330, 100)
        distinct = Candidate("fallback", 1, 1, 600, 900, 100)
        repeated = Candidate("fallback", 2, 1, 630, 930, 100)
        other_frame = Candidate("fallback", 3, 1, 631, 931, 100)
        self.assertEqual(list(distinct_supplements(
            [primary], [redundant, distinct, repeated, other_frame], 0.8,
        )), [distinct, other_frame])


@unittest.skipUnless(find_spec("Bio") and find_spec("pyrodigal"), "archive dependencies unavailable")
class WindowValidationTest(unittest.TestCase):
    def test_balanced_windows_keep_the_terminal_residues_on_both_strands(self):
        from experiments.paired_preprocessing.reporter_pilot import MAX_RESIDUES
        from experiments.preprocessing_exploration.validation import called_genes, grouped_windows, measure, select

        protein = "A" * MAX_RESIDUES + "W"
        coding_bases = len(protein) * 3
        for strand in (-1, 1):
            with self.subTest(strand=strand):
                gene = SimpleNamespace(
                    begin=1, end=coding_bases, strand=strand,
                    translate=lambda *, include_stop: protein,
                )
                finder = SimpleNamespace(find_genes=lambda sequence: [gene])
                target = {"strand": strand, "spans": [(0, coding_bases)], "codon_offset": 0,
                          "protein": protein}
                dna = "A" * coding_bases
                parents = called_genes(dna, finder)
                standard, = grouped_windows(parents, balanced=False)
                balanced, = grouped_windows(parents, balanced=True)
                self.assertEqual(len(standard), 1)
                self.assertEqual(len(balanced), 2)
                self.assertEqual(measure(standard, target)["peptide_available"], 0)
                self.assertEqual(measure(balanced, target)["long_complete_across_windows"], 1)
                chosen, cost = select(balanced, 2, MAX_RESIDUES)
                self.assertLessEqual(cost, MAX_RESIDUES)
                self.assertEqual(measure(chosen, target)["peptide_available"], 0)

    def test_mixed_guard_keeps_only_the_first_parents_windows_together(self):
        from experiments.preprocessing_exploration.validation import preserve_first_parent

        first = {"strand": -1, "spans": [[300, 600]]}
        sibling = {"strand": -1, "spans": [[0, 300]]}
        other = {"strand": 1, "spans": [[800, 1400]]}
        self.assertEqual(preserve_first_parent([first, other, sibling], [[first, sibling], [other]]),
                         [first, sibling, other])
