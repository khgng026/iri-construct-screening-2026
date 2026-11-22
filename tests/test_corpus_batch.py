"""Regression tests for the corpus audit coordinate handling.

The audit caches family evidence per protein. Two constructs that carry the same
protein at different coordinates must still receive their own DNA location: the
per-occurrence evidence must not be written back into the shared cache.
"""
import unittest

from experiments.paired_preprocessing.corpus_batch import (
    coding_dna,
    occurrence_anchor_evidence,
    protein_to_spans,
    transcript_positions,
)

COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def _revcomp(dna):
    return dna.translate(COMPLEMENT)[::-1]


class SharedProteinCoordinateTests(unittest.TestCase):
    def setUp(self):
        # One protein, encoded by a fixed codon string, deposited in two different records.
        self.protein = 'M' + 'A' * 49
        self.codons = 'ATG' + 'GCT' * 49
        # The anchor hit covers the whole protein. This dict stands for the protein-level
        # cache entry that both records share; the audit must never mutate it.
        self.cached = {'demo_family': {'anchor': 'DEMO', 'identity': 1.0, 'unit_coverage': 1.0,
                                       'anchor_coverage': 1.0, 'unit_region_aa': [0, 50]}}

    def _occurrence(self, dna, segments):
        positions = transcript_positions(segments, 1)
        coding = coding_dna(dna, segments, 1)
        return occurrence_anchor_evidence(self.cached, self.protein, positions, coding, 11)

    def test_shared_protein_keeps_distinct_per_record_coordinates(self):
        # Record A: forward strand near the start. Record B: minus strand further along.
        dna_a = 'C' * 100 + self.codons + 'C' * 40
        segments_a = [{'start0': 100, 'end0': 250, 'strand': 1}]
        dna_b = 'G' * 300 + _revcomp(self.codons) + 'G' * 20
        segments_b = [{'start0': 300, 'end0': 450, 'strand': -1}]

        evidence_a = self._occurrence(dna_a, segments_a)
        evidence_b = self._occurrence(dna_b, segments_b)

        # Each record translates to the protein, and each keeps its own DNA location.
        self.assertTrue(evidence_a['demo_family']['dna_region_retranslates'])
        self.assertTrue(evidence_b['demo_family']['dna_region_retranslates'])
        self.assertEqual(evidence_a['demo_family']['dna_spans'], [[100, 250]])
        self.assertEqual(evidence_b['demo_family']['dna_spans'], [[300, 450]])
        self.assertNotEqual(evidence_a['demo_family']['dna_spans'], evidence_b['demo_family']['dna_spans'])

        # The regression: computing record B must not reach back and change record A.
        self.assertEqual(evidence_a['demo_family']['dna_spans'], [[100, 250]])

    def test_cache_is_not_mutated(self):
        self._occurrence('C' * 100 + self.codons, [{'start0': 100, 'end0': 250, 'strand': 1}])
        self.assertNotIn('dna_spans', self.cached['demo_family'])
        self.assertNotIn('dna_region_retranslates', self.cached['demo_family'])
        self.assertEqual(self.cached['demo_family']['unit_region_aa'], [0, 50])

    def test_joined_minus_strand_region_spans_are_sorted_and_merged(self):
        # A split region on the minus strand still yields ascending, forward-sequence spans.
        positions = transcript_positions([{'start0': 50, 'end0': 56, 'strand': -1},
                                          {'start0': 10, 'end0': 16, 'strand': -1}], 1)
        self.assertEqual(protein_to_spans(positions, 0, 4), [[10, 16], [50, 56]])


if __name__ == '__main__':
    unittest.main()
