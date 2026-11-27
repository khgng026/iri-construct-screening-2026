import unittest

from experiments.preprocessing_exploration.frame_audit import boundary_metrics


class FrameAuditTests(unittest.TestCase):
    def test_exact_span_uses_coding_start_on_each_strand(self):
        for strand, coding in ((1, (9, 99)), (-1, (12, 102))):
            with self.subTest(strand=strand):
                result = boundary_metrics((strand, [(9, 102)]), [(strand, [coding])])
                self.assertTrue(all(result.values()))

    def test_wrong_frame_can_cover_the_locus(self):
        for strand in (1, -1):
            with self.subTest(strand=strand):
                result = boundary_metrics((strand, [(9, 102)]), [(strand, [(10, 100)])])
                self.assertTrue(result["overlap"])
                self.assertFalse(result["in_frame_overlap"])
                self.assertFalse(result["contains_coding_span"])

    def test_longer_in_frame_candidate_can_be_useful(self):
        result = boundary_metrics((1, [(9, 102)]), [(1, [(3, 111)])])
        self.assertTrue(result["contains_coding_span"])
        self.assertFalse(result["exact_coding_span"])

    def test_stop_normalized_exact_is_separate_from_untrimmed_overlap(self):
        result = boundary_metrics((-1, [(9, 15)]), [(-1, [(12, 15)])])
        self.assertFalse(result["overlap"])
        self.assertTrue(result["exact_coding_span"])
