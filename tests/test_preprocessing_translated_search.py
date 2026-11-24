import importlib.util
import unittest

from experiments.preprocessing_exploration.ranking import Candidate, rank_indices
from experiments.preprocessing_exploration.translated_search import Hit, Region, evaluate, parse_hits


class TranslatedSearchTest(unittest.TestCase):
    def test_mixed_ranking_keeps_first_and_prefers_length_in_remaining_slots(self):
        pool = [Candidate(1, ((0, length * 3),), length) for length in (800, 600, 250, 310)]
        self.assertEqual(rank_indices(pool, "keep_first_then_250", limit=2), [0, 2])

    def test_reverse_coordinate_mapping_preserves_frame(self):
        candidate = Region(100, 190, -1)
        hits = parse_hits("q1\tr1\t2\t9\t80\t1e-12", {"q1": ("p", 0, candidate)})
        self.assertEqual(hits["p"][0].region, Region(163, 187, -1))
        self.assertEqual(hits["p"][0].region.frame, candidate.frame)

    def test_hits_require_consistent_frame_and_selected_candidate(self):
        target = Region(10, 100, 1)
        candidates = {"p": [target, target]}
        hits = {"p": [
            Hit(Region(11, 101, 1), 150, "r1", 0),
            Hit(target, 88, "r2", 1),
        ]}
        truth = {"p": ("corroborated_gfp_like", target)}
        first = evaluate(truth, hits, candidates, 1)["reports"]["all"]
        second = evaluate(truth, hits, candidates, 2)["reports"]["all"]
        self.assertEqual(first["positives_recovered_at_80pct"], 0)
        self.assertEqual(second["positives_recovered_at_80pct"], 1)

    @unittest.skipUnless(
        importlib.util.find_spec("Bio") and importlib.util.find_spec("pyrodigal"),
        "archive coverage dependencies unavailable",
    )
    def test_resolved_segments_retain_coordinates_after_unknown(self):
        from experiments.preprocessing_exploration.archive_coverage import resolved_translation_segments

        dna = "NNN" + "GCT" * 30 + "TAA"
        candidates = resolved_translation_segments(dna)
        self.assertTrue(any(row["strand"] == 1 and row["spans"] == [[3, 93]]
                            and row["protein"] == "A" * 30 for row in candidates))
