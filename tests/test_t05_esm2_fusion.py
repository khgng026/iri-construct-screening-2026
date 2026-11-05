"""Model-free ESM-2 fusion input validation tests."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / "experiments" / "t05_esm2_fusion" / "compare_fusion.py"
SPEC = importlib.util.spec_from_file_location("compare_fusion", PATH)
fusion = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fusion)


class FusionTests(unittest.TestCase):
    def setUp(self):
        self.proteins = {
            "isolated": "ACDE", "fusion": "GGACDETT",
            "reference": "ACDF", "decoy": "WWWW"
        }
        self.row = {
            "pair_id": "fictional_test",
            "isolated_id": "isolated", "fusion_id": "fusion",
            "start_aa": "3", "end_aa": "6",
            "isolated_accession": "fictional-1", "fusion_accession": "fictional-2",
            "isolated_source_url": "https://example.org/1",
            "fusion_source_url": "https://example.org/2",
            "boundary_evidence": "fictional boundary",
            "review_status": "approved", "licensing_status": "approved",
            "safety_review": "approved",
            "reference_id": "reference", "decoy_id": "decoy"
        }

    def test_pair(self):
        p = fusion.validate_pair(self.row, self.proteins)
        self.assertEqual((p["start0"], p["end0"]), (2, 6))

    def test_bad_boundary(self):
        with self.assertRaises(ValueError):
            fusion.validate_pair(dict(self.row, start_aa="2"), self.proteins)

    def test_pending_review(self):
        with self.assertRaises(ValueError):
            fusion.validate_pair(dict(self.row, safety_review="pending"), self.proteins)

    def test_bad_reference(self):
        with self.assertRaises(ValueError):
            fusion.validate_pair(dict(self.row, reference_id="isolated"), self.proteins)

    def test_cosine(self):
        self.assertAlmostEqual(fusion.cosine([1, 0], [1, 0]), 1)
        self.assertAlmostEqual(fusion.cosine([1, 0], [0, 1]), 0)
        with self.assertRaises(ValueError):
            fusion.cosine([0, 0], [1, 0])

    def test_duplicate_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "input.fasta"
            p.write_text(">same\nACDE\n>same\nGGGG\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                fusion.parse_fasta(p)

    def test_invalid_residue(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "input.fasta"
            p.write_text(">test\nAC*E\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                fusion.parse_fasta(p)


if __name__ == "__main__":
    unittest.main()
