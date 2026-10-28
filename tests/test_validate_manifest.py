"""Unit tests for the metadata-only validation scaffold."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from iri_constructs.validate_manifest import REQUIRED_COLUMNS, validate_manifest

BASE = {
    "record_id": "test-one", "source": "placeholder", "source_url": "https://example.org/test",
    "accession": "FAKE", "accession_version": "1", "archetype": "single_cds",
    "sequence_length_nt": "1000", "feature_count": "2", "annotation_status": "pending",
    "review_status": "pending", "licensing_status": "pending", "safety_review": "pending",
    "retrieved_on": "2026-10-28", "checksum_sha256": "", "notes": "fake metadata",
}

class MetadataValidationTests(unittest.TestCase):
    def check_rows(self, rows, fields=REQUIRED_COLUMNS):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "records.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            return validate_manifest(path)

    def test_example_is_valid(self):
        self.assertEqual(self.check_rows([BASE]), [])

    def test_duplicate_record_id(self):
        self.assertTrue(any("repeated" in e for e in self.check_rows([BASE, BASE])))

    def test_reject_invalid_archetype(self):
        self.assertTrue(any("archetype" in e for e in self.check_rows([{**BASE, "archetype": "other"}])))

    def test_reject_invalid_number(self):
        self.assertTrue(any("feature_count" in e for e in self.check_rows([{**BASE, "feature_count": "-1"}])))

    def test_missing_column(self):
        fields = [f for f in REQUIRED_COLUMNS if f != "review_status"]
        self.assertTrue(any("Missing columns" in e for e in self.check_rows([], fields)))

    def test_invalid_date(self):
        self.assertTrue(any("retrieved_on" in e for e in self.check_rows([{**BASE, "retrieved_on": "next Tuesday"}])))

if __name__ == "__main__":
    unittest.main()
