"""Ensure hypothetical schema examples remain structurally valid."""
from pathlib import Path
import json
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.schema_check import check_schema_example

class ExampleSchemaTests(unittest.TestCase):
    def test_all_three_fictional_examples(self):
        examples = sorted((ROOT / "data" / "examples").glob("*.json"))
        self.assertEqual(len(examples), 3)
        for example in examples:
            with self.subTest(example=example.name):
                self.assertEqual(check_schema_example(example), [])

    def test_bad_reference_is_rejected(self):
        path = ROOT / "data" / "examples" / "example_single_cds.json"
        doc = json.loads(path.read_text(encoding="utf-8"))
        doc["units"][0]["feature_ids"] = ["missing-id"]
        with tempfile.TemporaryDirectory() as directory:
            f = Path(directory) / "bad.json"
            f.write_text(json.dumps(doc), encoding="utf-8")
            self.assertTrue(any("nonexistent" in err for err in check_schema_example(f)))

    def test_out_of_bounds_is_rejected(self):
        path = ROOT / "data" / "examples" / "example_single_cds.json"
        doc = json.loads(path.read_text(encoding="utf-8"))
        doc["features"][0]["location"]["end0"] = 99999
        with tempfile.TemporaryDirectory() as directory:
            f = Path(directory) / "bad.json"
            f.write_text(json.dumps(doc), encoding="utf-8")
            self.assertTrue(any("location" in err for err in check_schema_example(f)))
