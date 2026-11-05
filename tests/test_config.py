"""Verify reproducible standard-library configuration rules."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.config import load_config

class ConfigTests(unittest.TestCase):
    def test_default_configuration_loads(self):
        cfg = load_config(ROOT / "config" / "default.json")
        self.assertEqual(cfg["log_level"], "INFO")
        self.assertIn("schema_examples", cfg)

    def test_unrecognized_log_level_fails(self):
        cfg = load_config(ROOT / "config" / "default.json")
        cfg["log_level"] = "VERBOSE"
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "bad.json"
            p.write_text(json.dumps(cfg), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_config(p)

    def test_invalid_json_and_missing_settings(self):
        for text in ('[]', '{invalid-json', '{}', '{"project": false}'):
            with self.subTest(text=text), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'config.json'
                path.write_text(text, encoding='utf-8')
                with self.assertRaises(ValueError):
                    load_config(path)
