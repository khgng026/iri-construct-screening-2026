"""Cache integrity regressions; no pretrained execution in routine tests."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
spec = importlib.util.spec_from_file_location('cleared', Path(__file__).resolve().parents[1] / 'scripts/execute_cleared_retrieval.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class ClearedCacheTests(unittest.TestCase):

    def test_saved_score_or_provenance_tampering_blocks_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            names = ['report.json', 'execution-freeze.json', 'esm2-scores.json']
            for name in names:
                (root / name).write_text('{}')
            (root / 'cache-checksums.json').write_text(json.dumps({name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in names}))
            module.verify_score_cache(root)
            for name in names:
                (root / name).write_text('{"changed": true}')
                with self.assertRaises(module.InputError):
                    module.verify_score_cache(root)
                (root / name).write_text('{}')
