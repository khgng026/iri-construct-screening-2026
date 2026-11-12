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
            # Empty fictional cache must never pass the immutable experimental anchor.
            with self.assertRaises((module.InputError, KeyError)):
                module.verify_score_cache(root)
            for name in names:
                (root / name).write_text('{"changed": true}')
                with self.assertRaises(module.InputError):
                    module.verify_score_cache(root)
                (root / name).write_text('{}')

    def test_mandatory_gates_survive_optimized_python(self):
        import subprocess, sys
        script = str(Path(__file__).resolve().parents[1]/'scripts/execute_cleared_retrieval.py')
        code = "import runpy; m=runpy.run_path(%r); m['validate_authorization']({'review_comment_id':6800553713}, 'E3')" % script
        result = subprocess.run([sys.executable, '-O', '-c', code], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('InputError', result.stderr)
        with self.assertRaises(module.InputError):
            module.freeze_model({'resolved_revision':'wrong','weight_sha256':'wrong'},Path('unused'))
