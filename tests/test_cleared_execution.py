"""Cache integrity regressions; no pretrained execution in routine tests."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
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
            checksum = lambda name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            original = {name: checksum(name) for name in names}
            anchor = {'execution_freeze_sha256': original['execution-freeze.json'],
                      'local_report_sha256': original['report.json'],
                      'reproducibility': {'canonical_score_matrix_sha256': {
                          'esm2-scores.json': original['esm2-scores.json']}}}
            for name in names:
                with self.subTest(changed=name):
                    (root / name).write_text('{"changed": true}')
                    # Coordinated mutation also rewrites the adjacent checksums.
                    (root / 'cache-checksums.json').write_text(json.dumps({n: checksum(n) for n in names}))
                    with patch.object(module, 'immutable_e3', return_value=anchor):
                        with self.assertRaisesRegex(module.InputError, 'immutable'):
                            module.verify_score_cache(root)
                    (root / name).write_text('{}')

    def test_mandatory_gates_survive_optimized_python(self):
        import subprocess, sys
        script = str(Path(__file__).resolve().parents[1]/'scripts/execute_cleared_retrieval.py')
        checks = [
            "m['validate_authorization']({'review_comment_id':6800553713}, 'E3')",
            "m['freeze_model']({'resolved_revision':'wrong','weight_sha256':'wrong'}, __import__('pathlib').Path('unused'))",
            """import tempfile,json
from pathlib import Path
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp); (root/'docs/week2').mkdir(parents=True)
    (root/'source').write_text('original')
    f={'artifacts':{'source':{'sha256':m['digest_file'](root/'source')}}}
    f['freeze_sha256']=m['_hash'](f)
    (root/'docs/week2/review-freeze.json').write_text(json.dumps(f))
    m['verify_original_freeze'].__globals__['ROOT']=root
    m['verify_original_freeze']()
    (root/'source').write_text('changed')
    m['verify_original_freeze']()
""",
        ]
        for check in checks:
            with self.subTest(check=check.splitlines()[0]):
                code = "import runpy; m=runpy.run_path(%r); " % script + check
                result = subprocess.run([sys.executable, '-O', '-c', code], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('InputError', result.stderr)

    def test_exploratory_balanced_merge_retains_shared_budget_and_component_origin(self):
        spec = importlib.util.spec_from_file_location('analysis', Path(__file__).resolve().parents[1] / 'scripts/analyze_cleared_pilot.py')
        analysis = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(analysis)
        scores = {'a': {'r1': 0.8, 'r2': 0.1, 'r3': 0.2},
                  'b': {'r1': 90, 'r2': 100, 'r3': 99}}
        families = {'a':'A', 'b':'B'}
        refs = {'r1':'A', 'r2':'B', 'r3':'B'}
        hits = analysis.attributed_hits(scores, families, refs, k=2, balanced=True)
        self.assertEqual([(h['reference_id'],h['winning_component']) for h in hits], [('r1','a'),('r2','b')])
        self.assertTrue(all(h['component_attributed_correct'] for h in hits))
        self.assertEqual(hits, analysis.attributed_hits(dict(reversed(list(scores.items()))), families, refs, k=2, balanced=True))
        all_hits = analysis.attributed_hits(scores, families, refs, k=5, balanced=True)
        self.assertEqual(len(all_hits), 3)  # one shared catalogue, no duplicate budget credit
        with self.assertRaises(module.InputError):
            analysis.attributed_hits({'a':{'r1':float('nan')}}, families, refs)
