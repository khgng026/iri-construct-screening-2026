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
            "m['freeze_model']({'resolved_revision':m['ESM2_REVISION'],'weight_sha256':'wrong'}, __import__('pathlib').Path('unused'))",
            "m['validate_authorization']({'review_comment_id':6800553713,'review_comment_url':m['URL'],'review_freeze_sha256':'c90d063f7f17c9f929af2fddb4e3fe9bc719c9a4cae5c36f66b238f968e47c2e','decisions':{'E3':'HOLD'}}, 'E3')",
            """import tempfile,json
from pathlib import Path
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp); (root/'docs/week2').mkdir(parents=True)
    (root/'source').write_text('original')
    f={'artifacts':{'source':{'sha256':m['digest_file'](root/'source')}}}
    f['freeze_sha256']=m['_hash'](f)
    (root/'docs/week2/review-freeze.json').write_text(json.dumps(f))
    m['verify_original_freeze'].__globals__['ROOT']=root
    m['verify_original_freeze'].__globals__['EXPECTED_REVIEW_FREEZE_SHA256']=f['freeze_sha256']
    d={'authority_freeze_sha256':f['freeze_sha256'],'stages':{'E2':['source']},'artifacts':{'source':f['artifacts']['source']['sha256']}}
    (root/'dependencies.json').write_text(json.dumps(d))
    g=m['verify_stage_dependencies'].__globals__
    g['DEPENDENCIES_PATH']='dependencies.json'
    g['EXPECTED_DEPENDENCIES_SHA256']=m['digest_file'](root/'dependencies.json')
    m['verify_stage_dependencies']('E2')
    (root/'source').write_text('changed')
    m['verify_stage_dependencies']('E2')
""",
        ]
        for check in checks:
            with self.subTest(check=check.splitlines()[0]):
                code = "import runpy; m=runpy.run_path(%r); " % script + check
                result = subprocess.run([sys.executable, '-O', '-c', code], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('InputError', result.stderr)

        # Keep both acceptance and rejection paths covered under -O in routine CI.
        cases = [
            'ClearedCacheTests.test_stage_dependencies_ignore_held_assets_but_require_live_inputs',
            'ClearedCacheTests.test_normalized_alignment_mutation_rejected_against_verified_inputs',
            'ClearedCacheTests.test_saved_score_or_provenance_tampering_blocks_reuse',
        ]
        code = "import sys, unittest; sys.path.insert(0, %r); unittest.main(module='test_cleared_execution', argv=%r)" % (
            str(Path(__file__).resolve().parent), ['optimized-regressions', *cases])
        flags = ['-O', *(['-S'] if sys.flags.no_site else [])]
        result = subprocess.run([sys.executable, *flags, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

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

    def test_stage_dependencies_ignore_held_assets_but_require_live_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = ['source', 'reference', 'config', 'esm2-weight', 'E3-source']
            for name in paths:
                (root/name).write_text(name)
            artifacts = {name: {'sha256': module.digest_file(root/name)} for name in paths}
            artifacts['held-esmc-weight'] = {'sha256':'held'}
            artifacts['held-igem-original'] = {'sha256':'held'}
            freeze = {'artifacts':artifacts}
            freeze['freeze_sha256'] = module._hash(freeze)
            (root/'docs/week2').mkdir(parents=True)
            (root/'docs/week2/review-freeze.json').write_text(json.dumps(freeze))
            deps = {'authority_freeze_sha256':freeze['freeze_sha256'],
                    'stages':{'E2':['source','reference','config','esm2-weight'],
                              'E2-replay':['source','reference','config'],
                              'E3':paths},
                    'artifacts':{name:v['sha256'] for name,v in artifacts.items() if name in paths}}
            contract = root/'dependencies.json'
            contract.write_text(json.dumps(deps))
            with patch.multiple(module, ROOT=root, DEPENDENCIES_PATH='dependencies.json',
                                EXPECTED_REVIEW_FREEZE_SHA256=freeze['freeze_sha256'],
                                EXPECTED_DEPENDENCIES_SHA256=module.digest_file(contract)):
                for stage in ('E2','E3'):
                    module.verify_stage_dependencies(stage)  # held originals never existed
                (root/'E3-source').unlink()
                module.verify_stage_dependencies('E2')
                module.verify_stage_dependencies('E2', reuse_e3=True)
                with self.assertRaisesRegex(module.InputError, 'Missing stage dependency: E3-source'):
                    module.verify_stage_dependencies('E3')
                (root/'E3-source').write_text('E3-source')
                (root/'esm2-weight').unlink()
                module.verify_stage_dependencies('E2', reuse_e3=True)
                with self.assertRaisesRegex(module.InputError, 'Missing stage dependency'):
                    module.verify_stage_dependencies('E2')
                for name in ('source','reference','config'):
                    with self.subTest(required=name):
                        (root/name).write_text('changed')
                        with self.assertRaisesRegex(module.InputError, 'Changed stage dependency'):
                            module.verify_stage_dependencies('E2', reuse_e3=True)
                        (root/name).unlink()
                        with self.assertRaisesRegex(module.InputError, 'Missing stage dependency'):
                            module.verify_stage_dependencies('E2', reuse_e3=True)
                        (root/name).write_text(name)
                contract.write_text('{}')
                with self.assertRaisesRegex(module.InputError, 'Changed stage dependency manifest'):
                    module.verify_stage_dependencies('E2', reuse_e3=True)
        actual = json.loads((module.ROOT/module.DEPENDENCIES_PATH).read_text())
        for names in actual['stages'].values():
            self.assertFalse(any('esmc' in n.lower() or 'BBa_' in n for n in names))
        self.assertFalse(any('model-cache' in n for n in actual['stages']['E2-replay']))
        for stage, count in [('E2',1), ('E2-replay',1), ('E3',3)]:
            self.assertEqual(sum('/raw_sources/NCBI_' in n for n in actual['stages'][stage]), count)

    def test_normalized_alignment_mutation_rejected_against_verified_inputs(self):
        try:
            import Bio
        except ImportError:
            self.skipTest('optional Biopython alignment dependency')
        spec = importlib.util.spec_from_file_location('analysis', Path(__file__).resolve().parents[1] / 'scripts/analyze_cleared_pilot.py')
        analysis = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(analysis)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root/'R.json'
            ref.write_text(json.dumps({'sequence':{'value':'AA'}}))
            raw = {'q':{'uniprot-R':4.0}}
            anchor = {'reproducibility':{'canonical_score_matrix_sha256':{'alignment-scores.json':module._hash(raw)}},
                      'panel':{'input_fingerprints':{'uniprot-R':{'source_sha256':module.digest_file(ref),
                           'sequence_sha256':hashlib.sha256(b'AA').hexdigest()}}}}
            cache = {'raw_scores':raw, 'normalized_scores':{'q':{'uniprot-R':0.5}}}
            self.assertEqual(analysis.verified_normalized_alignment(cache, anchor, root), cache['normalized_scores'])
            cache['normalized_scores']['q']['uniprot-R'] = 0.9
            with self.assertRaisesRegex(module.InputError, 'normalized alignment scores'):
                analysis.verified_normalized_alignment(cache, anchor, root)
            cache['normalized_scores']['q']['uniprot-R'] = 0.5
            ref.write_text(json.dumps({'sequence':{'value':'CC'}}))
            with self.assertRaisesRegex(module.InputError, 'reference source'):
                analysis.verified_normalized_alignment(cache, anchor, root)
