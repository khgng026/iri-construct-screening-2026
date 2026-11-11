"""Fixture validation rejects unusable numerical output; no model download."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('validate_model_runtime',
    Path(__file__).resolve().parents[1] / 'scripts' / 'validate_model_runtime.py')
runtime = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runtime)
from iri_constructs.model_adapters import _validate_snapshot_revision


class ModelRuntimeValidationTests(unittest.TestCase):
    def test_pooling_is_per_residue_and_finite(self):
        result = runtime.representation_contract([[1.0, 4.0], [3.0, 8.0]], 2)
        self.assertEqual(result['residue_count'], 2)
        self.assertEqual(result['embedding_width'], 2)
        self.assertTrue(result['mean_pool_verified'])
        self.assertEqual(result, runtime.representation_contract([[1.0, 4.0], [3.0, 8.0]], 2))

    def test_bad_model_shapes_and_nonfinite_fail(self):
        for rows in [[], [[]], [[1.0], [2.0, 3.0]], [[float('nan')], [1.0]], [[float('inf')], [1.0]]]:
            with self.subTest(rows=rows), self.assertRaises(runtime.InputError):
                runtime.representation_contract(rows, 2)

    def test_cli_refuses_scientific_inputs_and_unsupported_models(self):
        with self.assertRaises(SystemExit):
            runtime.main(['--output', 'results/test-models', '--inputs', 'anything.fasta'])
        for name in ['esm3', 'esm2,esm2', 'anything']:
            with self.subTest(model=name), self.assertRaises(SystemExit):
                runtime.main(['--output', 'results/test-models', '--models', name])

    def test_missing_model_is_blocked_and_source_free(self):
        with patch.object(runtime, 'prepare_device', side_effect=runtime.InputError('CPU runtime unavailable')):
            result = runtime.check_model('esm2', runtime.ESM2_REVISION, True, 2, 0)
        self.assertEqual(result['status'], 'blocked')
        self.assertFalse(result['biological_inference_executed'])
        self.assertEqual(result['fixture_forward_calls'], 0)

    def test_pinned_revision_mismatch_and_mutable_snapshot_are_rejected(self):
        with self.assertRaises(runtime.InputError):
            _validate_snapshot_revision(Path('/cache/snapshots/' + 'b' * 40), 'a' * 40)
        for suffix in ['main', 'not-an-immutable-revision']:
            with self.subTest(suffix=suffix), self.assertRaises(runtime.InputError):
                _validate_snapshot_revision(Path('/cache/snapshots/' + suffix), 'main')
        _validate_snapshot_revision(Path('/cache/snapshots/' + 'a' * 40), 'main')


if __name__ == '__main__':
    unittest.main()
