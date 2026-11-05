"""Orchestration tests use numerical adapters, never pretrained-model measurements."""
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.experiment_cli import execute
from iri_constructs.runtime import InputError


class ComparisonRunnerTests(unittest.TestCase):
    def setUp(self):
        self.proteins = {k: ProteinInput(k,k,letter*3,'fictional') for k,letter in zip('abcd','ACDE')}
        self.config = {'families': [dict(input_id=k, family='one' if k in 'ab' else 'two',
            independence_group=k, evidence='Fictional unit test', annotation_status='known', review_status='approved') for k in 'abcd']}
        self.provenance = {'sources': [dict(record_id=k, archetype='single_cds') for k in 'abcd']}
        self.args = SimpleNamespace(inputs=Path('never-read'), experiment='comparison', validate_only=False,
            device='cpu', seed=0, models=['esm2','esmc'], esm2_revision='test', esmc_revision='test',
            max_residues=1022, allow_download=False)
        self.closed = []

    def run_mock(self, factory):
        report = {'limitations': []}
        with patch('iri_constructs.experiment_cli.load_proteins', return_value=(self.config,self.proteins,self.provenance,[])), \
             patch('iri_constructs.experiment_cli.prepare_device', return_value='cpu'):
            code = execute(self.args, report, factory)
        return code, report

    def factory(self, name, device, *args, **kwargs):
        closed = self.closed
        class Numbers:
            metadata = {'model': name, 'fictional_test_adapter': True}
            def residues(self, sequence):
                if name == 'esmc' and sequence.startswith('E'):
                    raise InputError('Fictional inference failure')
                value = [1., 0.] if sequence[0] in 'AC' else [0., 1.]
                return [value] * len(sequence)
            def close(self):
                closed.append(name)
        return Numbers()

    def test_common_panel_exact_ids_and_failure_coverage(self):
        code, report = self.run_mock(self.factory)
        self.assertEqual(code, 1)
        self.assertEqual(report['status'], 'partial')
        self.assertEqual(report['common_panel_ids'], ['a','b','c'])
        self.assertEqual(report['models']['esmc']['coverage']['embedded_units'], 3)
        self.assertEqual(report['models']['esmc']['retrieval_full_panel']['scored'], 0)
        panels = report['common_panel_retrieval']
        self.assertEqual(panels['esm2']['query_ids'], panels['esmc']['query_ids'])
        self.assertEqual(panels['esm2']['reference_ids'], panels['esmc']['reference_ids'])
        self.assertEqual(self.closed, ['esm2','esmc'])

    def test_model_load_failure_blocks_common_panel(self):
        def factory(name, *args, **kwargs):
            if name == 'esmc':
                raise ModuleNotFoundError('Secret or source content must not be saved')
            return self.factory(name,*args,**kwargs)
        code, report = self.run_mock(factory)
        self.assertEqual(report['common_panel_ids'], [])
        self.assertEqual(report['models']['esmc']['status'], 'blocked')
        self.assertNotIn('Secret', str(report))

    def test_validate_only_never_loads_models(self):
        self.args.validate_only = True
        def forbidden(*args, **kwargs):
            self.fail('Validation must not load a model')
        code, report = self.run_mock(forbidden)
        self.assertEqual(code, 0)
        self.assertEqual(report['status'], 'validated_only')
        self.assertFalse(report['inference_executed'])
