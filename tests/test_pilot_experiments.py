"""Numerical adapter fixtures exercise orchestration; no biological/model findings."""
import json
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.pilot_experiments import run_pilot_experiment
from iri_constructs.pilot_scoring import freeze_panel
from iri_constructs.runtime import InputError


class PilotExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'results')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.output = self.base / 'run'
        self.inputs_path = self.base / 'inputs.json'
        self.inputs_path.write_text(json.dumps({'records':[{'fictional_test':True}]}),encoding='utf-8')
        self.proteins = {k: ProteinInput(k,'parent-'+k,residue*3,'fictional-software-fixture')
                         for k,residue in zip('abcd','ACDE')}
        self.labels = [dict(input_id=k, family='one' if k in 'ab' else 'two',
            annotation_status='known',review_status='approved',evidence='Fictional software test only',
            independence_group='group-'+k) for k in 'abcd']
        self.inputs = {'families':self.labels}
        self.provenance = {'sources':[dict(record_id=p.construct_id,archetype='single_cds') for p in self.proteins.values()]}
        self.closed = []

    def protocol(self, experiment='D', top_k=1):
        config = json.loads((ROOT/'config'/'pilot'/('experiment_'+experiment+'.json')).read_text())
        config.update(top_k=top_k,retrieval_budget=top_k,bootstrap_replicates=100,min_independent_clusters=2)
        path = self.base/'protocol.json'
        path.write_text(json.dumps(config),encoding='utf-8')
        return path

    def factory(self, name, device, *args, **kwargs):
        output, closed = self.output, self.closed
        class NumericalAdapter:
            metadata = {'model':name,'fictional_numerical_adapter':True,'resolved_revision':'fixture'}
            def residues(self, sequence):
                assert (output/'protocol-freeze.json').exists()
                assert (output/'panel-freeze.json').exists()
                assert (output/('model-freeze-'+name+'.json')).exists()
                vector = [1.,0.] if sequence[0] in 'ACF' else [0.,1.]
                return [vector]*len(sequence)
            def close(self):
                closed.append(name)
        adapter = NumericalAdapter()
        adapter.device = device
        return adapter

    def run_numbers(self, experiment='D', factory=None, top_k=1):
        with patch('iri_constructs.pilot_experiments.load_proteins',return_value=(self.inputs,self.proteins,self.provenance,[])), \
             patch('iri_constructs.pilot_experiments.prepare_device',return_value='cpu'):
            return run_pilot_experiment(self.protocol(experiment,top_k),self.inputs_path,self.output,
                                        adapter_factory=factory or self.factory)

    def test_freeze_precedes_every_forward_and_same_panel_shared(self):
        result = self.run_numbers()
        self.assertEqual(result['status'],'ok')
        self.assertEqual(result['common_scored_query_ids'],list('abcd'))
        self.assertEqual(self.closed,['esm2','esmc'])
        self.assertTrue(all(m['scored_queries']==4 for m in result['methods'].values()))
        self.assertEqual(result['models']['esm2']['unique_forward_calls'],4)
        self.assertEqual(len({m['panel_sha256'] for m in result['methods'].values()}),1)

    def test_pending_empty_inputs_never_load_model(self):
        self.inputs_path.write_text(json.dumps({'schema_version':'t05-inputs-1','records':[]}),encoding='utf-8')
        def forbidden(*args,**kwargs):
            self.fail('No weights or records may load without admission')
        result = run_pilot_experiment(self.protocol(),self.inputs_path,self.output,adapter_factory=forbidden)
        self.assertEqual(result['status'],'blocked_no_reviewed_inputs')
        self.assertFalse(result['inference_executed'])
        self.assertFalse(result['scientific_scores_produced'])
        self.assertTrue((self.output/'protocol-freeze.json').exists())

    def test_insufficient_independent_panel_never_loads_model(self):
        for entry in self.labels:
            entry['family'] = 'one'
        def forbidden(*args,**kwargs):
            self.fail('Insufficient panel must not load a model')
        result = self.run_numbers(factory=forbidden)
        self.assertEqual(result['status'],'blocked_insufficient_panel')
        self.assertFalse(result['inference_executed'])

    def test_missing_model_is_blocked_and_sensitive_error_is_omitted(self):
        def factory(name,*args,**kwargs):
            if name=='esmc':
                raise ModuleNotFoundError('SENSITIVE_DETAIL')
            return self.factory(name,*args,**kwargs)
        result = self.run_numbers(factory=factory)
        self.assertEqual(result['status'],'partial')
        self.assertEqual(result['methods']['esmc']['status'],'blocked')
        self.assertEqual(result['common_scored_query_ids'],[])
        self.assertNotIn('SENSITIVE_DETAIL',json.dumps(result))

    def test_all_failed_inputs_do_not_claim_successful_inference(self):
        def factory(name,*args,**kwargs):
            adapter = self.factory(name,*args,**kwargs)
            def fail(sequence):
                raise InputError('Fictional cap incompatibility before a forward pass')
            adapter.residues = fail
            return adapter
        result = self.run_numbers(factory=factory)
        self.assertFalse(result['inference_executed'])
        self.assertFalse(result['scientific_scores_produced'])
        self.assertTrue(result['models']['esm2']['inference_attempted'])
        self.assertEqual(result['models']['esm2']['unique_forward_calls'],0)

    def test_budget_saturated_panel_is_descriptive_not_method_advantage(self):
        result = self.run_numbers(top_k=5)
        self.assertEqual(result['status'],'descriptive_only_budget_saturated_panel')
        self.assertEqual(result['paired_comparisons'][0]['paired_query_ids'],[])
        self.assertEqual(result['paired_comparisons'][0]['excluded_budget_saturated_query_ids'],list('abcd'))

    def test_annotation_free_parsing_is_not_scientific_success(self):
        with patch('iri_constructs.pilot.prepare_pilot',return_value={'status':'ok','scored_records':0}):
            result = run_pilot_experiment(self.protocol('A'),self.inputs_path,self.output)
        self.assertEqual(result['status'],'blocked_no_observable_annotations')
        self.assertFalse(result['scientific_scores_produced'])

    def fusion_fixture(self):
        self.proteins['f'] = ProteinInput('f','real-fusion-parent','G'+'A'*3+'C'*3+'T','fictional')
        self.proteins['positive'] = ProteinInput('positive','independent-positive','F'*3,'fictional')
        self.proteins['leak'] = ProteinInput('leak','different-parent',self.proteins['f'].sequence,'fictional')
        for iid in ('positive','leak'):
            self.labels.append(dict(input_id=iid,family='one',annotation_status='known',review_status='approved',
                evidence='Fictional software test',independence_group=iid))
        self.provenance['sources'].append(dict(record_id='real-fusion-parent',archetype='tagged_fusion'))
        self.provenance['sources'].extend([dict(record_id=self.proteins[k].construct_id,archetype='single_cds') for k in ('positive','leak')])
        self.inputs['fusion_cases'] = [dict(case_id='first',fusion_id='f',target_id='a',start0=1,end0=4,
            coordinate_evidence='Fictional software indices only',review_status='approved'),
            dict(case_id='second',fusion_id='f',target_id='b',start0=4,end0=7,
            coordinate_evidence='Fictional software indices only',review_status='approved')]

    def test_fusion_controls_exclude_full_fusion_copies_and_share_parent_cluster(self):
        self.fusion_fixture()
        result = self.run_numbers('B')
        self.assertEqual(result['status'],'ok')
        self.assertNotIn('leak',result['panel']['reference_ids'])
        excluded = next(x for x in result['fusion_reference_exclusions'] if x['input_id']=='leak')
        self.assertTrue(any('complete-fusion' in reason for reason in excluded['reasons']))
        fingerprints = [result['panel']['input_fingerprints'][q] for q in result['panel']['query_ids']]
        self.assertEqual({f['construct_id'] for f in fingerprints},{'real-fusion-parent'})
        self.assertEqual({f['archetype'] for f in fingerprints},{'tagged_fusion'})
        self.assertEqual(len({f['independence_group'] for f in fingerprints}),1)
        self.assertEqual(result['paired_comparisons'][0]['independent_clusters'],1)
        self.assertEqual(set(result['methods']['esm2/isolated_target']['per_construct']),{'real-fusion-parent'})

    def test_query_roles_are_lists_not_string_membership(self):
        self.fusion_fixture()
        self.inputs['query_ids'] = 'ab'
        result = self.run_numbers('B')
        self.assertEqual(result['status'],'blocked')
        self.assertIn('query_ids',result['reason'])

    def test_existing_output_refuses_overwrite(self):
        self.output.mkdir()
        (self.output/'report.json').write_text('{}',encoding='utf-8')
        with self.assertRaises(InputError):
            self.run_numbers()

    @unittest.skipUnless(importlib.util.find_spec('Bio'),'Optional Biopython alignment integration')
    def test_alignment_and_numerical_embeddings_share_panel_and_budget(self):
        result = self.run_numbers('C')
        self.assertEqual(result['status'],'ok')
        self.assertEqual(set(result['methods']),{'alignment_raw','alignment_reference_normalized','esm2'})
        self.assertTrue(all(m['overall']['returned_hits']==4 for m in result['methods'].values()))
        self.assertEqual(len({m['panel_sha256'] for m in result['methods'].values()}),1)
