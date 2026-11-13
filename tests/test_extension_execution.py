"""Cross-module v1.1 tests use fictional numerical fixtures, never biological findings."""
import copy
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.extension_experiments import main, run_reference_generalization
from iri_constructs.cost_quality import summarize_cost_quality
from iri_constructs.runtime import InputError


class ExtensionExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=ROOT/'results')
        self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name)
        self.output=self.base/'run'
        self.input_path=self.base/'inputs.json'
        self.input_path.write_text(json.dumps({'records':[{'fictional_test':True}],'evidence_kind':'fictional_software_fixture'}),encoding='utf-8')
        ids=['q1','q2']+['r'+str(i) for i in range(8)]
        alphabet='ACDEFGHIKL'
        self.proteins={i:ProteinInput(i,'parent-'+i,aa*(3+n),'f'*64) for n,(i,aa) in enumerate(zip(ids,alphabet))}
        self.families=[dict(input_id=i,family='one' if n%2==0 else 'two',independence_group='group-'+i,
                           annotation_status='known',review_status='approved',evidence='Fictional software metadata only') for n,i in enumerate(ids)]
        self.inputs={'families':self.families,'query_ids':['q1','q2'],'reference_ids':ids[2:],
                     'evidence_kind':'fictional_software_fixture'}
        self.provenance={'sources':[dict(record_id=p.construct_id,archetype='single_cds',review_status='approved',
            licensing_status='approved',safety_review='approved',annotation_status='approved') for p in self.proteins.values()],
            'input_config_sha256':'a'*64,'review_manifest_sha256':'b'*64}
        self.closed=[]

    def factory(self,name,device,*args,**kwargs):
        output,closed=self.output,self.closed
        proteins=self.proteins
        family_by_seq={proteins[e['input_id']].sequence:e['family'] for e in self.families}
        class NumericalAdapter:
            metadata={'model':name,'fictional_software_fixture':True,'resolved_revision':'fictional','weight_sha256':'c'*64}
            def residues(self,sequence):
                assert (output/'reference-tier-freeze.json').exists()
                assert (output/('model-freeze-'+name+'.json')).exists()
                vector=[1.,0.] if family_by_seq[sequence]=='one' else [0.,1.]
                return [vector]*len(sequence)
            def close(self):
                closed.append(name)
        adapter=NumericalAdapter();adapter.device=device
        return adapter

    def alignment(self,proteins,panel,**kwargs):
        self.assertTrue((self.output/'reference-tier-freeze.json').exists())
        labels={e['input_id']:e for e in self.families}
        scores={q:{r:1. if labels[q]['family']==labels[r]['family'] else 0. for r in panel['eligible_reference_ids'][q]} for q in panel['query_ids']}
        return {'status':'fixture','raw_scores':scores,'normalized_scores':copy.deepcopy(scores),
                'configuration':{'fictional_numerical_scores':True},'pairwise_score_calls':16,'reference_self_score_calls':8,
                'resources':{},'failures':[]}

    def run_fixture(self,factory=None):
        with patch('iri_constructs.extension_experiments.reviewed_sources',return_value=({},[SimpleNamespace(review={'annotation_status':'approved'})],{})), \
             patch('iri_constructs.extension_experiments.load_proteins',return_value=(self.inputs,self.proteins,self.provenance,[])), \
             patch('iri_constructs.extension_experiments.prepare_device',return_value='cpu'), \
             patch('iri_constructs.pilot_alignment.alignment_scores',side_effect=self.alignment):
            return run_reference_generalization(self.input_path,ROOT/'config/extensions/E3.json',self.output,
                                                 adapter_factory=factory or self.factory)

    def test_all_tiers_freeze_before_scores_and_missing_strict_evidence_blocks(self):
        result=self.run_fixture()
        self.assertEqual(result['status'],'completed_software_fixture')
        self.assertFalse(result['scientific_scores_produced'])
        self.assertEqual(self.closed,['esm2','esmc'])
        tiers=result['tiers']
        self.assertEqual(tiers['tier3_reviewed_relatedness']['status'],'blocked_relatedness_policy')
        for name in ('tier1_t07_invariants','tier2_additional_linkage'):
            self.assertEqual(tiers[name]['common_scored_query_ids'],['q1','q2'])
            self.assertEqual(len({m['panel_sha256'] for m in tiers[name]['methods'].values()}),1)
        records=json.loads((self.output/'cost-quality-records.json').read_text())['records']
        self.assertTrue(records)
        self.assertTrue(all(r['evidence_kind']=='fictional_software_fixture' for r in records))
        summary=summarize_cost_quality(records)
        self.assertEqual(summary['actual_admitted_research_records'],0)
        self.assertEqual(summary['status'],'descriptive_only_no_admitted_research')

    def test_reviewed_fixture_strict_policy_applies_before_all_method_scores(self):
        policy={'schema_version':'reviewed-relatedness-policy-1',
            'policy_kind':'reviewed_pair_exclusions','policy_scope':'sequence_relatedness',
            'review_status':'approved','declared_before_results':True,'policy_id':'fictional-policy',
            'reviewer':'fictional fixture','reviewed_on':'2026-11-05',
            'evidence':'Numerical software test only, no research admission',
            'excluded_pairs':[{'query_id':'q1','reference_id':'r0','reason':'fictional linkage',
                'evidence':'test metadata','query_sequence_sha256':self.proteins['q1'].sequence_sha256,
                'reference_sequence_sha256':self.proteins['r0'].sequence_sha256}]}
        (self.base/'policy.json').write_text(json.dumps(policy),encoding='utf-8')
        self.inputs['relatedness_policy']='policy.json'
        result=self.run_fixture()
        self.assertEqual(result['status'],'completed_software_fixture')
        self.assertFalse(result['scientific_scores_produced'])
        tier=result['tiers']['tier3_reviewed_relatedness']
        self.assertEqual(tier['common_scored_query_ids'],['q1','q2'])
        for method in tier['methods'].values():
            row=next(r for r in method['rows'] if r['query_id']=='q1')
            self.assertNotIn('r0',row['reference_ids'])

    def test_common_cli_keeps_all_four_pending_extensions_non_successful(self):
        for extension,expected in [('E1','no_inputs'),('E2','blocked_no_reviewed_inputs'),
                                   ('E3','blocked_no_reviewed_inputs'),('E4','blocked_no_measurements')]:
            output=self.base/extension
            args=['--protocol',str(ROOT/'config/extensions'/f'{extension}.json'),
                  '--inputs',str(ROOT/'config/extensions/inputs.pending.json'),'--output',str(output)]
            with self.subTest(extension=extension),redirect_stdout(io.StringIO()):
                self.assertEqual(main(args),1)
                report=json.loads((output/'report.json').read_text())
                self.assertEqual(report['status'],expected)
                self.assertFalse(report['inference_executed'])
                self.assertFalse(report['scientific_scores_produced'])

    def test_unavailable_model_does_not_silently_shrink_common_queries(self):
        def factory(name,*args,**kwargs):
            if name=='esmc':
                raise ModuleNotFoundError('PRIVATE_ERROR_DETAILS')
            return self.factory(name,*args,**kwargs)
        result=self.run_fixture(factory)
        self.assertEqual(result['tiers']['tier1_t07_invariants']['common_scored_query_ids'],[])
        self.assertEqual(result['models']['esmc']['status'],'blocked')
        self.assertNotIn('PRIVATE_ERROR_DETAILS',json.dumps(result))

    def test_same_parent_queries_share_uncertainty_and_e4_cluster_count(self):
        for qid in ('q1', 'q2'):
            old = self.proteins[qid]
            self.proteins[qid] = ProteinInput(qid, 'one-shared-parent', old.sequence, old.source_sha256)
        result = self.run_fixture()
        tier = result['reference_tiers']['tiers']['tier1_t07_invariants']
        self.assertEqual(tier['query_independence_group_count'], 2)
        self.assertEqual(tier['query_independence_cluster_count'], 1)
        comparisons = result['tiers']['tier1_t07_invariants']['paired_comparisons']
        self.assertTrue(comparisons)
        self.assertTrue(all(c['independent_clusters'] == 1 for c in comparisons))
        records = json.loads((self.output / 'cost-quality-records.json').read_text())['records']
        self.assertTrue(records)
        self.assertTrue(all(r['metrics']['independent_groups'] == 1 for r in records))

    def test_empty_reviewed_input_configuration_never_loads_sources_or_weights(self):
        self.input_path.write_text(json.dumps({'records':[]}),encoding='utf-8')
        with patch('iri_constructs.extension_experiments.reviewed_sources',side_effect=AssertionError('source read')):
            result=run_reference_generalization(self.input_path,ROOT/'config/extensions/E3.json',self.output)
        self.assertEqual(result['status'],'blocked_no_reviewed_inputs')
        self.assertFalse(result['inference_executed'])
        self.assertFalse(result['scientific_scores_produced'])

    def test_output_cannot_overwrite_a_previous_attempt(self):
        self.output.mkdir();(self.output/'report.json').write_text('{}')
        with self.assertRaises(InputError):
            run_reference_generalization(self.input_path,ROOT/'config/extensions/E3.json',self.output)
