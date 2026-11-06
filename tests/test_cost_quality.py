"""Fictional software metadata/numbers only; flags simulate gates, not source approvals.

The explicit fixture plot (when optional matplotlib exists) stays under an ignored
temporary test directory and is titled FICTIONAL SOFTWARE FIXTURE - NOT RESEARCH.
"""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.cost_quality import (plot_cost_quality,records_from_report,
                                         summarize_cost_quality,write_cost_quality_report)
from iri_constructs.runtime import InputError


def fictional_metadata_record(method='fixture-method-a',run_id='fixture-run'):
    # These fields deliberately simulate a declared admitted record to exercise
    # metadata validation; all identities, measurements and approvals are fictional.
    return dict(schema_version='cost-quality-record-1',run_id=run_id,method=method,
        evidence_kind='admitted_research',status='completed',actual_measurements=True,
        source_admission={'all_sources_approved':True,'review_manifest_sha256':'a'*64},
        panel_sha256='b'*64,protocol_sha256='c'*64,input_config_sha256='d'*64,budget=5,
        budget_scope='five total hits',comparison_task='fictional software test',reference_tier='fixture-tier',
        runtime_policy_sha256='e'*64,method_configuration_sha256='f'*64,
        cache_policy='exact-sequence within method',cache_state='new model and empty sequence cache; no warm-up',
        precision='float32',hardware={'platform':'Fictional OS','cpu':'Fictional CPU','device':'cpu','logical_cpu_count':2},
        software={'python':'fictional-test-version','packages':{'fictional-test':'1'}},
        scored_observation_ids=['fixture-q'+str(i) for i in range(8)],
        shared_execution_id=run_id+'/'+method,
        metrics={'coverage':1.0,'documented_family_recall_at_k':.5,'complete_family_recovery_rate':.25,
                 'independent_groups':8,'scored_observations':8,'requested_observations':8,
                 'budget_saturated_observations':0},
        costs={'cold_load_seconds':2.0,'inference_pooling_seconds':8.0,'total_wall_seconds':10.0,
               'model_calls':16,'unique_encoded_sequences':16,'processed_residues':160,
               'cpu_seconds':6.0,'gpu_seconds':None,'peak_process_rss_bytes':1024})


class CostQualityTests(unittest.TestCase):
    def setUp(self):
        self.a=fictional_metadata_record()
        self.b=fictional_metadata_record('fixture-method-b')
        self.b['method_configuration_sha256']='1'*64

    def test_identical_context_is_required_and_no_winner_is_selected(self):
        summary=summarize_cost_quality([self.a,self.b])
        self.assertEqual(summary['status'],'ok')
        self.assertEqual(summary['pair_comparability'][0]['status'],'matched_measured_context')
        self.assertEqual(len(summary['direct_comparison_groups']),1)
        self.assertEqual(summary['rows'][0]['scored_observations_per_inference_second'],1.0)
        self.assertIsNone(summary['automatic_winner'])

    def test_context_mismatches_are_never_direct_comparisons(self):
        changes={'panel_sha256':'2'*64,'protocol_sha256':'3'*64,'input_config_sha256':'4'*64,
            'runtime_policy_sha256':'5'*64,'budget':4,'cache_state':'warm cache',
            'cache_policy':'different policy','precision':'bfloat16','reference_tier':'other-tier',
            'hardware':{'platform':'Different','cpu':'different','device':'cpu'},
            'software':{'python':'other-version'}}
        for key,value in changes.items():
            with self.subTest(key=key):
                other={**self.b,key:value}
                result=summarize_cost_quality([self.a,other])
                pair=result['pair_comparability'][0]
                self.assertEqual(pair['status'],'rejected_direct_comparison')
                self.assertIn(key,pair['mismatched_fields'])

    def test_different_supported_query_identities_cannot_hide_behind_equal_counts(self):
        self.b['scored_observation_ids']=['different-fixture-q'+str(i) for i in range(8)]
        result=summarize_cost_quality([self.a,self.b])
        self.assertIn('scored_observation_ids',result['pair_comparability'][0]['mismatched_fields'])

    def test_same_method_name_cannot_mix_changed_model_configurations(self):
        self.b['run_id']='other-fixture-run'
        self.b['method']=self.a['method']
        pair=summarize_cost_quality([self.a,self.b])['pair_comparability'][0]
        self.assertEqual(pair['status'],'rejected_direct_comparison')
        self.assertTrue(any('different underlying' in key for key in pair['mismatched_fields']))

    def test_three_row_group_cannot_plot_ambiguous_same_method_variants(self):
        changed=copy.deepcopy(self.a)
        changed['run_id']='another-fictional-run'
        changed['method_configuration_sha256']='2'*64
        summary=summarize_cost_quality([self.a,changed,self.b])
        self.assertFalse(summary['rows'][0]['direct_comparison_eligible'])
        self.assertFalse(summary['rows'][1]['direct_comparison_eligible'])
        self.assertTrue(summary['rows'][2]['direct_comparison_eligible'])
        self.assertFalse(any(g['status']=='matched_methods_available' for g in summary['direct_comparison_groups']))
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            self.assertIsNone(plot_cost_quality(summary,Path(directory)/'forbidden.svg')['path'])

    def test_claimed_research_without_source_admission_is_not_counted_admitted(self):
        for admission in ({'all_sources_approved':False,'review_manifest_sha256':'a'*64},
                          {'all_sources_approved':True,'review_manifest_sha256':None}):
            record={**self.a,'source_admission':admission}
            result=summarize_cost_quality([record])
            self.assertEqual(result['actual_admitted_research_records'],0)
            self.assertEqual(result['status'],'descriptive_only_no_admitted_research')

    def test_missing_metrics_are_null_not_invented_zero(self):
        record={**self.a,'metrics':{},'costs':{},'scored_observation_ids':None,'evidence_kind':'unclassified'}
        result=summarize_cost_quality([record])
        row=result['rows'][0]
        self.assertIsNone(row['metrics']['coverage'])
        self.assertIsNone(row['costs']['cold_load_seconds'])
        self.assertIsNone(row['processed_residues_per_inference_second'])
        self.assertFalse(row['direct_comparison_eligible'])

    def test_false_measurement_or_source_approval_flags_exclude_comparisons(self):
        for change in ({'actual_measurements':False},{'source_admission':{'all_sources_approved':False}},
                       {'evidence_kind':'fictional_software_fixture'},{'status':'blocked'},
                       {'hardware':None},{'cache_state':None}):
            with self.subTest(change=change):
                summary=summarize_cost_quality([{**self.a,**change},self.b])
                self.assertFalse(summary['rows'][0]['direct_comparison_eligible'])
                self.assertEqual(summary['pair_comparability'][0]['status'],'rejected_direct_comparison')
        row=summarize_cost_quality([{**self.a,'actual_measurements':False}])['rows'][0]
        self.assertIsNone(row['scored_observations_per_inference_second'])

    def test_insufficient_groups_and_budget_saturation_are_explicit(self):
        for change in ({'independent_groups':2},{'budget_saturated_observations':8}):
            record=copy.deepcopy(self.a);record['metrics'].update(change)
            row=summarize_cost_quality([record])['rows'][0]
            self.assertFalse(row['direct_comparison_eligible'])
            self.assertEqual(row['metrics']['documented_family_recall_at_k'],.5)

    def test_nonfinite_impossible_counts_and_duplicate_record_identity_fail(self):
        for field,value in [('coverage',float('nan')),('documented_family_recall_at_k',1.1),
                            ('scored_observations',9),('scored_observations',True)]:
            record=copy.deepcopy(self.a);record['metrics'][field]=value
            with self.subTest(field=field,value=value),self.assertRaises(InputError):
                summarize_cost_quality([record])
        other=copy.deepcopy(self.a);other['costs']['cold_load_seconds']=3.0
        with self.assertRaises(InputError):
            summarize_cost_quality([self.a,other])

    def test_legacy_adapter_keeps_unknown_context_unclassified(self):
        report={'methods':{'fixture':{'status':'ok','overall':{'documented_family_recall_at_k':1.0}}},
                'sequence':'SHOULD NOT BE COPIED'}
        records=records_from_report(report)
        self.assertEqual(records[0]['evidence_kind'],'unclassified')
        self.assertFalse(records[0]['actual_measurements'])
        self.assertNotIn('sequence',str(records))
        self.assertIsNone(summarize_cost_quality(records)['rows'][0]['costs']['cold_load_seconds'])
        self.assertEqual(records_from_report({'records':[self.a]}),[self.a])

    def test_empty_report_and_fixture_records_never_write_research_plot(self):
        protocol=json.loads((ROOT/'config'/'extensions'/'E4.json').read_text())
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            output=Path(directory)/'run'
            result=write_cost_quality_report([],output,protocol=protocol)
            self.assertEqual(result['status'],'blocked_no_measurements')
            self.assertFalse((output/'cost-quality.svg').exists())
            self.assertTrue((output/'protocol-freeze.json').exists())
            self.assertFalse(result['inference_executed'])
            record={**self.a,'evidence_kind':'fictional_software_fixture'}
            fixture=write_cost_quality_report([record],Path(directory)/'fixture',protocol=protocol)
            self.assertEqual(fixture['status'],'descriptive_only_no_admitted_research')
            self.assertIsNone(fixture['plot']['path'])
            self.assertIn('NULL',(Path(directory)/'fixture'/'cost-quality.csv').read_text())

    def test_optional_plot_dependency_is_honestly_unavailable(self):
        summary=summarize_cost_quality([self.a,self.b])
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            path=Path(directory)/'plot.svg'
            with patch.dict(sys.modules,{'matplotlib':None}):
                result=plot_cost_quality(summary,path)
            self.assertEqual(result['status'],'unavailable_plot_dependency')
            self.assertFalse(path.exists())

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'Optional explicitly fictional plot fixture')
    def test_fixture_plot_is_labelled_and_never_claimed_as_research(self):
        records=[{**r,'evidence_kind':'fictional_software_fixture'} for r in (self.a,self.b)]
        summary=summarize_cost_quality(records)
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            path=Path(directory)/'fictional.svg'
            self.assertIsNone(plot_cost_quality(summary,path)['path'])
            result=plot_cost_quality(summary,path,_allow_fictional_test_fixture=True)
            self.assertEqual(result['status'],'written_fictional_software_fixture')
            self.assertIn('FICTIONAL SOFTWARE FIXTURE',path.read_text())
