"""Fictional numerical fixtures for E2; not biological or pretrained-model results."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.multi_cds import (freeze_multi_cds_panel, score_multi_cds, run_multi_cds,
                                    load_multi_cds_inputs, compare_multi_cds_conditions)
from iri_constructs.runtime import InputError, digest_file
from iri_constructs.validate_manifest import REQUIRED_COLUMNS


class MultiCdsFixture:
    def setUp(self):
        self.config=json.loads((ROOT/'config'/'extensions'/'E2.json').read_text())
        self.proteins={
            'u1':ProteinInput('u1','p','A'*4,'fixture','feature-1'),
            'u2':ProteinInput('u2','p','C'*3,'fixture','feature-2'),
        }
        for iid,aa in [('ra','A'),('rb','C'),('n1','D'),('n2','E'),('n3','F'),('n4','G')]:
            self.proteins[iid]=ProteinInput(iid,'parent-'+iid,aa*3+'W','fixture')
        self.labels={iid:dict(input_id=iid,family='one' if iid in ('u1','ra') else 'two' if iid in ('u2','rb') else iid,
            annotation_status='known',review_status='approved',evidence='Fictional software fixture only',
            independence_group='group-'+iid) for iid in self.proteins}
        self.parents={'p':dict(unit_ids=['u1','u2'],source_cds_count=2,archetype='multi_cds',
            source_sha256='fixture',documented_family_by_unit={'u1':'one','u2':'two'})}
        self.refs=['ra','rb','n1','n2','n3','n4']
        self.scores={'u1':dict(ra=.99,rb=.1,n1=.8,n2=.7,n3=.6,n4=.5),
                     'u2':dict(ra=.2,rb=.98,n1=.1,n2=.1,n3=.1,n4=.1)}

    def panel(self):
        return freeze_multi_cds_panel(self.config,self.proteins,self.labels,self.parents,self.refs)


class MultiCdsScoringTests(MultiCdsFixture,unittest.TestCase):
    def test_paired_uncertainty_retains_cluster_threshold_and_excludes_saturation(self):
        panel=self.panel()
        result=score_multi_cds(self.scores,panel,self.labels)
        comparison=compare_multi_cds_conditions(result,panel,self.config)
        estimate=comparison['endpoints']['matched']['documented_family_recall_at_k']
        self.assertEqual(estimate['paired_query_ids'],['p'])
        self.assertEqual(estimate['paired_cluster_mean_difference'],.5)
        self.assertEqual(estimate['status'],'insufficient_independent_clusters')
        self.assertIsNone(estimate['ci_lower'])
        # Only five eligible references: the whole catalogue is returned.
        self.refs.remove('n4')
        panel=self.panel()
        scores={u:{r:v for r,v in row.items() if r in self.refs} for u,row in self.scores.items()}
        comparison=compare_multi_cds_conditions(score_multi_cds(scores,panel,self.labels),panel,self.config)
        self.assertEqual(comparison['excluded_budget_saturated_parent_ids'],['p'])
        self.assertEqual(comparison['endpoints']['matched']['documented_family_recall_at_k']['paired_query_ids'],[])

    def test_parent_max_merge_returns_five_total_unique_references(self):
        panel=self.panel()
        result=score_multi_cds(self.scores,panel,self.labels)
        before=result['conditions']['longest_cds']['full']['rows'][0]
        after=result['conditions']['all_cds']['full']['rows'][0]
        self.assertEqual(before['documented_family_recall_at_k'],.5)
        self.assertEqual(after['documented_family_recall_at_k'],1.)
        self.assertEqual(after['reference_ids'],['ra','rb','n1','n2','n3'])
        self.assertEqual(len(set(after['reference_ids'])),5)
        self.assertEqual(result['matched_per_unit_control'][0]['unit_id'],'u1')

    def test_duplicate_documented_families_have_one_denominator(self):
        self.labels['u2']['family']='one'
        self.parents['p']['documented_family_by_unit']['u2']='one'
        result=score_multi_cds(self.scores,self.panel(),self.labels)
        self.assertEqual(result['conditions']['all_cds']['full']['rows'][0]['distinct_expected_families'],1)

    def test_equal_scores_use_reference_id_tie_break_after_max_merge(self):
        scores={u:{r:0. for r in reversed(self.refs)} for u in ('u1','u2')}
        result=score_multi_cds(scores,self.panel(),self.labels)
        self.assertEqual(result['conditions']['all_cds']['full']['rows'][0]['reference_ids'],sorted(self.refs)[:5])

    def test_no_independent_family_examples_are_abstentions_not_negatives(self):
        self.refs=['n1','n2','n3','n4']
        panel=self.panel()
        result=score_multi_cds({u:{r:self.scores[u][r] for r in self.refs} for u in ('u1','u2')},panel,self.labels)
        self.assertEqual(panel['status'],'blocked_insufficient_parent_panel')
        self.assertEqual(result['conditions']['all_cds']['matched']['scored_parents'],0)
        self.assertEqual(result['matched_per_unit_control'],[])
        self.assertTrue(result['matched_per_unit_abstentions'])

    def test_longest_tie_uses_stable_unit_id(self):
        self.proteins['u2']=ProteinInput('u2','p','C'*4,'fixture','feature-2')
        self.parents['p']['unit_ids'].reverse()
        self.assertEqual(self.panel()['parents']['p']['baseline_unit_id'],'u1')

    def test_longer_unlabelled_cds_is_withheld_not_favorable_silent_baseline(self):
        self.proteins['u2']=ProteinInput('u2','p','C'*8,'fixture','feature-2')
        del self.labels['u2']
        del self.parents['p']['documented_family_by_unit']['u2']
        parent=self.panel()['parents']['p']
        self.assertEqual(parent['baseline_unit_id'],'u1')
        self.assertEqual(parent['unknown_family_unit_ids'],['u2'])
        self.assertEqual(parent['parsed_cds_coverage'],1.)
        self.assertEqual(parent['model_label_eligible_cds_coverage'],.5)
        result=score_multi_cds({'u1':self.scores['u1']},self.panel(),self.labels)
        self.assertEqual(result['conditions']['all_cds']['full']['scored_parents'],0)
        self.assertEqual(result['conditions']['all_cds']['matched']['scored_parents'],1)

    def test_missing_cds_changes_full_and_matched_supported_denominators(self):
        del self.proteins['u2']
        del self.labels['u2']
        panel=self.panel()
        self.assertEqual(panel['parents']['p']['unsupported_documented_families'],['two'])
        result=score_multi_cds({'u1':self.scores['u1']},panel,self.labels)
        condition=result['conditions']['all_cds']
        self.assertEqual(condition['full']['documented_family_recall_at_k'],.5)
        self.assertEqual(condition['matched']['documented_family_recall_at_k'],1.)

    def test_empty_parent_and_empty_catalogue_are_blocked_not_perfect(self):
        self.parents={'empty':dict(unit_ids=[],source_cds_count=0,documented_family_by_unit={})}
        panel=self.panel()
        self.assertEqual(panel['status'],'blocked_insufficient_parent_panel')
        result=score_multi_cds({},panel,self.labels)
        self.assertIsNone(result['conditions']['all_cds']['full']['documented_family_recall_at_k'])
        self.assertEqual(result['common_scored_parent_ids'],[])

    def test_every_component_parent_copy_and_group_excluded_before_scores(self):
        for iid,protein,group in [
            ('copy',ProteinInput('copy','different','C'*3,'fixture'),'different'),
            ('parent',ProteinInput('parent','p','H'*3,'fixture'),'different'),
            ('linked',ProteinInput('linked','different','I'*3,'fixture'),'group-u2')]:
            self.proteins[iid]=protein
            self.labels[iid]=dict(self.labels['ra'],input_id=iid,independence_group=group)
            self.refs.append(iid)
        parent=self.panel()['parents']['p']
        self.assertEqual(parent['eligible_reference_ids'],sorted(self.refs[:6]))
        self.assertEqual({r['reference_id'] for r in parent['reference_exclusions']},{'copy','parent','linked'})

    def test_partial_or_nonfinite_unit_panels_cannot_favor_max_merge(self):
        panel=self.panel()
        for change in ('missing_reference','missing_unit','nonfinite'):
            scores=copy.deepcopy(self.scores)
            if change=='missing_reference': del scores['u2']['n4']
            elif change=='missing_unit': del scores['u2']
            else: scores['u2']['n4']=float('nan')
            result=score_multi_cds(scores,panel,self.labels)
            self.assertEqual(result['conditions']['longest_cds']['matched']['scored_parents'],1)
            self.assertEqual(result['conditions']['all_cds']['matched']['scored_parents'],0)

    def test_injected_reference_or_changed_freeze_and_label_rejected(self):
        panel=self.panel()
        scores=copy.deepcopy(self.scores)
        scores['u1']['self']=100.
        with self.assertRaises(InputError): score_multi_cds(scores,panel,self.labels)
        bad=copy.deepcopy(panel)
        bad['top_k']=6
        with self.assertRaises(InputError): score_multi_cds(self.scores,bad,self.labels)
        labels=copy.deepcopy(self.labels)
        labels['u2']['family']='after-score-change'
        with self.assertRaises(InputError): score_multi_cds(self.scores,panel,labels)

    def test_duplicate_units_bad_budget_and_posthoc_config_blocked(self):
        self.parents['p']['unit_ids']=['u1','u1']
        with self.assertRaises(InputError): self.panel()
        self.parents['p']['unit_ids']=['u1','u2']
        for change in ({'top_k':6},{'design_status':'edited-after-results'}):
            with self.assertRaises(InputError):
                freeze_multi_cds_panel({**self.config,**change},self.proteins,self.labels,self.parents,self.refs)


class MultiCdsExecutionTests(MultiCdsFixture,unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.temp=tempfile.TemporaryDirectory(dir=ROOT/'results')
        self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name)
        self.output=self.base/'run'
        self.inputs_path=self.base/'inputs.json'
        self.inputs_path.write_text(json.dumps({'records':[{'fixture':True}]}),encoding='utf-8')
        self.protocol=self.base/'protocol.json'
        self.protocol.write_text(json.dumps(self.config),encoding='utf-8')
        self.inputs={'reference_ids':self.refs,'evidence_kind':'fictional_software_fixture'}
        self.provenance={'review_manifest_sha256':'fixture-manifest','input_config_sha256':'fixture-input'}

    def factory(self,name,device,*args,**kwargs):
        output=self.output
        class Numerical:
            metadata={'model':name,'device':device,'fictional_software_fixture':True}
            def residues(self,sequence):
                assert (output/'protocol-freeze.json').exists()
                assert (output/'panel-freeze.json').exists()
                assert (output/('model-freeze-'+name+'.json')).exists()
                v=[1.,0.] if sequence[0]=='A' else [0.,1.] if sequence[0]=='C' else [.5,.5]
                return [v]*len(sequence)
            def close(self): pass
        model=Numerical()
        model.device=device
        return model

    def run_fixture(self,factory=None):
        with patch('iri_constructs.multi_cds.load_multi_cds_inputs',return_value=(
                self.inputs,self.proteins,self.labels,self.parents,self.provenance)), \
             patch('iri_constructs.multi_cds.prepare_device',return_value='cpu'):
            return run_multi_cds(self.inputs_path,self.protocol,self.output,adapter_factory=factory or self.factory)

    def test_freeze_precedes_inference_and_stage_counts_are_actual_increments(self):
        report=self.run_fixture()
        self.assertEqual(report['status'],'ok')
        self.assertFalse(report['scientific_scores_produced'])
        records=report['cost_quality_records']
        from iri_constructs.cost_quality import summarize_cost_quality
        summary=summarize_cost_quality(records)
        self.assertEqual(summary['actual_admitted_research_records'],0)
        self.assertEqual(summary['status'],'descriptive_only_no_admitted_research')
        self.assertEqual(records[0]['shared_execution_id'],records[1]['shared_execution_id'])
        self.assertEqual([r['logical_query_units'] for r in records],[1,2])
        self.assertEqual([r['costs']['model_calls'] for r in records],[7,8])
        self.assertIsNone(records[0]['costs']['total_wall_seconds'])
        self.assertEqual(records[0]['panel_sha256'],records[1]['panel_sha256'])
        self.assertEqual(json.loads((self.output/'cost-quality-records.json').read_text())['records'],records)

    def test_empty_pending_inputs_never_load_data_or_models(self):
        self.inputs_path.write_text(json.dumps({'records':[]}),encoding='utf-8')
        with patch('iri_constructs.multi_cds.load_multi_cds_inputs',side_effect=AssertionError('No source admission')):
            report=run_multi_cds(self.inputs_path,self.protocol,self.output)
        self.assertEqual(report['status'],'blocked_no_reviewed_inputs')
        self.assertFalse(report['inference_executed'])
        self.assertEqual(report['cost_quality_records'],[])

    def test_dependency_or_all_input_failure_never_success_or_fake_scores(self):
        def missing(*args,**kwargs): raise ModuleNotFoundError('PRIVATE_DETAIL')
        report=self.run_fixture(factory=missing)
        self.assertEqual(report['status'],'blocked_no_scores')
        self.assertFalse(report['inference_executed'])
        self.assertFalse(report['scores_produced'])
        self.assertNotIn('PRIVATE_DETAIL',json.dumps(report))

    def test_additional_cds_failure_keeps_baseline_and_marks_partial(self):
        def fail_more(*args,**kwargs):
            adapter=self.factory(*args,**kwargs)
            original=adapter.residues
            def residues(sequence):
                if sequence=='C'*3: raise InputError('Fictional extra CDS failure')
                return original(sequence)
            adapter.residues=residues
            return adapter
        report=self.run_fixture(factory=fail_more)
        self.assertEqual(report['status'],'partial')
        result=report['methods']['esm2']['results']['conditions']
        self.assertEqual(result['longest_cds']['matched']['scored_parents'],1)
        self.assertEqual(result['all_cds']['matched']['scored_parents'],0)


@unittest.skipUnless(importlib.util.find_spec('Bio'),'Optional Biopython source-admission integration')
class MultiCdsSourceTests(unittest.TestCase):
    def test_loader_retains_withheld_source_cds_and_rejects_pending_review(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from iri_constructs.annotations import normalize_record
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            base=Path(directory)
            record=SeqRecord(Seq('ATG'+'GCT'*2+'TAA'+'ATG'+'TGT'*2+'TAA'),id='FICTIONAL.1')
            record.annotations={'molecule_type':'DNA'}
            record.features=[SeqFeature(SimpleLocation(0,12,strand=1),type='CDS',qualifiers={'translation':['M'+'A'*2]}),
                             SeqFeature(SimpleLocation(12,24,strand=1),type='CDS',qualifiers={'note':['Fictional missing translation fixture']})]
            source=base/'fixture.gb'
            SeqIO.write(record,source,'genbank')
            ref=base/'ref.fasta'
            ref.write_text('>REF.1\n'+'W'*4+'\n',encoding='utf-8')
            row=dict(record_id='p',source='fictional-software-fixture',source_url='https://example.invalid/fixture',
                accession='FICTIONAL',accession_version='1',archetype='multi_cds',sequence_length_nt='24',feature_count='2',
                annotation_status='approved',review_status='approved',licensing_status='approved',safety_review='approved',
                retrieved_on='2026-11-05',checksum_sha256=digest_file(source),notes='Fictional software gate values, never a real source approval')
            doc,_,_=normalize_record(record,row,row['checksum_sha256'])
            unit=next(u for u in doc['units'] if u['model_input_status']=='pending_model_selection')
            iid='source/'+unit['unit_id']
            label=lambda i,f:dict(input_id=i,family=f,annotation_status='known',review_status='approved',
                independence_group=i,evidence='Fictional software family annotation only')
            config=dict(schema_version='t05-inputs-1',review_manifest='review.csv',parent_ids=['p'],reference_ids=['r'],
                records=[dict(input_id='source',record_id='p',path=source.name,format='genbank',molecule_type='dna',sequence_id=record.id),
                         dict(input_id='r',record_id='ref-parent',path=ref.name,format='fasta',molecule_type='protein',sequence_id='REF.1')],
                families=[label(iid,'fictional-one'),label('r','fictional-other')])
            path=base/'inputs.json'
            path.write_text(json.dumps(config),encoding='utf-8')
            refrow={**row,'record_id':'ref-parent','checksum_sha256':digest_file(ref),'source':'fictional-reference',
                    'accession':'REF','feature_count':'0'}
            def manifest(parent):
                with (base/'review.csv').open('w',encoding='utf-8',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=REQUIRED_COLUMNS)
                    writer.writeheader(); writer.writerows([parent,refrow])
            manifest(row)
            _,proteins,_,parents,_=load_multi_cds_inputs(path)
            self.assertEqual(parents['p']['source_cds_count'],2)
            self.assertEqual(len(parents['p']['unit_ids']),2)
            self.assertIn(iid,proteins)
            self.assertEqual(len(parents['p']['unit_exclusions']),1)
            manifest({**row,'annotation_status':'pending'})
            with patch('iri_constructs.multi_cds.read_source',side_effect=AssertionError('No pending source parsing')):
                with self.assertRaises(InputError): load_multi_cds_inputs(path)
