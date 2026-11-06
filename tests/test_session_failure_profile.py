"""Small fictional fixtures; no model/ranking work or real outcome files."""
import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import session_failure_profile as profile


def design():
    return {'version':'descriptive-failure-profile-v1','frozen_before_outcome_row_reads':True,
        'source_job':'fictional-comparison','producer_code_sha256_at_design':'fictional-producer-sha',
        'allowed_exports':['component_metrics.csv','cases.json'],
        'models':[{'job_id':f'model{i}','short_label':f'Fictional model {i}'} for i in range(5)],
        'population':{'component_observations_per_model':12,'joined_queries_per_model':4},
        'fraction_bins':['lt_0.25','0.25_to_lt_0.40','0.40_to_lt_0.60','0.60_to_lt_0.75','ge_0.75'],
        'positions':['N','internal','C'],'case_rules':['fictional fixed first2 worst and first1 disagreement']}


def metric(job,dimension,key,n,hits):
    return {'job_id':job,'label':job,'stratum_dimension':dimension,'stratum':key,
        'original_components':str(n),'scorable_components':str(n),'unscorable_components':'0',
        'parent_queries':str(min(n,4)),'distinct_groups':str(min(n,2)),
        'scorable_targets':str(min(n,3)),'hits_at_5':str(hits),'recall_at_5':str(hits/n) if n else '',**profile.FILTER}


def metrics():
    rows=[];d=design()
    for model in d['models']:
        job=model['job_id'];rows.append(metric(job,'overall','all',12,6))
        for key,n,hits in zip(d['fraction_bins'],[4,3,3,2,0],[2,1,2,1,0]):
            if n:rows.append(metric(job,'fraction',key,n,hits))
        rows.extend(metric(job,'position',key,4,2) for key in d['positions'])
    return rows


def observation(job,index,rank=15):
    return {'job_id':job,'mode':'model','arm':'whole','retrieval_budget':profile.FILTER['retrieval_budget'],
            'query_id':f'q{index}','group_id':f'g{index}','accession':'fictional-A',
            'component_count':3,'rank':rank}


def cases():
    values=[]
    for i in range(5):
        for j in range(3):values.append({'selection_rule':'worst_component_rank','case_type':'within_method',
                                        'observation':observation(f'model{i}',j,20-j)})
    for a in range(5):
        for b in range(a+1,5):
            for j in range(2):
                values.append({'selection_rule':'largest_absolute_rank_disagreement','case_type':'between_methods',
                               'left':observation(f'model{a}',j,2),'right':observation(f'model{b}',j,20-j),
                               'group_id':f'g{j}','right_minus_left_rank':18-j})
    return {'rules':dict(profile.SOURCE_RULES),'cases':values}


class ProfileTests(unittest.TestCase):
    def test_exact_common_panel_partitions_and_undefined_empty_bins(self):
        result=profile.extract_profile(metrics(),design())
        self.assertEqual(len(result),45)
        empty=[r for r in result if r['stratum']=='ge_0.75']
        self.assertEqual(len(empty),5)
        self.assertTrue(all(r['scorable_components']==0 and r['recall_at_5'] is None
                            and not r['source_row_present'] for r in empty))
        overall=[r for r in result if r['stratum_dimension']=='overall']
        self.assertTrue(all(r['hits_at_5']==6 and r['parent_queries']==4 for r in overall))

    def test_clean_full_coverage_oracle_and_duplicate_labels_cannot_contaminate(self):
        original=metrics();extras=[]
        for field,value in [('scope','clean'),('evaluation_population','full_job_coverage'),
                            ('arm','isolated_component'),('label_kind','duplicate_equivalent')]:
            row=copy.deepcopy(original[0]);row[field]=value;row['hits_at_5']='12';row['recall_at_5']='1.0';extras.append(row)
        self.assertEqual(profile.extract_profile(original+extras,design()),profile.extract_profile(original,design()))

    def test_duplicate_summary_rejected(self):
        rows=metrics();rows.append(copy.deepcopy(rows[0]))
        with self.assertRaisesRegex(ValueError,'Duplicate source'):profile.extract_profile(rows,design())

    def test_changed_panel_query_count_rejected(self):
        rows=metrics();rows[0]['parent_queries']='3'
        with self.assertRaisesRegex(ValueError,'joined panel'):profile.extract_profile(rows,design())

    def test_fraction_partition_numerator_must_match_original_summary(self):
        rows=metrics();rows[1]['hits_at_5']='3';rows[1]['recall_at_5']='0.75'
        with self.assertRaisesRegex(ValueError,'partition differs'):profile.extract_profile(rows,design())

    def test_original_recall_not_replaced_by_silent_recalculation(self):
        rows=metrics();rows[1]['recall_at_5']='0.6'
        with self.assertRaisesRegex(ValueError,'Saved recall disagrees'):profile.extract_profile(rows,design())

    def test_cross_model_stratum_support_must_match(self):
        rows=metrics()
        for row in rows:
            if row['job_id']=='model1' and row['stratum_dimension']=='position' and row['stratum']=='N':
                row['parent_queries']='3'
        with self.assertRaisesRegex(ValueError,'identical stratum denominators'):profile.extract_profile(rows,design())

    def test_case_selection_reuses_saved_order_with_fixed_maximum(self):
        source=cases();result=profile.select_cases(source,design())
        self.assertEqual(len(result['selected_cases']),20)
        selected=result['selected_cases']
        self.assertEqual([v['source_case_index'] for v in selected[:4]],[0,1,3,4])
        self.assertTrue(all(v['case']==source['cases'][v['source_case_index']] for v in selected))
        self.assertEqual(sum(x['selected'] for x in result['category_counts']),20)

    def test_case_rules_and_unmatched_disagreement_rejected(self):
        source=cases();source['rules']['worst']='a changed rule'
        with self.assertRaisesRegex(ValueError,'case rules changed'):profile.select_cases(source,design())
        source=cases();source['cases'][15]['right']['query_id']='different'
        with self.assertRaisesRegex(ValueError,'unmatched component'):profile.select_cases(source,design())

    def test_duplicate_group_in_source_case_category_rejected(self):
        source=cases();source['cases'][2]['observation']['group_id']='g0'
        with self.assertRaisesRegex(ValueError,'repeats an unordered group'):profile.select_cases(source,design())

    def test_source_alignment_examples_excluded_without_new_selection(self):
        source=cases();source['cases'].insert(0,{'selection_rule':'worst_component_rank','case_type':'within_method',
            'observation':{**observation('alignment-source',0,105),'mode':'alignment','arm':'local_alignment'}})
        result=profile.select_cases(source,design())
        self.assertEqual(len(result['selected_cases']),20)
        self.assertEqual(result['selected_cases'][0]['source_case_index'],1)

    def test_csv_roundtrip_preserves_empty_recall_and_support(self):
        result=profile.extract_profile(metrics(),design())
        parsed=list(csv.DictReader(io.StringIO(profile.csv_bytes(result).decode())))
        empty=[r for r in parsed if r['stratum']=='ge_0.75']
        self.assertTrue(all(r['recall_at_5']=='' and r['scorable_components']=='0' for r in empty))


class BindingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='iri-fictional-profile-');self.addCleanup(self.temp.cleanup)
        self.source=Path(self.temp.name).resolve()/'fictional-comparison';self.source.mkdir()
        self.design=design();(self.source/'exports').mkdir()
        # Binding checks intentionally do not parse these outcome bytes.
        (self.source/'exports/component_metrics.csv').write_bytes(b'fictional not-yet-parsed csv\n')
        (self.source/'exports/cases.json').write_bytes(b'fictional not-yet-parsed json\n')
        identity={'mode':'common_model_comparison','catalogue_accessions':[f'F{i:03}' for i in range(105)],
                  'sources':[{'job_id':x['job_id'],'mode':'model'} for x in self.design['models']],
                  'code_hashes':{'scripts/session_compare_models.py':self.design['producer_code_sha256_at_design']}}
        self.write(self.source/'identity.json',identity)
        ih=hashlib.sha256(profile.canonical(identity)).hexdigest();units={}
        exports={'path':str(self.source/'exports'),'artifacts':{p.name:{'sha256':profile.sha(p),'bytes':p.stat().st_size}
                 for p in (self.source/'exports').iterdir()}}
        for key,value in [('comparison-data',{'fictional':'body must not be read'}),('exports',exports)]:
            unit=self.source/'units'/key/'fictional-v1';unit.mkdir(parents=True);self.write(unit/'record.json',value)
            self.write(unit/'complete.json',{'key':key,'identity_hash':ih,'files':{'record.json':{
                'sha256':profile.sha(unit/'record.json'),'bytes':(unit/'record.json').stat().st_size}}})
            units[key]={'path':str(unit.relative_to(self.source)),'descriptor_sha256':profile.sha(unit/'complete.json')}
        completion=self.source/'completion-fixture.json'
        self.write(completion,{'identity_hash':ih,'expected_count':2,'units':units})
        self.design.update(identity_file_sha256=profile.sha(self.source/'identity.json'),
            completion_file=completion.name,completion_sha256=profile.sha(completion),
            exports_descriptor_sha256=units['exports']['descriptor_sha256'])

    @staticmethod
    def write(path,value):path.write_bytes(profile.canonical(value)+b'\n')

    def test_export_chain_binds_before_any_outcome_parse(self):
        paths,evidence=profile.bind_source(self.source,self.design)
        self.assertEqual(set(paths),{'component_metrics.csv','cases.json'})
        self.assertFalse(evidence['comparison_data_body_read']);self.assertFalse(evidence['raw_model_records_read'])
        self.assertEqual(evidence['exports']['cases.json']['sha256'],profile.sha(paths['cases.json']))

    def test_mutated_outcome_export_cannot_be_used(self):
        (self.source/'exports/cases.json').write_bytes(b'changed fictional bytes')
        with self.assertRaisesRegex(ValueError,'export bytes/hash changed'):
            profile.bind_source(self.source,self.design)

    def test_incomplete_or_modified_completion_rejected(self):
        path=self.source/self.design['completion_file'];value=json.loads(path.read_text())
        value['units'].pop('comparison-data');self.write(path,value)
        self.design['completion_sha256']=profile.sha(path)
        with self.assertRaisesRegex(ValueError,'completion is incomplete'):
            profile.bind_source(self.source,self.design)

    def test_exports_symlink_rejected_even_for_identical_bytes(self):
        path=self.source/'exports/cases.json';copy=self.source/'fictional-copy';copy.write_bytes(path.read_bytes())
        path.unlink();path.symlink_to(copy)
        with self.assertRaisesRegex(ValueError,'Symlink source'):
            profile.bind_source(self.source,self.design)


if __name__=='__main__':unittest.main()
