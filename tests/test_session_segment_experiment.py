"""Small source-integrity and label-separation tests; no model or real alignment."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
spec=importlib.util.spec_from_file_location('segment_experiment',ROOT/'scripts/session_segment_experiment.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
from iri_constructs.session_store import UnitStore,digest
from iri_constructs.session_benchmark import rank_scores,evaluate_query
from iri_constructs.session_segmentation import fixed_windows,interval_metrics,fuse_rankings


class SegmentExperimentTests(unittest.TestCase):
    def test_source_must_be_complete_and_checksum_valid(self):
        with tempfile.TemporaryDirectory() as d:
            session=Path(d);job=session/'job'
            identity={'mode':'alignment','query_ids':['q1']}
            store=UnitStore(job,identity)
            store.publish('query-q1',{'marker':'kept'})
            with self.assertRaises(ValueError):runner.open_store(job,session)
            store.publish('summary',{'marker':'summary'})
            store.complete(['query-q1','summary'])
            opened=runner.open_store(job,session)
            self.assertEqual(runner.completed_read(opened,'query-q1'),{'marker':'kept'})
            (opened.completed_paths['query-q1']/'record.json').write_text('{}')
            with self.assertRaises(ValueError):runner.open_store(job,session)

    def test_source_hashes_reject_same_accession_different_data(self):
        suite={'sequences':{'A':'numeric fixture one','B':'numeric fixture two'}}
        identity={'challenge_sha256':'challenge','catalogue_accessions':['A','B'],
                  'catalogue_sequence_sha256':{a:digest(v.encode()) for a,v in suite['sequences'].items()}}
        runner.validate_source_suite(identity,suite,'challenge')
        changed=copy.deepcopy(identity);changed['catalogue_sequence_sha256']['A']='wrong'
        with self.assertRaises(ValueError):runner.validate_source_suite(changed,suite,'challenge')
        with self.assertRaises(ValueError):runner.validate_source_suite(identity,suite,'another challenge')

    def fixture(self):
        query={'query_id':'q1','group_id':'A__C','suite':'pairs','component_count':2,'components':['A','C'],
               'length_aa':6,'sequence':'not used','spans':[{'accession':'A','start0':0,'end0':3,'fraction':.5},
                                                        {'accession':'C','start0':3,'end0':6,'fraction':.5}]}
        prediction_cfg={'segmenter':'fixed_windows','width':4,'stride':2}
        prediction=fixed_windows(6,4,2)
        proposal={'query':runner.public_query(query),'prediction_config':prediction_cfg,'prediction':prediction,
                  'boundary_evaluation':interval_metrics(prediction['segments'],query['spans']),
                  'logical_candidate_scan_residues':0,'logical_encoded_segment_residues':8,'logical_total_queried_residues':8}
        config={**prediction_cfg,'configuration_id':'window','fusion':'rrf','rrf_constant':60}
        names=['A','B','C'];eq={a:[a] for a in names}
        rankings=[rank_scores([.9,.8,.1],names,'cosine_similarity'),rank_scores([.1,.8,.9],names,'cosine_similarity')]
        fusion=fuse_rankings(rankings,'rrf')
        record={'query':runner.public_query(query),'configuration':config,'prediction':prediction,
                'boundary_evaluation':proposal['boundary_evaluation'],'logical_total_queried_residues':8,
                'segment_rankings':rankings,'ranking':fusion['ranking'],
                'evaluation':evaluate_query(fusion['ranking'],query,eq),'larger_budget':runner.union_results(rankings,query)}
        return query,prediction_cfg,proposal,config,names,eq,record

    def test_proposal_geometry_and_logical_work(self):
        q,c,p,*_=self.fixture();runner.validate_proposal(p,q,c)
        altered=copy.deepcopy(p);altered['logical_total_queried_residues']=6
        with self.assertRaises(ValueError):runner.validate_proposal(altered,q,c)
        altered=copy.deepcopy(p);altered['prediction']['segments'][0]['end0']=7
        with self.assertRaises(ValueError):runner.validate_proposal(altered,q,c)

    def test_metrics_keys_and_larger_budget_are_checked(self):
        q,_,p,c,names,eq,record=self.fixture()
        runner.validate_retrieval_record(record,q,c,p,names,eq)
        self.assertEqual(record['evaluation']['query_metrics']['exact']['component_recall_at_5'],1.)
        self.assertEqual(record['larger_budget']['1']['returned_hit_occurrences'],2)
        self.assertEqual(record['larger_budget']['1']['evaluation']['component_hits'],2)
        altered=copy.deepcopy(record);altered['larger_budget']['1']['evaluation']['component_hits']=1
        with self.assertRaises(ValueError):runner.validate_retrieval_record(altered,q,c,p,names,eq)
        altered=copy.deepcopy(record);altered['segment_rankings'][0].pop()
        with self.assertRaises(ValueError):runner.validate_retrieval_record(altered,q,c,p,names,eq)

    def test_expected_units_include_every_query_configuration(self):
        identity={'mode':'automatic_segment_retrieval_v1','query_ids':['q1','q2'],
                  'configs':[{'configuration_id':'a'},{'configuration_id':'b'}]}
        self.assertEqual(runner.expected_keys(identity),{'q1--a','q1--b','q2--a','q2--b','summary'})
        identity['configs']=[]
        with self.assertRaises(ValueError):runner.expected_keys(identity)

    def test_dev_representation_is_not_silently_reselected(self):
        expected={'model_id':'fixed model','revision':'fixed revision','device':'cpu','dtype':'float32','layer':'last'}
        method={**expected,'layer':6,'num_layers':6}
        runner.validate_dev_representation(method,{'selection_representation':expected})
        for key,value in [('model_id','another model'),('device','mps'),('layer',3)]:
            with self.assertRaises(ValueError):runner.validate_dev_representation({**method,key:value},{'selection_representation':expected})

    def test_setup_uses_dev_only_then_frozen_winners_only(self):
        with tempfile.TemporaryDirectory() as d:
            session=Path(d);inputs=session/'inputs';inputs.mkdir()
            challenge=inputs/'challenge_manifest.json';challenge.write_text('{}')
            configs=[{'configuration_id':'w1','segmenter':'fixed_windows'},
                     {'configuration_id':'w2','segmenter':'fixed_windows'},
                     {'configuration_id':'a1','segmenter':'alignment_guided'}]
            representation={'model_id':'fixed','revision':'fixed','device':'cpu','dtype':'float32','layer':'last'}
            design={'source_challenge_sha256':runner.file_hash(challenge),'implementation_sha256':{},
                    'selection_query_ids':['dev'],'test_query_ids':['test'],'clean_sanity_query_ids':['clean'],
                    'partition':{'queries':{'dev':'dev','test':'test','clean':'dev'}},
                    'candidate_configurations':configs,'selection_representation':representation}
            design_path=inputs/'design.json';design_path.write_text(json.dumps(design))
            suite={'queries':[{'query_id':'dev','component_count':2},{'query_id':'test','component_count':2},
                              {'query_id':'clean','component_count':1}]}
            args=SimpleNamespace(session=session,design=design_path,partition='dev',selection=None)
            with patch.object(runner,'load_suite',return_value=suite):
                result=runner.setup(args)
            self.assertEqual([q['query_id'] for q in result[3]],['dev'])
            self.assertEqual(len(result[4]),3)
            summary=session/'summary';summary.mkdir();(summary/'complete.json').write_text('{}')
            selection={'selected_by_segmenter':{'fixed_windows':'w2','alignment_guided':'a1'}}
            frozen={'design_sha256':runner.file_hash(design_path),'development_job':str(session/'development'),
                    'development_identity':'identity','development_summary_descriptor_sha256':runner.file_hash(summary/'complete.json'),
                    'selection':selection}
            selection_path=inputs/'selection.json';selection_path.write_text(json.dumps(frozen))
            dev=SimpleNamespace(identity={'mode':'automatic_segment_retrieval_v1','partition':'dev',
                'design_sha256':runner.file_hash(design_path),'query_ids':['dev'],'configs':configs,
                'model':{**representation,'layer':6,'num_layers':6}},identity_hash='identity',completed_paths={'summary':summary})
            args.partition='test';args.selection=selection_path
            with patch.object(runner,'load_suite',return_value=suite),patch.object(runner,'open_store',return_value=dev),patch.object(runner,'completed_read',return_value={'selection':selection}):
                result=runner.setup(args)
                self.assertEqual([q['query_id'] for q in result[3]],['test','clean'])
                self.assertEqual([c['configuration_id'] for c in result[4]],['w2','a1'])
                selection['selected_by_segmenter']['fixed_windows']='unknown'
                selection_path.write_text(json.dumps({**frozen,'selection':selection}))
                with self.assertRaises(ValueError):runner.setup(args)


if __name__=='__main__':unittest.main()
