"""Masked prediction and fusion checks using numeric text/metadata only."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from session_family_segmentation import (UnionAccumulator,expected_transfer_keys,family_union_results,
    mask_for_query,masked_alignment_ranker,masked_fusion,transfer_setup,validate_parent_development_method)
from iri_constructs.session_store import canonical,digest
from iri_constructs.session_benchmark import rank_scores
from iri_constructs.session_segmentation import SegmentationConfig,alignment_guided_segments,fuse_rankings


class FamilySegmentationTest(unittest.TestCase):
    def test_removed_reference_cannot_inform_predicted_boundary_even_from_full_cache(self):
        rows=rank_scores([200,100],['removed','survivor'])
        for row,span in zip(rows,([0,10],[10,30])):
            row.update(query_span=span,reference_coverage=1.,identity_fraction=1.)
        def lookup(sha,length):return {'sequence_sha256':sha,'length':length,'ranking':rows}
        def forbidden(*args):raise AssertionError('No fresh alignment expected')
        full={'removed':'old-reference','survivor':'surviving-reference'};remaining={'survivor':full['survivor']}
        ranker,counts=masked_alignment_ranker(remaining,full,lookup,forbidden)
        result=alignment_guided_segments('0123456789'*3,remaining,
            SegmentationConfig(minimum_span=5,maximum_segments=1),ranker=ranker)
        self.assertEqual(result['segments'][0]['seed_hit']['accession'],'survivor')
        self.assertEqual((result['segments'][0]['start0'],result['segments'][0]['end0']),(10,30))
        self.assertEqual(counts['source_full_alignment_hits'],1)

    def test_fresh_aligner_only_receives_remaining_catalogue_and_cache_is_exact(self):
        full={'A':'removed','B':'remaining'};remaining={'B':'remaining'};calls=[]
        def fresh(sequence,catalogue):
            calls.append((sequence,catalogue));return {'ranking':rank_scores([100],['B'])}
        ranker,counts=masked_alignment_ranker(remaining,full,lambda sha,n:None,fresh)
        first=ranker('0123',remaining);second=ranker('0123',remaining)
        self.assertEqual(first,second);self.assertEqual(len(calls),1)
        self.assertEqual(calls[0][1],remaining);self.assertEqual(counts['memory_cache_hits'],1)
        with self.assertRaises(ValueError):ranker('0123',full)

    def test_rrf_filters_each_segment_before_fusion(self):
        one=rank_scores([4,3,2,1],['A','B','C','D'],'cosine_similarity')
        two=rank_scores([2,1,4,3],['A','B','C','D'],'cosine_similarity')
        cfg={'fusion':'rrf','rrf_constant':60}
        masked,fused=masked_fusion([one,two],['B','C','D'],cfg)
        scores={r['accession']:r['score'] for r in fused['ranking']}
        self.assertAlmostEqual(scores['B'],1/61+1/63)
        wrong=fuse_rankings([one,two],'rrf')['ranking']
        wrong_b=next(r['score'] for r in wrong if r['accession']=='B')
        self.assertNotAlmostEqual(scores['B'],wrong_b)
        self.assertTrue(all({r['accession'] for r in rows}=={'B','C','D'} for rows in masked))

    def test_union_keeps_unscorable_family_and_absent_exact_distinct(self):
        rankings=[rank_scores([2,1],['B','C'],'cosine_similarity')]
        q={'components':['A','Z']};labels={'A':['A','B'],'Z':['Z']}
        result=family_union_results(rankings,q,labels,['B','C'])['1']
        e=result['evaluation'];self.assertEqual(e['scorable_components'],1)
        self.assertEqual(e['unscorable_components'],1);self.assertEqual(e['component_recall'],1)
        self.assertIsNone(e['all_components_hit']);self.assertIsNone(e['components'][0]['exact_hit'])
        acc=UnionAccumulator();acc.add('q',result);summary=acc.result()
        self.assertEqual(summary['all_components_scorable_queries'],0)
        self.assertIsNone(summary['all_components_hit'])
        self.assertEqual(summary['mean_returned_hit_occurrences'],1)

    def test_common_mask_must_include_every_component(self):
        q={'group_id':'AB','components':['A','B']};suite={'sequences':{'A':'a','B':'b','C':'c','D':'d'}}
        labels={'removal_policies':{'exact_removed':{'common_removed_accessions_by_group':{'AB':['A','B']},
                    'removed_accessions_by_target':{'A':['A'],'B':['B']}}}}
        self.assertEqual(mask_for_query(q,suite,labels,'exact_removed')['included_accessions'],['C','D'])
        labels['removal_policies']['exact_removed']['common_removed_accessions_by_group']['AB']=['A']
        with self.assertRaises(ValueError):mask_for_query(q,suite,labels,'exact_removed')

    def test_transfer_expected_units_include_every_policy_and_configuration(self):
        identity={'mode':'family_automatic_spans_v1','query_ids':['q1','q2'],
                  'configs':[{'configuration_id':'windows'},{'configuration_id':'alignment'}],
                  'policies':['exact_removed','near90_removed']}
        self.assertEqual(len(expected_transfer_keys(identity)),8)
        identity['mode']='family_automatic_retrieval_v1'
        self.assertEqual(len(expected_transfer_keys(identity)),9)
        self.assertIn('summary',expected_transfer_keys(identity))

    def test_missing_completed_primary_selection_blocks_before_any_prediction(self):
        with tempfile.TemporaryDirectory() as folder:
            session=Path(folder);inputs=session/'inputs';inputs.mkdir()
            primary=inputs/'primary.json';primary.write_text('{}')
            family=inputs/'family.json';family.write_text('{}')
            sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
            design={'version':'family-segmentation-transfer-v1','status':'prospective_frozen_before_transfer_scoring',
                    'primary_design_filename':'primary.json','family_labels_filename':'family.json',
                    'primary_selection_filename':'missing_selection.json','primary_design_sha256':sha(primary),
                    'family_labels_sha256':sha(family)}
            path=inputs/'transfer.json';path.write_text(json.dumps(design))
            with patch('session_family_segmentation.check_deadline'),patch('session_family_segmentation.primary_setup') as setup:
                with self.assertRaisesRegex(ValueError,'selection is not complete'):
                    transfer_setup(SimpleNamespace(session=session,design=path,include_clean=False))
                setup.assert_not_called()

    def test_transfer_requires_exact_selected_development_representation(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);development=root/'development';development.mkdir()
            method={'model_id':'same8M','revision':'same_revision','device':'cpu','dtype':'float32',
                    'layer':6,'num_layers':6,'attention_implementation':'sdpa',
                    'tokenizer_vocab_sha256':'original_vocab','packages':{'torch':'original_version'}}
            identity={'mode':'automatic_segment_retrieval_v1','partition':'dev','model':method}
            identity_path=development/'identity.json';identity_path.write_text(json.dumps(identity))
            selection=root/'selection.json'
            selection.write_text(json.dumps({'development_job':str(development),
                                             'development_identity':digest(canonical(identity))}))
            validate_parent_development_method(dict(method),selection)
            for change in ({'attention_implementation':'eager'}, {'tokenizer_vocab_sha256':'other_vocab'},
                           {'packages':{'torch':'different_version'}}):
                with self.subTest(change=change),self.assertRaisesRegex(ValueError,'exact selected development'):
                    validate_parent_development_method({**method,**change},selection)
            identity_path.write_text(json.dumps({**identity,'model':{**method,'attention_implementation':'eager'}}))
            with self.assertRaisesRegex(ValueError,'development identity changed'):
                validate_parent_development_method(method,selection)


if __name__=='__main__':unittest.main()
