"""Metadata admission gates and deterministic nesting; no network or proteins."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from session_expand_catalogue import (frozen_json,metadata_eligibility_errors,nested_prefixes,
    next_link,selection_key,sequence_candidates,source_family_statements)


class ExpandedCatalogueTest(unittest.TestCase):
    def source(self):
        return {'accession':'A001','reviewed':'reviewed','organism_id':'559292','fragment':'',
                'version':'4','sequence_version':'2','length':'42','protein_existence':'Evidence at protein level'}

    def fixture(self):
        design={'selection_seed':'fixture-seed','license_evidence_sha256':'license'}
        index={'file_sha256':'metadata','records':{a:{'metadata':{**self.source(),'accession':a},
                 'metadata_row_sha256':'row-'+a} for a in ('A001','A002','A003')}}
        review={'metadata_index_sha256':'metadata','records':{a:{'metadata_row_sha256':'row-'+a,
                 'review_status':'approved','safety_review':'approved','licensing_status':'approved',
                 'rationale':'Specific ordinary role reviewed for '+a,'source_annotation_basis':{'function':'fixture'},
                 'license_evidence_sha256':'license'} for a in index['records']}}
        return design,index,review

    def test_next_link_preserves_commas_and_escaped_relation(self):
        url='https://rest.uniprot.org/uniprotkb/search?fields=accession,id,reviewed&cursor=next'
        self.assertEqual(next_link({'Link':'<'+url+'>; rel="next", <https://example.invalid/last>; rel="last"'}),url)
        self.assertEqual(next_link({'link':'<'+url+'>; rel=\\"next\\"'}),url)
        self.assertIsNone(next_link({'Link':'<https://example.invalid/last>; rel="last"'}))

    def test_source_review_is_required_per_record_and_base_is_not_refetched(self):
        design,index,review=self.fixture();review['records']['A002']['safety_review']='pending'
        self.assertEqual(sequence_candidates(index,review,design,{'A001'}),['A003'])
        del review['records']['A003']['rationale']
        with self.assertRaisesRegex(ValueError,'individual evidence'):sequence_candidates(index,review,design,set())

    def test_review_source_and_licence_must_match(self):
        design,index,review=self.fixture()
        for key,value,pattern in [('metadata_row_sha256','changed','metadata row'),
                                  ('license_evidence_sha256','changed','licence')]:
            candidate=copy.deepcopy(review);candidate['records']['A002'][key]=value
            with self.assertRaisesRegex(ValueError,pattern):sequence_candidates(index,candidate,design,set())
        review['metadata_index_sha256']='different'
        with self.assertRaisesRegex(ValueError,'different source'):sequence_candidates(index,review,design,set())

    def test_approval_cannot_override_mechanical_source_gate(self):
        design,index,review=self.fixture();index['records']['A002']['metadata']['protein_existence']='Uncertain'
        with self.assertRaisesRegex(ValueError,'mechanical source gate'):sequence_candidates(index,review,design,set())
        for change in ({'reviewed':'unreviewed'},{'organism_id':'other'},{'accession':'A001-2'},
                       {'fragment':'Fragment'},{'version':''},{'sequence_version':'unknown'},
                       {'length':'1023'},{'length':'0'}):
            with self.subTest(change=change):self.assertTrue(metadata_eligibility_errors({**self.source(),**change}))
        self.assertEqual(metadata_eligibility_errors(self.source()),[])

    def test_hash_order_and_nested_achieved_endpoints(self):
        design,index,review=self.fixture();ordered=sequence_candidates(index,review,design,set())
        self.assertEqual(ordered,sorted(index['records'],key=lambda a:selection_key(a,design)))
        panels=nested_prefixes(['base1','base2'],ordered,[3,4,8])
        self.assertEqual(panels['3']['accessions'],panels['4']['accessions'][:3])
        self.assertEqual(panels['8']['achieved_unique_sequences'],5)
        self.assertFalse(panels['8']['endpoint_reached'])
        with self.assertRaises(ValueError):nested_prefixes(['base'],['base'],[2])

    def test_family_characterization_uses_explicit_statement_only(self):
        self.assertEqual(source_family_statements('SIMILARITY: Belongs to the Example family. A subfamily. {ECO:0000}'),
                         ['Belongs to the Example family'])
        self.assertEqual(source_family_statements('Contains a shared domain. Similar ordinary function.'),[])

    def test_freeze_refuses_silent_replacement(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'frozen.json';frozen_json(path,{'value':1});frozen_json(path,{'value':1})
            with self.assertRaisesRegex(ValueError,'different contents'):frozen_json(path,{'value':2})
            self.assertEqual(json.loads(path.read_text()),{'value':1})


if __name__=='__main__':unittest.main()
