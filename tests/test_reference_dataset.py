"""Fictional reference-evidence guards, never scientific family approvals."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from iri_constructs.reference_dataset import (build_source_queries, canonical_hash, freeze_proposal,
                                              proposed_eligibility, summarize_reference)
from iri_constructs.runtime import InputError


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


class ReferenceDatasetTests(unittest.TestCase):
    def setUp(self):
        self.spec = {'accession':'Q00000', 'taxon_id':123, 'proposed_family':'gfp_like_homologs',
                     'lineage_group':'fictional-natural-origin',
                     'independence_caveat':'Fictional metadata for software validation only'}
        self.document = {
            'primaryAccession':'Q00000', 'entryType':'UniProtKB reviewed (Swiss-Prot)',
            'uniProtkbId':'FICTIONAL', 'sequence':{'value':'ACDEFGHIKL', 'length':10},
            'entryAudit':{'entryVersion':1, 'sequenceVersion':1},
            'organism':{'taxonId':123, 'scientificName':'Fictional reference'},
            'proteinDescription':{'recommendedName':{'fullName':{'value':'Fictional homolog'}}},
            'uniProtKBCrossReferences':[{'database':'Pfam','id':'PF01353'}],
            'comments':[{'commentType':'SIMILARITY','texts':[{'value':'Belongs to the GFP family'}]}],
            'references':[{'citation':{'title':'Fictional citation', 'citationCrossReferences':[
                {'database':'PubMed','id':'fictional'}]}}],
        }

    def summary(self, doc=None, spec=None):
        return summarize_reference(doc or self.document, spec or self.spec,
                                   raw_sha256='a'*64, retrieved_at_utc='2026-11-11T00:00:00Z')

    def test_exact_reviewed_provenance_and_no_sequences(self):
        summary = self.summary()
        self.assertEqual(summary['protein_sequence_sha256'], digest('ACDEFGHIKL'))
        self.assertEqual(summary['entry_version'], 1)
        self.assertNotIn('ACDEFGHIKL', json.dumps(summary))
        self.assertFalse(summary['family_panel_approved'])
        self.assertFalse(summary['run_authorized'])
        self.assertFalse(summary['source_admitted_for_scoring'])

    def test_identity_version_sequence_taxon_or_review_mismatch_rejected(self):
        for changes in ({'primaryAccession':'Q00001'}, {'entryType':'Unreviewed'},
                        {'entryAudit':{}}, {'entryAudit':{'entryVersion':True,'sequenceVersion':1}},
                        {'sequence':{'value':'ACD*', 'length':4}},
                        {'sequence':{'value':'ACDE', 'length':5}},
                        {'organism':{'taxonId':999}}):
            with self.subTest(changes=changes), self.assertRaises(InputError):
                self.summary({**self.document, **changes})

    def test_shared_function_or_name_cannot_replace_family_evidence(self):
        for changes in ({'uniProtKBCrossReferences':[]}, {'comments':[]}, {'references':[]}):
            with self.subTest(changes=changes), self.assertRaises(InputError):
                self.summary({**self.document, **changes})
        spec = {**self.spec,'proposed_family':'coleopteran_luciferase_homologs'}
        with self.assertRaises(InputError):
            self.summary(spec=spec)

    def query(self):
        return {'query_id':'q1', 'parent_id':'parent1', 'protein_sequence_sha256':digest('query'),
                'proposed_family':'gfp_like_homologs', 'lineage_group':'query-origin',
                'family_evidence':'Fictional evidence', 'known_ancestor_accessions':['A00000'],
                'all_parent_component_sha256':[digest('query'),digest('unreviewed-other-component')]}

    def refs(self):
        return [{'reference_id':'r'+str(i), 'accession':'A0000'+str(i),
                 'protein_sequence_sha256':digest('reference'+str(i)),
                 'lineage_group':'origin'+str(i), 'proposed_family':('gfp_like_homologs' if i<5
                   else 'coleopteran_luciferase_homologs')} for i in range(9)]

    def test_all_component_copies_parent_lineages_and_ancestors_excluded(self):
        refs = self.refs()
        refs[1]['protein_sequence_sha256'] = digest('unreviewed-other-component')
        refs[2]['parent_id'] = 'parent1'
        refs[3]['lineage_group'] = 'query-origin'
        result = proposed_eligibility(self.query(), refs)
        self.assertEqual(result['eligible_count'], 5)
        self.assertEqual(set(result['removed_references']), {'r0','r1','r2','r3'})
        self.assertIn('exact_copy_of_any_documented_parent_component',result['removed_references']['r1'])
        self.assertIn('same_source_parent',result['removed_references']['r2'])
        self.assertIn('shared_origin_or_vector_lineage',result['removed_references']['r3'])
        self.assertTrue(result['budget_saturated'])
        self.assertFalse(result['scientific_clearance'])

    def test_missing_family_or_lineage_cannot_become_independence(self):
        for key in ['family_evidence','lineage_group','proposed_family','all_parent_component_sha256']:
            query = self.query()
            query.pop(key)
            with self.subTest(key=key),self.assertRaises(InputError):
                proposed_eligibility(query,self.refs())

    def test_insufficient_positive_explicit_hold(self):
        refs = self.refs()[5:]
        result = proposed_eligibility(self.query(),refs)
        self.assertEqual(result['positive_count'],0)
        self.assertEqual(result['support_status'],'HOLD_NO_INDEPENDENT_POSITIVE_OR_COMPETITOR')
        self.assertFalse(result['scoring_performed'])

    def test_freeze_is_deterministic_and_fingerprint_bound(self):
        first = freeze_proposal(self.refs(),[self.query()],source_manifest_sha256='a'*64,
                                catalogue_spec_sha256='b'*64)
        reverse = freeze_proposal(list(reversed(self.refs())),[self.query()],source_manifest_sha256='a'*64,
                                  catalogue_spec_sha256='b'*64)
        self.assertEqual(first,reverse)
        payload = copy.deepcopy(first)
        declared = payload.pop('proposal_sha256')
        self.assertEqual(canonical_hash(payload),declared)
        changed = freeze_proposal(self.refs(),[self.query()],source_manifest_sha256='c'*64,
                                  catalogue_spec_sha256='b'*64)
        self.assertNotEqual(first['proposal_sha256'], changed['proposal_sha256'])
        self.assertFalse(first['scientific_clearance'])

    def test_duplicate_roles_fail(self):
        refs = self.refs()
        with self.assertRaises(InputError):
            freeze_proposal(refs+refs[:1],[self.query()],source_manifest_sha256='a'*64,catalogue_spec_sha256='b'*64)
        with self.assertRaises(InputError):
            freeze_proposal(refs,[self.query(),self.query()],source_manifest_sha256='a'*64,catalogue_spec_sha256='b'*64)
        refs[1]['protein_sequence_sha256'] = refs[0]['protein_sequence_sha256']
        with self.assertRaisesRegex(InputError, 'duplicate reference proteins'):
            freeze_proposal(refs,[self.query()],source_manifest_sha256='a'*64,catalogue_spec_sha256='b'*64)

    def source_bundle(self):
        components = [
            {'protein_accession_version':'FICTIONAL.1','translation_state':'known',
             'translation_sha256':digest('query'),'translation_length_aa':10,
             'segments':[{'start0':0,'end0':33,'strand':1}],'product':'Fictional annotated component'},
            {'protein_accession_version':'FICTIONAL.2','translation_state':'known',
             'translation_sha256':digest('unsupported'),'translation_length_aa':20,
             'segments':[{'start0':40,'end0':103,'strand':1}],'product':'Fictional unreviewed marker'},
        ]
        source = {'source_id':'parent1','status':'READY_FOR_REVIEW','eligible_for_preparatory_QA':True,
                  'source_url':'https://example.org/fictional','source_retrieval':{'sha256':'a'*64},
                  'cds_evidence':components}
        unit = {'input_id':'parent1/unit1','source_id':'parent1','source_sha256':'a'*64,
                'sequence_sha256':digest('query'),'length_aa':10,'feature_ids':['feature1'],
                'source_segments':components[0]['segments']}
        spec = {'source_id':'parent1','protein_accession_version':'FICTIONAL.1',
                'proposed_family':'gfp_like_homologs','lineage_group':'fictional-origin',
                'known_ancestor_accessions':['Q00000'],'family_evidence':'Fictional evidence only'}
        return {'records':[source]},{'units':[unit]},{'queries':[spec]}

    def test_source_bound_queries_retain_unsupported_full_denominator(self):
        result = build_source_queries(*self.source_bundle())
        query = result['queries'][0]
        self.assertEqual(query['query_id'],'parent1/unit1')
        self.assertEqual(set(query['all_parent_component_sha256']),{digest('query'),digest('unsupported')})
        coverage = result['parent_component_coverage']['parent1']
        self.assertEqual(coverage['source_documented_cds_count'],2)
        self.assertEqual(coverage['proposed_supported_component_count'],1)
        self.assertEqual(coverage['source_longest_translation_aa'],20)
        self.assertEqual(coverage['proposed_eligible_longest_translation_aa'],10)
        self.assertEqual(len(coverage['unsupported_components']),1)
        self.assertFalse(query['family_panel_approved'])
        self.assertFalse(query['run_authorized'])

    def test_held_source_or_changed_provenance_blocks_query(self):
        source,units,specs = self.source_bundle()
        source['records'][0]['eligible_for_preparatory_QA'] = False
        with self.assertRaisesRegex(InputError,'held'):
            build_source_queries(source,units,specs)
        source,units,specs = self.source_bundle()
        units['units'][0]['source_sha256'] = 'b'*64
        with self.assertRaisesRegex(InputError,'original source evidence'):
            build_source_queries(source,units,specs)

    def test_query_requires_exact_unit_and_source_translation(self):
        source,units,specs = copy.deepcopy(self.source_bundle())
        units['units'][0]['source_segments'] = [{'start0':1,'end0':33,'strand':1}]
        with self.assertRaisesRegex(InputError,'deterministic normalized unit'):
            build_source_queries(source,units,specs)


if __name__ == '__main__':
    unittest.main()
