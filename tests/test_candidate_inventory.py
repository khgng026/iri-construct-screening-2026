"""Candidate discovery must never silently become biological-data admission."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.candidate_inventory import (GATES, build_pilot_review_queue, load_candidate_inventory,
    review_queue, validate_candidate_records, validate_pilot_review_queue, load_pilot_shortlist)
from iri_constructs.runtime import InputError


class CandidateInventoryTests(unittest.TestCase):
    def setUp(self):
        self.path = ROOT / 'data/manifests/t06_candidate_inventory.json'
        self.rows = load_candidate_inventory(self.path)

    def test_twenty_candidates_ten_proposed_and_zero_admitted(self):
        queue = review_queue(self.rows)
        self.assertEqual(queue['candidate_count'], 20)
        self.assertEqual(queue['proposed_pilot_count'], 10)
        self.assertEqual(queue['approved_for_scoring'], 0)
        self.assertEqual(set(queue['primary_archetype_counts']),
                         {'single_cds', 'multi_cds', 'tagged_fusion', 'backbone_heavy'})
        self.assertTrue(all(e['status'] == 'excluded_from_scoring' and e['blockers'] for e in queue['entries']))

    def test_approval_cannot_be_entered_in_discovery_schema(self):
        for gate in GATES:
            rows = deepcopy(self.rows)
            rows[0][gate] = 'approved'
            with self.subTest(gate=gate), self.assertRaises(InputError):
                validate_candidate_records(rows)

    def test_pilot_queue_binds_inventory_and_cannot_be_promoted_to_approval(self):
        queue = build_pilot_review_queue(self.path)
        self.assertEqual(queue['record_count'],10)
        self.assertEqual(queue['approved_record_count'],0)
        self.assertFalse(queue['experimental_input_ready'])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'queue.json'
            path.write_text(json.dumps(queue),encoding='utf-8')
            validate_pilot_review_queue(path,self.path)
            queue['approved_record_count']=1
            path.write_text(json.dumps(queue),encoding='utf-8')
            with self.assertRaises(InputError):
                validate_pilot_review_queue(path,self.path)

    def test_duplicate_record_or_source_is_rejected(self):
        rows = deepcopy(self.rows)
        rows.append(deepcopy(rows[0]))
        with self.assertRaisesRegex(InputError, 'Duplicate candidate record_id'):
            validate_candidate_records(rows)
        rows[-1]['record_id'] = 'another-local-id'
        with self.assertRaisesRegex(InputError, 'source identity'):
            validate_candidate_records(rows)

    def test_missing_source_and_credentialed_urls_rejected(self):
        for change in ({'source_url': ''}, {'terms_url': 'file:///local'},
                       {'evidence_urls': []}, {'source_url': 'https://user:password@example.invalid/'}):
            rows = deepcopy(self.rows)
            rows[0].update(change)
            with self.subTest(change=change), self.assertRaises(InputError):
                validate_candidate_records(rows)

    def test_unknown_is_explicit_and_zero_length_is_invalid(self):
        self.assertTrue(any(r['accession_version'] == 'unknown' for r in self.rows))
        for change in ({'sequence_length_nt': 0}, {'feature_count': True}, {'accession_version': ''},
                       {'pilot_priority': 0}, {'checksum_sha256': 'invented'}):
            rows = deepcopy(self.rows)
            rows[0].update(change)
            with self.subTest(change=change), self.assertRaises(InputError):
                validate_candidate_records(rows)
        for field in ('sequence_length_nt', 'feature_count', 'checksum_sha256'):
            rows = deepcopy(self.rows)
            del rows[0][field]
            with self.subTest(missing=field), self.assertRaisesRegex(InputError, 'must be explicit'):
                validate_candidate_records(rows)

    def test_sequence_payload_is_forbidden(self):
        for key in ('sequence', 'translation', 'qualifiers', 'embeddings'):
            rows = deepcopy(self.rows)
            rows[0][key] = 'fictional-payload'
            with self.subTest(key=key), self.assertRaisesRegex(InputError, 'metadata-only'):
                validate_candidate_records(rows)

    def test_wrong_version_empty_inventory_and_invalid_date_fail(self):
        doc = json.loads(self.path.read_text(encoding='utf-8'))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'candidate.json'
            for change in ({'schema_version': '0.1'}, {'records': []}, {'metadata_retrieved_on': 'not-a-date'}):
                path.write_text(json.dumps({**doc, **change}), encoding='utf-8')
                with self.subTest(change=change), self.assertRaises(InputError):
                    load_candidate_inventory(path)

    def test_queue_is_deterministic_and_does_not_mutate_rows(self):
        before = deepcopy(self.rows)
        self.assertEqual(review_queue(self.rows), review_queue(list(reversed(self.rows))))
        self.assertEqual(self.rows, before)

    def test_approved_selection_has_exact_order_but_no_source_admissions(self):
        shortlist = ROOT / 'data/manifests/pilot_shortlist_v1_1.json'
        selection = load_pilot_shortlist(shortlist, self.path)
        queue = build_pilot_review_queue(self.path, shortlist)
        expected = ['ncbi-J01749.1','addgene-212936','addgene-18964','addgene-87121','addgene-18673',
                    'addgene-64034','addgene-108228','addgene-21073','addgene-51970','addgene-255006']
        self.assertEqual([r['record_id'] for r in queue['records']], expected)
        self.assertEqual(queue['approved_record_count'], 0)
        self.assertFalse(queue['experimental_input_ready'])
        self.assertTrue(all(r['status'] == 'excluded_from_scoring' for r in queue['records']))
        self.assertEqual(selection['optional_controls'][0]['record_id'], 'ncbi-L09137.2')
        validate_pilot_review_queue(ROOT / 'data/manifests/t08_pilot_review_queue.json', self.path, shortlist)
        validate_pilot_review_queue(ROOT / 'data/manifests/t08_pilot_review_queue_v1_0.json', self.path)

    def test_selection_manifest_cannot_inject_source_approvals(self):
        original = json.loads((ROOT / 'data/manifests/pilot_shortlist_v1_1.json').read_text(encoding='utf-8'))
        for where in ('top','record','control'):
            selection = deepcopy(original)
            target = selection if where == 'top' else selection['records'][0] if where == 'record' else selection['optional_controls'][0]
            target['review_status'] = 'approved'
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'selection.json'
                path.write_text(json.dumps(selection),encoding='utf-8')
                with self.subTest(where=where), self.assertRaises(InputError):
                    load_pilot_shortlist(path,self.path)

    def test_duplicate_or_unknown_selected_identity_is_rejected(self):
        original = json.loads((ROOT / 'data/manifests/pilot_shortlist_v1_1.json').read_text(encoding='utf-8'))
        for replacement in ('not-a-candidate', original['records'][0]['record_id']):
            selection = deepcopy(original)
            selection['records'][1]['record_id'] = replacement
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'selection.json'
                path.write_text(json.dumps(selection),encoding='utf-8')
                with self.assertRaises(InputError):
                    load_pilot_shortlist(path,self.path)
