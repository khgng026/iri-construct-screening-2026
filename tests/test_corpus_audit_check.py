"""Check saved audit metadata using fictional records and coordinates."""
from copy import deepcopy
import unittest

from experiments.paired_preprocessing.check_corpus_audit import check_audit


def _fixture():
    versions = ['record-a', 'record-b']
    acquisition = {'schema': 'paired-corpus-acquisition-1', 'samples': {'sample': {'files': [{
        'file': 'batch-01', 'byte_identical': True, 'sha256': 'batch-hash', 'frozen_sha256': 'batch-hash',
        'frozen_versions': versions, 'returned_versions': list(versions),
    }]}}}
    records = []
    for accession, offset in zip(versions, (0, 40), strict=True):
        region = [[offset + 3, offset + 6], [offset + 9, offset + 15]]
        records.append({
            'sample': 'sample', 'accession': accession, 'status': 'ok', 'batch_file': 'batch-01', 'length': 100,
            'named_subregions': [{'location_status': 'exact', 'strand': '+', 'segments': region,
                                  'named_families': ['family-a']}] if offset == 0 else [],
            'cds': [{
                'extraction_state': 'exported', 'protein_sha256': 'shared-protein', 'frozen_unit_match': True,
                'frozen_label_match': True, 'strand': 1, 'segments': [[offset, offset + 6], [offset + 9, offset + 24]],
                'length_aa': 7, 'anchor_evidence': {'family-a': {'unit_region_aa': [1, 4], 'dna_spans': region,
                                                                'dna_region_retranslates': True}},
            }],
        })
    audit = {'schema': 'paired-corpus-audit-1', 'constructs': records,
             'checks': {key: 2 for key in ('records', 'exported_units', 'frozen_unit_match', 'frozen_label_match',
                                           'anchor_regions', 'anchor_regions_retranslate')}}
    return acquisition, audit


class AuditMetadataTests(unittest.TestCase):
    def test_accepts_valid_joined_regions_and_separates_source_support(self):
        acquisition, audit = _fixture()
        report = check_audit(acquisition, audit)
        self.assertEqual(report['findings'], {})
        self.assertEqual(report['saved_counts']['anchor_regions_matching_a_named_source_feature'], 1)
        self.assertEqual(report['shared_protein_group_diagnostics'],
                         {'groups': 1, 'largest_group': 2, 'distinct_protein_hashes': 1})

    def test_catches_overwritten_coordinates_despite_passing_counters(self):
        acquisition, audit = _fixture()
        left, right = [record['cds'][0] for record in audit['constructs']]
        left['anchor_evidence'] = deepcopy(right['anchor_evidence'])
        report = check_audit(acquisition, audit)
        self.assertEqual(report['findings'], {'regions_outside_parent_cds': 1})
        self.assertEqual(report['saved_counts']['anchor_regions_retranslate'], 2)

    def test_rejects_a_region_that_covers_the_gap_in_a_joined_cds(self):
        acquisition, audit = _fixture()
        audit['constructs'][0]['cds'][0]['anchor_evidence']['family-a']['dna_spans'] = [[6, 15]]
        self.assertEqual(check_audit(acquisition, audit)['findings'], {'regions_outside_parent_cds': 1})

    def test_equal_total_record_count_can_still_hide_a_missing_record(self):
        acquisition, audit = _fixture()
        audit['constructs'][1] = deepcopy(audit['constructs'][0])
        self.assertEqual(check_audit(acquisition, audit)['findings'],
                         {'missing_audit_records': 1, 'duplicate_audit_records': 1})

    def test_detects_nonidentical_reported_snapshot(self):
        acquisition, audit = _fixture()
        batch = acquisition['samples']['sample']['files'][0]
        batch['sha256'] = 'different'
        batch['returned_versions'][1] = 'different-version'
        self.assertEqual(check_audit(acquisition, audit)['findings'],
                         {'nonidentical_batch_files': 1, 'batch_version_list_mismatches': 1})


if __name__ == '__main__':
    unittest.main()
