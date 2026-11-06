"""E1 uses fictional software fixtures, never admitted biological research records."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.annotation_quality import (audit_quality, classify_quality,
                                               run_annotation_quality, summarize_quality)
from iri_constructs.annotations import normalize_record
from iri_constructs.runtime import InputError
from iri_constructs.runtime import digest_file
from iri_constructs.validate_manifest import REQUIRED_COLUMNS

REVIEW = dict(record_id='fictional-quality', source='fictional-software-test',
              source_url='https://example.invalid/quality', accession='FICTIONAL', accession_version='1',
              archetype='single_cds', notes='Fictional software fixture, not biological source approval')


class QualityClassificationTests(unittest.TestCase):
    def test_fixed_strata_and_bad_counts(self):
        self.assertEqual(classify_quality(dict(source_features=1, source_cds=0, source_issue_features=0)), 'annotation_limited')
        self.assertEqual(classify_quality(dict(source_features=1, source_cds=1, source_issue_features=0)), 'well_annotated')
        self.assertEqual(classify_quality(dict(source_features=2, source_cds=1, source_issue_features=1)), 'partially_annotated')
        for item in ({}, dict(source_features=0, source_cds=1, source_issue_features=0),
                     dict(source_features=True, source_cds=0, source_issue_features=0)):
            with self.assertRaises(InputError):
                classify_quality(item)

    def test_empty_inputs_produce_no_quality_results(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'results') as directory:
            base = Path(directory)
            config = base / 'inputs.json'
            config.write_text(json.dumps(dict(schema_version='t05-inputs-1', records=[])), encoding='utf-8')
            with patch('iri_constructs.annotation_quality.reviewed_sources', side_effect=AssertionError('No admission')):
                report = run_annotation_quality(config, base / 'out')
            self.assertEqual(report['status'], 'no_inputs')
            self.assertEqual(report['quality_records'], [])
            self.assertFalse(report['scientific_scores_produced'])
            self.assertFalse(report['inference_executed'])

    def test_candidate_inventory_cannot_be_used_as_source_admission(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'results') as directory:
            report = run_annotation_quality(ROOT / 'data/manifests/t06_candidate_inventory.json', Path(directory) / 'out')
            self.assertEqual(report['status'], 'blocked')
            self.assertFalse(report['scientific_scores_produced'])

    def test_empty_aggregate_has_null_denominators(self):
        summary = summarize_quality([])
        self.assertIsNone(summary['metrics'])
        self.assertEqual(summary['record_usable_model_input_rate']['denominator'], 0)
        self.assertIsNone(summary['record_usable_model_input_rate']['value'])


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython source quality audits')
class QualityAuditTests(unittest.TestCase):
    def setUp(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        self.Seq, self.Feature, self.Location = Seq, SeqFeature, SimpleLocation
        self.record = SeqRecord(Seq('ATG' + 'GCT'*2 + 'TAA'), id='FICTIONAL.1')
        self.record.annotations = dict(molecule_type='DNA')
        self.record.features = [SeqFeature(SimpleLocation(0,12,strand=1), type='CDS',
            qualifiers={'translation':['M'+'A'*2], 'product':['fictional fixture']})]

    def audit(self, **kwargs):
        doc, sidecar, proteins = normalize_record(self.record, REVIEW, 'a'*64)
        return audit_quality(self.record, doc, sidecar, proteins, source_sha256='a'*64, **kwargs)

    def test_complete_annotation_and_traceability(self):
        result = self.audit()
        self.assertEqual(result['quality_stratum'], 'well_annotated')
        self.assertEqual(result['counts']['erroneous_invention_count'], 0)
        self.assertEqual(result['metrics']['useful_model_input_rate'], dict(numerator=1,denominator=1,value=1.0))
        self.assertTrue(result['source_checksum_binding'])
        self.assertEqual(result['source_preservation_audit_status'], 'agrees_with_source')

    def test_source_only_and_empty_features_are_not_perfect_preservation(self):
        for features in ([], [self.Feature(self.Location(0,12,strand=1), type='source')]):
            self.record.features = features
            result = self.audit()
            self.assertEqual(result['quality_stratum'], 'annotation_limited')
            self.assertEqual(result['status'], 'no_observable_annotations')
            self.assertIsNone(result['metrics']['observable_feature_preservation']['value'])
            self.assertIsNone(result['metrics']['useful_model_input_rate']['value'])
            self.assertTrue(result['record_abstention_justified'])

    def test_missing_and_conflicting_annotations_are_preserved_abstentions(self):
        self.record.features[0].qualifiers.pop('translation')
        self.record.features.append(self.Feature(self.Location(0,12,strand=1), type='CDS',
            qualifiers={'translation':['M'+'C'*2]}))
        result = self.audit()
        self.assertEqual(result['status'], 'audited')
        self.assertEqual(result['quality_stratum'], 'partially_annotated')
        self.assertEqual(result['counts']['missing_translation_cds'], 1)
        self.assertEqual(result['counts']['translation_conflict_cds'], 1)
        self.assertEqual(result['metrics']['justified_model_abstention_rate'], dict(numerator=2,denominator=2,value=1.0))
        self.assertEqual(result['counts']['erroneous_invention_count'], 0)

    def test_fuzzy_and_unknown_strands_are_not_silently_exact(self):
        from Bio.SeqFeature import BeforePosition
        for location, counter in [(self.Location(BeforePosition(0),12,strand=1),'fuzzy_location_features'),
                                  (self.Location(0,12),'unknown_strand_features')]:
            self.record.features[0].location = location
            result = self.audit()
            self.assertEqual(result['counts'][counter], 1)
            self.assertEqual(result['quality_stratum'], 'partially_annotated')
            self.assertEqual(result['counts']['usable_model_input_cds'], 0)
            self.assertTrue(result['record_abstention_justified'])

    def test_reverse_multipart_can_be_well_annotated(self):
        from Bio.SeqFeature import CompoundLocation
        sequence = str(self.record.seq)
        self.record.seq = self.Seq(sequence[6:]+sequence[:6]).reverse_complement()
        self.record.features[0].location = CompoundLocation([
            self.Location(0,6,strand=-1), self.Location(6,12,strand=-1)])
        result = self.audit()
        self.assertEqual(result['quality_stratum'], 'well_annotated')
        self.assertEqual(result['counts']['multipart_features'], 1)
        self.assertEqual(result['counts']['reverse_strand_features'], 1)
        self.assertEqual(result['counts']['usable_model_input_cds'], 1)

    def test_cap_is_model_compatibility_not_source_quality(self):
        result = self.audit(max_residues=2)
        self.assertEqual(result['quality_stratum'], 'well_annotated')
        self.assertEqual(result['counts']['unsupported_length_cds'], 1)
        self.assertEqual(result['counts']['usable_model_input_cds'], 0)
        self.assertEqual(result['metrics']['justified_model_abstention_rate']['value'], 1.0)

    def test_feature_and_protein_invention_detected(self):
        doc, sidecar, proteins = normalize_record(self.record, REVIEW, 'a'*64)
        invented = copy.deepcopy(doc['features'][0])
        invented['feature_id'] = 'fabricated-software-error'
        doc['features'].append(invented)
        proteins['untraceable-export'] = 'M'+'A'*2
        result = audit_quality(self.record,doc,sidecar,proteins)
        self.assertEqual(result['counts']['invented_feature_count'], 1)
        self.assertEqual(result['counts']['unsupported_protein_export_count'], 1)
        self.assertEqual(result['status'], 'failed')

    def test_parser_loss_does_not_reclassify_good_source_as_bad(self):
        doc, sidecar, proteins = normalize_record(self.record, REVIEW, 'a'*64)
        result = audit_quality(self.record, doc, sidecar, {})
        self.assertEqual(result['quality_stratum'], 'well_annotated')
        self.assertEqual(result['status'], 'failed')
        self.assertIn('valid_source_translation_lost', result['failed_checks'])

    def test_bad_provenance_and_uncertainty_claims_are_failed(self):
        from Bio.SeqFeature import BeforePosition
        self.record.features[0].location = self.Location(BeforePosition(0),12,strand=1)
        doc, sidecar, proteins = normalize_record(self.record, REVIEW, 'a'*64)
        doc['features'][0]['annotation_status'] = 'known'
        result = audit_quality(self.record, doc, sidecar, proteins, source_sha256='b'*64)
        self.assertIn('source_uncertainty_not_preserved', result['failed_checks'])
        self.assertIn('source_checksum_binding_failed', result['failed_checks'])

    def run_fixture_pipeline(self, *, max_residues=1022, failed_audit=False):
        from Bio import SeqIO
        with tempfile.TemporaryDirectory(dir=ROOT / 'results') as directory:
            base = Path(directory)
            source = base / 'fictional.gb'
            SeqIO.write(self.record,source,'genbank')
            review = {**REVIEW, 'sequence_length_nt':str(len(self.record)), 'feature_count':str(len(self.record.features)),
                'annotation_status':'approved', 'review_status':'approved', 'licensing_status':'approved',
                'safety_review':'approved', 'retrieved_on':'2026-11-05', 'checksum_sha256':digest_file(source)}
            review['notes'] = 'Fictional unit-test gates only, never a real-source permission or annotation approval'
            with (base / 'review.csv').open('w',encoding='utf-8',newline='') as handle:
                writer = csv.DictWriter(handle,fieldnames=REQUIRED_COLUMNS)
                writer.writeheader(); writer.writerow(review)
            config = base / 'inputs.json'
            config.write_text(json.dumps({'schema_version':'t05-inputs-1','review_manifest':'review.csv',
                'records':[dict(input_id='fixture',record_id=review['record_id'],path=source.name,
                                format='genbank',molecule_type='dna',sequence_id=self.record.id)]}),encoding='utf-8')
            if failed_audit:
                original = audit_quality
                def deliberately_failed(*args,**kwargs):
                    item = original(*args,**kwargs)
                    item.update(status='failed',failed_checks=['fictional_preservation_failure'])
                    item['counts']['erroneous_invention_count'] = 1
                    return item
                with patch('iri_constructs.annotation_quality.audit_quality',side_effect=deliberately_failed):
                    report = run_annotation_quality(config,base/'out',max_residues=max_residues)
            else:
                report = run_annotation_quality(config,base/'out',max_residues=max_residues)
            self.assertNotIn(str(self.record.seq),json.dumps(report))
            return report

    def test_full_gated_pipeline_and_declared_cap(self):
        report = self.run_fixture_pipeline(max_residues=2)
        self.assertEqual(report['status'],'ok')
        self.assertEqual(report['configuration']['max_residues'],2)
        self.assertEqual(report['quality_records'][0]['counts']['unsupported_length_cds'],1)
        self.assertEqual(report['overall']['metrics']['useful_model_input_rate']['value'],0.0)
        self.assertFalse(report['inference_executed'])
        self.assertFalse(report['scientific_scores_produced'])
        self.assertTrue(report['annotation_measurements_produced'])


    def test_source_only_pipeline_never_claims_observable_measurement(self):
        self.record.features=[self.Feature(self.Location(0,12,strand=1),type='source')]
        report=self.run_fixture_pipeline()
        self.assertEqual(report['status'],'no_observable_annotations')
        self.assertFalse(report['scientific_scores_produced'])
        self.assertIsNone(report['overall']['metrics']['observable_feature_preservation']['value'])

    def test_failed_audits_are_not_scientific_successes(self):
        report=self.run_fixture_pipeline(failed_audit=True)
        self.assertEqual(report['status'],'partial')
        self.assertFalse(report['scientific_scores_produced'])
        self.assertEqual(report['overall']['records'],0)
        self.assertEqual(report['parser_failure_summary']['erroneous_invention_count'],1)
