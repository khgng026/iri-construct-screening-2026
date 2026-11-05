"""Generated fictional software fixtures, never biological source approvals/results."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.annotations import normalize_record
from iri_constructs.pilot import audit_annotation_preservation, prepare_pilot
from iri_constructs.runtime import digest_file
from iri_constructs.validate_manifest import REQUIRED_COLUMNS


class PilotFixture:
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'results')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.path = self.base / 'fictional.fasta'
        self.path.write_text('>FICTIONAL.1\n' + 'A'*12 + '\n', encoding='utf-8')
        self.row = dict(record_id='fictional-one', source='fictional-software-fixture',
            source_url='https://example.invalid/fixture', accession='FICTIONAL', accession_version='1',
            archetype='single_cds', sequence_length_nt='12', feature_count='0', annotation_status='approved',
            review_status='approved', licensing_status='approved', safety_review='approved',
            retrieved_on='2026-11-05', checksum_sha256=digest_file(self.path),
            notes='Fictional unit-test statuses, never a real biological source approval')
        self.entries = [dict(input_id='fictional', record_id='fictional-one', path=self.path.name,
                             format='fasta', molecule_type='dna', sequence_id='FICTIONAL.1')]

    def config(self, entries=None, rows=None):
        with (self.base / 'review.csv').open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
            writer.writeheader()
            writer.writerows(rows if rows is not None else [self.row])
        path = self.base / 'inputs.json'
        path.write_text(json.dumps(dict(schema_version='t05-inputs-1', review_manifest='review.csv',
            records=self.entries if entries is None else entries)), encoding='utf-8')
        return path

    def run_pipeline(self, **kwargs):
        return prepare_pilot(self.config(**kwargs), self.base / 'output')


class PilotGateTests(PilotFixture, unittest.TestCase):
    def test_empty_panel_is_no_inputs_without_reading_any_source(self):
        config = self.base / 'empty.json'
        config.write_text(json.dumps(dict(schema_version='t05-inputs-1', records=[])), encoding='utf-8')
        with patch('iri_constructs.pilot.reviewed_sources', side_effect=AssertionError('Must not load sources')):
            report = prepare_pilot(config, self.base / 'output')
        self.assertEqual(report['status'], 'no_inputs')
        self.assertEqual(report['scored_records'], 0)
        self.assertEqual(sorted(p.name for p in (self.base / 'output').iterdir()), ['report.json'])

    def test_pending_uncertain_rejected_gates_block_before_parsing(self):
        for gate in ('review_status','licensing_status','safety_review','annotation_status'):
            for state in ('pending','uncertain','rejected'):
                with self.subTest(gate=gate,state=state), tempfile.TemporaryDirectory(dir=self.base) as output:
                    row = {**self.row, gate: state}
                    config = self.config(rows=[row])
                    with patch('iri_constructs.pilot.read_source', side_effect=AssertionError('Must not parse unapproved source')):
                        report = prepare_pilot(config, Path(output))
                    self.assertEqual(report['status'], 'blocked')
                    self.assertEqual(report['scored_records'], 0)
                    self.assertEqual(sorted(p.name for p in Path(output).iterdir()), ['report.json'])

    def test_incomplete_provenance_and_bad_checksum_block(self):
        for change in ({'source_url': ''}, {'notes': ''}, {'accession_version': ''},
                       {'checksum_sha256': '0'*64}, {'retrieved_on': 'unknown'}):
            with self.subTest(change=change), tempfile.TemporaryDirectory(dir=self.base) as output:
                with patch('iri_constructs.pilot.read_source', side_effect=AssertionError('Must not parse')):
                    report = prepare_pilot(self.config(rows=[{**self.row, **change}]), Path(output))
                self.assertEqual(report['status'], 'blocked')

    def test_duplicate_parent_and_path_record_aliases_block(self):
        variants = [
            [self.entries[0], {**self.entries[0], 'input_id': 'alias'}],
            [self.entries[0], {**self.entries[0], 'input_id': 'alias', 'record_id': 'other-parent'}],
        ]
        for entries in variants:
            with tempfile.TemporaryDirectory(dir=self.base) as output:
                report = prepare_pilot(self.config(entries=entries), Path(output))
                self.assertEqual(report['status'], 'blocked')
                self.assertIn('Duplicate', report['reason'])

    def test_identical_file_copies_cannot_inflate_parent_counts(self):
        copy_path = self.base / 'copy.fasta'
        copy_path.write_bytes(self.path.read_bytes())
        alias = {**self.entries[0], 'input_id': 'alias', 'record_id': 'other-parent', 'path': copy_path.name}
        rows = [self.row, {**self.row, 'record_id': 'other-parent'}]
        with patch('iri_constructs.pilot.read_source', side_effect=AssertionError('Must not parse duplicated source')):
            report = self.run_pipeline(entries=self.entries + [alias], rows=rows)
        self.assertEqual(report['status'], 'blocked')
        self.assertIn('Duplicate checksummed', report['reason'])

    def test_protein_only_records_are_not_construct_annotation_inputs(self):
        report = self.run_pipeline(entries=[{**self.entries[0], 'molecule_type': 'protein'}])
        self.assertEqual(report['status'], 'blocked')
        self.assertIn('parent DNA', report['reason'])

    def test_invalid_identifier_cannot_inject_export_headers(self):
        report = self.run_pipeline(entries=[{**self.entries[0], 'input_id': 'bad\nheader'}])
        self.assertEqual(report['status'], 'blocked')
        self.assertIn('metadata tokens', report['reason'])

    def test_cli_no_inputs_is_nonzero_and_path_independent(self):
        path = self.config(entries=[])
        output = self.base / 'output'
        proc = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'prepare_pilot.py'),
            '--inputs', str(path), '--output', str(output)], cwd=self.base, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(json.loads((output / 'report.json').read_text())['status'], 'no_inputs')


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython parser and preservation tests')
class PilotParsingTests(PilotFixture, unittest.TestCase):
    def setUp(self):
        super().setUp()
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        self.Seq, self.SeqIO = Seq, SeqIO
        self.Feature, self.Location = SeqFeature, SimpleLocation
        self.record = SeqRecord(Seq('ATG' + 'GCT'*2 + 'TAA'), id='FICTIONAL.1', description='Fictional software fixture')
        self.record.annotations = {'molecule_type': 'DNA', 'topology': 'linear'}
        self.record.features = [SeqFeature(SimpleLocation(0,12,strand=1), type='CDS',
            qualifiers={'translation': ['M'+'A'*2], 'product': ['fictional fixture']}),
            SeqFeature(SimpleLocation(0,3,strand=1), type='promoter', qualifiers={'note': ['fictional noncoding fixture']})]

    def prepare_record(self):
        self.path = self.base / 'fictional.gb'
        self.SeqIO.write(self.record, self.path, 'genbank')
        self.entries[0].update(path=self.path.name, format='genbank')
        self.row.update(sequence_length_nt=str(len(self.record)), feature_count=str(len(self.record.features)),
                        checksum_sha256=digest_file(self.path))

    def test_preservation_and_translation_exports_are_auditable(self):
        self.prepare_record()
        report = self.run_pipeline()
        self.assertEqual(report['status'], 'ok')
        self.assertEqual(report['scored_records'], 1)
        counts = report['experiment_a']['overall']['counts']
        self.assertEqual(counts['source_features'], 2)
        self.assertEqual(counts['original_qualifiers_preserved'], 2)
        self.assertEqual(counts['coordinate_agreement_features'], 2)
        self.assertEqual(counts['annotated_translations_verified'], 1)
        self.assertEqual(counts['retained_noncoding_features'], 1)
        output = self.base / 'output'
        with (output/'units.fasta').open(encoding='utf-8') as stream:
            self.assertEqual(len(list(self.SeqIO.parse(stream,'fasta'))), 1)
        text = (output/'report.json').read_text()
        self.assertNotIn(str(self.record.seq), text)
        self.assertFalse(report['inference_executed'])

    def test_reverse_joined_coordinates_preserved(self):
        from Bio.SeqFeature import CompoundLocation
        coding = str(self.record.seq)
        rotated = coding[6:] + coding[:6]
        self.record.seq = self.Seq(rotated).reverse_complement()
        self.record.features = [self.Feature(CompoundLocation([
            self.Location(0,6,strand=-1), self.Location(6,12,strand=-1)]), type='CDS',
            qualifiers={'translation': ['M'+'A'*2]})]
        self.record.annotations['topology'] = 'circular'
        self.prepare_record()
        report = self.run_pipeline()
        self.assertEqual(report['status'], 'ok')
        self.assertEqual(report['records'][0]['counts']['strand_agreement_features'], 1)
        doc = json.loads((self.base/'output'/'normalized.json').read_text())[0]
        self.assertIsNone(doc['features'][0]['location']['start0'])
        self.assertEqual([s['strand'] for s in doc['features'][0]['segments']], [-1,-1])

    def test_missing_and_conflicting_translation_rates(self):
        self.record.features[0].qualifiers.pop('translation')
        self.record.features.append(self.Feature(self.Location(0,12,strand=1), type='CDS',
            qualifiers={'translation': ['M'+'C'*2]}))
        self.prepare_record()
        report = self.run_pipeline()
        counts = report['experiment_a']['overall']['counts']
        self.assertEqual(counts['source_cds'], 2)
        self.assertEqual(counts['missing_translation_cds'], 1)
        self.assertEqual(counts['ambiguous_cds'], 1)
        self.assertEqual(counts['withheld_cds'], 2)
        self.assertEqual(counts['annotated_proteins_emitted'], 0)
        self.assertIsNone(report['experiment_a']['overall']['rates']['exported_translation_agreement'])

    def test_corrupt_normalized_coordinate_or_qualifier_is_detected(self):
        document, audit, extracted = normalize_record(self.record, self.row, self.row['checksum_sha256'])
        cds = next(f for f in document['features'] if f['feature_label']=='CDS')
        broken = copy.deepcopy(document)
        next(f for f in broken['features'] if f['feature_label']=='CDS')['location']['start0'] = 1
        result = audit_annotation_preservation(self.record, broken, audit, extracted)
        self.assertIn('coordinate_agreement_failed', result['failed_checks'])
        broken_audit = copy.deepcopy(audit)
        broken_audit['features'][cds['feature_id']]['qualifiers']['product'] = ['unsupported invented label']
        self.assertIn('qualifier_preservation_failed',
            audit_annotation_preservation(self.record, document, broken_audit, extracted)['failed_checks'])

    def test_invalid_coordinates_are_failed_not_successful_results(self):
        self.record.features[0].location = self.Location(0,99,strand=1)
        self.prepare_record()
        report = self.run_pipeline()
        self.assertEqual(report['status'], 'blocked')
        self.assertEqual(report['parsed_records'], 0)
        self.assertEqual(report['scored_records'], 0)
        self.assertEqual(report['experiment_a']['overall']['rates']['feature_retention'], None)
        self.assertEqual(report['experiment_a']['failures_by_archetype'], {'single_cds': 1})

    def test_fuzzy_annotation_cannot_be_silently_exported_as_exact(self):
        from Bio.SeqFeature import BeforePosition
        self.record.features[0].location = self.Location(BeforePosition(0),12,strand=1)
        document, audit, extracted = normalize_record(self.record, self.row, self.row['checksum_sha256'])
        cds = next(f for f in document['features'] if f['feature_label']=='CDS')
        unit = next(u for u in document['units'] if u['feature_ids']==[cds['feature_id']])
        self.assertEqual(extracted, {})
        # A broken exporter must not turn an uncertain source span into a scored protein.
        cds['annotation_status'] = 'known'
        unit['model_input_status'] = 'pending_model_selection'
        result = audit_annotation_preservation(self.record, document, audit, {unit['unit_id']: 'M'+'A'*2})
        self.assertEqual(result['status'], 'failed')
        self.assertIn('unresolved_cds_exported', result['failed_checks'])

    def test_fasta_is_unobservable_not_perfect_annotation_evidence(self):
        report = self.run_pipeline()
        self.assertEqual(report['status'], 'ok')
        self.assertEqual(report['scored_records'], 0)
        self.assertEqual(report['experiment_a']['status'], 'no_observable_annotations')
        self.assertIsNone(report['experiment_a']['overall']['rates']['feature_retention'])
